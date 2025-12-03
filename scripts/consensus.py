#!/usr/bin/env python3
"""
BELT Consensus Module
Generate consensus sequences from grouped BAM files
"""

import os
import re
import sys
import subprocess
import tempfile
import logging
import shlex
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import pysam

# Set up logging
logger = logging.getLogger(__name__)


# Safe regex for sanitizing UMI names
SAFE = re.compile(r"[^A-Za-z0-9._-]")

def sanitize_umi_name(umi: str) -> str:
    """Sanitize UMI name for safe file naming."""
    return SAFE.sub("_", umi or "NA")

def write_group_fastq(umi: str, group_iter, fq_path: Path):
    """Write a group of reads to FASTQ format."""
    with open(fq_path, "w", buffering=1024*1024) as fq:
        for read in group_iter:
            read_name = read.query_name
            sequence = read.query_sequence or ""
            quality = read.qual if read.qual is not None else ("I" * len(sequence))
            fq.write(f"@{read_name}\n{sequence}\n+\n{quality}\n")


def umi_consensus_bcftools(umi: str, bam: Path, ref: Path, threads: int = 1, skip_extraction: bool = False, 
                          indel_sensitivity: float = 0.01) -> str:
    """
    Generate consensus sequence for one UMI group using bcftools.
    
    Args:
        umi: UMI identifier
        bam: Path to grouped BAM file
        ref: Path to reference FASTA
        threads: Number of threads for bcftools operations (default: 1)
        skip_extraction: If True, bam is already pre-extracted for this UMI (default: False)
        indel_sensitivity: Controls bcftools' indel calling heuristics (default: 0.01)
                          - Maps to `--gap-frac` parameter (minimum fraction of reads with gap to call indel)
                          - Smaller values (e.g., 0.05, 0.01) increase recall for large deletions
                          - Lower values may also increase spurious indel calls
                          - Recommended: 0.05-0.1 for balanced calling, 0.01-0.03 for maximum sensitivity
    """
    import sys
    
    # define a function to run a command with error capture
    def run(cmd, description=""):
        """Run a command and capture stderr for better error messages."""
        try:
            # Run command capturing output as bytes first to avoid Unicode issues
            result = subprocess.run(
                cmd, 
                shell=True, 
                check=True,
                capture_output=True,
                text=False  # Capture as bytes first
            )
            return result
        except subprocess.CalledProcessError as e:
            error_msg = f"Error running: {description or cmd}\n"
            error_msg += f"Exit code: {e.returncode}\n"
            # Try to decode output safely
            try:
                if e.stdout:
                    stdout_str = e.stdout.decode('utf-8', errors='replace')
                    if stdout_str.strip():  # Only include if non-empty
                        error_msg += f"STDOUT: {stdout_str}\n"
                if e.stderr:
                    stderr_str = e.stderr.decode('utf-8', errors='replace')
                    if stderr_str.strip():  # Only include if non-empty
                        error_msg += f"STDERR: {stderr_str}\n"
            except Exception as decode_err:
                # If decoding fails, just report the error without output
                error_msg += f"Could not decode command output (may contain binary data): {decode_err}\n"
            raise RuntimeError(error_msg) from e
    
    # Extract BAM for this UMI (skip if bam is already a pre-extracted UMI BAM)
    out_bam = f"umi.{umi}.bam"
    
    if skip_extraction and Path(bam).name.endswith('.bam'):
        # BAM is already pre-extracted for this UMI, just copy it and its index
        import shutil
        shutil.copy2(str(bam), out_bam)
        # Copy index if it exists
        bam_index = Path(str(bam) + '.bai')
        if bam_index.exists():
            shutil.copy2(str(bam_index), f"{out_bam}.bai")
        # Check read count
        count_result = subprocess.run(
            f"samtools view -c {shlex.quote(out_bam)}",
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        read_count = int(count_result.stdout.strip())
    else:
        # Extract BAM for this UMI from grouped BAM
        # Use samtools view with tag filtering - much faster than iterating entire BAM
        try:
            # Method 1: Use samtools view with grep (fastest for large BAMs)
            # Extract header and reads with matching tag
            # Pattern matches: header lines (^@) or lines with ZU:Z:UMI followed by tab/space/end
            bam_quoted = shlex.quote(str(bam))
            out_bam_quoted = shlex.quote(out_bam)
            # Escape special regex characters in UMI for grep -E
            import re
            umi_escaped = re.escape(umi)
            # More robust pattern: match ZU:Z:UMI followed by tab, space, or end of line
            # This handles cases where the tag is at the end of the line or followed by whitespace
            grep_pattern = f"^@|ZU:Z:{umi_escaped}([[:space:]]|$)"
            result = subprocess.run(
                f"samtools view -h {bam_quoted} | grep -E {shlex.quote(grep_pattern)} | samtools view -Sb -o {out_bam_quoted}",
                shell=True,
                check=True,
                capture_output=True,
                text=False
            )
            
            # Check if BAM was created and has reads
            if not Path(out_bam).exists():
                raise FileNotFoundError(f"BAM file not created: {out_bam}")
            
            # Quick check: count reads
            count_result = subprocess.run(
                f"samtools view -c {out_bam_quoted}",
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            read_count = int(count_result.stdout.strip())
            
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as e:
            # Fallback to pysam if samtools method fails
            logger.debug(f"samtools extraction failed, falling back to pysam: {e}")
            with pysam.AlignmentFile(bam, "rb") as bam_in, \
                 pysam.AlignmentFile(out_bam, "wb", header=bam_in.header) as bam_out:
                read_count = 0
                all_umis = set()  # Track all UMIs found for debugging
                for read in bam_in:
                    # Try multiple ways to get the UMI tag
                    read_umi = None
                    if read.has_tag("ZU"):
                        read_umi = read.get_tag("ZU")
                    elif "ZU" in dict(read.tags):
                        read_umi = dict(read.tags)["ZU"]
                    
                    if read_umi:
                        all_umis.add(read_umi)
                        if read_umi == umi:
                            bam_out.write(read)
                            read_count += 1
                
                # If no reads found, provide helpful error message
                if read_count == 0:
                    if all_umis:
                        raise ValueError(
                            f"No reads found for UMI '{umi}'. "
                            f"Found UMIs in BAM: {sorted(all_umis)[:10]}{'...' if len(all_umis) > 10 else ''}"
                        )
                    else:
                        raise ValueError(
                            f"No reads found for UMI '{umi}' and no UMI tags found in BAM file. "
                            f"BAM file may be empty or missing ZU tags."
                        )
    
    # Check if we found any reads for this UMI
    if read_count == 0:
        # Additional check: try to see what UMIs are in the BAM
        try:
            with pysam.AlignmentFile(bam, "rb") as bam_in:
                all_umis = set()
                for read in bam_in.fetch(until_eof=True):
                    if read.has_tag("ZU"):
                        all_umis.add(read.get_tag("ZU"))
                if all_umis:
                    raise ValueError(
                        f"No reads found for UMI '{umi}'. "
                        f"Found UMIs in BAM: {sorted(all_umis)[:10]}{'...' if len(all_umis) > 10 else ''}"
                    )
                else:
                    raise ValueError(
                        f"No reads found for UMI '{umi}' and no UMI tags found in BAM file. "
                        f"BAM file may be empty or missing ZU tags."
                    )
        except Exception:
            # If we can't inspect the BAM, just raise the original error
            raise ValueError(f"No reads found for UMI {umi}")
    
    # Index BAM (use threads if available, skip if already indexed)
    bai_file = Path(f"{out_bam}.bai")
    if not bai_file.exists():
        out_bam_quoted = shlex.quote(out_bam)
        run(f"samtools index -@ {threads} {out_bam_quoted}",
            f"indexing BAM for UMI {umi}")
    
    # Variant call
    # Note: bcftools threading support varies by version, so we don't use it here
    # The main speedup comes from parallelizing multiple UMIs, not threading individual commands
    # Options for better large deletion detection:
    # - -B: disable BAQ (Base Alignment Quality) which can filter out large indels
    # - -A: count anomalous read pairs (helps with large deletions)
    # - -a FORMAT/AD,FORMAT/DP: output AD and DP fields for better filtering
    # - -m: multiallelic caller (better for indels)
    # - -v: output variants only
    # - --gap-frac: minimum fraction of gapped reads supporting indel (configurable, lower = more sensitive to large dels)
    #   Note: If --gap-frac is not supported, bcftools will use default (typically 0.2)
    ref_quoted = shlex.quote(str(ref))
    out_bam_quoted = shlex.quote(out_bam)
    
    # Variant call with improved parameters for indel detection
    # Indel-related tuning based on indel_sensitivity
    # indel_sensitivity ~ minimum fraction of reads with a gap to consider an indel.
    # Example: 0.1 -> 10% of reads gapped; 0.05 -> 5%; 0.01 -> very permissive.
    gap_frac = max(1e-4, min(indel_sensitivity, 0.5))  # clamp to [1e-4, 0.5]
    
    # For per-UMI calling with relatively low depth, we want:
    # - relatively permissive gap fraction (from indel_sensitivity)
    # Note: Using well-supported bcftools parameters
    # -F/--gap-frac: minimum fraction of reads with gap (primary parameter for indel sensitivity)
    # This is the most important parameter for controlling indel calling sensitivity
    # Additional parameters may vary by bcftools version, so we use only the core parameter
    # Use --gap-frac (long form) as it's more explicit and less likely to conflict
    indel_opts = f"--gap-frac {gap_frac} "  # minimum fraction of reads with gap
    
    logger.debug(
        f"UMI {umi}: indel_sensitivity={indel_sensitivity}, "
        f"gap_frac={gap_frac}, indel_opts='{indel_opts.strip()}'"
    )
    
    vcf_output = f"umi.{umi}.vcf.gz"
    vcf_output_quoted = shlex.quote(vcf_output)
    
    # Build mpileup command with indel options
    mpileup_cmd = (
        f"bcftools mpileup -Ou "
        f"-B -A "                       # BAQ off, count anomalous pairs
        f"-d 10000 "                    # max depth
        f"-Q 0 -q 0 "                   # Remove quality filters
        f"{indel_opts}"
        f"-a FORMAT/AD,FORMAT/DP "      # per-sample depths
        f"-f {ref_quoted} {out_bam_quoted}"
    )
    
    # Build call command
    call_cmd = (
        f"bcftools call -mv "
        f"--ploidy 1 --keep-alts "
        f"-Oz -o {vcf_output_quoted}"
    )
    
    run(f"{mpileup_cmd} | {call_cmd}",
        f"variant calling for UMI {umi}")
    
    # Check if VCF was created and has content
    vcf_file = f"umi.{umi}.vcf.gz"
    if not Path(vcf_file).exists():
        raise FileNotFoundError(f"VCF file not created: {vcf_file}")
    
    # Check VCF file size (should be > 0)
    if Path(vcf_file).stat().st_size == 0:
        raise ValueError(f"VCF file is empty: {vcf_file}")
    
    # Index VCF (bcftools index -t creates tabix index, no threading option)
    run(f"bcftools index -t umi.{umi}.vcf.gz",
        f"indexing VCF for UMI {umi}")
    
    # Normalize VCF (bcftools norm doesn't support threading in most versions)
    # - -m -any: split multiallelic sites into multiple records (important for indels)
    # - -N: do not recalculate INFO fields (keeps original calling info)
    # Note: Normalization is important for proper consensus generation
    ref_quoted = shlex.quote(str(ref))
    run(f"bcftools norm -f {ref_quoted} -m -any -N umi.{umi}.vcf.gz -Oz -o umi.{umi}.norm.vcf.gz",
        f"normalizing VCF for UMI {umi}")
    
    # Check if normalized VCF was created
    norm_vcf = f"umi.{umi}.norm.vcf.gz"
    if not Path(norm_vcf).exists():
        raise FileNotFoundError(f"Normalized VCF file not created: {norm_vcf}")
    
    # Index normalized VCF
    run(f"bcftools index -t umi.{umi}.norm.vcf.gz",
        f"indexing normalized VCF for UMI {umi}")
    
    # Get sample name
    try:
        sample_name = subprocess.check_output(
            f"bcftools query -l umi.{umi}.norm.vcf.gz",
            shell=True,
            text=True
        ).strip()
        if not sample_name:
            raise ValueError(f"No sample name found in normalized VCF for UMI {umi}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to get sample name from VCF for UMI {umi}: {e}")
    
    # Generate consensus (bcftools consensus doesn't support threading)
    # - -H 1: use the first (and only) haplotype
    # Note: bcftools consensus applies variants from VCF to reference
    # Indels are included by default in consensus generation
    ref_quoted = shlex.quote(str(ref))
    sample_name_quoted = shlex.quote(sample_name)
    consensus_file = f"umi.{umi}.consensus.fa"
    run(f"bcftools consensus -f {ref_quoted} -s {sample_name_quoted} -H 1 -o {consensus_file} umi.{umi}.norm.vcf.gz",
        f"generating consensus for UMI {umi}")
    
    # Apply CIGAR-based large deletion patching as fallback
    # This ensures deletions visible in CIGAR strings are applied even if bcftools missed them
    # First, check what deletions bcftools already called to avoid double-application
    try:
        # Check normalized VCF for deletions that were already called
        # Use min_len=1 to catch all deletions, not just large ones
        vcf_deletions = summarize_large_deletions(Path(norm_vcf), min_len=1)
        vcf_del_positions = {(pos, del_size) for pos, _, _, del_size in vcf_deletions}
        logger.debug(
            f"VCF deletions found for UMI {umi}: {len(vcf_deletions)} deletions "
            f"at positions {[(pos, size) for pos, _, _, size in vcf_deletions[:5]]}"
        )
        
        # Use the extracted BAM if it exists, otherwise use the input BAM
        if skip_extraction:
            extracted_bam_path = Path(bam)
        else:
            extracted_bam_path = Path(out_bam)
        
        if extracted_bam_path.exists():
            large_dels = detect_large_deletions_from_cigar(
                extracted_bam_path,
                min_len=20,
                min_support_frac=0.5
            )
            
            # Filter out deletions that bcftools already called
            # Compare by position and length (within 5bp tolerance for coordinate differences)
            uncalled_dels = []
            for ev in large_dels:
                pos = ev['pos']
                length = ev['length']
                # Check if this deletion was already called by bcftools
                already_called = False
                for vcf_pos, vcf_del_size in vcf_del_positions:
                    # Allow 5bp tolerance for position and length differences
                    if abs(pos - vcf_pos) <= 5 and abs(length - vcf_del_size) <= 5:
                        already_called = True
                        break
                
                if not already_called:
                    uncalled_dels.append(ev)
            
            if uncalled_dels:
                logger.debug(
                    f"CIGAR-based large deletion detection for UMI {umi}: "
                    f"{len(uncalled_dels)} uncalled events (out of {len(large_dels)} total)"
                )
                for ev in uncalled_dels:
                    logger.debug(
                        f"  {ev['contig']}:{ev['pos']} len={ev['length']} "
                        f"support={ev['support']}/{ev['depth']} "
                        f"({ev['support_frac']:.2f})"
                    )
                apply_large_deletions_to_consensus(Path(consensus_file), uncalled_dels)
                logger.debug(f"Applied {len(uncalled_dels)} CIGAR-based large deletions to {consensus_file}")
            elif large_dels:
                logger.debug(
                    f"All {len(large_dels)} CIGAR-detected deletions for UMI {umi} "
                    f"were already called by bcftools"
                )
            
            # Also check for insertions that bcftools might have missed
            # First, check what insertions bcftools already called to avoid double-application
            vcf_insertions = summarize_large_insertions(Path(norm_vcf), min_len=1)
            vcf_ins_positions = {(pos, ins_size) for pos, _, _, ins_size in vcf_insertions}
            logger.debug(
                f"VCF insertions found for UMI {umi}: {len(vcf_insertions)} insertions "
                f"at positions {[(pos, size) for pos, _, _, size in vcf_insertions[:5]]}"
            )
            
            large_ins = detect_large_insertions_from_cigar(
                extracted_bam_path,
                min_len=5,  # Smaller threshold for insertions
                min_support_frac=0.5
            )
            
            # Filter out insertions that bcftools already called
            uncalled_ins = []
            for ev in large_ins:
                pos = ev['pos']
                length = ev['length']
                # Check if this insertion was already called by bcftools
                already_called = False
                for vcf_pos, vcf_ins_size in vcf_ins_positions:
                    # Allow 5bp tolerance for position and length differences
                    if abs(pos - vcf_pos) <= 5 and abs(length - vcf_ins_size) <= 5:
                        already_called = True
                        break
                
                if not already_called:
                    uncalled_ins.append(ev)
            
            if uncalled_ins:
                logger.debug(
                    f"CIGAR-based large insertion detection for UMI {umi}: "
                    f"{len(uncalled_ins)} uncalled events (out of {len(large_ins)} total)"
                )
                for ev in uncalled_ins:
                    logger.debug(
                        f"  {ev['contig']}:{ev['pos']} len={ev['length']} "
                        f"support={ev['support']}/{ev['depth']} "
                        f"({ev['support_frac']:.2f})"
                    )
                apply_large_insertions_to_consensus(Path(consensus_file), uncalled_ins)
                logger.debug(f"Applied {len(uncalled_ins)} CIGAR-based large insertions to {consensus_file}")
            elif large_ins:
                logger.debug(
                    f"All {len(large_ins)} CIGAR-detected insertions for UMI {umi} "
                    f"were already called by bcftools"
                )
    except Exception as e:
        logger.warning(
            f"Failed to apply CIGAR-based large indels for UMI {umi}: {e}"
        )
        # Don't fail the function if CIGAR patching fails
    
    return consensus_file    


def debug_mpileup_for_umi(umi_bam: Path, ref: Path, region: str, indel_sensitivity: float = 0.01) -> Path:
    """
    Debug helper to dump raw mpileup output for a specific region.
    
    This function helps inspect whether large deletions are present in the mpileup
    stream before bcftools call processes them.
    
    Args:
        umi_bam: Path to UMI-specific BAM file
        ref: Path to reference FASTA
        region: Genomic region (e.g., "chr1:1000-2000")
        indel_sensitivity: Indel sensitivity parameter (default: 0.05)
    
    Returns:
        Path to temporary text file with mpileup output
    """
    import tempfile
    
    # Compute indel options (same as umi_consensus_bcftools)
    gap_frac = max(1e-4, min(indel_sensitivity, 0.5))
    indel_opts = (
        f"-F {gap_frac} "
        f"-m 1 "
        f"-L 100000 "
        f"--indel-size 1000 "
        f"--indel-bias 2.0 "
    )
    
    # Create temporary output file
    tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
    tmp_path = Path(tmp_file.name)
    tmp_file.close()
    
    ref_quoted = shlex.quote(str(ref))
    bam_quoted = shlex.quote(str(umi_bam))
    
    mpileup_cmd = (
        f"bcftools mpileup -Ou "
        f"-B -A -d 10000 "
        f"-Q 0 -q 0 "
        f"{indel_opts}"
        f"-a AD,DP,IDV,IMF "
        f"-f {ref_quoted} -r {region} {bam_quoted}"
    )
    
    # Pipe through bcftools view to get text output
    subprocess.run(
        f"{mpileup_cmd} | bcftools view -Ov > {tmp_path}",
        shell=True,
        check=True,
        capture_output=True
    )
    
    return tmp_path


def summarize_large_deletions(vcf_path: Path, min_len: int = 20) -> list:
    """
    Summarize large deletions found in a VCF file.
    
    Args:
        vcf_path: Path to VCF file (can be .vcf, .vcf.gz, or .bcf)
        min_len: Minimum deletion length to report (default: 20bp)
    
    Returns:
        List of tuples: (position, ref_length, alt_length, deletion_size)
    """
    import pysam
    
    deletions = []
    
    try:
        with pysam.VariantFile(str(vcf_path)) as vcf:
            for record in vcf:
                # Check if it's a deletion
                # Simple deletion: ref longer than alt
                if len(record.ref) > len(record.alts[0]):
                    del_size = len(record.ref) - len(record.alts[0])
                    if del_size >= min_len:
                        deletions.append((
                            record.pos,
                            len(record.ref),
                            len(record.alts[0]),
                            del_size
                        ))
                # Structural variant deletion
                elif '<DEL>' in str(record.alts[0]).upper():
                    # Try to get size from INFO field
                    if 'SVLEN' in record.info:
                        svlen = abs(record.info['SVLEN'])
                        if isinstance(svlen, (list, tuple)):
                            svlen = svlen[0]
                        if svlen >= min_len:
                            deletions.append((
                                record.pos,
                                len(record.ref),
                                0,  # Unknown alt length
                                svlen
                            ))
    except Exception as e:
        logger.warning(f"Could not read VCF file {vcf_path}: {e}")
    
    return deletions


def summarize_large_insertions(vcf_path: Path, min_len: int = 5) -> list:
    """
    Summarize large insertions found in a VCF file.
    
    Args:
        vcf_path: Path to VCF file (can be .vcf, .vcf.gz, or .bcf)
        min_len: Minimum insertion length to report (default: 5bp)
    
    Returns:
        List of tuples: (position, ref_length, alt_length, insertion_size)
    """
    import pysam
    
    insertions = []
    
    try:
        with pysam.VariantFile(str(vcf_path)) as vcf:
            for record in vcf:
                # Check if it's an insertion
                # Simple insertion: alt longer than ref
                if len(record.alts[0]) > len(record.ref):
                    ins_size = len(record.alts[0]) - len(record.ref)
                    if ins_size >= min_len:
                        insertions.append((
                            record.pos,
                            len(record.ref),
                            len(record.alts[0]),
                            ins_size
                        ))
                # Structural variant insertion
                elif '<INS>' in str(record.alts[0]).upper() or '<INS:ME>' in str(record.alts[0]).upper():
                    # Try to get size from INFO field
                    if 'SVLEN' in record.info:
                        svlen = abs(record.info['SVLEN'])
                        if isinstance(svlen, (list, tuple)):
                            svlen = svlen[0]
                        if svlen >= min_len:
                            insertions.append((
                                record.pos,
                                len(record.ref),
                                0,  # Unknown alt length
                                svlen
                            ))
    except Exception as e:
        logger.warning(f"Could not read VCF file {vcf_path}: {e}")
    
    return insertions


def detect_large_deletions_from_cigar(
    umi_bam: Path,
    min_len: int = 20,
    min_support_frac: float = 0.5,
) -> list:
    """
    Detect large deletions directly from CIGAR strings.

    Args:
        umi_bam: Path to BAM for a single UMI.
        min_len: Minimum deletion length to consider (bp).
        min_support_frac: Minimum fraction of reads supporting the deletion.

    Returns:
        List of events: dicts with keys:
            {
                "contig": str,
                "pos": int,        # 1-based reference coordinate of deletion start
                "length": int,     # deletion length in bp
                "support": int,    # number of supporting reads
                "depth": int,      # total reads used
                "support_frac": float
            }
    """
    import collections

    events = collections.Counter()
    contig_for_event = {}
    total_reads = 0

    with pysam.AlignmentFile(str(umi_bam), "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_reads += 1
            ref_pos = read.reference_start  # 0-based
            contig = read.reference_name

            if read.cigartuples is None:
                continue

            for op, length in read.cigartuples:
                # pysam CIGAR op codes:
                # 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, ...
                if op == 2 and length >= min_len:
                    # record deletion event at this reference position
                    key = (contig, ref_pos, length)
                    events[key] += 1
                    contig_for_event[key] = contig
                # advance reference position for ops that consume reference
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                    ref_pos += length

    if total_reads == 0:
        return []

    large_dels = []
    for (contig, ref_start0, length), support in events.items():
        support_frac = support / total_reads
        if support_frac >= min_support_frac:
            large_dels.append({
                "contig": contig,
                "pos": ref_start0 + 1,  # convert to 1-based
                "length": length,
                "support": support,
                "depth": total_reads,
                "support_frac": support_frac,
            })

    # sort by position for reproducibility
    large_dels.sort(key=lambda e: (e["contig"], e["pos"], e["length"]))

    logger.debug(
        f"CIGAR-based large deletion detection in {umi_bam}: "
        f"{len(large_dels)} events (min_len={min_len}, "
        f"min_support_frac={min_support_frac}, total_reads={total_reads})"
    )

    return large_dels


def detect_large_insertions_from_cigar(
    umi_bam: Path,
    min_len: int = 5,
    min_support_frac: float = 0.5,
) -> list:
    """
    Detect large insertions directly from CIGAR strings.

    Args:
        umi_bam: Path to BAM for a single UMI.
        min_len: Minimum insertion length to consider (bp).
        min_support_frac: Minimum fraction of reads supporting the insertion.

    Returns:
        List of events: dicts with keys:
            {
                "contig": str,
                "pos": int,        # 1-based reference coordinate after insertion
                "length": int,     # insertion length in bp
                "inserted_seq": str, # inserted sequence (if available)
                "support": int,    # number of supporting reads
                "depth": int,      # total reads used
                "support_frac": float
            }
    """
    import collections

    events = collections.Counter()
    contig_for_event = {}
    total_reads = 0

    with pysam.AlignmentFile(str(umi_bam), "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_reads += 1
            ref_pos = read.reference_start  # 0-based
            contig = read.reference_name

            if read.cigartuples is None:
                continue

            for op, length in read.cigartuples:
                # pysam CIGAR op codes:
                # 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, ...
                if op == 1 and length >= min_len:  # Insertion
                    # record insertion event at this reference position
                    # Insertion is after the current reference position
                    key = (contig, ref_pos, length)
                    events[key] += 1
                    contig_for_event[key] = contig
                # advance reference position for ops that consume reference
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                    ref_pos += length

    if total_reads == 0:
        return []

    large_ins = []
    for (contig, ref_start0, length), support in events.items():
        support_frac = support / total_reads
        if support_frac >= min_support_frac:
            large_ins.append({
                "contig": contig,
                "pos": ref_start0 + 1,  # convert to 1-based (insertion is after this position)
                "length": length,
                "support": support,
                "depth": total_reads,
                "support_frac": support_frac,
            })

    # sort by position for reproducibility
    large_ins.sort(key=lambda e: (e["contig"], e["pos"], e["length"]))

    logger.debug(
        f"CIGAR-based large insertion detection in {umi_bam}: "
        f"{len(large_ins)} events (min_len={min_len}, "
        f"min_support_frac={min_support_frac}, total_reads={total_reads})"
    )

    return large_ins


def apply_large_deletions_to_consensus(
    consensus_fasta: Path,
    large_dels: list,
) -> None:
    """
    Apply large deletions to a consensus FASTA sequence in-place.

    Assumes:
        - consensus FASTA contains one or more contigs,
        - coordinates in large_dels are 1-based on the corresponding contig,
        - consensus currently has no large deletions (i.e., matches reference length).

    Args:
        consensus_fasta: Path to consensus FASTA file (will be modified in-place).
        large_dels: List of events from detect_large_deletions_from_cigar.
    """
    if not large_dels:
        return

    # Load FASTA into memory (small per-UMI consensus, so this is safe)
    contigs = {}
    order = []
    with open(consensus_fasta, "r") as fh:
        name = None
        seq_chunks = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(seq_chunks)
                    order.append(name)
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            contigs[name] = "".join(seq_chunks)
            order.append(name)

    if not contigs:
        logger.warning(f"No contigs found in consensus FASTA {consensus_fasta}")
        return

    # Group deletions by contig
    dels_by_contig = {}
    for ev in large_dels:
        c = ev["contig"]
        dels_by_contig.setdefault(c, []).append(ev)

    # Apply deletions for each contig independently
    for contig, events in dels_by_contig.items():
        if contig not in contigs:
            logger.warning(
                f"Contig {contig} from large deletion events not found "
                f"in consensus FASTA {consensus_fasta}"
            )
            continue

        seq = list(contigs[contig])
        # Sort deletions descending by start so earlier indices are unaffected
        events_sorted = sorted(events, key=lambda e: e["pos"], reverse=True)

        for ev in events_sorted:
            start = ev["pos"]          # 1-based
            length = ev["length"]
            start_idx = start - 1      # 0-based
            end_idx = start_idx + length

            if start_idx < 0 or start_idx >= len(seq):
                logger.warning(
                    f"Large deletion start {ev} out of bounds for contig {contig} "
                    f"of length {len(seq)}"
                )
                continue

            end_idx = min(end_idx, len(seq))
            logger.debug(
                f"Applying large deletion to {contig}: "
                f"pos={ev['pos']}, length={length}, "
                f"support_frac={ev['support_frac']:.3f}"
            )
            del seq[start_idx:end_idx]

        contigs[contig] = "".join(seq)

    # Write updated FASTA back to disk
    with open(consensus_fasta, "w") as fh:
        for name in order:
            fh.write(f">{name}\n")
            seq = contigs[name]
            # wrap at 80 columns for readability
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    logger.debug(
        f"Applied CIGAR-based large deletions to consensus FASTA {consensus_fasta}"
    )


def apply_large_insertions_to_consensus(
    consensus_fasta: Path,
    large_ins: list,
) -> None:
    """
    Apply large insertions to a consensus FASTA sequence in-place.

    Args:
        consensus_fasta: Path to consensus FASTA file (will be modified in-place).
        large_ins: List of events from detect_large_insertions_from_cigar.
    """
    if not large_ins:
        return

    # Load FASTA into memory
    contigs = {}
    order = []
    with open(consensus_fasta, "r") as fh:
        name = None
        seq_chunks = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(seq_chunks)
                    order.append(name)
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            contigs[name] = "".join(seq_chunks)
            order.append(name)

    if not contigs:
        logger.warning(f"No contigs found in consensus FASTA {consensus_fasta}")
        return

    # Group insertions by contig
    ins_by_contig = {}
    for ev in large_ins:
        c = ev["contig"]
        ins_by_contig.setdefault(c, []).append(ev)

    # Apply insertions for each contig independently
    # Sort ascending by position so later insertions don't affect earlier coordinates
    for contig, events in ins_by_contig.items():
        if contig not in contigs:
            logger.warning(
                f"Contig {contig} from large insertion events not found "
                f"in consensus FASTA {consensus_fasta}"
            )
            continue

        seq = list(contigs[contig])
        # Sort insertions ascending by position
        events_sorted = sorted(events, key=lambda e: e["pos"])

        for ev in events_sorted:
            pos = ev["pos"]          # 1-based, insertion goes after this position
            length = ev["length"]
            insert_idx = pos          # 0-based index where to insert (after this position)

            if insert_idx < 0 or insert_idx > len(seq):
                logger.warning(
                    f"Large insertion position {ev} out of bounds for contig {contig} "
                    f"of length {len(seq)}"
                )
                continue

            # Insert placeholder sequence (N's) - we don't have the actual inserted sequence
            # This is a limitation, but for length checking it's sufficient
            inserted_seq = "N" * length
            logger.debug(
                f"Applying large insertion to {contig}: "
                f"pos={ev['pos']}, length={length}, "
                f"support_frac={ev['support_frac']:.3f}"
            )
            seq.insert(insert_idx, inserted_seq)

        contigs[contig] = "".join(seq)

    # Write updated FASTA back to disk
    with open(consensus_fasta, "w") as fh:
        for name in order:
            fh.write(f">{name}\n")
            seq = contigs[name]
            # wrap at 80 columns for readability
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    logger.debug(
        f"Applied CIGAR-based large insertions to consensus FASTA {consensus_fasta}"
    )





def run_consensus_for_umi_bcftools(grouped_bam: Path, ref_fa: Path, meta_dir: Path,
                                   con_dir: Path, bam_dir: Path, umi_threads: int,
                                   model: str, log_file: Path, umi_safe: str, umi: str,
                                   pre_extracted_bam: Path = None, indel_sensitivity: float = 0.1) -> str:
    """
    Generate consensus sequence for one UMI group using bcftools.
    
    This function uses umi_consensus_bcftools to process a single UMI from a grouped BAM.
    
    Args:
        grouped_bam: Path to grouped BAM file with ZU tags (used if pre_extracted_bam is None)
        ref_fa: Path to reference FASTA
        meta_dir: Directory for intermediate files
        con_dir: Directory for final consensus output
        bam_dir: Directory for BAM files
        umi_threads: Number of threads to use (not used by bcftools directly)
        model: Model name (not used by bcftools)
        log_file: Path to log file
        umi_safe: Sanitized UMI name for output files
        umi: Original UMI string
        pre_extracted_bam: Optional pre-extracted BAM file for this UMI (much faster)
        indel_sensitivity: Minimum fraction of reads supporting an indel (default: 0.1)
                          Lower values (e.g., 0.05) help detect large deletions
    
    Returns:
        Path to generated consensus FASTA file
    """
    import os
    import shutil
    
    logfh = open(log_file, "a") if log_file else subprocess.DEVNULL
    
    # Create directories
    meta_dir.mkdir(parents=True, exist_ok=True)
    con_dir.mkdir(parents=True, exist_ok=True)
    bam_dir.mkdir(parents=True, exist_ok=True)
    
    original_cwd = os.getcwd()
    consensus_fasta = (con_dir / f"{umi_safe}.fasta").resolve()
    
    # Convert paths to absolute paths before changing directories
    grouped_bam = Path(grouped_bam).resolve()
    ref_fa = Path(ref_fa).resolve()
    con_dir = Path(con_dir).resolve()
    bam_dir = Path(bam_dir).resolve()
    meta_dir = Path(meta_dir).resolve()
    
    try:
        if logfh is not subprocess.DEVNULL:
            logfh.write(f"✅ Processing UMI {umi_safe} with bcftools...\n")
            logfh.write(f"   Using grouped BAM: {grouped_bam}\n")
            logfh.write(f"   Using reference: {ref_fa}\n")
            logfh.flush()
        
        # Verify BAM file exists
        if not grouped_bam.exists():
            raise FileNotFoundError(f"Grouped BAM file not found: {grouped_bam}")
        
        # Change to meta_dir for bcftools processing (it uses relative paths)
        os.chdir(meta_dir)
        
        # Use pre-extracted BAM if available, otherwise use grouped BAM
        # Convert pre_extracted_bam to Path if it's not None
        if pre_extracted_bam is not None:
            pre_extracted_bam = Path(pre_extracted_bam).resolve()
            bam_to_use = pre_extracted_bam if pre_extracted_bam.exists() else grouped_bam
            skip_extraction = pre_extracted_bam.exists()
        else:
            bam_to_use = grouped_bam
            skip_extraction = False
        
        # Call umi_consensus_bcftools with absolute paths and threading
        result_path = umi_consensus_bcftools(umi, bam_to_use, ref_fa, threads=umi_threads, 
                                            skip_extraction=skip_extraction, indel_sensitivity=indel_sensitivity)
        
        # Convert result_path to absolute path (it's relative to meta_dir)
        result_path_abs = (Path(meta_dir) / result_path).resolve()
        
        # Ensure destination directory exists
        consensus_fasta.parent.mkdir(parents=True, exist_ok=True)
        
        # Move consensus file to final location
        if result_path_abs.exists():
            shutil.move(str(result_path_abs), str(consensus_fasta))
            
            # Update FASTA header with sanitized UMI name
            with open(consensus_fasta, 'r') as f:
                content = f.read()
            
            # Replace header with sanitized UMI name
            lines = content.split('\n')
            if lines and lines[0].startswith('>'):
                lines[0] = f">{umi_safe}"
            
            with open(consensus_fasta, 'w') as f:
                f.write('\n'.join(lines))
            
            # Copy extracted BAM to final location if it exists
            extracted_bam = meta_dir / f"umi.{umi}.bam"
            if extracted_bam.exists():
                final_bam = bam_dir / f"{umi_safe}.bam"
                shutil.copy2(extracted_bam, final_bam)
                # Copy index if it exists
                bai_file = extracted_bam.with_suffix('.bam.bai')
                if bai_file.exists():
                    shutil.copy2(bai_file, final_bam.with_suffix('.bam.bai'))
            
            # CIGAR-based large deletion patching on the consensus
            try:
                if extracted_bam.exists():
                    large_dels = detect_large_deletions_from_cigar(
                        extracted_bam,
                        min_len=20,          # adjust as needed
                        min_support_frac=0.5 # e.g. >=50% of reads
                    )
                    if large_dels:
                        if logfh is not subprocess.DEVNULL:
                            logfh.write(
                                f"   CIGAR-based large deletions for {umi_safe}: "
                                f"{len(large_dels)} events\n"
                            )
                            for ev in large_dels:
                                logfh.write(
                                    f"      {ev['contig']}:{ev['pos']} "
                                    f"len={ev['length']} "
                                    f"support={ev['support']}/{ev['depth']} "
                                    f"({ev['support_frac']:.2f})\n"
                                )
                            logfh.flush()
                        apply_large_deletions_to_consensus(consensus_fasta, large_dels)
                else:
                    logger.debug(
                        f"No extracted BAM found for {umi_safe}; "
                        f"skipping CIGAR-based large deletion patching."
                    )
            except Exception as e:
                logger.warning(
                    f"Failed to apply CIGAR-based large deletions for {umi_safe}: {e}"
                )
        
        if logfh is not subprocess.DEVNULL:
            logfh.write(f"🎉 Consensus FASTA created: {consensus_fasta}\n")
            logfh.flush()
            
    except Exception as e:
        if logfh is not subprocess.DEVNULL:
            logfh.write(f"❌ Error processing {umi_safe}: {e}\n")
            logfh.flush()
        raise
    finally:
        os.chdir(original_cwd)
        if logfh is not subprocess.DEVNULL:
            logfh.close()
    
    return str(consensus_fasta)

def run_consensus_for_umi(reads_fq: Path, ref_fa: Path, meta_dir: Path, 
                         con_dir: Path, bam_dir: Path, umi_threads: int, 
                         model: str, log_file: Path, umi_safe: str) -> str:
    """
    Generate consensus sequence for one UMI group (default method - uses bcftools).
    
    This function is a wrapper that uses bcftools for consensus calling.
    It requires the grouped BAM file to be passed through the workflow.
    
    Args:
        reads_fq: Path to FASTQ file with reads for this UMI (not used for bcftools)
        ref_fa: Path to reference FASTA
        meta_dir: Directory for intermediate files
        con_dir: Directory for final consensus output
        bam_dir: Directory for BAM files
        umi_threads: Number of threads to use
        model: Model name (not used by bcftools)
        log_file: Path to log file
        umi_safe: Sanitized UMI name for output files
    
    Returns:
        Path to generated consensus FASTA file
    """
    raise NotImplementedError(
        "Default consensus method not implemented. Use --consensus-method bcftools"
    ) 

def check_dependencies():
    """Check if required external tools are available."""
    import shutil
    required_tools = ['samtools', 'bcftools']
    missing_tools = []
    
    for tool in required_tools:
        # Use shutil.which to check if tool is in PATH
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    
    if missing_tools:
        raise RuntimeError(
            f"Missing required tools: {', '.join(missing_tools)}\n"
            "Please install: samtools, bcftools"
        )

def process_grouped_bam(grouped_bam: str, reference: str, output_dir: str, 
                       sample: str, tag: str = "ZU", threads: int = 2, 
                       procs: int = None, umi_threads: int = 2, 
                       model: str = "r1041_e82_400bps_sup_variant_v5.0.0",
                       tmp_root: str = None, manifest: str = None, 
                       log_file: str = None, consensus_method: str = "default",
                       indel_sensitivity: float = 0.1) -> list:
    """
    Process grouped BAM file to generate consensus sequences.
    
    Args:
        grouped_bam: Path to grouped BAM file
        reference: Path to reference genome FASTA
        output_dir: Output directory
        sample: Sample name
        tag: UMI tag name (default: "ZU")
        threads: Total threads (not used directly)
        procs: Number of parallel UMI processes
        umi_threads: Threads per UMI for consensus calling
        model: Model name (not used by bcftools, kept for compatibility)
        tmp_root: Temporary directory root
        manifest: Output manifest file path
        log_file: Log file path
        consensus_method: Consensus calling method ("bcftools" or "default")
        indel_sensitivity: Minimum fraction of reads supporting an indel (default: 0.1)
                          Lower values (e.g., 0.05) improve detection of large deletions
                          that may have fewer supporting reads
    
    Returns:
        List of produced consensus FASTA file paths
    """
    # Check dependencies based on method
    if consensus_method == "bcftools":
        # Bcftools method only requires bcftools and samtools
        import shutil
        required_tools = ['bcftools', 'samtools']
        missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing_tools:
            raise RuntimeError(
                f"Missing required tools for bcftools method: {', '.join(missing_tools)}\n"
                "Please install: bcftools, samtools"
            )
    else:
        # Default method - check all dependencies
        check_dependencies()
    
    # Auto-detect number of processes if not specified
    if procs is None:
        import os
        cpu_count = os.cpu_count() or 4
        # Use most of available CPUs, but leave some headroom
        # Reserve 1-2 CPUs for system, divide rest between parallel UMIs and per-UMI threads
        procs = max(1, cpu_count // max(umi_threads, 2))
        logger.info(f"Auto-detected {procs} parallel processes (CPU count: {cpu_count}, umi_threads: {umi_threads})")
    
    # Setup directories - Option 4 structure: bam/ and fasta/
    out_root = Path(output_dir) / sample
    con_dir = out_root / "fasta"          # Consensus FASTA files
    bam_dir = out_root / "bam"            # Aligned BAM files
    tmp_root = Path(tmp_root) if tmp_root else (out_root / "tmp")
    
    con_dir.mkdir(parents=True, exist_ok=True)
    bam_dir.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)
    
    log_file = Path(log_file) if log_file else (out_root / "runner.log")
    
    produced = []
    futures = []
    
    def get_umi_paths(umi_safe: str):
        """Get paths for UMI-specific files."""
        umi_root = tmp_root / umi_safe
        reads_fq = umi_root / f"{umi_safe}_reads.fastq"  # matches consensus.sh UMI parsing
        meta_dir = umi_root / "meta"                     # safe for cleanup inside consensus.sh
        meta_dir.mkdir(parents=True, exist_ok=True)
        return reads_fq, meta_dir
    
    # Pre-extract all UMI BAMs in one pass (MAJOR SPEEDUP for bcftools method)
    # This eliminates the N×BAM_reads overhead of extracting each UMI separately
    umi_list = []
    umi_bam_files = {}  # Map UMI -> pre-extracted BAM path
    
    # Ensure grouped_bam is an absolute path
    grouped_bam = Path(grouped_bam).resolve()
    if not grouped_bam.exists():
        raise FileNotFoundError(f"Grouped BAM file not found: {grouped_bam}")
    
    if consensus_method == "bcftools":
        # Pre-extract all UMI BAMs in a single pass through the grouped BAM
        logger.info(f"Pre-extracting UMI BAMs from {grouped_bam}...")
        
        # Create a temporary directory for pre-extracted BAMs
        pre_extract_dir = tmp_root / "pre_extracted_bams"
        pre_extract_dir.mkdir(parents=True, exist_ok=True)
        
        # Two-pass approach to avoid "too many open files" error
        # First pass: identify all UMIs
        # Second pass: write reads to BAM files in batches
        umi_counts = {}   # Map UMI -> read count
        all_umis = set()
        
        try:
            # First pass: identify all unique UMIs
            logger.info("Identifying UMIs...")
            with pysam.AlignmentFile(grouped_bam, "rb") as bam_in:
                header = bam_in.header
                
                for read in bam_in:
                    try:
                        # Safely get UMI tag
                        umi = read.get_tag(tag) if read.has_tag(tag) else None
                        if umi:
                            all_umis.add(umi)
                            umi_counts[umi] = umi_counts.get(umi, 0) + 1
                    except (KeyError, ValueError) as e:
                        # Skip reads without the tag or with invalid tags
                        logger.debug(f"Skipping read without {tag} tag: {e}")
                        continue
            
            # Second pass: write reads to BAM files in batches to avoid too many open files
            logger.info(f"Writing {len(all_umis)} UMI BAMs in batches (max {50} files open at once)...")
            max_open_files = 50  # Limit concurrent open files
            umi_list_sorted = sorted(all_umis)  # Process in consistent order
            
            # Process in batches
            for batch_start in range(0, len(umi_list_sorted), max_open_files):
                batch_end = min(batch_start + max_open_files, len(umi_list_sorted))
                batch_umis = umi_list_sorted[batch_start:batch_end]
                
                # Open writers for this batch
                batch_writers = {}
                for umi in batch_umis:
                    umi_safe = sanitize_umi_name(umi)
                    umi_bam = pre_extract_dir / f"{umi_safe}.bam"
                    batch_writers[umi] = pysam.AlignmentFile(str(umi_bam), "wb", header=header)
                
                # Second pass through BAM: write reads for this batch
                with pysam.AlignmentFile(grouped_bam, "rb") as bam_in:
                    for read in bam_in:
                        try:
                            umi = read.get_tag(tag) if read.has_tag(tag) else None
                            if umi and umi in batch_writers:
                                batch_writers[umi].write(read)
                        except (KeyError, ValueError):
                            continue
                
                # Close all writers in this batch
                for writer in batch_writers.values():
                    writer.close()
                
                logger.debug(f"Wrote batch {batch_start//max_open_files + 1} ({len(batch_umis)} UMIs)")
            
            # Build umi_list and index pre-extracted BAMs in parallel
            logger.info(f"Indexing {len(umi_counts)} pre-extracted UMI BAMs...")
            
            # Index BAMs in parallel using ThreadPoolExecutor (I/O-bound operation)
            def index_bam(umi_bam_path):
                subprocess.run(
                    ['samtools', 'index', str(umi_bam_path)],
                    check=False,
                    capture_output=True
                )
            
            indexing_tasks = []
            for umi, count in umi_counts.items():
                if count > 0:  # Only process UMIs with reads
                    umi_safe = sanitize_umi_name(umi)
                    umi_bam = pre_extract_dir / f"{umi_safe}.bam"
                    umi_bam_files[umi] = umi_bam
                    indexing_tasks.append(umi_bam)
                    reads_fq, meta_dir = get_umi_paths(umi_safe)
                    reads_fq.parent.mkdir(parents=True, exist_ok=True)
                    umi_list.append((umi, umi_safe, reads_fq, meta_dir, None))
            
            # Index all BAMs in parallel
            with ThreadPoolExecutor(max_workers=min(procs * 2, len(indexing_tasks))) as index_pool:
                list(index_pool.map(index_bam, indexing_tasks))
            
            logger.info(f"Pre-extracted and indexed {len(umi_list)} UMI BAMs")
            
            if len(umi_list) == 0:
                logger.warning(f"No UMIs found in grouped BAM {grouped_bam}")
                return []
            
        except Exception as e:
            logger.warning(f"Pre-extraction failed, falling back to per-UMI extraction: {e}")
            import traceback
            logger.debug(f"Pre-extraction traceback: {traceback.format_exc()}")
            # Fallback: collect UMI list without pre-extraction
            try:
                with pysam.AlignmentFile(grouped_bam, "rb") as bam:
                    seen_umis = set()
                    for read in bam.fetch(until_eof=True):
                        try:
                            if read.has_tag(tag):
                                umi = read.get_tag(tag)
                                if umi and umi not in seen_umis:
                                    seen_umis.add(umi)
                                    umi_safe = sanitize_umi_name(umi)
                                    reads_fq, meta_dir = get_umi_paths(umi_safe)
                                    reads_fq.parent.mkdir(parents=True, exist_ok=True)
                                    umi_list.append((umi, umi_safe, reads_fq, meta_dir, None))
                        except (KeyError, ValueError):
                            continue
            except Exception as fallback_error:
                logger.error(f"Fallback UMI collection also failed: {fallback_error}")
                raise
    
    # Check if we have any UMIs to process
    if len(umi_list) == 0:
        logger.warning(f"No UMIs found in grouped BAM {grouped_bam}")
        return []
    
    logger.info(f"Processing {len(umi_list)} UMIs with {procs} parallel processes")
    
    # Process grouped BAM file with parallel execution
    # Increase max_workers for better parallelization
    # Use more aggressive parallelism: allow more queued tasks
    max_queue_size = max(procs * 8, 50)  # Allow more queued tasks
    
    with ProcessPoolExecutor(max_workers=procs) as pool:
        for umi, umi_safe, reads_fq, meta_dir, group_iter in umi_list:
            # Use bcftools method (only supported method)
            # For bcftools, use pre-extracted BAM if available
            pre_extracted = umi_bam_files.get(umi) if umi_bam_files else None
            # Convert Path to absolute path string for ProcessPoolExecutor serialization
            if pre_extracted is not None:
                pre_extracted = Path(pre_extracted).resolve()
            consensus_func = run_consensus_for_umi_bcftools
            fut = pool.submit(
                consensus_func,
                Path(grouped_bam).resolve(), Path(reference).resolve(), 
                Path(meta_dir).resolve(), Path(con_dir).resolve(),
                Path(bam_dir).resolve(), umi_threads, model, 
                Path(log_file).resolve() if log_file else None, umi_safe, umi,
                pre_extracted, indel_sensitivity
            )
            
            futures.append(fut)
            
            # Throttle queued tasks to prevent memory issues (more aggressive)
            # Process completed tasks as they come in
            if len(futures) >= max_queue_size:
                # Process all completed tasks
                completed = []
                remaining = []
                for fut in futures:
                    if fut.done():
                        try:
                            completed.append(fut.result())
                        except Exception as e:
                            logger.error(f"Error processing UMI: {e}")
                    else:
                        remaining.append(fut)
                produced.extend(completed)
                futures = remaining
        
        # Process remaining futures
        for done in as_completed(futures):
            try:
                produced.append(done.result())
            except Exception as e:
                logger.error(f"Error processing UMI: {e}")
    
    # Write manifest file if requested
    if manifest:
        with open(manifest, "w") as fh:
            for path in produced:
                fh.write(path + "\n")
    
    return produced

def main():
    """Main function for standalone execution."""
    import argparse
    
    ap = argparse.ArgumentParser(description="Generate consensus sequences from grouped BAM file")
    ap.add_argument("--grouped-bam", required=True, help="Input grouped BAM file")
    ap.add_argument("--tag", default="ZU", help="UMI tag name")
    ap.add_argument("--reference", required=True, help="Reference genome fasta file")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--sample", required=True, help="Sample name")
    ap.add_argument("--threads", type=int, default=2, help="Total threads (not used directly)")
    ap.add_argument("--procs", type=int, default=None, 
                    help="Parallel UMI processes (default: auto-detect from CPU count)")
    ap.add_argument("--umi-threads", type=int, default=2, 
                    help="Threads per UMI for consensus calling (default: 2)")
    ap.add_argument("--model", default="r1041_e82_400bps_sup_variant_v5.0.0", 
                    help="Model name (not used by bcftools, kept for compatibility)")
    ap.add_argument("--tmp-root", help="Temporary directory root")
    ap.add_argument("--manifest", help="Output manifest file")
    ap.add_argument("--log-file", help="Log file path")
    ap.add_argument("--consensus-method", default="bcftools", choices=["bcftools"],
                    help="Consensus calling method (default: 'bcftools')")
    ap.add_argument("--indel-sensitivity", type=float, default=0.1,
                    help="Minimum fraction of reads supporting an indel (default: 0.1). "
                         "Lower values (e.g., 0.05) improve detection of large deletions")
    ap.add_argument("--debug-umi-bam", help="Optional: path to a single UMI BAM to debug large deletions")
    
    args = ap.parse_args()
    
    if args.debug_umi_bam:
        bam_path = Path(args.debug_umi_bam)
        large_dels = detect_large_deletions_from_cigar(
            bam_path,
            min_len=20,
            min_support_frac=0.5,
        )
        print(f"Detected {len(large_dels)} large deletions in {bam_path}")
        for ev in large_dels:
            print(
                f"{ev['contig']}:{ev['pos']} len={ev['length']} "
                f"support={ev['support']}/{ev['depth']} "
                f"({ev['support_frac']:.2f})"
            )
        sys.exit(0)
    
    try:
        produced = process_grouped_bam(
            grouped_bam=args.grouped_bam,
            reference=args.reference,
            output_dir=args.outdir,
            sample=args.sample,
            tag=args.tag,
            threads=args.threads,
            procs=args.procs,
            umi_threads=args.umi_threads,
            model=args.model,
            tmp_root=args.tmp_root,
            manifest=args.manifest,
            log_file=args.log_file,
            consensus_method=args.consensus_method,
            indel_sensitivity=args.indel_sensitivity
        )
        
        # Print produced files
        for path in produced:
            print(path)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
