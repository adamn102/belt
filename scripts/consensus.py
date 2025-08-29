#!/usr/bin/env python3
"""
BELT Consensus Module
Generate consensus sequences from grouped BAM files
"""

import os
import re
import sys
import itertools
import subprocess
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam


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

def run_consensus_for_umi(reads_fq: Path, ref_fa: Path, meta_dir: Path, 
                         con_dir: Path, bam_dir: Path, umi_threads: int, 
                         model: str, log_file: Path, umi_safe: str) -> str:
    """Run consensus.sh for one UMI and return the consensus FASTA path."""
    logfh = open(log_file, "a") if log_file else subprocess.DEVNULL
    
    try:
        # Run the consensus.sh script
        subprocess.run(
            [
                "./scripts/consensus.sh",
                str(reads_fq),
                str(ref_fa),
                str(meta_dir),
                str(con_dir),
                str(umi_threads),
                str(model),
                str(bam_dir),  # optional 7th arg handled by consensus.sh
            ],
            stdout=logfh,
            stderr=logfh,
            check=True,
        )
    finally:
        if logfh is not subprocess.DEVNULL:
            logfh.close()
    
    # Clean up temporary FASTQ
    try:
        os.remove(reads_fq)
    except OSError:
        pass
    
    # consensus.sh writes <UMI>_consensus.fasta
    return str(Path(con_dir) / f"{umi_safe}_consensus.fasta")

def check_dependencies():
    """Check if required external tools are available."""
    required_tools = ['minimap2', 'samtools', 'medaka_variant', 'bcftools']
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)
    
    if missing_tools:
        raise RuntimeError(
            f"Missing required tools: {', '.join(missing_tools)}\n"
            "Please install: minimap2, samtools, medaka, bcftools"
        )

def process_grouped_bam(grouped_bam: str, reference: str, output_dir: str, 
                       sample: str, tag: str = "ZU", threads: int = 2, 
                       procs: int = 1, umi_threads: int = 1, 
                       model: str = "r1041_e82_400bps_sup_variant_v5.0.0",
                       tmp_root: str = None, manifest: str = None, 
                       log_file: str = None) -> list:
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
        model: Medaka model name
        tmp_root: Temporary directory root
        manifest: Output manifest file path
        log_file: Log file path
    
    Returns:
        List of produced consensus FASTA file paths
    """
    # Check dependencies
    check_dependencies()
    
    # Setup directories
    out_root = Path(output_dir) / sample
    con_dir = out_root / "consensus"
    bam_dir = out_root / "bam"
    tmp_root = Path(tmp_root) if tmp_root else (out_root / "tmp")
    
    con_dir.mkdir(parents=True, exist_ok=True)
    bam_dir.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)
    
    log_file = Path(log_file) if log_file else (con_dir / "runner.log")
    
    produced = []
    futures = []
    
    def get_umi_paths(umi_safe: str):
        """Get paths for UMI-specific files."""
        umi_root = tmp_root / umi_safe
        reads_fq = umi_root / f"{umi_safe}_reads.fastq"  # matches consensus.sh UMI parsing
        meta_dir = umi_root / "meta"                     # safe for cleanup inside consensus.sh
        meta_dir.mkdir(parents=True, exist_ok=True)
        return reads_fq, meta_dir
    
    # Process grouped BAM file
    with pysam.AlignmentFile(grouped_bam, "rb") as bam, \
         ProcessPoolExecutor(max_workers=procs) as pool:
        
        # Group reads by UMI tag
        key = lambda r: r.get_tag(tag)
        for umi, group_iter in itertools.groupby(bam.fetch(until_eof=True), key=key):
            umi_safe = sanitize_umi_name(umi)
            reads_fq, meta_dir = get_umi_paths(umi_safe)
            reads_fq.parent.mkdir(parents=True, exist_ok=True)
            
            # Write FASTQ for this UMI group
            write_group_fastq(umi, group_iter, reads_fq)
            
            # Submit consensus job
            fut = pool.submit(
                run_consensus_for_umi,
                reads_fq, Path(reference), meta_dir, con_dir,
                bam_dir, umi_threads, model, log_file, umi_safe
            )
            futures.append(fut)
            
            # Throttle queued tasks to prevent memory issues
            if len(futures) >= procs * 4:
                for done in as_completed(futures[:procs]):
                    produced.append(done.result())
                futures = futures[procs:]
        
        # Process remaining futures
        for done in as_completed(futures):
            produced.append(done.result())
    
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
    ap.add_argument("--procs", type=int, default=1, help="Parallel UMI processes")
    ap.add_argument("--umi-threads", type=int, default=1, help="Threads per UMI for consensus calling")
    ap.add_argument("--model", default="r1041_e82_400bps_sup_variant_v5.0.0", help="Medaka model")
    ap.add_argument("--tmp-root", help="Temporary directory root")
    ap.add_argument("--manifest", help="Output manifest file")
    ap.add_argument("--log-file", help="Log file path")
    
    args = ap.parse_args()
    
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
            log_file=args.log_file
        )
        
        # Print produced files
        for path in produced:
            print(path)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
