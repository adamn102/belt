import pysam
import pandas as pd
import sys
import os
from collections import Counter

# Add the current directory to the path so we can import local modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from utils import log
except ImportError:
    # Create a simple mock log function if utils is not available
    class MockLog:
        def log_message(self, message):
            print(f"[LOG] {message}")
    
    log = MockLog


def long_read_umi_grouping(reference, alignedReads, correctedUMIs, outputDir, grouping="UMI1", threshold=2):
    """
    Generate grouped FASTQ files from long reads using corrected UMIs.
    Supports variable numbers of UMIs (UMI1, UMI2, UMI3, etc.) or legacy format (5pUMI, 3pUMI).
    
    :param reference: Path to the reference genome.
    :param alignedReads: Path to the BAM file with aligned reads.
    :param correctedUMIs: Path to the CSV file with corrected UMIs (must have 'Read Name' column and UMI columns).
    :param outputDir: Directory where the output files will be saved.
    :param grouping: Which UMI column to group by (e.g., 'UMI1', '5pUMI'). Default: 'UMI1'.
    :param threshold: Minimum number of reads per UMI group to include. Default: 2.
    """
    os.makedirs(outputDir, exist_ok=True)

    sampleName = os.path.basename(alignedReads).split(".")[0]
    logger = log()
    logger.log_message("Starting long read UMI grouping...")

    # Read in corrected UMIs CSV
    df = pd.read_csv(
        correctedUMIs,
        low_memory=False,
        keep_default_na=False  # so "NA" stays "NA" instead of NaN, optional
    )
    
    # Verify required columns exist
    if "Read Name" not in df.columns:
        raise ValueError(f"'Read Name' column not found in {correctedUMIs}")
    
    if grouping not in df.columns:
        raise ValueError(f"Grouping column '{grouping}' not found in {correctedUMIs}. Available columns: {list(df.columns)}")
    
    # Clean up UMI grouping
    df[grouping] = df[grouping].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    # Group by selected UMI grouping
    grouped = df.groupby(grouping).agg({
        "Read Name": ["count", lambda x: list(x)]
    })


    grouped.columns = ["ReadCount", "ReadNames"]

    # Filter groups by threshold
    filtered = grouped[grouped["ReadCount"] >= threshold]

    filtered_csv = os.path.join(outputDir, f"{sampleName}_filtered_umi_groups.csv")
    filtered.reset_index().to_csv(filtered_csv, index=False)
    logger.log_message(f"Filtered UMI groups saved to {filtered_csv}")

    # Create output dir per sample
    sample_output_dir = os.path.join(outputDir, sampleName)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Read name -> UMI lookup
    read_to_umis = df.set_index("Read Name")[df.columns[1:]].to_dict(orient="index")

    # print(read_to_umis)

    # Iterate over UMI groups
    for umi, row in filtered.iterrows():
        out_fastq = os.path.join(sample_output_dir, f"{umi}_reads.fastq")

        extract_reads_to_fastq(
            read_names=row["ReadNames"],
            bam_path=alignedReads,
            output_fastq_path=out_fastq,
            read_umi_lookup=read_to_umis
        )

    logger.log_message(f"Grouped FASTQ files written to {sample_output_dir}")


def extract_reads_to_fastq(
        read_names,
        bam_path,
        output_fastq_path,
        read_umi_lookup
        ):
    """
    Extract reads from BAM and write FASTQ with UMI annotations.
    Handles variable numbers of UMIs (UMI1, UMI2, etc.) or legacy format (5pUMI, 3pUMI).

    Parameters:
        read_names (Iterable[str]): Read names to extract.
        bam_path (str): Path to BAM file.
        output_fastq_path (str): Output FASTQ file path.
        read_umi_lookup (dict): Read name -> UMI info dictionary.
    """
    read_names_set = set(read_names)
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    n_written = 0

    with open(output_fastq_path, "w") as fastq_out:
        for read in bamfile.fetch(until_eof=True):
            if read.query_name in read_names_set and not read.is_secondary and not read.is_supplementary:
                seq = read.query_sequence
                qual = read.qual

                if seq and qual:
                    umi_info = read_umi_lookup.get(read.query_name, {})
                    
                    # Build UMI string dynamically based on available UMI columns
                    umi_parts = []
                    if umi_info:
                        # Sort keys to ensure consistent ordering (UMI1, UMI2, etc.)
                        for umi_key in sorted(umi_info.keys()):
                            umi_value = umi_info.get(umi_key, "NA")
                            umi_parts.append(f"{umi_key}:{umi_value}")
                    
                    # Build header with all UMIs
                    if umi_parts:
                        umi_string = "|".join(umi_parts)
                        header = f"{read.query_name}|{umi_string}"
                    else:
                        header = f"{read.query_name}|UMI:NA"
                    
                    # Write FASTQ format
                    fastq_out.write(f"@{header}\n{seq}\n+\n{qual}\n")
                    n_written += 1

    bamfile.close()

    if n_written == 0:
        print(f"Warning: No reads written for {output_fastq_path}")

import os, subprocess, pysam, pandas as pd

def long_read_umi_grouping_fast(
        reference, alignedReads, correctedUMIs, outputDir,
        grouping="UMI1", threshold=5, threads=8, write_fastq=False
        ):
    """
    Fast UMI grouping using BAM tagging and sorting.
    Supports variable numbers of UMIs (UMI1, UMI2, etc.) or legacy format (5pUMI, 3pUMI).
    
    :param reference: Path to the reference genome.
    :param alignedReads: Path to the BAM file with aligned reads.
    :param correctedUMIs: Path to the CSV file with corrected UMIs.
    :param outputDir: Directory where the output files will be saved.
    :param grouping: Which UMI column to group by (e.g., 'UMI1', '5pUMI'). Default: 'UMI1'.
    :param threshold: Minimum number of reads per UMI group. Default: 5.
    :param threads: Number of threads for samtools sort. Default: 8.
    :param write_fastq: If True, also write FASTQ files per UMI group with full UMI annotations. Default: False.
    """
    os.makedirs(outputDir, exist_ok=True)
    sample = os.path.basename(alignedReads).split(".")[0]
    logger = log()
    logger.log_message("UMI grouping via tag+sort...")

    # 1) Load CSV and determine available UMI columns
    df = pd.read_csv(
        correctedUMIs,
        low_memory=False,
        keep_default_na=False
    )
    
    # Verify required columns exist
    if "Read Name" not in df.columns:
        raise ValueError(f"'Read Name' column not found in {correctedUMIs}")
    
    if grouping not in df.columns:
        available_umis = [col for col in df.columns if col.startswith('UMI') or col.endswith('UMI')]
        raise ValueError(f"Grouping column '{grouping}' not found. Available UMI columns: {available_umis}")
    
    # Get all UMI columns for dtype specification
    umi_columns = [col for col in df.columns if col != "Read Name"]
    dtype_dict = {"Read Name": "string"}
    dtype_dict.update({col: "string" for col in umi_columns})
    
    # Reload with proper dtypes (only Read Name and UMI columns)
    df = pd.read_csv(
        correctedUMIs,
        usecols=["Read Name"] + umi_columns,
        dtype=dtype_dict,
        low_memory=False,
        keep_default_na=False
    )
    
    gcol = grouping
    df[gcol] = df[gcol].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    filt_csv = os.path.join(outputDir, f"{sample}_filtered_umi_groups.csv")
    keepers.rename("ReadCount").rename_axis(gcol).reset_index().to_csv(filt_csv, index=False)
    if keepers.empty:
        logger.log_message("No UMIs pass threshold; done.");
        return

    # map read name -> UMI for kept groups
    name_to_umi = dict(df.loc[df[gcol].isin(keepers.index), ["Read Name", gcol]].itertuples(index=False, name=None))
    del df

    # 2) Tag BAM with ZU:<string> = UMI (only for kept reads)
    tagged_bam = os.path.join(outputDir, f"{sample}.tagged.bam")
    with pysam.AlignmentFile(alignedReads, "rb") as bam_in, \
            pysam.AlignmentFile(tagged_bam, "wb", template=bam_in) as bam_out:
        for r in bam_in.fetch(until_eof=True):
            umi = name_to_umi.get(r.query_name)
            if umi is None:
                continue
            r.set_tag("ZU", str(umi), value_type='Z')
            bam_out.write(r)
    # index optional, but sort will read sequentially
    # 3) Sort by tag
    grouped_bam = os.path.join(outputDir, f"{sample}.grouped.bam")
    subprocess.run(
        ["samtools", "sort", "-@",
         str(threads), "-t", "ZU", "-o", grouped_bam, tagged_bam],
        check=True
    )
    # Remove tagged bam (optional cleanup)
    os.remove(tagged_bam)
    
    logger.log_message(f"Grouped BAM created: {grouped_bam}")
    
    # Optional: Write FASTQ files per UMI group with full UMI annotations
    if write_fastq:
        logger.log_message("Writing FASTQ files with full UMI annotations...")
        
        # Create read name -> all UMIs lookup
        read_to_all_umis = df.set_index("Read Name")[umi_columns].to_dict(orient="index")
        
        out_dir = os.path.join(outputDir, sample)
        os.makedirs(out_dir, exist_ok=True)
        
        curr_umi, fh = None, None
        written = 0
        
        with pysam.AlignmentFile(grouped_bam, "rb") as bam:
            for r in bam.fetch(until_eof=True):
                umi = r.get_tag("ZU")
                
                # Open new file when UMI changes
                if umi != curr_umi:
                    if fh:
                        fh.close()
                    curr_umi = umi
                    # Sanitize UMI for filename
                    safe_umi = umi.replace("/", "_").replace("\\", "_")
                    fastq_path = os.path.join(out_dir, f"{safe_umi}_reads.fastq")
                    fh = open(fastq_path, "w", buffering=1024 * 1024)
                
                # Get read info
                rn = r.query_name
                seq = r.query_sequence or ""
                qual = r.qual if r.qual is not None else ("I" * len(seq))
                
                # Build header with all UMI info
                umi_info = read_to_all_umis.get(rn, {})
                umi_parts = []
                if umi_info:
                    for umi_key in sorted(umi_info.keys()):
                        umi_value = umi_info.get(umi_key, "NA")
                        umi_parts.append(f"{umi_key}:{umi_value}")
                
                if umi_parts:
                    umi_string = "|".join(umi_parts)
                    header = f"{rn}|{umi_string}"
                else:
                    header = f"{rn}|UMI:NA"
                
                # Write FASTQ entry
                fh.write(f"@{header}\n{seq}\n+\n{qual}\n")
                written += 1
        
        if fh:
            fh.close()
        
        logger.log_message(f"Wrote {written} reads across {len(keepers)} UMI groups to {out_dir}")
    
    logger.log_message(f"Done. Processed {len(keepers)} UMI groups passing threshold.")
