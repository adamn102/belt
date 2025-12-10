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

import os, subprocess, pysam, pandas as pd

def long_read_umi_grouping_fast(
        alignedReads, 
        correctedUMIs, 
        outputDir,
        grouping="UMI1", 
        threshold=5, 
        threads=8, 
        write_fastq=False, 
        log_dir=None,
        umi2_majority_threshold=0.5
        ):
    """
    Fast UMI grouping using BAM tagging and sorting.
    Supports variable numbers of UMIs (UMI1, UMI2, etc.) or legacy format (5pUMI, 3pUMI).
    
    :param alignedReads: Path to the BAM file with aligned reads.
    :param correctedUMIs: Path to the CSV file with corrected UMIs.
    :param outputDir: Directory where the output files will be saved.
    :param grouping: Which UMI column to group by (e.g., 'UMI1', '5pUMI'). Default: 'UMI1'.
    :param threshold: Minimum number of reads per UMI group. Default: 5.
    :param threads: Number of threads for samtools sort. Default: 8.
    :param write_fastq: If True, also write FASTQ files per UMI group with full UMI annotations. Default: False.
    :param log_dir: Directory where log files will be saved. If None, uses outputDir.
    :param umi2_majority_threshold: When grouping by UMI1 and UMI2 exists, the fraction (0-1) 
                                     of reads that must have the majority UMI2 species for that 
                                     UMI1 group to be kept. Only reads matching the majority UMI2 
                                     are retained. Default: 0.5 (50%).
    
    Additional filtering when grouping by UMI1 and UMI2 exists:
    - Only reads with both UMI1 and UMI2 (non-empty) are considered for grouping
    - For each UMI1 group passing the read count threshold, identify the most common UMI2 species
      among reads that have both UMIs
    - If the majority UMI2 represents >= umi2_majority_threshold fraction of reads with both UMIs, keep the group
    - Within kept groups, only retain reads that have both UMIs and match the majority UMI2 species
    - Metrics for this filter are recorded in the JSON log.
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

    # Apply UMI2-majority filter only when grouping by UMI1 and UMI2 exists
    umi2_filter_applied = False
    groups_passing_umi2 = None
    fraction_passing_umi2 = None
    reads_before_umi2_filter = 0
    reads_after_umi2_filter = 0
    umi1_to_majority_umi2 = {}  # Track majority UMI2 for each UMI1
    
    if gcol == "UMI1" and ("UMI2" in df.columns):
        umi2_filter_applied = True
        passing_umi1_values = []
        umi1_majority_counts = {}  # Track majority counts for each passing UMI1
        # Ensure UMI2 invalid placeholders are NA for counting
        df["UMI2"] = df["UMI2"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
        
        for umi1_value in keepers.index:
            # Only consider reads with both UMI1 and UMI2 (non-empty)
            sub_all = df.loc[df[gcol] == umi1_value, [gcol, "UMI2"]]
            total_in_umi1 = int(sub_all.shape[0])
            reads_before_umi2_filter += total_in_umi1
            
            # Filter to only reads with both UMIs (non-NA UMI2)
            sub = sub_all[sub_all["UMI2"].notna()]
            reads_with_both_umis = int(sub.shape[0])
            
            if reads_with_both_umis == 0:
                continue
                
            vc2 = sub["UMI2"].value_counts(dropna=True)
            if vc2.empty:
                continue
                
            majority_umi2 = vc2.index[0]
            majority_count = int(vc2.iloc[0])
            # Calculate fraction based on reads with both UMIs, not all reads
            majority_fraction = majority_count / reads_with_both_umis
            
            # Check if majority UMI2 meets threshold AND has enough reads after filtering
            if majority_fraction >= umi2_majority_threshold and majority_count >= threshold:
                passing_umi1_values.append(umi1_value)
                umi1_to_majority_umi2[umi1_value] = majority_umi2
                umi1_majority_counts[umi1_value] = majority_count
                reads_after_umi2_filter += majority_count
        
        # Update keepers to use majority counts instead of original counts
        if umi1_majority_counts:
            # Create new Series with majority counts for passing groups
            keepers = pd.Series(umi1_majority_counts, name=keepers.name)
        else:
            # No groups passed, filter keepers to empty
            keepers = keepers.loc[keepers.index.isin(passing_umi1_values)]
        groups_passing_umi2 = int(len(keepers))
        groups_passing_threshold = int(counts[counts >= threshold].shape[0])
        fraction_passing_umi2 = (
            (groups_passing_umi2 / groups_passing_threshold)
            if groups_passing_threshold > 0 else 0.0
        )

    # Prepare log and filtered CSV after all filters
    logs_dir = log_dir if log_dir else outputDir
    os.makedirs(logs_dir, exist_ok=True)
    log_json_path = os.path.join(logs_dir, f"{sample}_umi_group.log.json")
    filt_csv = os.path.join(outputDir, f"{sample}_filtered_umi_groups.csv")
    keepers.rename("ReadCount").rename_axis(gcol).reset_index().to_csv(filt_csv, index=False)

    if keepers.empty:
        try:
            import json
            from datetime import datetime, UTC
            total_unique = int(counts.shape[0])
            passing_threshold = int(counts[counts >= threshold].shape[0])
            fraction_passing_threshold = (passing_threshold / total_unique) if total_unique > 0 else 0.0
            log_payload = {
                "sample": sample,
                "grouping": gcol,
                "threshold": int(threshold),
                "umi2_majority_threshold": float(umi2_majority_threshold),
                "total_unique_umis": total_unique,
                "groups_passing_threshold": passing_threshold,
                "fraction_passing_threshold": fraction_passing_threshold,
                "umi2_majority_filter_applied": umi2_filter_applied,
                "groups_passing_umi2_majority": int(groups_passing_umi2 or 0),
                "fraction_passing_umi2_majority": float(fraction_passing_umi2 or 0.0),
                "reads_before_umi2_filter": int(reads_before_umi2_filter),
                "reads_after_umi2_filter": int(reads_after_umi2_filter),
                "reads_filtered_by_umi2": int(reads_before_umi2_filter - reads_after_umi2_filter),
                "aligned_bam": alignedReads,
                "corrected_umis_csv": correctedUMIs,
                "output_dir": outputDir,
                "created_at": datetime.now(UTC).isoformat().replace('+00:00', 'Z')
            }
            with open(log_json_path, "w") as jf:
                json.dump(log_payload, jf)
            logger.log_message(f"UMI grouping log written: {log_json_path}")
        except Exception as e:
            print(f"Warning: failed to write JSON log: {e}")
        logger.log_message("No UMIs pass filters; done.");
        return log_json_path

    # Write JSON log with requested metrics when there are keepers
    try:
        import json
        from datetime import datetime, UTC
        total_unique = int(counts.shape[0])
        passing_threshold = int(counts[counts >= threshold].shape[0])
        fraction_passing_threshold = (passing_threshold / total_unique) if total_unique > 0 else 0.0
        if umi2_filter_applied:
            groups_passing_umi2 = int(len(keepers))
            fraction_passing_umi2 = (groups_passing_umi2 / passing_threshold) if passing_threshold > 0 else 0.0
        log_payload = {
            "sample": sample,
            "grouping": gcol,
            "threshold": int(threshold),
            "umi2_majority_threshold": float(umi2_majority_threshold),
            "total_unique_umis": total_unique,
            "groups_passing_threshold": passing_threshold,
            "fraction_passing_threshold": fraction_passing_threshold,
            "umi2_majority_filter_applied": umi2_filter_applied,
            "groups_passing_umi2_majority": int(groups_passing_umi2 or 0),
            "fraction_passing_umi2_majority": float(fraction_passing_umi2 or 0.0),
            "reads_before_umi2_filter": int(reads_before_umi2_filter),
            "reads_after_umi2_filter": int(reads_after_umi2_filter),
            "reads_filtered_by_umi2": int(reads_before_umi2_filter - reads_after_umi2_filter),
            "aligned_bam": alignedReads,
            "corrected_umis_csv": correctedUMIs,
            "output_dir": outputDir,
            "created_at": datetime.now(UTC).isoformat().replace('+00:00', 'Z')
        }
        with open(log_json_path, "w") as jf:
            json.dump(log_payload, jf)
        logger.log_message(f"UMI grouping log written: {log_json_path}")
    except Exception as e:
        print(f"Warning: failed to write JSON log: {e}")

    # map read name -> UMI for kept groups
    # When UMI2 filtering is applied, only keep reads with both UMIs and matching majority UMI2
    if umi2_filter_applied and umi1_to_majority_umi2:
        # Filter reads to only include those with both UMIs and matching majority UMI2
        kept_rows = []
        for umi1_value in keepers.index:
            majority_umi2 = umi1_to_majority_umi2.get(umi1_value)
            if majority_umi2 is not None:
                # Only keep reads with both UMIs (UMI2 is not NA) and matching majority UMI2
                mask = (df[gcol] == umi1_value) & (df["UMI2"].notna()) & (df["UMI2"] == majority_umi2)
                kept_rows.append(df.loc[mask, ["Read Name", gcol]])
        if kept_rows:
            filtered_df = pd.concat(kept_rows, ignore_index=True)
            name_to_umi = dict(filtered_df.itertuples(index=False, name=None))
        else:
            name_to_umi = {}
    else:
        # No UMI2 filtering - keep all reads in passing groups
        name_to_umi = dict(df.loc[df[gcol].isin(keepers.index), ["Read Name", gcol]].itertuples(index=False, name=None))

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
    return log_json_path
