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
        outputDir,
        correctedUMIs=None, 
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
    :param correctedUMIs: Optional path to CSV with corrected UMIs. If None or empty, UMIs are read from BAM tags U1, U2, ...
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

    # 1) Load UMI data: from CSV or from BAM tags (U1, U2, ...)
    if correctedUMIs and str(correctedUMIs).strip():
        df = pd.read_csv(
            correctedUMIs,
            low_memory=False,
            keep_default_na=False
        )
        if "Read Name" not in df.columns:
            raise ValueError(f"'Read Name' column not found in {correctedUMIs}")
        if grouping not in df.columns:
            available_umis = [col for col in df.columns if col.startswith('UMI') or col.endswith('UMI')]
            raise ValueError(f"Grouping column '{grouping}' not found. Available UMI columns: {available_umis}")
        umi_columns = [col for col in df.columns if col != "Read Name"]
        dtype_dict = {"Read Name": "string"}
        dtype_dict.update({col: "string" for col in umi_columns})
        df = pd.read_csv(
            correctedUMIs,
            usecols=["Read Name"] + umi_columns,
            dtype=dtype_dict,
            low_memory=False,
            keep_default_na=False
        )
    else:
        # Build df from BAM tags U1, U2, ... (U1=UMI1, U2=UMI2, ...)
        logger.log_message("Reading UMI data from BAM tags (U1, U2, ...)...")
        # First pass: find max UMI tag index present in the BAM
        max_ui = 0
        with pysam.AlignmentFile(alignedReads, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                i = 1
                while read.has_tag(f"U{i}"):
                    i += 1
                max_ui = max(max_ui, i - 1)
        if max_ui == 0:
            raise ValueError(f"No reads with UMI tags (U1, U2, ...) found in BAM: {alignedReads}")
        umi_names_bam = [f"UMI{j+1}" for j in range(max_ui)]
        rows = []
        with pysam.AlignmentFile(alignedReads, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                qname = read.query_name
                row = {"Read Name": qname}
                for j in range(max_ui):
                    tag = f"U{j+1}"
                    if read.has_tag(tag):
                        val = read.get_tag(tag)
                        row[umi_names_bam[j]] = val if val else pd.NA
                    else:
                        row[umi_names_bam[j]] = pd.NA
                rows.append(row)
        df = pd.DataFrame(rows)
        umi_columns = umi_names_bam
        if grouping not in df.columns:
            raise ValueError(f"Grouping column '{grouping}' not found. BAM has UMI columns: {umi_columns}")
    
    gcol = grouping
    df[gcol] = df[gcol].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]

    # Performance optimization: Limit each UMI group to first 1000 reads early
    # This reduces memory usage and speeds up all subsequent operations
    # Note: For UMI2 filtering, we'll limit after UMI2 calculations to ensure
    # accurate majority detection. For non-UMI2 cases, early limiting is safe.
    MAX_READS_PER_UMI_GROUP = 1000
    will_apply_umi2_filter = (gcol == "UMI1" and "UMI2" in df.columns)
    
    # Only apply early limiting if we're NOT doing UMI2 filtering
    # (UMI2 filtering needs to see all reads to calculate accurate majority)
    if not keepers.empty and not will_apply_umi2_filter:
        # Filter DataFrame to only include passing UMI groups
        passing_mask = df[gcol].isin(keepers.index) & (df[gcol].notna())
        passing_df = df.loc[passing_mask].copy()
        
        # Limit each UMI group to first 1000 reads
        limited_rows = []
        for umi_value in keepers.index:
            umi_reads = passing_df[passing_df[gcol] == umi_value]
            if len(umi_reads) > MAX_READS_PER_UMI_GROUP:
                umi_reads = umi_reads.head(MAX_READS_PER_UMI_GROUP)
            limited_rows.append(umi_reads)
        
        if limited_rows:
            # Replace df with limited version for all subsequent operations
            limited_df = pd.concat(limited_rows, ignore_index=True)
            # Update keepers to reflect actual limited counts
            limited_counts = limited_df[gcol].value_counts(dropna=True)
            keepers = limited_counts[limited_counts.index.isin(keepers.index)]
            # Use limited_df for all subsequent DataFrame operations
            df = limited_df

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
        
        # Performance optimization: Limit DataFrame after UMI2 filtering is complete
        # This reduces memory and speeds up subsequent operations
        # IMPORTANT: Only limit reads that match the majority UMI2 for each group
        if not keepers.empty:
            limited_rows = []
            for umi_value in keepers.index:
                majority_umi2 = umi1_to_majority_umi2.get(umi_value)
                if majority_umi2 is not None:
                    # Only get reads with both UMIs matching the majority
                    mask = (df[gcol] == umi_value) & (df["UMI2"].notna()) & (df["UMI2"] == majority_umi2)
                    majority_reads = df.loc[mask]
                    
                    # Limit to first 1000 reads per UMI group
                    if len(majority_reads) > MAX_READS_PER_UMI_GROUP:
                        majority_reads = majority_reads.head(MAX_READS_PER_UMI_GROUP)
                    
                    limited_rows.append(majority_reads)
            
            if limited_rows:
                limited_df = pd.concat(limited_rows, ignore_index=True)
                # Update keepers to reflect actual limited counts (but only for groups that still pass threshold)
                limited_counts = limited_df[gcol].value_counts(dropna=True)
                # Only keep groups that still meet threshold after limiting
                keepers = limited_counts[limited_counts >= threshold]
                df = limited_df

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
    # Ensure perfect UMI concordance: each read must have the exact correct UMI1 and UMI2
    # Note: DataFrame is already limited to MAX_READS_PER_UMI_GROUP per group (if applicable)
    if umi2_filter_applied and umi1_to_majority_umi2:
        # Filter reads to only include those with both UMIs and matching majority UMI2
        # Create lookup that includes both UMI1 and UMI2 for validation
        kept_rows = []
        read_to_umi1_umi2 = {}  # For validation: read_name -> (UMI1, UMI2)
        
        for umi1_value in keepers.index:
            majority_umi2 = umi1_to_majority_umi2.get(umi1_value)
            if majority_umi2 is not None:
                # Only keep reads with both UMIs (UMI2 is not NA) and matching majority UMI2
                # Perfect concordance: UMI1 must exactly match umi1_value AND UMI2 must exactly match majority_umi2
                mask = (df[gcol] == umi1_value) & (df["UMI2"].notna()) & (df["UMI2"] == majority_umi2)
                matching_reads = df.loc[mask, ["Read Name", gcol, "UMI2"]]
                
                # Limit to first 1000 reads per UMI group
                if len(matching_reads) > MAX_READS_PER_UMI_GROUP:
                    matching_reads = matching_reads.head(MAX_READS_PER_UMI_GROUP)
                
                # Validate perfect concordance: ensure all reads have exact correct UMIs
                for _, row in matching_reads.iterrows():
                    read_name = row["Read Name"]
                    read_umi1 = row[gcol]
                    read_umi2 = row["UMI2"]
                    
                    # Only include if UMI1 and UMI2 match exactly
                    if read_umi1 == umi1_value and read_umi2 == majority_umi2:
                        kept_rows.append({"Read Name": read_name, gcol: read_umi1})
                        read_to_umi1_umi2[read_name] = (read_umi1, read_umi2)
        
        if kept_rows:
            filtered_df = pd.DataFrame(kept_rows)
            name_to_umi = dict(filtered_df[["Read Name", gcol]].itertuples(index=False, name=None))
        else:
            name_to_umi = {}
            read_to_umi1_umi2 = {}
    else:
        # No UMI2 filtering - keep all reads in passing groups
        # Still ensure UMI1 concordance (exact match required)
        passing_mask = df[gcol].isin(keepers.index) & (df[gcol].notna())
        passing_reads = df.loc[passing_mask, ["Read Name", gcol]]
        # Validate that UMI1 matches exactly
        validated_passing = passing_reads[passing_reads[gcol].isin(keepers.index)]
        
        # Limit each UMI group to first 1000 reads
        if len(validated_passing) > 0:
            limited_reads = []
            for umi_value in validated_passing[gcol].unique():
                umi_reads = validated_passing[validated_passing[gcol] == umi_value]
                if len(umi_reads) > MAX_READS_PER_UMI_GROUP:
                    umi_reads = umi_reads.head(MAX_READS_PER_UMI_GROUP)
                limited_reads.append(umi_reads)
            
            if limited_reads:
                validated_passing = pd.concat(limited_reads, ignore_index=True)
        
        name_to_umi = dict(validated_passing.itertuples(index=False, name=None)) if len(validated_passing) > 0 else {}
        read_to_umi1_umi2 = None  # Not needed when UMI2 filtering is not applied

    # 2) Tag BAM with ZU:<string> = UMI (only for kept reads)
    # Ensure perfect UMI concordance: validate each read has exact correct UMIs before writing
    tagged_bam = os.path.join(outputDir, f"{sample}.tagged.bam")
    
    # Create validation lookup for perfect concordance checking when UMI2 filtering is applied
    # Only include reads that are in name_to_umi (already limited to 1000 per group)
    if umi2_filter_applied and umi1_to_majority_umi2:
        # Create lookup: read_name -> (expected_UMI1, expected_UMI2) for final validation
        # Only include reads that are in name_to_umi to maintain consistency with the limit
        validation_lookup = {}
        for read_name in name_to_umi.keys():
            read_row = df[df["Read Name"] == read_name]
            if len(read_row) > 0:
                read_umi1 = read_row[gcol].iloc[0]
                read_umi2 = read_row["UMI2"].iloc[0] if "UMI2" in read_row.columns else None
                if read_umi2 is not None:
                    validation_lookup[read_name] = (read_umi1, read_umi2)
    else:
        validation_lookup = None
    
    with pysam.AlignmentFile(alignedReads, "rb") as bam_in, \
            pysam.AlignmentFile(tagged_bam, "wb", template=bam_in) as bam_out:
        reads_written = 0
        reads_skipped = 0
        
        for r in bam_in.fetch(until_eof=True):
            umi = name_to_umi.get(r.query_name)
            if umi is None:
                continue
            
            # Final validation for perfect UMI concordance when UMI2 filtering is applied
            if validation_lookup is not None:
                expected_umis = validation_lookup.get(r.query_name)
                if expected_umis is None:
                    # Read not in validation lookup - skip to ensure perfect concordance
                    reads_skipped += 1
                    continue
                
                # Verify against original dataframe to ensure perfect match
                read_row = df[df["Read Name"] == r.query_name]
                if len(read_row) == 0:
                    reads_skipped += 1
                    continue
                
                actual_umi1 = read_row[gcol].iloc[0]
                actual_umi2 = read_row["UMI2"].iloc[0] if "UMI2" in read_row.columns else None
                expected_umi1, expected_umi2 = expected_umis
                
                # Only write if UMIs match exactly (perfect concordance)
                if actual_umi1 != expected_umi1 or (actual_umi2 is not None and actual_umi2 != expected_umi2):
                    reads_skipped += 1
                    continue
            
            r.set_tag("ZU", str(umi), value_type='Z')
            bam_out.write(r)
            reads_written += 1
        
        if reads_skipped > 0:
            logger.log_message(f"Perfect concordance validation: {reads_skipped} reads skipped, {reads_written} reads written")
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
