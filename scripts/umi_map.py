from utils import log
from extract import *
from utils import umi_clustering, correct_UMI_dict
from pybktree import BKTree
from umi_tools import UMIClusterer
from collections import Counter
import pandas as pd
import os
import numpy as np
import pysam


def umi_extraction_map(
        reference,
        alignedReads,
        minReadLength,
        barcodeMotifs,
        whitelist,
        outputDir,
        max_distance=3,
        logOutputDir=None,
        write_csv=False):
    """
    Extract and correct UMIs, then write a tagged BAM (U1, U2, ... tags).
    Optionally write a CSV of corrected UMIs (write_csv=True).
    """
    sampleName = os.path.basename(alignedReads).split(".")[0]

    # Initialize output directories
    os.makedirs(outputDir, exist_ok=True)
    if logOutputDir is None:
        logOutputDir = outputDir
    os.makedirs(logOutputDir, exist_ok=True)

    # Set up logging
    logger = log()
    from datetime import datetime, UTC
    import json
    log_entries = []
    start_time = datetime.now(UTC).isoformat().replace('+00:00', 'Z')

    def add_log(message, level="INFO"):
        logger.log_message(message)
        log_entries.append({
            "timestamp": datetime.now(UTC).isoformat().replace('+00:00', 'Z'),
            "level": level,
            "message": message
        })

    add_log("Extracting UMIs from aligned reads...")

    # confirm barcodeMotifs is a list
    if isinstance(barcodeMotifs, str):
        barcodeMotifs = [barcodeMotifs]
    # count the total number of motifs
    total_motifs = len(barcodeMotifs)

    # find the raw UMIs and return the motif indices and the raw UMI dictionary
    motif_indices, raw_UMI_dict = findUMI(reference, alignedReads, barcodeMotifs, minReadLength)

    if motif_indices is None or raw_UMI_dict is None:
        add_log("UMI extraction failed. Please check the input files and parameters.", level="ERROR")
        # Write minimal failure log and exit
        end_time = datetime.now(UTC).isoformat().replace('+00:00', 'Z')
        out_path = os.path.join(logOutputDir, f"{sampleName}_umi_extraction_map.log.json")
        with open(out_path, "w") as jf:
            json.dump({
                "sample": sampleName,
                "status": "failed",
                "start_time": start_time,
                "end_time": end_time,
                "logs": log_entries
            }, jf, indent=2)
        return
    add_log(f"Found motifs at indices: {motif_indices}")
    add_log(f"UMI extraction of {sampleName} with {len(barcodeMotifs)} barcode motifs finished.")

    # convert raw_UMI_dict into a datafame with a row for each read and a column for each motif
    df = pd.DataFrame.from_dict(raw_UMI_dict, orient='index', columns=barcodeMotifs)
    df.index.name = 'Read Name'
    df.reset_index(inplace=True)

    # update the column names to be the motif names, e.g. UMI1, UMI2, UMI3, etc., skip the first column
    umi_names = [f"UMI{i+1}" for i in range(total_motifs)]
    df.columns = ["Read Name"] + umi_names
    
    # Log statistics about missing UMIs and compute presence combinations
    total_reads = len(df)
    missing_stats = {}
    # Compute per-column "only this UMI present" counts and "all UMIs present" count
    umi_only_counts = {}
    if umi_names:
        all_present_mask = df[umi_names].notna().all(axis=1)
        all_umis_present_count = int(all_present_mask.sum())
        # Compute per-UMI unique counts restricted to rows with all UMIs present
        unique_counts_all_present = {
            umi_col: int(df.loc[all_present_mask, umi_col].dropna().nunique())
            for umi_col in umi_names
        }
    else:
        all_umis_present_count = 0
        unique_counts_all_present = {}
    for umi_col in umi_names:
        missing_count = df[umi_col].isna().sum()
        if missing_count > 0:
            add_log(f"{umi_col}: {missing_count}/{total_reads} reads have missing UMI ({missing_count/total_reads*100:.1f}%)")
        else:
            add_log(f"{umi_col}: All reads have valid UMI")
        missing_stats[umi_col] = {
            "missing": int(missing_count),
            "total_reads": int(total_reads)
        }
        # Only-this-UMI-present mask: this UMI present, all other UMIs missing
        other_cols = [c for c in umi_names if c != umi_col]
        if other_cols:
            only_mask = df[umi_col].notna() & df[other_cols].isna().all(axis=1)
        else:
            # If there is only one UMI column, "only" equals present
            only_mask = df[umi_col].notna()
        umi_only_counts[umi_col] = int(only_mask.sum())
    add_log(f"Total reads: {total_reads}")
    add_log(f"Reads with all UMIs present: {all_umis_present_count}")

    # Save a copy of raw UMIs before processing (for Hamming distance analysis)
    raw_df = df.copy()

    # Handle different whitelist modes
    wl_mode = str(whitelist).lower()

    if wl_mode == "false":
        add_log("Whitelist set to 'false'. Using raw UMIs without mapping.")
        # No correction needed - df already contains raw UMIs
        # This will be handled in the output section where raw_df == df

    elif wl_mode == "cluster":
        add_log("Whitelist mode set to 'cluster'. Performing UMI-tools clustering.")
        # Cluster and correct UMIs for each column
        corrected_df = df[["Read Name"]].copy()
        clustering_stats = {
            "per_column": {},
            "overall": {
                "total_raw": 0,
                "total_corrected": 0,
                "total_clusters": 0,
                "total_singletons": 0,
                "mean_cluster_size": None
            }
        }
        overall_cluster_sizes = []
        
        for column in umi_names:
            add_log(f"Clustering {column}...")
            # Get valid UMIs (non-NaN)
            valid_mask = df[column].notna()
            valid_umis = df.loc[valid_mask, column].tolist()
            total_valid = len(valid_umis)
            
            if total_valid == 0:
                corrected_df[column] = df[column]
                clustering_stats["per_column"][column] = {
                    "total_raw": 0,
                    "unique_raw": 0,
                    "total_corrected": 0,
                    "reduction_in_unique_umis": 0,
                    "reduction_percentage": 0.0,
                    "total_clusters": 0,
                    "total_singletons": 0,
                    "mean_cluster_size": None,
                    "cluster_size_distribution": {},
                    "raw_umi_counts": {},
                    "corrected_umi_counts": {}
                }
                add_log(f"{column}: No valid UMIs to cluster.")
                continue
            
            # Count raw UMIs before clustering
            raw_umi_counts = Counter(valid_umis)
            unique_raw = len(raw_umi_counts)
            
            # Log raw UMI counts (top 20 most frequent)
            add_log(f"{column}: Raw UMIs - {unique_raw} unique UMIs from {total_valid} total reads")
            if unique_raw > 0:
                top_raw = raw_umi_counts.most_common(20)
                add_log(f"{column}: Top raw UMIs by count:")
                for umi, count in top_raw:
                    add_log(f"{column}:   {umi}: {count} reads")
                if unique_raw > 20:
                    add_log(f"{column}:   ... and {unique_raw - 20} more unique UMIs")
            
            # Perform clustering using UMI-tools API
            correction_map, cluster_info = cluster_umis_umitools(
                valid_umis, 
                threshold=max_distance,
                cluster_method="directional"
            )
            
            # Apply correction mapping
            corrected_umis = []
            for i, umi in enumerate(df[column]):
                if pd.isna(umi):
                    corrected_umis.append(None)
                else:
                    corrected_umis.append(correction_map.get(umi, umi))
            
            corrected_df[column] = corrected_umis
            
            # Count corrected UMIs after clustering
            corrected_umi_counts = Counter(c for c in corrected_umis if c is not None)
            unique_corrected = len(corrected_umi_counts)
            
            # Log corrected UMI counts (top 20 most frequent)
            reduction = unique_raw - unique_corrected
            reduction_pct = (reduction / unique_raw * 100) if unique_raw > 0 else 0.0
            add_log(f"{column}: Corrected UMIs - {unique_corrected} unique UMIs ({reduction} reduction, {reduction_pct:.1f}%)")
            if unique_corrected > 0:
                top_corrected = corrected_umi_counts.most_common(20)
                add_log(f"{column}: Top corrected UMIs by count:")
                for umi, count in top_corrected:
                    add_log(f"{column}:   {umi}: {count} reads")
                if unique_corrected > 20:
                    add_log(f"{column}:   ... and {unique_corrected - 20} more unique UMIs")
            
            # Calculate statistics for this column
            num_clusters = cluster_info["num_clusters"]
            num_singletons = cluster_info["num_singletons"]
            cluster_sizes = cluster_info["cluster_sizes"]
            overall_cluster_sizes.extend(cluster_sizes)
            mean_cluster_size = float(np.mean(cluster_sizes)) if cluster_sizes else None
            
            clustering_stats["per_column"][column] = {
                "total_raw": int(total_valid),
                "unique_raw": int(unique_raw),
                "total_corrected": int(unique_corrected),
                "reduction_in_unique_umis": int(reduction),
                "reduction_percentage": float(reduction_pct),
                "total_clusters": int(num_clusters),
                "total_singletons": int(num_singletons),
                "mean_cluster_size": mean_cluster_size,
                "cluster_size_distribution": dict(Counter(cluster_sizes)),
                "raw_umi_counts": {umi: int(count) for umi, count in raw_umi_counts.items()},
                "corrected_umi_counts": {umi: int(count) for umi, count in corrected_umi_counts.items()}
            }
            
            clustering_stats["overall"]["total_raw"] += int(total_valid)
            clustering_stats["overall"]["total_corrected"] += int(unique_corrected)
            clustering_stats["overall"]["total_clusters"] += int(num_clusters)
            clustering_stats["overall"]["total_singletons"] += int(num_singletons)
            
            add_log(
                f"{column}: {total_valid} UMIs clustered into {num_clusters} clusters "
                f"({num_singletons} singletons); mean cluster size: {mean_cluster_size:.2f}"
            )
        
        # Finalize overall statistics
        if overall_cluster_sizes:
            clustering_stats["overall"]["mean_cluster_size"] = float(np.mean(overall_cluster_sizes))
        
        df = corrected_df
        add_log(f"Overall clustering: {clustering_stats['overall']['total_clusters']} clusters "
                f"from {clustering_stats['overall']['total_raw']} raw UMIs.")

    else:
        # TODO: Add more whitelist strict mode where non-matching UMIs are not allowed 
        # TODO: Add output data of matching and non-matching UMIs (e.g. rates or counts)
        add_log("Mapping raw UMIs to whitelist...")
        whitelist_df = load_whitelist(whitelist)
        
        # for each column in the raw UMI dataframe, use BK_scan to find the best match in the whitelist
        # perform this operation for each column in the raw UMI dataframe, are return a dataframe with the best matches
        # TODO: Test this is increasingly divergent UMI sequences

        # create a new dataframe with the best matches for each read, initiate the df with the read names

        best_matches_df = df[["Read Name"]].copy()
        # Collect per-column and overall mapping/divergence metrics for logging and JSON summary
        mapping_stats = {
            "per_column": {},
            "overall": {
                "total_raw": 0,
                "total_mapped": 0,
                "fraction_mapped": 0.0,
                "exact_matches": 0,
                "distance_histogram": {},
                "mean_distance": None,
                "median_distance": None,
                "length_mismatches": 0,
                "search_errors": 0,
                "no_matches": 0
            }
        }
        overall_distances = []
        # Track rows that contain any non-mapping UMI (for strict drop)
        nonmapping_row_mask = np.zeros(len(df), dtype=bool)
        for column in df.columns[1:]:
            best_matches = BK_scan(df[column], whitelist_df[column], max_distance)
            # update the best_matches_df with the best matches for each read
            # Handle NaN values - preserve them in the output
            mapped_umis = []
            for i, umi in enumerate(df[column]):
                if pd.isna(umi):
                    # If original UMI was NaN, keep it as NaN
                    mapped_umis.append(None)
                else:
                    # Get the best match (first element of tuple)
                    mapped, dist = best_matches[umi]
                    mapped_umis.append(mapped)
                    # Mark row for drop if this present UMI did not map within allowed distance
                    if not isinstance(dist, int):
                        nonmapping_row_mask[i] = True
            best_matches_df[column] = mapped_umis
            # Compute stats for this column
            valid_raw = [umi for umi in df[column] if pd.notna(umi)]
            total_valid = len(valid_raw)
            mapped_count = 0
            exact_count = 0
            length_mismatch = 0
            search_error = 0
            no_match = 0
            distances = []
            for umi in valid_raw:
                mapped, dist = best_matches[umi]
                if isinstance(dist, int):
                    mapped_count += 1
                    distances.append(dist)
                    if dist == 0:
                        exact_count += 1
                else:
                    if dist == "LengthMismatch":
                        length_mismatch += 1
                    elif dist == "SearchError":
                        search_error += 1
                    elif dist == "NoMatch":
                        no_match += 1
            frac_mapped = (mapped_count / total_valid) if total_valid else 0.0
            distance_hist = dict(Counter(distances)) if distances else {}
            mean_dist = float(np.mean(distances)) if distances else None
            median_dist = float(np.median(distances)) if distances else None
            mapping_stats["per_column"][column] = {
                "total_raw": int(total_valid),
                "mapped": int(mapped_count),
                "fraction_mapped": float(frac_mapped),
                "exact_matches": int(exact_count),
                "distance_histogram": distance_hist,
                "mean_distance": mean_dist,
                "median_distance": median_dist,
                "length_mismatches": int(length_mismatch),
                "search_errors": int(search_error),
                "no_matches": int(no_match)
            }
            # Update overall tallies
            mapping_stats["overall"]["total_raw"] += int(total_valid)
            mapping_stats["overall"]["total_mapped"] += int(mapped_count)
            mapping_stats["overall"]["exact_matches"] += int(exact_count)
            mapping_stats["overall"]["length_mismatches"] += int(length_mismatch)
            mapping_stats["overall"]["search_errors"] += int(search_error)
            mapping_stats["overall"]["no_matches"] += int(no_match)
            overall_distances.extend(distances)
            # Log concise per-column summary
            if total_valid:
                if distances:
                    add_log(f"{column}: {mapped_count}/{total_valid} ({frac_mapped*100:.1f}%) mapped; mean={mean_dist:.2f}, median={median_dist:.2f}, exact={exact_count}")
                else:
                    add_log(f"{column}: {mapped_count}/{total_valid} ({frac_mapped*100:.1f}%) mapped; no distance stats (no mapped UMIs)")
            else:
                add_log(f"{column}: no valid raw UMIs")
        # Apply strict filtering: drop any read where a present UMI failed to map
        total_rows_before = len(best_matches_df)
        if nonmapping_row_mask.any():
            rows_dropped_non_mapping = int(nonmapping_row_mask.sum())
            best_matches_df = best_matches_df.loc[~nonmapping_row_mask].reset_index(drop=True)
        else:
            rows_dropped_non_mapping = 0
        rows_kept_after_mapping = int(len(best_matches_df))
        add_log(f"Dropping {rows_dropped_non_mapping} reads due to non-mapping UMIs; kept {rows_kept_after_mapping}/{total_rows_before}.")
        # Debug print retained intentionally
        print(best_matches_df)

        df = best_matches_df
        # Finalize overall stats and log
        if mapping_stats["overall"]["total_raw"]:
            overall_frac = mapping_stats["overall"]["total_mapped"] / mapping_stats["overall"]["total_raw"]
        else:
            overall_frac = 0.0
        mapping_stats["overall"]["fraction_mapped"] = float(overall_frac)
        if overall_distances:
            mapping_stats["overall"]["distance_histogram"] = dict(Counter(overall_distances))
            mapping_stats["overall"]["mean_distance"] = float(np.mean(overall_distances))
            mapping_stats["overall"]["median_distance"] = float(np.median(overall_distances))
            add_log(
                f"Overall mapping: {mapping_stats['overall']['total_mapped']}/" \
                f"{mapping_stats['overall']['total_raw']} ({overall_frac*100:.1f}%) mapped; " \
                f"mean={mapping_stats['overall']['mean_distance']:.2f}, " \
                f"median={mapping_stats['overall']['median_distance']:.2f}, " \
                f"exact={mapping_stats['overall']['exact_matches']}"
            )
        else:
            add_log(
                f"Overall mapping: {mapping_stats['overall']['total_mapped']}/" \
                f"{mapping_stats['overall']['total_raw']} ({overall_frac*100:.1f}%) mapped; " \
                f"no distance stats (no mapped UMIs)"
            )
        # Mapping stats will be embedded in the main JSON log summary
        add_log("Whitelist mapping stats embedded in main log JSON.")



    # --- Standardized output: Include both raw and corrected UMIs for Hamming distance analysis
    # Format: Read Name, UMI1_Raw, UMI1, UMI2_Raw, UMI2, etc.
    # Where UMI1, UMI2, etc. are the corrected/mapped UMIs (backward compatible)
    
    # Merge raw and corrected UMIs into output dataframe
    output_df = df[["Read Name"]].copy()
    
    # Add raw and corrected UMI columns
    for col in umi_names:
        # Add raw UMI column (from raw_df)
        raw_col = f"{col}_Raw"
        if col in raw_df.columns:
            output_df[raw_col] = raw_df[col]
        else:
            # If raw column doesn't exist (shouldn't happen), use corrected as fallback
            output_df[raw_col] = df[col] if col in df.columns else None
        
        # Add corrected UMI column (from df - this is backward compatible)
        if col in df.columns:
            output_df[col] = df[col]
        else:
            output_df[col] = None
    
    # Ensure consistent data types for CSV output (compatible with downstream steps)
    # Convert UMI columns to string type, ensuring proper handling of NaN values
    for col in output_df.columns:
        if col != "Read Name":
            output_df[col] = output_df[col].apply(lambda x: str(x) if pd.notna(x) else pd.NA)
    
    # Ensure Read Name is string type
    output_df["Read Name"] = output_df["Read Name"].astype(str)

    # Build read_name -> corrected UMI values for BAM tagging (U1=UMI1, U2=UMI2, ...)
    # and read_name -> raw UMI values (R1=UMI1_Raw, R2=UMI2_Raw, ...)
    read_to_umis = {}
    read_to_raw_umis = {}

    # Corrected UMIs from df (current corrected/mapped table)
    for _, row in df[["Read Name"] + umi_names].iterrows():
        rn = str(row["Read Name"])
        vals = []
        for col in umi_names:
            v = row[col] if col in row else None
            if pd.isna(v) or v == "":
                vals.append(None)
            else:
                vals.append(str(v))
        read_to_umis[rn] = vals

    # Raw UMIs from raw_df (pre-correction snapshot)
    for _, row in raw_df[["Read Name"] + umi_names].iterrows():
        rn = str(row["Read Name"])
        vals = []
        for col in umi_names:
            v = row[col] if col in row else None
            if pd.isna(v) or v == "":
                vals.append(None)
            else:
                vals.append(str(v))
        read_to_raw_umis[rn] = vals

    # Write tagged BAM: same alignments with U1, U2, ... tags (primary output)
    tagged_bam_path = os.path.join(outputDir, f"{sampleName}_umi_tagged.bam")
    add_log(f"Writing UMI-tagged BAM: {tagged_bam_path}")
    with pysam.AlignmentFile(alignedReads, "rb") as bam_in, \
         pysam.AlignmentFile(tagged_bam_path, "wb", template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            qname = read.query_name
            umi_vals = read_to_umis.get(qname)
            raw_vals = read_to_raw_umis.get(qname)

            # Corrected UMIs: U1 = UMI1, U2 = UMI2, ...
            if umi_vals is not None:
                for i, val in enumerate(umi_vals):
                    if val is not None:
                        tag_name = f"U{i + 1}"
                        read.set_tag(tag_name, val, value_type="Z")

            # Raw UMIs: R1 = raw UMI1, R2 = raw UMI2, ...
            if raw_vals is not None:
                for i, val in enumerate(raw_vals):
                    if val is not None:
                        tag_name = f"R{i + 1}"
                        read.set_tag(tag_name, val, value_type="Z")
            bam_out.write(read)
    add_log("Finished writing UMI-tagged BAM.")

    if write_csv:
        csv_path = os.path.join(outputDir, f"{sampleName}_Paired_UMIs_read_names.csv")
        output_df.to_csv(csv_path, index=False, na_rep="")
        add_log(f"Wrote optional CSV: {csv_path}")

    add_log("Finished UMI extraction and mapping.")
    add_log("Begin Analysis...")

    # Write JSON log summary
    end_time = datetime.now(UTC).isoformat().replace('+00:00', 'Z')
    out_path = os.path.join(logOutputDir, f"{sampleName}_umi_extraction_map.log.json")
    # Include mapping/divergence stats if available (only in whitelist mapping mode)
    summary = {
        "sample": sampleName,
        "status": "success",
        "start_time": start_time,
        "end_time": end_time,
        "reference": reference,
        "alignment": alignedReads,
        "min_read_length": int(minReadLength),
        "barcode_motifs": list(barcodeMotifs),
        "whitelist_mode": wl_mode,
        "tagged_bam": tagged_bam_path,
        "max_distance": int(max_distance),
        "total_reads": int(total_reads),
        "reads_with_all_umis": int(all_umis_present_count),
        "unique_umis_when_all_present": unique_counts_all_present,
        "reads_with_only": umi_only_counts,
        "missing_umi_stats": missing_stats,
        "logs": log_entries
    }
    if write_csv:
        summary["output_csv"] = os.path.join(outputDir, f"{sampleName}_Paired_UMIs_read_names.csv")
    # Attach clustering or mapping stats depending on mode
    if wl_mode == "cluster":
        try:
            summary["clustering_stats"] = clustering_stats
        except NameError:
            summary["clustering_stats"] = {"error": "clustering stats unavailable"}
    elif wl_mode not in ("false",):
        try:
            summary["whitelist_mapping_stats"] = mapping_stats
            # No separate stats file; all stats are embedded in this summary
            # Include strict-filtering row counts if available
            try:
                summary["rows_dropped_non_mapping"] = int(rows_dropped_non_mapping)
                summary["rows_kept_after_mapping"] = int(rows_kept_after_mapping)
            except NameError:
                pass
        except NameError:
            # mapping_stats not defined if something failed prior to computation
            summary["whitelist_mapping_stats"] = {"error": "mapping stats unavailable"}
    with open(out_path, "w") as jf:
        json.dump(summary, jf, indent=2)


def export_umi_csv_from_bam(bam_path, output_csv_path):
    """
    Optional: Export a CSV of UMI data from a BAM file that has UMI tags (U1, U2, ...).
    Writes Read Name, UMI1, UMI2, ... in the same format as the legacy Paired_UMIs CSV.

    :param bam_path: Path to BAM file with U1, U2, ... tags (e.g. from umi_extraction_map).
    :param output_csv_path: Path for the output CSV file.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        max_ui = 0
        for read in bam.fetch(until_eof=True):
            i = 1
            while read.has_tag(f"U{i}"):
                i += 1
            max_ui = max(max_ui, i - 1)
    if max_ui == 0:
        raise ValueError(f"No reads with UMI tags (U1, U2, ...) found in BAM: {bam_path}")
    umi_cols = [f"UMI{j+1}" for j in range(max_ui)]
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qname = read.query_name
            if not read.has_tag("U1"):
                continue
            row = {"Read Name": qname}
            for j in range(max_ui):
                tag = f"U{j+1}"
                if read.has_tag(tag):
                    val = read.get_tag(tag)
                    row[umi_cols[j]] = val if val else ""
                else:
                    row[umi_cols[j]] = ""
            rows.append(row)
    if not rows:
        raise ValueError(f"No reads with U1 tag found in BAM: {bam_path}")
    df = pd.DataFrame(rows)
    out_dir = os.path.dirname(os.path.abspath(output_csv_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    df.to_csv(output_csv_path, index=False, na_rep="")
    return output_csv_path


def cluster_umis_umitools(umi_list, threshold=1, cluster_method="directional"):
    """
    Cluster UMIs using UMI-tools UMIClusterer API.
    
    Args:
        umi_list: List of UMI sequences (strings)
        threshold: Edit distance threshold for clustering (default: 1)
        cluster_method: Clustering method ("directional", "adjacency", etc.)
        
    Returns:
        tuple: (correction_map, cluster_info)
            - correction_map: dict mapping original UMI to corrected (primary) UMI
            - cluster_info: dict with cluster statistics
    """
    if not umi_list:
        return {}, {
            "num_clusters": 0,
            "num_singletons": 0,
            "cluster_sizes": []
        }
    
    # Create UMIClusterer instance
    umi_cluster = UMIClusterer(cluster_method=cluster_method)
    
    # Count UMIs
    umi_counter = Counter(umi_list)
    
    # Convert to bytes-encoded dict for UMI-tools (required by API)
    umi_count_dict = {}
    for umi, count in umi_counter.items():
        umi_count_dict[umi.encode("utf-8")] = count
    
    # Cluster UMIs using UMI-tools API
    clustered_umis = umi_cluster(umi_count_dict, threshold=threshold)
    
    # Create mapping: original UMI -> corrected (primary) UMI
    # The first UMI in each cluster is the primary (corrected) UMI
    correction_map = {}
    cluster_sizes = []
    num_singletons = 0
    
    for cluster in clustered_umis:
        primary_umi = cluster[0].decode('utf-8')
        cluster_size = len(cluster)
        cluster_sizes.append(cluster_size)
        
        if cluster_size == 1:
            num_singletons += 1
        
        # Map all UMIs in cluster to the primary UMI
        for umi_bytes in cluster:
            original_umi = umi_bytes.decode('utf-8')
            correction_map[original_umi] = primary_umi
    
    cluster_info = {
        "num_clusters": len(clustered_umis),
        "num_singletons": num_singletons,
        "cluster_sizes": cluster_sizes
    }
    
    return correction_map, cluster_info


def string_hamming(s1, s2):
    """Compute Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be the same length for Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def BK_scan(experimental_umis, whitelist, max_distance=3):
    """Find best approximate match in whitelist for each UMI (with dedup + exact-match fast path).
    
    Handles missing UMIs (NaN/None values) by preserving them in the output.
    """
    import pandas as pd
    
    # Ensure max_distance is an integer
    max_distance = int(max_distance)
    
    # Keep original input to preserve NaN values and their positions
    original_umis = list(experimental_umis)
    
    # Filter out NaN values from whitelist before building tree
    wl_unique = [w for w in dict.fromkeys(whitelist) if pd.notna(w)]
    if not wl_unique:
        # If whitelist is empty, return mapping that preserves NaN
        return {umi: (umi, "NoMatch") if pd.notna(umi) else (None, None) 
                for umi in original_umis}
    
    # Make sure all strings in whitelist are the same length
    if not all(len(umi) == len(wl_unique[0]) for umi in wl_unique):
        raise ValueError("All strings in the whitelist must be the same length.")
    
    # Filter valid (non-NaN) UMIs from experimental data for processing
    valid_umis = [umi for umi in original_umis if pd.notna(umi)]
    
    # Get expected UMI length from whitelist
    expected_length = len(wl_unique[0])
    
    wl_set = set(wl_unique)
    tree = BKTree(string_hamming, wl_unique)

    # Dedup queries but preserve order
    dedup = list(dict.fromkeys(valid_umis))

    # Scan only once per unique UMI (only valid UMIs, no NaN here)
    uniq_results = {}
    for umi in dedup:
        # Check if UMI has the correct length - if not, mark as NoMatch
        if len(umi) != expected_length:
            uniq_results[umi] = (umi, "LengthMismatch")
            continue
            
        # O(1) exact match
        if umi in wl_set:
            uniq_results[umi] = (umi, 0)
            continue

        try:
            matches = tree.find(umi, max_distance)
        except (ValueError, TypeError) as e:
            # Handle edge cases where hamming distance fails
            print(f"Warning: Could not find matches for UMI '{umi}' (length {len(umi)}): {e}")
            print(f"Expected length: {expected_length}, Whitelist sample: {wl_unique[0]}")
            uniq_results[umi] = (umi, "SearchError")
            continue
            
        if not matches:
            uniq_results[umi] = (umi, "NoMatch")
        else:
            # deterministic tie-break: smallest distance, then lexicographically smallest UMI
            min_dist = min(d for d, _m in matches)
            best = min(m for d, m in matches if d == min_dist)
            uniq_results[umi] = (best, min_dist)

    # Expand back to all input UMIs (including NaN values from original input)
    # For NaN keys, we need special handling since NaN != NaN in Python
    result = {}
    for umi in original_umis:
        if pd.isna(umi):
            result[umi] = (None, None)
        else:
            result[umi] = uniq_results[umi]
    
    return result

def load_whitelist(whitelist_file):
    """Load whitelist CSV into two lists of UMIs."""
    try:
        whitelist = pd.read_csv(whitelist_file)
        if whitelist.empty:
            raise ValueError("Whitelist file missing required columns or empty.")
        return whitelist
    except FileNotFoundError:
        raise ValueError(f"Whitelist file {whitelist_file} not found.")
