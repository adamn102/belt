from belt.scripts.utils import log
from belt.scripts.extract import *
from belt.scripts.utils import umi_clustering, correct_UMI_dict
from pybktree import BKTree
from umi_tools import UMIClusterer
from collections import Counter
import pandas as pd
import os
import numpy as np


def umi_extraction_map(
        reference,
        alignedReads,
        minReadLength,
        barcodeMotif,
        whitelist,
        outputDir,
        analysisOut,
        max_distance=1
):
    sampleName = os.path.basename(alignedReads).split(".")[0]

    # Initialize output directories
    os.makedirs(outputDir, exist_ok=True)
    os.makedirs(analysisOut, exist_ok=True)

    # Set up logging
    logger = log()
    logger.log_message("Extracting UMIs from aligned reads...")

    # Extract raw UMIs
    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(
        reference, alignedReads, barcodeMotif, minReadLength
    )
    if raw_umi5_list is None or raw_umi3_list is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return

    raw_umi5_set = list(set(raw_umi5_list))
    raw_umi3_set = list(set(raw_umi3_list))

    corrected_UMI_dict = None  # will become {read_name: (5pUMI, 3pUMI)}

    # Handle different whitelist modes
    wl_mode = str(whitelist).lower()
    if wl_mode == "false":
        logger.log_message("Whitelist set to 'false'. Using raw UMIs without mapping.")
        # Pass through raw UMIs per read
        corrected_UMI_dict = {rn: (umi5, umi3) for rn, (umi5, umi3) in raw_UMI_dict.items()}

    elif wl_mode == "cluster":
        logger.log_message("Whitelist mode set to 'cluster'. Performing UMI-tools clustering.")
        umi5p_clusters = umi_clustering(raw_umi5_list)
        umi3p_clusters = umi_clustering(raw_umi3_list)

        logger.log_message(f"UMI correction on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")
        corrected_UMI_dict = correct_UMI_dict(raw_UMI_dict, umi5p_clusters, umi3p_clusters)
        logger.log_message(
            f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

    else:
        logger.log_message("Mapping raw UMIs to whitelist...")
        wl_5, wl_3 = load_whitelist(whitelist)

        results5p = BK_scan(raw_umi5_set, wl_5, max_distance)
        results3p = BK_scan(raw_umi3_set, wl_3, max_distance)

        # Save whitelist match tables (optional / for QC)
        out_df5p = pd.DataFrame.from_dict(results5p, orient="index", columns=["Best Match", "Distance"])
        out_df5p.index.name = "5pUMI"
        out_df5p.reset_index(inplace=True)
        out_df5p.to_csv(f"{outputDir}/{sampleName}_5pUMI_whitelist_matches.csv", index=False)

        out_df3p = pd.DataFrame.from_dict(results3p, orient="index", columns=["Best Match", "Distance"])
        out_df3p.index.name = "3pUMI"
        out_df3p.reset_index(inplace=True)
        out_df3p.to_csv(f"{outputDir}/{sampleName}_3pUMI_whitelist_matches.csv", index=False)

        logger.log_message("Updating UMIs with best matches (ignoring distance in final output)...")
        # Build corrected dict that mirrors cluster output (no distance columns)
        corrected_UMI_dict = {}
        for read_name, (umi5, umi3) in raw_UMI_dict.items():
            best5, _d5 = results5p.get(umi5, (umi5, "NA"))
            best3, _d3 = results3p.get(umi3, (umi3, "NA"))
            corrected_UMI_dict[read_name] = (best5, best3)

    # --- Standardized output (always like 'cluster'): Read Name, 5pUMI, 3pUMI ---
    df = pd.DataFrame.from_dict(
        corrected_UMI_dict,
        orient="index",
        columns=["5pUMI", "3pUMI"]
    )
    df.index.name = "Read Name"
    df.reset_index(inplace=True)

    # Write final paired UMI results (uniform schema)
    df.to_csv(f"{outputDir}/{sampleName}_Paired_UMIs_read_names.csv", index=False)

    logger.log_message("✅ Finished UMI extraction and mapping.")
    logger.log_message("Begin Analysis...")

    # Grouped UMI Analysis (unchanged, now safely uses corrected_UMI_dict)
    grouped_umi_dict = groupByReference(corrected_UMI_dict, reference="5pUMI")

    # Flatten grouped dict to rows
    rows = [(outer_key, inner_key, count)
            for outer_key, inner_dict in grouped_umi_dict.items()
            for inner_key, count in inner_dict.items()]
    gdf = pd.DataFrame(rows, columns=["5pUMI", "3pUMI", "count"])

    # Clean + add pair id
    gdf = gdf[(gdf["5pUMI"] != "") & (gdf["3pUMI"] != "")]
    gdf["UMI_pair"] = gdf["5pUMI"] + "_" + gdf["3pUMI"]
    gdf = gdf[gdf["UMI_pair"] != ""]

    gdf.to_csv(f"{analysisOut}/{sampleName}_groupedUMIs.csv", index=False)

    filtered_umi_dict, scrapped_dict = filterGroups(grouped_umi_dict, COUNT_THRESHOLD=10)
    logger.log_message(
        f"Grouping by reference UMI on the lists of {sampleName} with barcode {barcodeMotif} motif finished."
    )


def groupByReference(umi_dict, reference="5pUMI"):
    """
    Groups UMIs by their reference sequence.
    umi_dict: {read_name: (5pUMI, 3pUMI)}
    Returns: {keyUMI: {satelliteUMI: count, ...}, ...}
    """
    grouped_umis = {}
    for read_name, umi_pair in umi_dict.items():
        corrected_5p, corrected_3p = umi_pair
        if reference == "5pUMI":
            key, satellite = corrected_5p, corrected_3p
        elif reference == "3pUMI":
            key, satellite = corrected_3p, corrected_5p
        else:
            print("Invalid reference UMI. Choose '5pUMI' or '3pUMI'. Defaulting to '5pUMI'.")
            key, satellite = corrected_5p, corrected_3p
        if key not in grouped_umis:
            grouped_umis[key] = {}
        grouped_umis[key][satellite] = grouped_umis[key].get(satellite, 0) + 1
    return grouped_umis


def filterGroups(grouped_dict, COUNT_THRESHOLD=5):
    """
    Filters groups based on a threshold. Also calculates the normalized entropy for each group.
    Returns (filtered_dict, scrapped_dict) where each value has {"satellites": {...}, "entropy": float}.
    """
    filtered_dict = {}
    scrapped_dict = {}

    for key, satellites in grouped_dict.items():
        counts = np.array(list(satellites.values()))
        if counts.sum() > 0:
            probabilities = counts / counts.sum()
            # normalized entropy (0..1)
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10)) / np.log(len(counts))
            if entropy < 1e-10 or len(counts) == 1:
                entropy = 0.0
        else:
            entropy = 0.0

        entry = {"satellites": satellites, "entropy": float(entropy)}
        if sum(satellites.values()) >= COUNT_THRESHOLD and key != "" and entropy < 0.66:
            filtered_dict[key] = entry
        else:
            scrapped_dict[key] = entry

    return filtered_dict, scrapped_dict


def string_hamming(s1, s2):
    """Compute Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be the same length for Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


# def BK_scan(experimental_umis, whitelist, max_distance=1):
#     """Find best approximate match in whitelist for each UMI."""
#     tree = BKTree(string_hamming, whitelist)
#     results = dict()
#     for umi in experimental_umis:
#         matches = tree.find(umi, max_distance)
#         if not matches:
#             results[umi] = (umi, "NoMatch")
#         else:
#             best_match = min(matches, key=lambda x: x[0])
#             distance, matched_umi = best_match
#             results[umi] = (matched_umi, distance)
#     return results

def BK_scan(experimental_umis, whitelist, max_distance=1):
    """Find best approximate match in whitelist for each UMI (with dedup + exact-match fast path)."""
    # Dedup whitelist too (just in case) and build helpers
    wl_unique = list(dict.fromkeys(whitelist))
    wl_set = set(wl_unique)
    print(len(wl_set), "unique whitelist UMIs")
    tree = BKTree(string_hamming, wl_unique)

    # Dedup queries but preserve order
    dedup = list(dict.fromkeys(experimental_umis))
    print(len(dedup), "unique experimental UMIs")

    # Scan only once per unique UMI
    uniq_results = {}
    for umi in dedup:
        # O(1) exact match
        if umi in wl_set:
            uniq_results[umi] = (umi, 0)
            continue

        matches = tree.find(umi, max_distance)
        if not matches:
            uniq_results[umi] = (umi, "NoMatch")
        else:
            # deterministic tie-break: smallest distance, then lexicographically smallest UMI
            min_dist = min(d for d, _m in matches)
            best = min(m for d, m in matches if d == min_dist)
            uniq_results[umi] = (best, min_dist)

    # Expand back to all input UMIs (duplicates included)
    return {umi: uniq_results[umi] for umi in experimental_umis}

def load_whitelist(whitelist_file):
    """Load whitelist CSV into two lists of UMIs."""
    try:
        whitelist = pd.read_csv(whitelist_file)
        whitelist_5pUMI = whitelist["5pUMI"].dropna()
        whitelist_3pUMI = whitelist["3pUMI"].dropna()
        if whitelist_5pUMI.empty or whitelist_3pUMI.empty:
            raise ValueError("Whitelist file missing required columns or empty.")
        return whitelist_5pUMI, whitelist_3pUMI
    except FileNotFoundError:
        raise ValueError(f"Whitelist file {whitelist_file} not found.")
