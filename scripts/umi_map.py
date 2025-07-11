from belt.scripts.utils import log
from belt.scripts.extract import *
from pybktree import BKTree
import pandas as pd
import os

def umi_extraction_map(
    reference,
    alignedReads,
    minReadLength,
    barcodeMotif,
    whitelist,
    outputDir,
    max_distance=2
):
    sampleName = os.path.basename(alignedReads).split(".")[0]

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

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

    # Handle different whitelist modes
    if whitelist.lower() == "false":
        logger.log_message("Whitelist set to 'false'. Using raw UMIs without mapping.")
        # Here, just assign raw UMI as 'Best Match' with Distance=0
        results5p = {umi: (umi, 0) for umi in raw_umi5_set}
        results3p = {umi: (umi, 0) for umi in raw_umi3_set}

    elif whitelist.lower() == "cluster":
        logger.log_message("Whitelist mode set to 'cluster'. Performing UMI-tools clustering.")
        # **Placeholder for clustering**
        # You would replace this block with real clustering logic (e.g., umi_tools)
        # For now, treat as "no clustering" to avoid breaking the pipeline.
        results5p = {umi: (umi, "Clustered") for umi in raw_umi5_set}
        results3p = {umi: (umi, "Clustered") for umi in raw_umi3_set}

    else:
        logger.log_message("Mapping raw UMIs to whitelist...")
        umi5_list, umi3_list = load_whitelist(whitelist)

        results5p = BK_scan(raw_umi5_set, umi5_list, max_distance)
        results3p = BK_scan(raw_umi3_set, umi3_list, max_distance)

        # Save whitelist match tables
        out_df5p = pd.DataFrame.from_dict(results5p, orient="index", columns=["Best Match", "Distance"])
        out_df5p.index.name = "5pUMI"
        out_df5p.reset_index(inplace=True)
        out_df5p.to_csv(f"{outputDir}/{sampleName}_5pUMI_whitelist_matches.csv", index=False)

        out_df3p = pd.DataFrame.from_dict(results3p, orient="index", columns=["Best Match", "Distance"])
        out_df3p.index.name = "3pUMI"
        out_df3p.reset_index(inplace=True)
        out_df3p.to_csv(f"{outputDir}/{sampleName}_3pUMI_whitelist_matches.csv", index=False)

    logger.log_message("Updating UMIs with best matches...")

    # Update final dictionary of read -> (best matches)
    for read_name, (umi5, umi3) in raw_UMI_dict.items():
        best_match_5p, distance_5p = results5p.get(umi5, (umi5, "NA"))
        best_match_3p, distance_3p = results3p.get(umi3, (umi3, "NA"))
        raw_UMI_dict[read_name] = (best_match_5p, distance_5p, best_match_3p, distance_3p)

    df = pd.DataFrame.from_dict(
        raw_UMI_dict,
        orient="index",
        columns=["5pUMI", "5pDistance", "3pUMI", "3pDistance"]
    )
    df.index.name = "Read Name"
    df.reset_index(inplace=True)

    df.to_csv(f"{outputDir}/{sampleName}_Paired_UMIs_read_names.csv", index=False)

    logger.log_message("✅ Finished UMI extraction and mapping.")


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


def string_hamming(s1, s2):
    """Compute Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be the same length for Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def BK_scan(experimental_umis, whitelist, max_distance=2):
    """Find best approximate match in whitelist for each UMI."""
    tree = BKTree(string_hamming, whitelist)
    results = dict()
    for umi in experimental_umis:
        matches = tree.find(umi, max_distance)
        if not matches:
            results[umi] = (umi, "NoMatch")
        else:
            best_match = min(matches, key=lambda x: x[0])
            distance, matched_umi = best_match
            results[umi] = (matched_umi, distance)
    return results
