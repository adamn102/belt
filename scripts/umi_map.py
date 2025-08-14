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
    max_distance=2
):
    sampleName = os.path.basename(alignedReads).split(".")[0]

    # Initialize output directory
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Initialize analysis directory
    if not os.path.exists(analysisOut):
        os.makedirs(analysisOut)

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

    # Handle different whitelist modes
    if whitelist.lower() == "false":
        logger.log_message("Whitelist set to 'false'. Using raw UMIs without mapping.")
        # Here, just assign raw UMI as 'Best Match' with Distance=0
        results5p = {umi: (umi, 0) for umi in raw_umi5_set}
        results3p = {umi: (umi, 0) for umi in raw_umi3_set}

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

    ## TODO - implement clustering logic
    ## TODO - update command to be more specific
    elif whitelist.lower() == "cluster":
        logger.log_message("Whitelist mode set to 'cluster'. Performing UMI-tools clustering.")
        umi5p_clusters = umi_clustering(raw_umi5_list)
        umi3p_clusters = umi_clustering(raw_umi3_list)

        logger.log_message(f"UMI correction on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")
        corrected_UMI_dict = correct_UMI_dict(raw_UMI_dict, umi5p_clusters, umi3p_clusters)

        logger.log_message(
            f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

        df = pd.DataFrame.from_dict(
            corrected_UMI_dict,
            orient="index",
            columns=["5pUMI", "3pUMI"]
        )
        df.index.name = "Read Name"
        df.reset_index(inplace=True)



        # Create results dictionaries with clustered UMIs

        results5p = {umi: (umi, "Clustered") for umi in umi5p_clusters}
        results3p = {umi: (umi, "Clustered") for umi in umi3p_clusters}

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

    # Write final paired UMI results
    df.to_csv(f"{outputDir}/{sampleName}_Paired_UMIs_read_names.csv", index=False)

    logger.log_message("✅ Finished UMI extraction and mapping.")
    logger.log_message("Begin Analysis...")

    # Grouped UMI Analysis
    grouped_umi_dict = groupByReference(corrected_UMI_dict, reference="5pUMI")

    # flatten the list of tuples
    rows = [(outer_key, inner_key, count)
            for outer_key, inner_dict in grouped_umi_dict.items()
            for inner_key, count in inner_dict.items()]

    df = pd.DataFrame(rows, columns=["5pUMI", "3pUMI", "count"])

    # drop rows with a blank value
    df = df[df["5pUMI"] != ""]
    df = df[df["3pUMI"] != ""]

    # concatentate the 5pUMI and 3pUMI columns to create a unique identifier
    df["UMI_pair"] = df["5pUMI"] + "_" + df["3pUMI"]

    # drop rows with a blank UMI_pair
    df = df[df["UMI_pair"] != ""]

    # write the DataFrame to a CSV file
    df.to_csv(f"{analysisOut}/{sampleName}_groupedUMIs.csv", index=False)

    filtered_umi_dict, scrapped_dict = filterGroups(grouped_umi_dict, COUNT_THRESHOLD=10)
    logger.log_message(
        f"Grouping by reference UMI on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")


    # # Write csv files with the filtered UMIs
    # with open(f"{analysisOut}/{sampleName}_filteredUMIs.tsv", "w") as f:
    #     f.write("Reference UMI\tSatellite UMIs\tCounts\tEntropy\n")
    #     for key, entry in filtered_umi_dict.items():
    #         satellites = entry["satellites"]
    #         entropy = entry["entropy"]
    #         satellite_umis = repr(list(satellites.keys()))
    #         counts = repr(list(satellites.values()))
    #         f.write(f"{key}\t{satellite_umis}\t{counts}\t{entropy}\n")
    #
    # # Write csv files with the scrapped UMIs
    # with open(f"{analysisOut}/{sampleName}_scrappedUMIs.tsv", "w") as f:
    #     f.write("Reference UMI\tSatellite UMIs\tCounts\tEntropy\n")
    #     for key, entry in scrapped_dict.items():
    #         satellites = entry["satellites"]
    #         entropy = entry["entropy"]
    #         satellite_umis = repr(list(satellites.keys()))
    #         counts = repr(list(satellites.values()))
    #         f.write(f"{key}\t{satellite_umis}\t{counts}\t{entropy}\n")
    # logger.log_message(f"Grouping by reference UMI on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")
    #
    # stats = {
    # "Sample": sampleName,
    # "Total input reads": len(raw_UMI_dict),
    # "Passed filter (with both UMIs)": len(corrected_UMI_dict),
    # "Unique raw 5' UMIs": len(set(raw_umi5_list)),
    # "Unique raw 3' UMIs": len(set(raw_umi3_list)),
    # "5' UMI clusters": len(umi5p_clusters),
    # "3' UMI clusters": len(umi3p_clusters),
    # "UMI groups formed (5' as ref)": len(grouped_umi_dict),
    # "Filtered UMI groups (kept)": len(filtered_umi_dict),
    # "Scrapped UMI groups": len(scrapped_dict),
    # "Whitelist written": f"{analysisOut}/{sampleName}_whitelistedUMIs.csv"
    # }
    #
    # logger.whitelist_summary_statistics(stats)




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


def groupByReference(umi_dict, reference="5pUMI"):
    """
    Groups UMIs by their reference sequence.
    - filtered_umi_dict: dict of structure {
          "reference_UMI": {
              "satellites": {"sat1": count1, "sat2": count2, ...},
              "entropy": float
          },
          ...
      }
    :param umi_dict: Dictionary with read names as keys and UMIs as values.
    :param reference: Reference sequence to group by (default is "5pUMI").
    :return: Dictionary with grouped UMIs.
    """
    grouped_umis = {}
    for read_name, umi_pair in umi_dict.items():
        # umi_pair is a list
        corrected_5p = umi_pair[0]
        corrected_3p = umi_pair[1]
        if reference == "5pUMI":
            key = corrected_5p
            satellite = corrected_3p

        elif reference == "3pUMI":
            key = corrected_3p
            satellite = corrected_5p

        else:
            print("Invalid reference UMI. Choose '5pUMI' or '3pUMI'. Defaulting to '5pUMI'.")
            key = corrected_5p
            satellite = corrected_3p
        if key not in grouped_umis:
            grouped_umis[key] = {}
        if satellite not in grouped_umis[key]:
            grouped_umis[key][satellite] = 0
        grouped_umis[key][satellite] += 1
    return grouped_umis


def filterGroups(grouped_dict, COUNT_THRESHOLD=5):
    """
    Filters groups based on a threshold. Also calculates the entropy for each group.
    Adds the entropy value to each key in the returned dictionaries.
    :param grouped_dict: Dictionary with grouped UMIs.
    :return: Filtered dictionary and scrapped_dict.
    """
    filtered_dict = {}
    scrapped_dict = {}
    entropy_list = []

    for key, satellites in grouped_dict.items():
        counts = np.array(list(satellites.values()))
        if counts.sum() > 0:
            probabilities = counts / counts.sum()
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))/np.log(len(counts))
            if entropy < 1e-10 or len(counts) == 1:
                entropy = 0
            else:
                entropy_list.append(entropy)
        else:
            entropy = 0

        entry = {"satellites": satellites, "entropy": entropy}
        if sum(satellites.values()) >= COUNT_THRESHOLD and key != "" and satellites != "" and entropy < 0.66:
            filtered_dict[key] = entry
        else:
            scrapped_dict[key] = entry
        return filtered_dict, scrapped_dict
