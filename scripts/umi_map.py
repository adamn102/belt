from belt.scripts.extract import *
from pybktree import BKTree
import pandas as pd


def umi_extraction_map(reference, alignedReads, minReadLength, barcodeMotif, whitelist, outputDir):
    # Parse the input parameters
    sampleName = alignedReads.split("/")[-1].split(".")[0]

    # Check if the output directory exists, if not create it
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    logger = log()
    logger.log_message("Extracting UMIs from aligned reads...")

    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(reference, alignedReads,
                                                                                     barcodeMotif, minReadLength)
    if raw_umi5_list is None or raw_umi3_list is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return

    logger.log_message("Mapping raw UMIs to whitelist...")
    umi5_list, umi3_list = load_whitelist(whitelist)

    # create a mapping dictionary using the raw umis and the white list
    # return the minimal set of raw umis
    raw_umi5_set = list(set(raw_umi5_list))
    raw_umi3_set = list(set(raw_umi3_list))

    # Create BK-trees for fast approximate matching
    results5p = BK_scan(raw_umi5_set, umi5_list, max_distance=3)
    results3p = BK_scan(raw_umi3_set, umi3_list, max_distance=3)

    # write the updated UMIs to a csv file
    out_df5p = pd.DataFrame.from_dict(results5p, orient='index', columns=['Best Match', 'Distance'])
    out_df5p.index.name = '5pUMI'
    out_df5p.reset_index(inplace=True)
    # fill blank values with 0
    out_df5p['Distance'] = out_df5p['Distance'].fillna("NA")
    out_df5p['Best Match'] = out_df5p['Best Match'].fillna("NA")
    out_df5p.to_csv(f"{outputDir}/{sampleName}_5pUMI_whitelist_matches.csv", index=True)

    out_df3p = pd.DataFrame.from_dict(results3p, orient='index', columns=['Best Match', 'Distance'])
    out_df3p.index.name = '3pUMI'
    out_df3p.reset_index(inplace=True)
    out_df3p['Distance'] = out_df3p['Distance'].fillna("NA")
    out_df3p['Best Match'] = out_df3p['Best Match'].fillna("NA")
    out_df3p.to_csv(f"{outputDir}/{sampleName}_3pUMI_whitelist_matches.csv", index=True)

    logger.log_message("Updating UMIs with best matches from whitelist...")
    # Update UMIs in raw_UMI_dict with the best matches from the whitelist
    # iterate though the dictionary and update the 5pUMI and 3pUMI with the best match
    for read_name, (umi5, umi3) in raw_UMI_dict.items():
        best_match_5p, distance_5p = results5p.get(umi5, (umi5, "NoBarcode"))
        best_match_3p, distance_3p = results3p.get(umi3, (umi3, "NoBarcode"))
        raw_UMI_dict[read_name] = (best_match_5p, distance_5p, best_match_3p, distance_3p)

    df = pd.DataFrame.from_dict(
        raw_UMI_dict,
        orient='index',
        columns=['5pUMI', '5pDistance', '3pUMI', '3pDistance']
    )
    df.index.name = 'Read Name'
    df.reset_index(inplace=True)

    df.to_csv(f"{outputDir}/{sampleName}_Paired_UMIs_read_names.csv", index=False)


def load_whitelist(whitelist_file):
    """
    Load whitelist from a file.
    return a list for each UMI column in the whitelist file.
    """
    try:
        whitelist = pd.read_csv(whitelist_file)
        whitelist_5pUMI = whitelist['5pUMI'].dropna()
        whitelist_3pUMI = whitelist['3pUMI'].dropna()
        if whitelist_5pUMI.empty or whitelist_3pUMI.empty:
            raise ValueError("Whitelist file is empty or missing required columns.")
        return whitelist_5pUMI, whitelist_3pUMI
    except FileNotFoundError:
        raise ValueError(f"Whitelist file {whitelist_file} not found.")


def string_hamming(s1, s2):
    """Compute Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be the same length for Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def BK_scan(experimental_umis, whitelist, max_distance=2):
    """
    Perform a BK-tree scan to find approximate matches for a given query.
    Returns a list of matches within the specified distance.
    """
    # Create BK-trees for fast approximate matching
    # return a dictionary with the raw read as the key
    # and the matched_umi and distance as values
    tree = BKTree(string_hamming, whitelist)

    results = dict()
    for umi in experimental_umis:
        matches = tree.find(umi, max_distance)
        if not matches:
            results[umi] = (umi, "NoMatch")
        else:
            # Pick the match with minimum distance
            best_match = min(matches, key=lambda x: x[0])
            distance, matched_umi = best_match
            results[umi] = (matched_umi, distance)

    return results
