from belt.scripts.utils import log
from belt.scripts.utils import calculator
from belt.scripts.extract import *
from umi_tools import UMIClusterer
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

def umi_whitelist(reference, alignedReads, minReadLength, barcodeMotif, outputDir, clusterMethod ="default"):


#This portion generates the output directory structure and checks if the necessary directories exist.
    # Initialize logging
    logger = log()
    logger.log_message(f"UMI extraction of {alignedReads} with barcode {barcodeMotif} motif started.")
#Parse the input parameters
    sampleName = alignedReads.split("/")[-1].split(".")[0]
    outputDir = outputDir + "/" + sampleName
    # Check if the output directory exists, if not create it
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    # Check if the html subdirectory exists, if not create it
    html_dir = os.path.join(outputDir, "html")
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)

    logger = log()
    calc = calculator()
    logger.log_message("Starting UMI whitelist creation...")

#Extracting UMIs from the aligned reads based on the provided barcode motif.
    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(reference, alignedReads,barcodeMotif, minReadLength)
    logger.log_message("Clustering UMIs...")

# Clustering the UMIs and creating a corrected list
    # TODO: Add more clustering methods
    if clusterMethod == "default":
        umi5p_clusters = umi_clustering(raw_umi5_list)
        umi3p_clusters = umi_clustering(raw_umi3_list)
    logger.log_message(f"UMI correction on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")
    corrected_UMI_dict = correct_UMI_dict(raw_UMI_dict, umi5p_clusters, umi3p_clusters)
    logger.log_message(f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

    write_csv(["Read Name", "5pUMI", "3pUMI"], [[read_name for read_name in corrected_UMI_dict.keys()],
                                                  [list(umi_pair.keys())[0] for umi_pair in
                                                   corrected_UMI_dict.values()],
                                                  [list(umi_pair.values())[0] for umi_pair in
                                                   corrected_UMI_dict.values()]],
              f"{outputDir}/{sampleName}_correctedUMIs_read_names.csv")
    
    groouped_umi_dict = groupByReference(corrected_UMI_dict, reference="5pUMI")
    write_csv(["Reference UMI", "Satelitte UMI", "Count"], [[key for key in groouped_umi_dict.keys()],
                                              [satellite for satellite in
                                               [list(sat_dict.keys()) for sat_dict in groouped_umi_dict.values()]],
                                              [count for count in
                                               [list(sat_dict.values()) for sat_dict in groouped_umi_dict.values()]]],
              f"{outputDir}/{sampleName}_groupedUMIs_read_names.csv")
    logger.log_message(f"Grouping by reference UMI on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")




def umi_clustering(umi_list):
        umi_cluster = UMIClusterer(cluster_method="adjacency")
        umi_counter = Counter(umi_list)
        values = list(umi_counter.values())
        keys = list(umi_counter.keys())
        bit_keys = [umi.encode("utf-8") for umi in keys]
        umi_count_dict = dict(zip(bit_keys, values))
        clustered_umis = umi_cluster(umi_count_dict, 1)

        # Create a dictionary with corrected UMIs as keys and uncorrected UMIs as values (list)
        umi_corrected_dict = {}
        for cluster in clustered_umis:
            primary_umi = cluster[0].decode('utf-8')  # The primary (corrected) UMI
            umi_corrected_dict[primary_umi] = [umi.decode('utf-8') for umi in cluster]  # All UMIs in the cluster

        return umi_corrected_dict

def correct_UMI_dict(raw_UMI_dict, correction_dict_5p, correction_dict_3p):
    corrected_dict = {}
    for read_name, umi_pair in raw_UMI_dict.items():
        umi5p, umi3p = umi_pair
        # Correct 5p UMI
        corrected_5p = next(
            (corr for corr, uncorrected_list in correction_dict_5p.items() if umi5p in uncorrected_list),
            umi5p
        )
        # Correct 3p UMI
        corrected_3p = next(
            (corr for corr, uncorrected_list in correction_dict_3p.items() if umi3p in uncorrected_list),
            umi3p
        )
        corrected_dict[read_name] = {corrected_5p: corrected_3p}
    return corrected_dict

def groupByReference(umi_dict, reference="5pUMI"):
    """
    Groups UMIs by their reference sequence.
    :param umi_dict: Dictionary with read names as keys and UMIs as values.
    :param reference: Reference sequence to group by (default is "5pUMI").
    :return: Dictionary with grouped UMIs.
    """
    grouped_umis = {}
    for read_name, umi_pair in umi_dict.items():
        # umi_pair is a dict: {corrected_5p: corrected_3p}
        corrected_5p = list(umi_pair.keys())[0]
        corrected_3p = list(umi_pair.values())[0]
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