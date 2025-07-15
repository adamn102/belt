from belt.scripts.utils import log
from belt.scripts.utils import calculator
from belt.scripts.extract import *
from umi_tools import UMIClusterer
from collections import Counter


def umi_whitelist(reference, alignedReads, minReadLength, barcodeMotif, clusterMethod, outputDir):


    ## TODO: Flag outliers and filter them out (write them out in a report/log).


    # Parse the input parameters
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
    logger.log_message("Extracting UMIs from aligned reads...")
    # find the raw UMIs

    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(reference, alignedReads,
                                                                                     barcodeMotif, minReadLength)

    logger.log_message(
        f"Motif5 found at index {motif5_index}, Motif3 found at index {motif3_index} in reference sequence.")
    if raw_umi5_list is None or raw_umi3_list is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return
    logger.log_message(f"UMI extraction of {sampleName} with barcode {barcodeMotif} motif finished.")

    logger.log_message("Clustering UMIs...")

    # Clustering the UMIs and creating a corrected list
    # TODO: Add more clustering methods
    if clusterMethod == "default":
        umi5p_clusters = umi_clustering(raw_umi5_list)
        umi3p_clusters = umi_clustering(raw_umi3_list)

    logger.log_message(f"UMI clustering of {sampleName} with barcode {barcodeMotif} motif finished.")
    corrected_umi5_list = correctUMIS(raw_umi5_list, umi5p_clusters)
    corrected_umi3_list = correctUMIS(raw_umi3_list, umi3p_clusters)
    logger.log_message(f"UMI correction on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")
    corrected_UMI_dict = correct_UMI_dict(raw_UMI_dict, umi5p_clusters, umi3p_clusters)
    logger.log_message(
        f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

    # report the UMI Statistics
    umi5p_stats = calculate_umi_stats(umi5p_clusters, corrected_umi5_list, raw_umi5_list)
    umi3p_stats = calculate_umi_stats(umi3p_clusters, corrected_umi3_list, raw_umi3_list)
    logger.log_message(f"UMI statistics calculation of {sampleName} with barcode {barcodeMotif} motif finished.")

    # write csv files with the clustered UMIs
    write_csv(["5pUMI", "3pUMI"], [corrected_umi5_list, corrected_umi3_list],
              f"{outputDir}/{sampleName}_correctedUMIs.csv")

    write_csv(["Read Name", "5pUMI", "3pUMI"], [[read_name for read_name in corrected_UMI_dict.keys()],
                                                  [list(umi_pair.keys())[0] for umi_pair in
                                                   corrected_UMI_dict.values()],
                                                  [list(umi_pair.values())[0] for umi_pair in
                                                   corrected_UMI_dict.values()]],
              f"{outputDir}/{sampleName}_correctedUMIs_read_names.csv")

    logger.log_message(f"UMI CSV files of {sampleName} with barcode {barcodeMotif} motif created.")
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi5p_clusters).keys()), sample_size=400),
                                f"{outputDir}/html/hamming_distances_5p.html", 400)
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi3p_clusters).keys()), sample_size=400),
                                f"{outputDir}/html/hamming_distances_3p.html", 400)
    calc.plot_umi_frequencies(corrected_umi5_list, f"{outputDir}/html/umi_freq_5p.html")
    calc.plot_umi_frequencies(corrected_umi3_list, f"{outputDir}/html/umi_freq_3p.html")
    calc.hamming_distance_matrix(list(umi5p_clusters.keys()), f"{outputDir}/html/umi_freq_matrix_5p.html",
                                 sample_size=400)
    calc.hamming_distance_matrix(list(umi3p_clusters.keys()), f"{outputDir}/html/umi_freq_matrix_3p.html",
                                 sample_size=400)
    logger.log_message(f"UMI plots of {sampleName} with barcode {barcodeMotif} motif created.")

    createHTML(umi5p_stats, umi3p_stats, f"{outputDir}", sampleName, pi_photo_path="../alex.png")
    logger.log_message(f"HTML report created at {outputDir}/{sampleName}_Report.html")
    logger.write_log(f"{outputDir}/{sampleName}_log.txt")


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


def correctUMIS(umi_list, correction_dict):
    corrected_umi_list = [
        next(
            (corrected for corrected, uncorrected_list in correction_dict.items() if umi in uncorrected_list),
            umi  # Default to the original UMI if no correction is found
        )
        for umi in umi_list
    ]
    return corrected_umi_list


def calculate_umi_stats(umi_dict, correctedList, rawList):
    """
    Calculate statistics for the UMIs.
    :param correctedList: List of corrected UMIs
    :param umi_dict: Dictionary of UMIs where keys are corrected UMIs and values are lists of uncorrected UMIs
    :param rawList: List of raw UMIs
    :return: List of statistics [raw_count, clusters, chao1, chao1_updated, entropy]
    """
    calc = calculator()

    raw_count = len(rawList)
    clusters = len(list(umi_dict.keys()))
    chao1 = calc.chao1(correctedList)
    chao1_updated = calc.balazs_chao1(correctedList)
    shannon_entropy = calc.shannon_entropy(correctedList)

    return raw_count, clusters, chao1, chao1_updated, shannon_entropy


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


def createHTML(umi5p_stats, umi3p_stats, output_dir, sample_name, pi_photo_path):
    """
    Creates a styled HTML summary report.
    """
    import os
    import shutil

    # Copy alex.png to the html directory
    dest_path = os.path.join(output_dir, "html", "alex.png")
    shutil.copy2(pi_photo_path, dest_path)

    # Update photo path to be relative to the HTML file
    pi_photo_path = "html/alex.png"

    # Read HTML template from file
    template_path = os.path.join(os.path.dirname(__file__), 'template.html')
    with open(template_path, 'r') as template_file:
        template = template_file.read()
    # Round umi5p_stats and umi3p_stats to 2 digits
    umi5p_stats = [round(x, 2) if isinstance(x, float) else x for x in umi5p_stats]
    umi3p_stats = [round(x, 2) if isinstance(x, float) else x for x in umi3p_stats]
    # Format template with values
    html_content = template.format(
        sample_name=sample_name,
        pi_photo_path=pi_photo_path,
        umi5p_raw_count=umi5p_stats[0],
        umi5p_clusters=umi5p_stats[1],
        umi5p_chao1=umi5p_stats[2],
        umi5p_chao1_updated=umi5p_stats[3],
        umi5p_entropy=umi5p_stats[4],
        umi3p_raw_count=umi3p_stats[0],
        umi3p_clusters=umi3p_stats[1],
        umi3p_chao1=umi3p_stats[2],
        umi3p_chao1_updated=umi3p_stats[3],
        umi3p_entropy=umi3p_stats[4],
        output_dir=output_dir
    )

    # Write formatted HTML to output file
    with open(f"{output_dir}/{sample_name}_Report.html", "w") as html_file:
        html_file.write(html_content)
