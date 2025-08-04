from belt.scripts.utils import log
from belt.scripts.utils import calculator
from belt.scripts.extract import *
from umi_tools import UMIClusterer
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

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

    logger.log_message(f"UMI extraction is being filtered for outlier barcodes...")
    # Filter out outlier barcodes
    raw_umi5_list, raw_umi3_list, raw_UMI_dict = filter_outlier_barcodes(raw_umi5_list, raw_umi3_list, raw_UMI_dict)

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
    logger.log_message(f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

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

    createHTML(umi5p_stats, umi3p_stats, f"{outputDir}", sampleName, pi_photo_path=os.path.join(os.path.dirname(__file__), "alex.png"))
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
    from scipy.stats import binom

    # Copy alex.png to the html directory
    dest_path = os.path.join(output_dir, "html", "alex.png")
    shutil.copy2(pi_photo_path, dest_path)

    # Update photo path to be relative to the HTML file
    pi_photo_path = os.path.join(os.path.dirname(__file__), "alex.png")

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


def filter_outlier_barcodes(raw_umi5_list, raw_umi3_list, raw_UMI_dict):
    """
    Filters out outlier barcodes based on the frequency of nucleotide sequences in the raw UMIs.
    Prints diagnostics of how much is filtered out and writes filtered UMIs to a CSV file.
    :param raw_umi5_list: List of raw 5' UMIs
    :param raw_umi3_list: List of raw 3' UMIs
    :param raw_UMI_dict: Dictionary of raw UMIs with read names as keys
    :return: Filtered lists and dictionary
    """

    std_dev_limit = 0.80  # Standard deviation threshold for outlier detection
    print("Filtering out outlier barcodes baseed on std:", std_dev_limit)
    # Count frequencies of each UMI
    umi5_counts = Counter(raw_umi5_list)
    umi3_counts = Counter(raw_umi3_list)

    # Define outlier threshold: remove UMIs with counts below (mean - 3*std)
    def get_threshold(counts):
        if not counts:
            return 0
        values = np.array(list(counts.values()))
        mean = values.mean()
        std = values.std()
        threshold = max(1, int(mean - std_dev_limit * std))
        return threshold

    threshold5 = get_threshold(umi5_counts)
    threshold3 = get_threshold(umi3_counts)

    # Identify UMIs to keep
    umi5_keep = {umi for umi, count in umi5_counts.items() if count >= threshold5 and count > 1}
    umi3_keep = {umi for umi, count in umi3_counts.items() if count >= threshold3 and count > 1}

    # Filter lists and dict
    filtered_umi5_list = [umi for umi in raw_umi5_list if umi in umi5_keep]
    filtered_umi3_list = [umi for umi in raw_umi3_list if umi in umi3_keep]

    filtered_UMI_dict = {}
    for read_name, umi_pair in raw_UMI_dict.items():
        umi5, umi3 = umi_pair
        if umi5 in umi5_keep and umi3 in umi3_keep:
            filtered_UMI_dict[read_name] = umi_pair

    # Diagnostics
    n5_before = len(raw_umi5_list)
    n5_after = len(filtered_umi5_list)
    n3_before = len(raw_umi3_list)
    n3_after = len(filtered_umi3_list)
    print(f"5' UMI: {n5_before-n5_after} filtered out of {n5_before} ({100*(n5_before-n5_after)/n5_before:.2f}%)")
    print(f"3' UMI: {n3_before-n3_after} filtered out of {n3_before} ({100*(n3_before-n3_after)/n3_before:.2f}%)")

    # Write filtered UMIs to CSV file in current working directory
    output_csv = os.path.join(os.getcwd(), "filtered_umis.csv")
    import plotly.graph_objs as go
    import plotly.offline as pyo

    # Visualize 5' UMI counts
    from scipy.stats import binom
    umi5_values = np.array(list(umi5_counts.values()))
    umi5_x = np.arange(umi5_values.min(), umi5_values.max() + 1)
    umi5_hist = np.histogram(umi5_values, bins=umi5_x)[0]

    # Binomial fit (mean and n, p estimation)
    n5 = umi5_values.max()
    p5 = umi5_values.mean() / n5 if n5 > 0 else 0
    binom_y5 = binom.pmf(umi5_x, n5, p5) * len(umi5_values)

    # Threshold line
    threshold5_line = go.Scatter(
        x=[threshold5, threshold5],
        y=[0, max(umi5_hist.max(), binom_y5.max())],
        mode='lines',
        name=f"Threshold ({threshold5})",
        line=dict(color='red', dash='dash')
    )

    fig5 = go.Figure([
        go.Bar(x=umi5_x, y=umi5_hist, name="5' UMI Counts"),
        go.Scatter(x=umi5_x, y=binom_y5, mode='lines', name='Binomial Fit'),
        threshold5_line
    ])
    fig5.update_layout(title="5' UMI Count Distribution with Binomial Overlay",
                      xaxis_title="UMI Count",
                      yaxis_title="Frequency")
    pyo.plot(fig5, filename="umi5_count_distribution.html", auto_open=False)

    # Visualize 3' UMI counts
    umi3_values = np.array(list(umi3_counts.values()))
    umi3_x = np.arange(umi3_values.min(), umi3_values.max() + 1)
    umi3_hist = np.histogram(umi3_values, bins=umi3_x)[0]

    n3 = umi3_values.max()
    p3 = umi3_values.mean() / n3 if n3 > 0 else 0
    binom_y3 = binom.pmf(umi3_x, n3, p3) * len(umi3_values)

    threshold3_line = go.Scatter(
        x=[threshold3, threshold3],
        y=[0, max(umi3_hist.max(), binom_y3.max())],
        mode='lines',
        name=f"Threshold ({threshold3})",
        line=dict(color='red', dash='dash')
    )

    fig3 = go.Figure([
        go.Bar(x=umi3_x, y=umi3_hist, name="3' UMI Counts"),
        go.Scatter(x=umi3_x, y=binom_y3, mode='lines', name='Binomial Fit'),
        threshold3_line
    ])
    fig3.update_layout(title="3' UMI Count Distribution with Binomial Overlay",
                      xaxis_title="UMI Count",
                      yaxis_title="Frequency")
    pyo.plot(fig3, filename="umi3_count_distribution.html", auto_open=False)
    with open(output_csv, "w") as f:
        f.write("5pUMI,3pUMI\n")
        for read_name, umi_pair in filtered_UMI_dict.items():
            f.write(f"{umi_pair[0]},{umi_pair[1]}\n")
    print(f"Filtered UMIs written to {output_csv}")
""" 
    plt.hist(list(umi5_counts.values()), bins=30, alpha=0.7, label="5' UMI Counts")
    plt.axvline(threshold5, color='red', linestyle='dashed', linewidth=2, label='Threshold')
    mean_5 = umi5_values.mean()
    std_5 = umi5_values.std()
    plt.axvline(mean_5, color='green', linestyle='solid', linewidth=2, label='Mean')
    plt.axvline(mean_5 - std_5, color='purple', linestyle='dotted', linewidth=2, label='Mean - 1 Std')
    plt.axvline(mean_5 + std_5, color='orange', linestyle='dotted', linewidth=2, label='Mean + 1 Std')
    plt.xlabel("Count")
    plt.ylabel("Frequency")
    plt.title("5' UMI Count Distribution")
    plt.legend()

    plt.savefig("/Users/nomur/Desktop/test-data/output/LP250/html")
    plt.close() """

    return filtered_umi5_list, filtered_umi3_list, filtered_UMI_dict