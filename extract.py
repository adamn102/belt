import argparse
from utils import log
from utils import calculator
import pysam
from Bio.Seq import Seq
from umi_tools import UMIClusterer
from Bio import SeqIO
from collections import Counter
import os


def extract_UMIs(reference, alignedReads, minReadLength ,barcodeMotif,outputDir):

    sampleName = alignedReads.split("/")[-1].split(".")[0]
    outputDir = outputDir + "/" + sampleName
    # Check if the output directory exists, if not create it
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    # Check if the html subdirectory exists, if not create it
    html_dir = os.path.join(outputDir, "html")
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)

    #Initializing the tools
    logger = log()
    calc = calculator()
    logger.log_message(f"UMI extraction of {sampleName} with barcode {barcodeMotif} motif started.")

    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(reference, alignedReads, barcodeMotif, minReadLength)

    

    logger.log_message(f"Motif5 found at index {motif5_index}, Motif3 found at index {motif3_index} in reference sequence.")
    if raw_umi5_list is None or raw_umi3_list is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return
    logger.log_message(f"UMI extraction of {sampleName} with barcode {barcodeMotif} motif finished.")
    #write csv files with the extracted UMIs
    write_csv(["5' UMI Raw","3' UMI Raw"], [raw_umi5_list,raw_umi3_list], f"{outputDir}/{sampleName}_rawUMIs.csv")

    #Clustering the UMIs and creating a corrected list
    umi5p_clusters = umi_clustering(raw_umi5_list)
    umi3p_clusters = umi_clustering(raw_umi3_list)
    logger.log_message(f"UMI clustering of {sampleName} with barcode {barcodeMotif} motif finished.")
    corrected_umi5_list = correctUMIS(raw_umi5_list, umi5p_clusters)
    corrected_umi3_list = correctUMIS(raw_umi3_list, umi3p_clusters)
    logger.log_message(f"UMI correction on the lists of {sampleName} with barcode {barcodeMotif} motif finished.")
    corrected_UMI_dict = correct_UMI_dict(raw_UMI_dict, umi5p_clusters, umi3p_clusters)
    logger.log_message(f"UMI correction on the dictionaries of {sampleName} with barcode {barcodeMotif} motif finished.")

    umi5p_stats = calculate_umi_stats(umi5p_clusters,corrected_umi5_list,raw_umi5_list)
    umi3p_stats = calculate_umi_stats(umi3p_clusters,corrected_umi3_list,raw_umi3_list)
    logger.log_message(f"UMI statistics calculation of {sampleName} with barcode {barcodeMotif} motif finished.")
    
    #write csv files with the clustered UMIs
    write_csv(["5' UMI Corrected","3' UMI Corrected"], [corrected_umi5_list,corrected_umi3_list], f"{outputDir}/{sampleName}_correctedUMIs.csv")
    write_csv(["Read Name", "5' UMI", "3' UMI"], [[read_name for read_name in corrected_UMI_dict.keys()],
                                                    [list(umi_pair.keys())[0] for umi_pair in corrected_UMI_dict.values()],
                                                    [list(umi_pair.values())[0] for umi_pair in corrected_UMI_dict.values()]], 
              f"{outputDir}/{sampleName}_correctedUMIs_read_names.csv")
    logger.log_message(f"UMI CSV files of {sampleName} with barcode {barcodeMotif} motif created.")
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi5p_clusters).keys()), sample_size=400),f"{outputDir}/html/hamming_distances_5p.html"400)
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi3p_clusters).keys()), sample_size=400),f"{outputDir}/html/hamming_distances_3p.html",400)
    calc.plot_umi_frequencies(corrected_umi5_list, f"{outputDir}/html/umi_freq_5p.html")
    calc.plot_umi_frequencies(corrected_umi3_list, f"{outputDir}/html/umi_freq_3p.html")
    calc.hamming_distance_matrix(list(umi5p_clusters.keys()), f"{outputDir}/html/umi_freq_matrix_5p.html",sample_size=400)
    calc.hamming_distance_matrix(list(umi3p_clusters.keys()), f"{outputDir}/html/umi_freq_matrix_3p.html",sample_size=400)
    logger.log_message(f"UMI plots of {sampleName} with barcode {barcodeMotif} motif created.")

    createHTML(umi5p_stats, umi3p_stats,f"{outputDir}", sampleName, pi_photo_path="alex.png")
    logger.log_message(f"HTML report created at {outputDir}/{sampleName}_Report.html")
    logger.write_log(f"{outputDir}/{sampleName}_log.txt")

def findUMI(reference, alignedReads, barcodeMotif, minReadLength):
    """
    Extracts UMI from the aligned reads based on the provided barcode motif.
    """
    # Intializing variables and finding the barcode motif in the reference sequence.
    raw_umi5_list = []
    raw_umi3_list = []
    raw_UMI_dict = {}
    barcodeLength = len(barcodeMotif)
    motif5_index = None
    motif3_index = None
    input_bam = alignedReads

    if not input_bam.endswith('.bam'):
        print(f"Input file {input_bam} is not a BAM file.")
        return
    if not reference.endswith('.fasta') and not reference.endswith('.fa'):
        print(f"Reference file {reference} is not a FASTA file.")
        return
    
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    motif5_index = ref_seq.find(barcodeMotif)
    if motif5_index == -1:
        print(f"First barcode region not found in reference sequence.")
        return
    motif3_index = ref_seq.rfind(barcodeMotif)
    if motif3_index == -1:
        print(f"Second barcode region not found in reference sequence.")
        return

    #The big for loop that goes through the aligned reads and extracts the UMI.
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    umi5_start = motif5_index + 1  # 1-based reference coordinate
    umi5_end = motif5_index + barcodeLength  # exclusive
    umi3_start = motif3_index +1  # 1-based reference coordinate
    umi3_end = motif3_index + barcodeLength

    for read in bamfile.fetch(start=umi5_start - 1, end=umi3_end):# This should popuulate the raw UMI lists with the 5p and 3p UMI sequences.
        if read.query_sequence is None or len(read.query_sequence)< minReadLength:
            continue  # Skip reads without sequence
        ref_positions = read.get_reference_positions(full_length=True)
        # Find query indices corresponding to the 5p UMI region
        query_start_5p = None
        query_end_5p = None
        for i, ref_pos in enumerate(ref_positions):
            if ref_pos == umi5_start - 1 and query_start_5p is None:
                query_start_5p = i
            if ref_pos == umi5_end - 1:
                query_end_5p = i + 1  # end is exclusive
                break
    
        if query_start_5p is not None and query_end_5p is not None:
            seq5p = read.query_sequence[query_start_5p:query_end_5p]
            if len(seq5p)== barcodeLength:
                raw_umi5_list.append(seq5p)
                raw_UMI_dict[read.query_name] = [seq5p,""]

            

        # Find query indices corresponding to the 3p UMI region
        query_start_3p = None
        query_end_3p = None
        for i, ref_pos in enumerate(ref_positions):
            if ref_pos == umi3_start - 1 and query_start_3p is None:
                query_start_3p = i
            if ref_pos == umi3_end - 1:
                query_end_3p = i + 1  # end is exclusive
                break
        if query_start_3p is not None and query_end_3p is not None:
            seq3p = read.query_sequence[query_start_3p:query_end_3p]
            if len(seq3p)== barcodeLength:
                raw_umi3_list.append(seq3p)
                if read.query_name in raw_UMI_dict:
                    raw_UMI_dict[read.query_name] = [raw_UMI_dict[read.query_name][0], seq3p]
                else:
                    raw_UMI_dict[read.query_name] = [None, seq3p]

    bamfile.close()
    return raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict      

def write_csv(headers, data, output_file):
    """
    Writes headers and data to a CSV file.
    Handles columns of unequal length by filling missing values with empty strings.

    :param headers: List of column headers.
    :param data: List of lists, where each sublist is a column.
    :param output_file: Path to the output CSV file.
    """
    import csv
    # Transpose columns to rows, filling missing values with ''
    max_len = max(len(col) for col in data)
    rows = []
    for i in range(max_len):
        row = []
        for col in data:
            row.append(col[i] if i < len(col) else '')
        rows.append(row)
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(rows)

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

def umi_clustering(umi_list):
    umi_clusterer = UMIClusterer(cluster_method="directional")
    umi_counter = Counter(umi_list)
    values = list(umi_counter.values())
    keys = list(umi_counter.keys())
    bit_keys = [umi.encode("utf-8") for umi in keys]
    umi_count_dict = dict(zip(bit_keys, values))
    clustered_umis = umi_clusterer(umi_count_dict, 1)

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

def calculate_umi_stats(dict, correctedList ,rawList):
    """
    Calculate statistics for the UMIs.
    :param correctedList: List of corrected UMIs
    :param rawList: List of raw UMIs
    :return: List of statistics [raw_count, clusters, chao1, chao1_updated, entropy]
    """
    calc = calculator()

    raw_count = len(rawList)
    clusters = len(list(dict.keys()))
    chao1 = calc.chao1(correctedList)
    chao1_updated = calc.balazs_chao1(correctedList)
    shannon_entropy = calc.shannon_entropy(correctedList)

    return raw_count, clusters, chao1, chao1_updated, shannon_entropy

