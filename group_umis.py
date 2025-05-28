import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from utils import log
from utils import calculator
from umi_tools import UMIClusterer
from collections import Counter



def group_and_extract(input_bam, reference, barcodeMotif, doubleBarcode=True, output_dir=None):
    logger = log()
    calc = calculator()
    barcodeLength = len(barcodeMotif)
    umi_pair_dict = {}
    
    #Find the index of the motif5 and motif3 in the refernece sequence
    motif5_index = None
    motif3_index = None
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    logger.log_message(f"Reference sequence loaded from {reference}")
    logger.log_message(f"Reference sequence length: {len(ref_seq)}")

    motif5_index = ref_seq.find(barcodeMotif)
    if motif5_index == -1:
        logger.log_message(f"First barcode region not found in reference sequence.")
        return
    motif3_index = ref_seq.rfind(barcodeMotif)
    if motif3_index == -1:
        logger.log_message(f"Second barcode region not found in reference sequence.")
        return
    logger.log_message(f"Motif5 found at index {motif5_index}, Motif3 found at index {motif3_index} in reference sequence.")

    # Open the BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    logger.log_message(f"BAM file opened: {input_bam}")
    #Extract 5' UMI and 3' UMI from the reads
    raw_umi5_list = []
    raw_umi3_list = []

    # Prepare to write to a new BAM file with updated query names
    output_bam_path = f"{output_dir}/grouped_reads.bam"
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bamfile)

    umi5_start = motif5_index + 1  # 1-based reference coordinate
    umi5_end = motif5_index + barcodeLength  # exclusive
    umi3_start = motif3_index +1  # 1-based reference coordinate
    umi3_end = motif3_index + barcodeLength  # exclusive

    umi5_extracted_count = 0
    umi3_extracted_count = 0
    both_umis_extracted_count = 0

    for read in bamfile.fetch(start=umi5_start - 1, end=umi3_end):
        # Map reference positions to query positions
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
        umi5_found = False
        if query_start_5p is not None and query_end_5p is not None:
            seq5p = read.query_sequence[query_start_5p:query_end_5p]
            raw_umi5_list.append(seq5p)
            read.query_name = f"{read.query_name}_5p_umi:{seq5p}"
            umi5_extracted_count += 1
            umi5_found = True

        # Find query indices corresponding to the 3p UMI region
        query_start_3p = None
        query_end_3p = None
        for i, ref_pos in enumerate(ref_positions):
            if ref_pos == umi3_start - 1 and query_start_3p is None:
                query_start_3p = i
            if ref_pos == umi3_end - 1:
                query_end_3p = i + 1  # end is exclusive
                break
        umi3_found = False
        if query_start_3p is not None and query_end_3p is not None:
            seq3p = read.query_sequence[query_start_3p:query_end_3p]
            raw_umi3_list.append(seq3p)
            read.query_name = f"{read.query_name}_3p_umi:{seq3p}"
            umi3_extracted_count += 1
            umi3_found = True

        if umi5_found and umi3_found:
            both_umis_extracted_count += 1
            if seq5p not in umi_pair_dict:
                umi_pair_dict[seq5p] = []
            umi_pair_dict[seq5p].append(seq3p)

        #Update the query name with the UMI
        output_bam.write(read)

    logger.log_message(f"5' UMI extracted from {umi5_extracted_count} reads.")
    logger.log_message(f"3' UMI extracted from {umi3_extracted_count} reads.")
    logger.log_message(f"Both UMIs extracted from {both_umis_extracted_count} reads.")
        
    output_bam.write(read)
    umi_5p_correction_dict = umi_correction(raw_umi5_list)
    umi_3p_correction_dict = umi_correction(raw_umi3_list)
    
    corrected_umi_pair_dict = {}
    for umi_5p, umi_3p_array in umi_pair_dict.items():
        # Correct the 5' UMI
        corrected_umi_5p = next(
            (corrected for corrected, uncorrected_list in umi_5p_correction_dict.items() if umi_5p in uncorrected_list),
            umi_5p  # Default to the original UMI if no correction is found
        )

        # Correct the 3' UMIs
        corrected_umi_3p_list = [
            next(
                (corrected for corrected, uncorrected_list in umi_3p_correction_dict.items() if umi_3p in uncorrected_list),
                umi_3p  # Default to the original UMI if no correction is found
            )
            for umi_3p in umi_3p_array
        ]

        # Add the corrected 5' UMI and its associated corrected 3' UMIs to the new dictionary
        if corrected_umi_5p not in corrected_umi_pair_dict:
            corrected_umi_pair_dict[corrected_umi_5p] = []
        corrected_umi_pair_dict[corrected_umi_5p].extend(corrected_umi_3p_list)

        # Update umi_5p_correction_dict by clustering 5' UMIs that share the same 3' UMI in corrected_umi_pair_dict
        # For each unique 3' UMI, find all 5' UMIs associated with it, and cluster them if more than one exists

    # Export UMI data as CSV
    csv_filename = "output_files/umi_extraction.csv"
    import csv
    
    with open(csv_filename, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['umi_5p', 'corrected_umi_5p', 'umi_3p', 'corrected_umi_3p'])
        for umi5, umi3 in zip( raw_umi5_list, raw_umi3_list):
            # Find corrected UMIs
            corrected_umi5 = next(
                (corrected for corrected, uncorrected_list in umi_5p_correction_dict.items() if umi5 in uncorrected_list),
                umi5
            )
            corrected_umi3 = next(
                (corrected for corrected, uncorrected_list in umi_3p_correction_dict.items() if umi3 in uncorrected_list),
                umi3
            )
            csv_writer.writerow([ umi5, corrected_umi5 ,umi3 ,corrected_umi3])
    print(f"UMI extraction data exported to {csv_filename}")

    corrected_umi_3p = correctUMIS(raw_umi3_list, umi_3p_correction_dict)
    corrected_umi_5p = correctUMIS(raw_umi5_list, umi_5p_correction_dict)
    corrected_umi_pair_dict = correctUMI_pair_dict(corrected_umi_pair_dict, umi_5p_correction_dict, umi_3p_correction_dict)

    calc.hamming_distance_matrix(list(umi_5p_correction_dict.keys()), "output_files/5p_umi_hamming_distances2.png")
    calc.hamming_distance_matrix(list(umi_3p_correction_dict.keys()), "output_files/3p_umi_hamming_distances2.png")
    bunch_of_Stats(raw_umi5_list, raw_umi3_list)
    calc.plot_umi_frequencies(corrected_umi_5p, "output_files/5p_umi_frequencies.png")
    calc.plot_umi_frequencies(corrected_umi_3p, "output_files/3p_umi_frequencies.png")
    calc.concordance(corrected_umi_pair_dict, "output_files/umi_concordance_plot.png")
   
   
    bamfile.close()
    output_bam.close()

    logger.log_message(f"Written updated BAM file with new query names to {output_bam_path}")


def umi_correction(umi_list):
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

def bunch_of_Stats(umi_5p_list,umi_3p_list):
    # Calculate percentages for 3' UMIs
    num_start_with_a_3p = sum(1 for umi in umi_3p_list if umi.startswith('A'))
    num_start_with_g_3p = sum(1 for umi in umi_3p_list if umi.startswith('G'))
    num_start_with_t_3p = sum(1 for umi in umi_3p_list if umi.startswith('T'))
    num_start_with_c_3p = sum(1 for umi in umi_3p_list if umi.startswith('C'))

    percentage_start_with_a_3p = (num_start_with_a_3p / len(umi_3p_list)) * 100 if umi_3p_list else 0
    percentage_start_with_g_3p = (num_start_with_g_3p / len(umi_3p_list)) * 100 if umi_3p_list else 0
    percentage_start_with_t_3p = (num_start_with_t_3p / len(umi_3p_list)) * 100 if umi_3p_list else 0
    percentage_start_with_c_3p = (num_start_with_c_3p / len(umi_3p_list)) * 100 if umi_3p_list else 0

    average_length_3p = sum(len(umi) for umi in umi_3p_list) / len(umi_3p_list) if umi_3p_list else 0

    # Calculate percentages for 5' UMIs
    num_start_with_a_5p = sum(1 for umi in umi_5p_list if umi.startswith('A'))
    num_start_with_g_5p = sum(1 for umi in umi_5p_list if umi.startswith('G'))
    num_start_with_t_5p = sum(1 for umi in umi_5p_list if umi.startswith('T'))
    num_start_with_c_5p = sum(1 for umi in umi_5p_list if umi.startswith('C'))

    percentage_start_with_a_5p = (num_start_with_a_5p / len(umi_5p_list)) * 100 if umi_5p_list else 0
    percentage_start_with_g_5p = (num_start_with_g_5p / len(umi_5p_list)) * 100 if umi_5p_list else 0
    percentage_start_with_t_5p = (num_start_with_t_5p / len(umi_5p_list)) * 100 if umi_5p_list else 0
    percentage_start_with_c_5p = (num_start_with_c_5p / len(umi_5p_list)) * 100 if umi_5p_list else 0

    average_length_5p = sum(len(umi) for umi in umi_5p_list) / len(umi_5p_list) if umi_5p_list else 0

    # Print the results
    print("\n5' UMI Statistics:")
    print(f"Percentage of 5' UMIs starting with 'A': {percentage_start_with_a_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'G': {percentage_start_with_g_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'T': {percentage_start_with_t_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'C': {percentage_start_with_c_5p:.2f}%")
    print(f"Average length of 5' UMIs: {average_length_5p:.2f}")

    print("3' UMI Statistics:")
    print(f"Percentage of 3' UMIs starting with 'A': {percentage_start_with_a_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'G': {percentage_start_with_g_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'T': {percentage_start_with_t_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'C': {percentage_start_with_c_3p:.2f}%")
    print(f"Average length of 3' UMIs: {average_length_3p:.2f}")

def correctUMIS(umi_list, correction_dict):

    corrected_umi_list = [
        next(
            (corrected for corrected, uncorrected_list in correction_dict.items() if umi in uncorrected_list),
            umi  # Default to the original UMI if no correction is found
        )
        for umi in umi_list
    ]
    return corrected_umi_list

def correctUMI_pair_dict(dict_umi_pair, umi_5p_correction_dict, umi_3p_correction_dict):
    corrected_umi_pair_dict = {}
    for umi_5p, umi_3p_array in dict_umi_pair.items():
        # Correct the 5' UMI
        corrected_umi_5p = next(
            (corrected for corrected, uncorrected_list in umi_5p_correction_dict.items() if umi_5p in uncorrected_list),
            umi_5p  # Default to the original UMI if no correction is found
        )

        # Correct the 3' UMIs
        corrected_umi_3p_list = [
            next(
                (corrected for corrected, uncorrected_list in umi_3p_correction_dict.items() if umi_3p in uncorrected_list),
                umi_3p  # Default to the original UMI if no correction is found
            )
            for umi_3p in umi_3p_array
        ]

        # Add the corrected 5' UMI and its associated corrected 3' UMIs to the new dictionary
        if corrected_umi_5p not in corrected_umi_pair_dict:
            corrected_umi_pair_dict[corrected_umi_5p] = []
        corrected_umi_pair_dict[corrected_umi_5p].extend(corrected_umi_3p_list)

    return corrected_umi_pair_dict

if __name__ == "__main__":
    input_bam = "/Users/nartsavomur/Desktop/Remote_Work/test_data_belt/test_data/NL43-Test-Aligned.bam"
    reference = "/Users/nartsavomur/Desktop/Remote_Work/test_data_belt/test_data/NL43-Test-Aligned.fasta"  # Example reference sequence
    barcodeMotif = "NNNNNNNNNNNN"
    doubleBarcode = True
    output_dir = "/Users/nartsavomur/belt/output_files"

    group_and_extract(input_bam, reference, barcodeMotif, doubleBarcode, output_dir)
