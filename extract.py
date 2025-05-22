from Bio.Seq import Seq
from Bio import SeqIO
from difflib import SequenceMatcher
from utils import log, calculator
from umi_tools import UMIClusterer
from collections import Counter
import math
import os
from collections import Counter
import Levenshtein

def extract_reads(input_fastq, min_read_length, output_fastq, motif_5p, motif_3p, fuzzy, window):
    logger = log()
    logger.log_message(f"Extracting reads from {input_fastq}")
    logger.log_message(f"Outputting to {output_fastq}")
    written = 0
    skipped = 0 # This is the number of reads that are skipped 1) if the read is shorter than the minimum read length, 2) if the read does not contain both UMIs
    no_umi5p = 0
    no_umi3p = 0
    no_umi = 0

    # Advanced UMI extraction, and corrections metrics 
    #TODO: add these counters to the correct places
    umi5p_num_before = 0
    umi3p_num_before = 0
    umi5p_num_after = 0
    umi3p_num_after = 0

    # create UMI lists
    read = []
    read_name = []
    umi_5p_list = [] # keeps track of the 5' UMI each instance seperately
    umi_3p_list = []# keeps track of the 3' UMI each instance seperately
    umi_pair_dict = {}  # Dictionary to store 5' UMI as key and associated 3' UMIs as values



    with open(input_fastq, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            seq = str(record.seq)
            rc_seq = str(record.seq.reverse_complement())
            if len(seq) <= min_read_length:
                skipped += 1
            else:
                # Check for motifs in the forward sequence
                #print(i, is_fuzzy_match(seq[:window], motif_5p, fuzzy),seq )
                match_5p, matched_5p = is_fuzzy_match(seq[:window], motif_5p, fuzzy)
                match_3p, matched_3p = is_fuzzy_match(seq[-window:], motif_3p, fuzzy)

                if match_5p and match_3p:
                    umi_5p_start = seq.find(matched_5p) + len(matched_5p)
                    umi_3p_end = seq.rfind(matched_3p)
                    umi_5p = seq[umi_5p_start:umi_5p_start + 18] if umi_3p_end - umi_5p_start >= 18 else None
                    umi_3p = seq[umi_3p_end - 18:umi_3p_end] if umi_3p_end - umi_5p_start >= 18 else None
                    new_record = record
                else:
                    # Check for motifs in the reverse complement sequence
                    match_5p, matched_5p = is_fuzzy_match(rc_seq[:window], motif_5p, fuzzy)
                    match_3p, matched_3p = is_fuzzy_match(rc_seq[-window:], motif_3p, fuzzy)

                    if match_5p and match_3p:
                        record.seq = Seq(rc_seq)
                        umi_5p_start = rc_seq.find(matched_5p) + len(matched_5p)
                        umi_3p_end = rc_seq.rfind(matched_3p)
                        umi_5p = rc_seq[umi_5p_start:umi_5p_start + 32] if umi_3p_end - umi_5p_start >= 32 else None
                        umi_3p = rc_seq[umi_3p_end - 32:umi_3p_end] if umi_3p_end - umi_5p_start >= 32 else None
                        new_record = record
                    else:
                        skipped += 1
                        continue

            # Skip reads without proper UMIs
            if not (umi_5p and is_fuzzy_match(umi_5p[:10],"ACCGTCGGAT",3) and is_fuzzy_match(umi_5p[-10:],"CTACGCAGGT",3)):
                no_umi5p += 1
                continue
            elif not (umi_3p and is_fuzzy_match(umi_3p[:10],"CAATGCTCGA",3) and is_fuzzy_match(umi_3p[-10:],"TGTTGGCCTC",3)):
                no_umi3p += 1
                continue
            elif not umi_5p and not umi_3p:
                no_umi += 1
                continue
            else:
                umi_5p = umi_5p[10:-10]
                umi_3p = umi_3p[10:-10]

                # only keep reads with both UMIs
                if umi_5p and umi_3p:
                    read.append(record)
                    read_name.append(record.id)
                    umi_5p_list.append(umi_5p)
                    umi_3p_list.append(umi_3p)
                    # Track UMI pairs
                    if umi_5p not in umi_pair_dict:
                        umi_pair_dict[umi_5p] = []
                    umi_pair_dict[umi_5p].append(umi_3p)

        # count number of each 5p_umi
        umi_5p_correction_dict = strict_direct_clustering(umi_5p_list)
        umi_3p_correction_dict = strict_direct_clustering(umi_3p_list)

    # Export UMI data as CSV
    csv_filename = "output_files/umi_extraction.csv"
    import csv
    with open(csv_filename, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['read_name', 'umi_5p', 'umi_3p', 'corrected_umi_5p', 'corrected_umi_3p'])
        for name, umi5, umi3 in zip(read_name, umi_5p_list, umi_3p_list):
            # Find corrected UMIs
            corrected_umi5 = next(
                (corrected for corrected, uncorrected_list in umi_5p_correction_dict.items() if umi5 in uncorrected_list),
                umi5
            )
            corrected_umi3 = next(
                (corrected for corrected, uncorrected_list in umi_3p_correction_dict.items() if umi3 in uncorrected_list),
                umi3
            )
            csv_writer.writerow([name, umi5, umi3, corrected_umi5, corrected_umi3])
    print(f"UMI extraction data exported to {csv_filename}")
        
    with open(output_fastq, 'w') as out_handle:
        for i in range(len(read_name)):
            # Map uncorrected UMIs to their corrected counterparts
            corrected_umi_5p = next(
                (corrected for corrected, uncorrected_list in umi_5p_correction_dict.items() if umi_5p_list[i] in uncorrected_list),
                None
            )
            corrected_umi_3p = next(
                (corrected for corrected, uncorrected_list in umi_3p_correction_dict.items() if umi_3p_list[i] in uncorrected_list),
                None
            )

            # Handle cases where no corrected UMI is found
            if corrected_umi_5p is None or corrected_umi_3p is None:
                skipped += 1
                continue

            # Append corrected UMIs to the read name
            new_record = read[i]
            new_record.id += f":UMI5-{corrected_umi_5p}-UMI3-{corrected_umi_3p}"
            new_record.description = new_record.id

            # Write the read to the general output file
            SeqIO.write(new_record, out_handle, "fastq")
            written += 1

    



    logger.extraction_log(input_fastq, output_fastq, written, skipped, no_umi5p, no_umi3p, no_umi)
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
    print("3' UMI Statistics:")
    print(f"Percentage of 3' UMIs starting with 'A': {percentage_start_with_a_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'G': {percentage_start_with_g_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'T': {percentage_start_with_t_3p:.2f}%")
    print(f"Percentage of 3' UMIs starting with 'C': {percentage_start_with_c_3p:.2f}%")
    print(f"Average length of 3' UMIs: {average_length_3p:.2f}")

    print("\n5' UMI Statistics:")
    print(f"Percentage of 5' UMIs starting with 'A': {percentage_start_with_a_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'G': {percentage_start_with_g_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'T': {percentage_start_with_t_5p:.2f}%")
    print(f"Percentage of 5' UMIs starting with 'C': {percentage_start_with_c_5p:.2f}%")
    print(f"Average length of 5' UMIs: {average_length_5p:.2f}")

    calc = calculator()
    clustered_umi_3p = correctUMIS(umi_3p_list, umi_3p_correction_dict)
    clustered_umi_5p = correctUMIS(umi_5p_list, umi_5p_correction_dict)

    
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

    # Calculate hamming distances
    #print(umi_5p_list)
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi_5p_list).keys()), sample_size=30),"output_files/5p_preCpluster_hamming_distances.png")
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi_5p_correction_dict).keys()), sample_size=10),"output_files/5p_postCluster_hamming_distances.png")
    calc.plot_umi_frequencies(clustered_umi_5p, "output_files/5p_UMI_freq.png")
 
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi_3p_list).keys()), sample_size=30),"output_files/3p_preCpluster_hamming_distances.png")
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(Counter(umi_3p_correction_dict).keys()), sample_size=10),"output_files/3p_postCluster_hamming_distances.png")
    calc.plot_umi_frequencies(clustered_umi_3p, "output_files/3p_UMI_freq.png")
    
    corrected_umi_pairings = calc.plot_most_frequent_3p_percentage_and_count_per_5p(output_fastq, "output_files/concordance.png")
    separate_fastq_by_umi5("output_files/reads_after_extract.fastq", "output_files/UMI5_separated",corrected_umi_pairings)
   
    logger.write_log("output_files/log.log")
    calc.export_5p_3p_UMIs(output_fastq, "output_files/5p_3p_UMI_pairs.csv")



def umi_correction(umi_list,threshold):
    umi_clusterer = UMIClusterer(cluster_method="adjacency")
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

def is_fuzzy_match(read_seq, motif, max_mismatches):
    best_match = None
    best_mismatches = max_mismatches + 1  # Initialize with a value greater than max_mismatches
    #print(read_seq)

    for i in range(len(read_seq) - len(motif) + 1):
        substring = read_seq[i:i + len(motif)]
        #print(substring)
        mismatches = sum(1 for a, b in zip(substring, motif) if a != b)
        #print(mismatches)

        if mismatches <= max_mismatches and mismatches < best_mismatches:
            best_match = substring
            best_mismatches = mismatches
            #print(best_match)

    if best_match is not None:
        return True, best_match
    else:
        return False, None
    
def lengthCheck(read, min_length):# boolean operator to check is the individual read is equal or greater than the minimum read length
    if min_length >= len(str(read)):
        True
    else:
        False

def correctUMIS(umi_list, correction_dict):

    corrected_umi_list = [
        next(
            (corrected for corrected, uncorrected_list in correction_dict.items() if umi in uncorrected_list),
            umi  # Default to the original UMI if no correction is found
        )
        for umi in umi_list
    ]
    return corrected_umi_list

def separate_fastq_by_umi5(input_fastq, output_dir, corrected_umi_pairings):
    """
    Separates a FASTQ file into multiple FASTQ files based on the corrected 5' and 3' UMI pairings.

    :param input_fastq: Path to the input FASTQ file.
    :param output_dir: Directory where the separated FASTQ files will be saved.
    :param corrected_umi_pairings: Dictionary mapping corrected 5' UMIs to their most frequent 3' UMIs.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    umi5_files = {}

    # Parse the input FASTQ file
    for record in SeqIO.parse(input_fastq, "fastq"):
        # Extract the corrected 5' and 3' UMIs from the read name
        umi5 = None
        umi3 = None
        if "UMI5-" in record.id and "-UMI3-" in record.id:
            umi5 = record.id.split("UMI5-")[1].split("-UMI3")[0]
            umi3 = record.id.split("-UMI3-")[1]

        # Check if the 5' UMI exists in the corrected pairings and matches the most frequent 3' UMI
        if umi5 in corrected_umi_pairings and corrected_umi_pairings[umi5] == umi3:
            # Create a new file for this UMI5 if it doesn't exist
            if umi5 not in umi5_files:
                umi5_filename = os.path.join(output_dir, f"UMI5_{umi5}_3p-{corrected_umi_pairings[umi5]}.fastq")
                umi5_files[umi5] = open(umi5_filename, 'w')

            # Add the most frequent 3' UMI to the read name
            record.id += f":MostFrequent3p-{corrected_umi_pairings[umi5]}"
            record.description = record.id

            # Write the record to the corresponding UMI5 file
            SeqIO.write(record, umi5_files[umi5], "fastq")

    # Close all open files
    for file in umi5_files.values():
        file.close()

def recursive_umi_clustering(umi_list, first_threshold=1, subcluster_threshold=1, max_cluster_size=10):
    umi_clusterer = UMIClusterer(cluster_method="adjacency")
    umi_counter = Counter(umi_list)
    keys = list(umi_counter.keys())
    bit_keys = [umi.encode("utf-8") for umi in keys]
    umi_count_dict = dict(zip(bit_keys, umi_counter.values()))
    
    clustered_umis = umi_clusterer(umi_count_dict, first_threshold)

    final_clusters = []

    for cluster in clustered_umis:
        decoded_cluster = [umi.decode("utf-8") for umi in cluster]
        if len(decoded_cluster) > max_cluster_size:
            # Subcluster large clusters using edit distance directly
            subclusters = subcluster_by_edit_distance(decoded_cluster, threshold=subcluster_threshold)
            final_clusters.extend(subclusters)
        else:
            final_clusters.append(decoded_cluster)

    # Build result dict
    umi_corrected_dict = {}
    for subcluster in final_clusters:
        primary_umi = subcluster[0]
        umi_corrected_dict[primary_umi] = subcluster

    return umi_corrected_dict

def subcluster_by_edit_distance(umis, threshold=1):
    """
    Clusters a list of UMIs using Levenshtein edit distance.
    """
    clusters = []
    unassigned = set(umis)

    while unassigned:
        seed = unassigned.pop()
        current_cluster = [seed]
        to_check = set(unassigned)
        for umi in to_check:
            if Levenshtein.distance(seed, umi) <= threshold:
                current_cluster.append(umi)
                unassigned.remove(umi)
        clusters.append(current_cluster)

    return clusters


def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

def strict_direct_clustering(umi_list, max_distance=1):
    umi_counter = Counter(umi_list)
    unclustered = set(umi_counter.keys())
    cluster_dict = {}

    while unclustered:
        # Select most abundant UMI as the seed
        seed = max(unclustered, key=lambda x: umi_counter[x])
        cluster = [seed]
        unclustered.remove(seed)

        to_remove = []
        for umi in unclustered:
            if hamming_distance(seed, umi) <= max_distance:
                cluster.append(umi)
                to_remove.append(umi)
        for umi in to_remove:
            unclustered.remove(umi)

        # Store as seed: [all cluster UMIs]
        cluster_dict[seed] = cluster

    return cluster_dict