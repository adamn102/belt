from Bio.Seq import Seq
from Bio import SeqIO
from difflib import SequenceMatcher
from utils import log, calculator
from umi_tools import UMIClusterer
from collections import Counter

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
    umi5p_num_before = 0
    umi3p_num_before = 0
    umi5p_num_after = 0
    umi3p_num_after = 0

    # create UMI lists
    read = []
    read_name = []
    umi_5p_list = []
    umi_3p_list = []

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
                    umi_5p = seq[umi_5p_start:umi_5p_start + 12] if umi_3p_end - umi_5p_start >= 12 else None
                    umi_3p = seq[umi_3p_end - 13:umi_3p_end+1] if umi_3p_end - umi_5p_start >= 12 else None
                    new_record = record
                else:
                    # Check for motifs in the reverse complement sequence
                    match_5p, matched_5p = is_fuzzy_match(rc_seq[:window], motif_5p, fuzzy)
                    match_3p, matched_3p = is_fuzzy_match(rc_seq[-window:], motif_3p, fuzzy)

                    if match_5p and match_3p:
                        record.seq = Seq(rc_seq)
                        umi_5p_start = rc_seq.find(matched_5p) + len(matched_5p)
                        umi_3p_end = rc_seq.rfind(matched_3p)
                        umi_5p = rc_seq[umi_5p_start:umi_5p_start + 12] if umi_3p_end - umi_5p_start >= 12 else None
                        umi_3p = rc_seq[umi_3p_end - 13:umi_3p_end+1] if umi_3p_end - umi_5p_start >= 12 else None
                        new_record = record
                    else:
                        skipped += 1
                        continue

                # Skip reads without UMIs, and log the item
                if not umi_5p:
                    no_umi5p += 1
                    continue
                elif not umi_3p:
                    no_umi3p += 1
                    continue
                elif not umi_5p and not umi_3p:
                    no_umi += 1
                    continue



                # only keep reads with both UMIs
                if umi_5p and umi_3p:
                    read.append(record)
                    read_name.append(record.id)
                    umi_5p_list.append(umi_5p)
                    umi_3p_list.append(umi_3p)
                    #print(f"Read {record.id} has 5' UMI: {umi_5p} and 3' UMI: {umi_3p}")

        # temporarily store data as dictionary
        read_dict = dict(zip(read_name, zip(read, umi_5p_list, umi_3p_list)))

        # count number of each 5p_umi
        umi_5p_correction_dict = umi_correction(umi_5p_list)
        umi_3p_correction_dict = umi_correction(umi_3p_list)



    with open(output_fastq, 'w') as out_handle:
        for i in range(len(read_name)):

            # Correct UMIs
            umi_5p = umi_5p_correction_dict[umi_5p_list[i]]
            umi_3p = umi_3p_correction_dict[umi_3p_list[i]]

            # Append corrected UMIs to the read name
            new_record = read[i]
            new_record.id += f":UMI5-{umi_5p}-UMI3-{umi_3p}"
            new_record.description = new_record.id

            SeqIO.write(new_record, out_handle, "fastq")
            written += 1
    logger.extraction_log(input_fastq, output_fastq, written, skipped, no_umi5p, no_umi3p, no_umi)
    calc = calculator()
    # Calculate hamming distances
    print(umi_5p_list)
    calc.plot_hamming_distances(calc.calculate_hamming_distances(umi_5p_list, sample_size=100))
    #calc.plot_hamming_distances(calc.calculate_hamming_distances(umi_3p_list, sample_size=100))
    
    print(list(umi_5p_correction_dict.keys()))
    calc.plot_hamming_distances(calc.calculate_hamming_distances(list(umi_5p_correction_dict.keys()), sample_size=100))
    #calc.plot_hamming_distances(calc.calculate_hamming_distances(umi_3p, sample_size=100))
    



def umi_correction(umi_list):
    umi_clusterer = UMIClusterer(cluster_method="adjacency")
    umi_counter = Counter(umi_list)
    values = list(umi_counter.values())
    keys = list(umi_counter.keys())
    bit_keys = [umi.encode("utf-8") for umi in keys]
    umi_count_dict = dict(zip(bit_keys, values))
    clustered_umis = umi_clusterer(umi_count_dict, threshold=3)

    umi_corrected_dict = {}
    for cluster in clustered_umis:
        primary_umi = cluster[0].decode('utf-8')
        for umi in cluster:
            umi_corrected_dict[umi.decode('utf-8')] = primary_umi
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

