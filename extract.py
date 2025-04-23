import click
from Bio.Seq import Seq
from Bio import SeqIO
from difflib import SequenceMatcher
from utils import log
from umi_tools import UMIClusterer
from collections import Counter

@click.command()
@click.argument('input_fastq')
@click.option('--output', required=True, help="Output FASTQ file")
@click.option('--min_read_length', default=3900, help="Minimum expected read length")
@click.option('--motif5', required=True, help="5' motif sequence")
@click.option('--motif3', required=True, help="3' motif sequence")
@click.option('--fuzzy', default=2, help="Number of fuzzy mismatches allowed")
@click.option('--window', default=100, help="Motif search window size")
def cli(input_fastq, min_read_length, output, motif5, motif3, fuzzy, window):
    """Extract UMIs from reads based on 5' and 3' motifs."""
    extract_reads(input_fastq, min_read_length, output, motif5, motif3, fuzzy, window)


def extract_reads(input_fastq, min_read_length, output_fastq, motif_5p, motif_3p, fuzzy, window):
    log(f"Extracting UMIs from {input_fastq}")
    written = 0
    skipped = 0


    # create UMI lists
    read = []
    read_name = []
    umi_5p_list = []
    umi_3p_list = []

    with open(input_fastq, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, "fastq"):#looping through the reads the first time, check size first extract UMI is size is good.
            if not lengthCheck(record.seq, min_read_length):
                continue
            seq = str(record.seq)
            rc_seq = str(record.seq.reverse_complement())

            # Check for motifs in the forward sequence
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

            # Skip reads without both UMIs
            if not umi_5p or not umi_3p:
                skipped += 1
                continue

            # only keep reads with both UMIs
            if umi_5p and umi_3p:
                read.append(record)
                read_name.append(record.id)
                umi_5p_list.append(umi_5p)
                umi_3p_list.append(umi_3p)

        # temporarily store data as dictionary
        read_dict = dict(zip(read_name, zip(read, umi_5p_list, umi_3p_list)))
        #print(read_dict)

        # count number of each 5p_umi
        umi_5p_correction_dict = umi_correction(umi_5p_list)
        umi_3p_correction_dict = umi_correction(umi_3p_list)

        print(umi_5p_correction_dict)
        print(umi_3p_correction_dict)

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

    for i in range(len(read_seq) - len(motif) + 1):
        substring = read_seq[i:i + len(motif)]
        mismatches = sum(1 for a, b in zip(substring, motif) if a != b)

        if mismatches <= max_mismatches and mismatches < best_mismatches:
            best_match = substring
            best_mismatches = mismatches

    if best_match is not None:
        return True, best_match
    return False, None

def lengthCheck(read, min_length):# boolean operator to check is the individual read is equal or greater than the minimum read length
    if min_length >= len(str(read)):
        True
    else:
        False
