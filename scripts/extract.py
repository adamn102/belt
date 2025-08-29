import sys
import os

# Add the current directory to the path so we can import local modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from utils import log
except ImportError:
    # Create a simple mock log function if utils is not available
    class MockLog:
        def log_message(self, message):
            print(f"[LOG] {message}")
    
    log = MockLog
import pysam
from Bio import SeqIO
import os
import pandas as pd


def extract_UMIs(reference, alignedReads, minReadLength, barcodeMotif, outputDir):
    sampleName: object = alignedReads.split("/")[-1].split(".")[0]

    # Check if the output directory exists, if not create it
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Initializing ttools
    logger = log()
    logger.log_message(f"UMI extraction of {sampleName} with barcode {barcodeMotif} motif started.")
    raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = findUMI(reference, alignedReads,
                                                                                     barcodeMotif, minReadLength)

    logger.log_message(
        f"Motif5 found at index {motif5_index}, Motif3 found at index {motif3_index} in reference sequence.")
    if raw_umi5_list is None or raw_umi3_list is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return
    logger.log_message(f"UMI extraction of {sampleName} with barcode {barcodeMotif} motif finished.")
    # Write csv files with the extracted UMIs
    # write_csv(["5' UMI Raw", "3' UMI Raw"], [raw_umi5_list, raw_umi3_list], f"{outputDir}/{sampleName}_rawUMIs.csv")

    # write raw_UMI_dict as a .csv with the read names as the first column,
    # 5' UMI as the second column and 3' UMI as the third column
    # do this using pandas
    df = pd.DataFrame.from_dict(raw_UMI_dict, orient='index', columns=['5pUMI', '3pUMI'])
    df.index.name = 'Read Name'
    df.reset_index(inplace=True)
    df.to_csv(f"{outputDir}/{sampleName}_rawUMIs_read_names.csv", index=False)


def findUMI(reference, alignedReads, barcodeMotif, minReadLength):
    """
    Extracts UMI from the aligned reads based on the provided barcode motif.
    """
    # Initializing variables and finding the barcode motif in the reference sequence.
    raw_umi5_list = []
    raw_umi3_list = []
    raw_UMI_dict = {}
    barcodeLength = len(barcodeMotif)
    input_bam = alignedReads

    if not input_bam.endswith('.bam'):
        print(f"Input file {input_bam} is not a BAM file.")
        return
    if not reference.endswith('.fasta') and not reference.endswith('.fa'):
        print(f"Reference file {reference} is not a FASTA file.")
        return

    # TODO: Add flexibility to use a variable number of barcodes (e.g, 1-3 barcodes)
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

    # The big for loop that iterates through the aligned reads and extracts the UMIs.
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    umi5_start = motif5_index + 1  # 1-based reference coordinate
    umi5_end = motif5_index + barcodeLength  # exclusive
    umi3_start = motif3_index + 1  # 1-based reference coordinate
    umi3_end = motif3_index + barcodeLength

    # This should populate the raw UMI lists with the 5p and 3p UMI sequences.
    for read in bamfile.fetch(start=umi5_start - 1,
                              end=umi3_end):
        if read.query_sequence is None or len(read.query_sequence) < minReadLength:
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
            if len(seq5p) == barcodeLength:
                raw_umi5_list.append(seq5p)
                raw_UMI_dict[read.query_name] = [seq5p, ""]

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
            if len(seq3p) == barcodeLength:
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
