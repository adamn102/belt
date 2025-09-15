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


def extract_UMIs(reference, alignedReads, minReadLength, barcodeMotifs, outputDir):
    sampleName: object = alignedReads.split("/")[-1].split(".")[0]

    # Check if the output directory exists, if not create it
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Initializing ttools
    logger = log()
    logger.log_message(f"UMI extraction of {sampleName} with {len(barcodeMotifs)} barcode motifs started.")
    motif_indices, raw_UMI_dict = findUMI(reference, alignedReads, barcodeMotifs, minReadLength)
                                                                       
    if motif_indices is None or raw_UMI_dict is None:
        logger.log_message("UMI extraction failed. Please check the input files and parameters.")
        return
    
    logger.log_message(f"Found motifs at indices: {motif_indices}")
    logger.log_message(f"UMI extraction of {sampleName} with {len(barcodeMotifs)} barcode motifs finished.")
    # Write csv files with the extracted UMIs
    # write_csv(["5' UMI Raw", "3' UMI Raw"], [raw_umi5_list, raw_umi3_list], f"{outputDir}/{sampleName}_rawUMIs.csv")

    # write raw_UMI_dict as a .csv with the read names as the first column,
    # and UMI columns for each motif found
    # do this using pandas
    umi_columns = [f"UMI{i+1}" for i in range(len(barcodeMotifs))]
    df = pd.DataFrame.from_dict(raw_UMI_dict, orient='index', columns=umi_columns)
    df.index.name = 'Read Name'
    df.reset_index(inplace=True)
    df.to_csv(f"{outputDir}/{sampleName}_rawUMIs_read_names.csv", index=False)


def findUMI(reference, alignedReads, barcodeMotifs, minReadLength):
    """
    Extracts UMI from the aligned reads based on the provided barcode motifs.
    Finds all motifs in one direction and extracts UMIs from all positions.
    """
    # Convert single motif to list for backward compatibility
    if isinstance(barcodeMotifs, str):
        barcodeMotifs = [barcodeMotifs]
    
    # Initializing variables
    raw_UMI_dict = {}
    input_bam = alignedReads

    if not input_bam.endswith('.bam'):
        print(f"Input file {input_bam} is not a BAM file.")
        return None, None
    if not reference.endswith('.fasta') and not reference.endswith('.fa'):
        print(f"Reference file {reference} is not a FASTA file.")
        return None, None

    # Read reference sequence and find all motif positions
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    
    motif_indices = []
    motif_lengths = []
    
    # Find all motif positions (going in one direction only)
    for motif in barcodeMotifs:
        motif_length = len(motif)
        motif_lengths.append(motif_length)
        
        # Find all occurrences of this motif
        motif_positions = []
        start = 0
        while True:
            pos = ref_seq.find(motif, start)
            if pos == -1:
                break
            motif_positions.append(pos)
            start = pos + 1
        
        if not motif_positions:
            print(f"Motif '{motif}' not found in reference sequence.")
            return None, None
        
        motif_indices.append(motif_positions)
    
    # Get the range to fetch reads (from first motif to last motif)
    all_positions = []
    for positions in motif_indices:
        all_positions.extend(positions)
    
    if not all_positions:
        print("No motif positions found.")
        return None, None
    
    min_pos = min(all_positions)
    max_pos = max(all_positions) + max(motif_lengths)
    
    # Extract UMIs from all motif positions
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch(start=min_pos, end=max_pos):
        if read.query_sequence is None or len(read.query_sequence) < minReadLength:
            continue  # Skip reads without sequence
        
        ref_positions = read.get_reference_positions(full_length=True)
        read_umis = [None] * len(barcodeMotifs)  # Initialize with None for each motif
        
        # Extract UMI for each motif at each position
        for motif_idx, (motif, positions, motif_length) in enumerate(zip(barcodeMotifs, motif_indices, motif_lengths)):
            for pos in positions:
                umi_start = pos + 1  # 1-based reference coordinate
                umi_end = pos + motif_length  # exclusive
                
                # Find query indices corresponding to this UMI region
                query_start = None
                query_end = None
                
                for i, ref_pos in enumerate(ref_positions):
                    if ref_pos == umi_start - 1 and query_start is None:
                        query_start = i
                    if ref_pos == umi_end - 1:
                        query_end = i + 1  # end is exclusive
                        break
                
                if query_start is not None and query_end is not None:
                    seq = read.query_sequence[query_start:query_end]
                    if len(seq) == motif_length:
                        read_umis[motif_idx] = seq
                        break  # Use first valid UMI for this motif
        
        # Only add to dictionary if at least one UMI was found
        if any(umi is not None for umi in read_umis):
            raw_UMI_dict[read.query_name] = read_umis

    bamfile.close()
    
    # Flatten motif_indices for return (just the positions, not grouped by motif)
    flat_indices = []
    for positions in motif_indices:
        flat_indices.extend(positions)
    
    return flat_indices, raw_UMI_dict


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
