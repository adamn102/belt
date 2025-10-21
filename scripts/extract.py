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

def extract_umis(reference, alignedReads, minReadLength, barcodeMotifs, outputPath):
    sampleName: object = alignedReads.split("/")[-1].split(".")[0]

    # Check if outputPath is a directory or a file
    if outputPath.endswith('.csv'):
        # Direct file path - create parent directory if needed
        output_dir = os.path.dirname(outputPath)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        csv_file = outputPath
    else:
        # Directory path - create directory and use default filename
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
        csv_file = os.path.join(outputPath, f"{sampleName}.belt_extract.csv")

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
    df.to_csv(csv_file, index=False)


def findUMI_legacy(reference, alignedReads, barcodeMotif, minReadLength):
    """
    Legacy version of findUMI that extracts 5' and 3' UMIs separately.
    Finds first occurrence (5' motif) and last occurrence (3' motif) in reference sequence.
    Returns: (raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict)
    """
    # Initializing variables
    raw_umi5_list = []
    raw_umi3_list = []
    raw_UMI_dict = {}
    input_bam = alignedReads

    if not input_bam.endswith('.bam'):
        print(f"Input file {input_bam} is not a BAM file.")
        return None, None, None, None, None
    if not reference.endswith('.fasta') and not reference.endswith('.fa'):
        print(f"Reference file {reference} is not a FASTA file.")
        return None, None, None, None, None

    # Read reference sequence and find motif positions
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    
    # Find 5' motif (first occurrence) and 3' motif (last occurrence)
    motif5_index = ref_seq.find(barcodeMotif)
    motif3_index = ref_seq.rfind(barcodeMotif)
    
    if motif5_index == -1:
        print(f"Motif '{barcodeMotif}' not found in reference sequence.")
        return None, None, None, None, None
    
    motif_length = len(barcodeMotif)
    
    # Calculate UMI regions
    umi5_start = motif5_index + 1  # Skip the first base (motif marker)
    umi5_end = motif5_index + motif_length
    umi3_start = motif3_index + 1
    umi3_end = motif3_index + motif_length
    
    # Extract UMIs from both regions
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch(start=min(umi5_start, umi3_start) - 1, end=max(umi5_end, umi3_end)):
        if read.query_sequence is None or len(read.query_sequence) < minReadLength:
            continue
        
        ref_positions = read.get_reference_positions(full_length=True)
        
        # Extract 5' UMI
        umi5 = None
        query_start = None
        query_end = None
        
        for i, ref_pos in enumerate(ref_positions):
            if ref_pos == umi5_start - 1 and query_start is None:
                query_start = i
            if ref_pos == umi5_end - 1:
                query_end = i + 1
                break
        
        if query_start is not None and query_end is not None:
            seq = read.query_sequence[query_start:query_end]
            if len(seq) == motif_length:
                umi5 = seq
        
        # Extract 3' UMI (if different from 5')
        umi3 = None
        if motif3_index != motif5_index:
            query_start = None
            query_end = None
            
            for i, ref_pos in enumerate(ref_positions):
                if ref_pos == umi3_start - 1 and query_start is None:
                    query_start = i
                if ref_pos == umi3_end - 1:
                    query_end = i + 1
                    break
            
            if query_start is not None and query_end is not None:
                seq = read.query_sequence[query_start:query_end]
                if len(seq) == motif_length:
                    umi3 = seq
        else:
            umi3 = umi5
        
        # Add to lists and dictionary if both UMIs found
        if umi5 is not None and umi3 is not None:
            raw_umi5_list.append(umi5)
            raw_umi3_list.append(umi3)
            raw_UMI_dict[read.query_name] = (umi5, umi3)
    
    bamfile.close()
    
    return raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict


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
    
    # Find motif positions sequentially (each motif gets the next available position)
    start_search = 0
    for i, motif in enumerate(barcodeMotifs):
        motif_length = len(motif)
        motif_lengths.append(motif_length)
        
        # Find the next occurrence of this motif starting from start_search
        pos = ref_seq.find(motif, start_search)
        if pos == -1:
            print(f"Motif '{motif}' not found in reference sequence starting from position {start_search}.")
            return None, None
        
        # Add this position to the motif indices
        motif_indices.append([pos])
        
        # Update start_search to begin after this motif for the next search
        start_search = pos + 1
    
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
