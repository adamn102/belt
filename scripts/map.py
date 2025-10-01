#!/usr/bin/env python3
"""
BELT Mapping Pipeline - Refactored

This module provides streamlined functionality for:
1. Mapping FASTQ files to reference sequences
2. Demultiplexing FASTQ files and then mapping

Core functionality only - focused and clean.
"""

# Standard library imports
import csv
import os
import re
import subprocess
from pathlib import Path

# Third-party imports
try:
    import pysam
    import yaml
    from Bio import SeqIO
    THIRD_PARTY_AVAILABLE = True
except ImportError:
    THIRD_PARTY_AVAILABLE = False

# Local imports
from utils import log

# Create a log instance
logger = log()

# Constants
DEFAULT_MIN_LENGTH = 1000
DEFAULT_MIN_QUALITY = 7
DEFAULT_THREADS = 4
DEFAULT_SEARCH_WINDOW = 200
DEFAULT_MISMATCH_THRESHOLD = 7

# Output directory names
DEMULTIPLEXED_DIR = "demultiplexed_fastq"
FILTERED_FASTQ_NAME = "filtered_input.fastq"
UNASSIGNED_NAME = "unassigned"
MOTIF_CSV_NAME = "motif_extraction_results.csv"


# =============================================================================
# CORE MAPPING FUNCTIONS
# =============================================================================

def filter_fastq(input_fastq, output_fastq, min_length=DEFAULT_MIN_LENGTH, min_quality=DEFAULT_MIN_QUALITY, verbose=True):
    """
    Filter FASTQ file by length and quality.
    
    Args:
        input_fastq (str): Path to input FASTQ file
        output_fastq (str): Path to output filtered FASTQ file
        min_length (int): Minimum read length
        min_quality (int): Minimum mean quality score
        verbose (bool): Print progress messages
    """
    if not THIRD_PARTY_AVAILABLE:
        raise ImportError("BioPython is required for FASTQ filtering")
    
    if verbose:
        logger.log_message(f"Filtering FASTQ: {input_fastq}")
        logger.log_message(f"Min length: {min_length}, Min quality: {min_quality}")
    
    filtered_count = 0
    total_count = 0
    
    with open(output_fastq, 'w') as outfile:
        for record in SeqIO.parse(input_fastq, "fastq"):
            total_count += 1
            
            # Check length
            if len(record.seq) < min_length:
                continue
                
            # Check quality
            if record.letter_annotations["phred_quality"]:
                mean_quality = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
                if mean_quality < min_quality:
                    continue
            
            # Write filtered read
            SeqIO.write(record, outfile, "fastq")
            filtered_count += 1
    
    if verbose:
        logger.log_message(f"Filtered {filtered_count}/{total_count} reads ({filtered_count/total_count*100:.1f}%)")


def align_and_sort_reads(reference, fastq, output_sorted_bam, threads=DEFAULT_THREADS, verbose=True):
    """
    Align reads using minimap2 and sort the output BAM file.
    
    Args:
        reference (str): Path to reference FASTA file
        fastq (str): Path to input FASTQ file
        output_sorted_bam (str): Path to output sorted BAM file
        threads (int): Number of threads to use
        verbose (bool): Print progress messages
    """
    if verbose:
        logger.log_message(f"Aligning reads: {fastq} -> {output_sorted_bam}")
    
    # Create temporary BAM file
    temp_bam = output_sorted_bam.replace('.bam', '_temp.bam')
    
    try:
        # Run minimap2 alignment
        minimap2_cmd = [
            'minimap2', '-ax', 'map-ont', '--secondary=no',
            '-t', str(threads), reference, fastq
        ]
        
        # Align with minimap2 (quiet mode)
        with open(temp_bam, 'wb') as bam_file:
            with open(os.devnull, 'w') as devnull:
                subprocess.run(minimap2_cmd, stdout=bam_file, stderr=devnull, check=True)
        
        # Sort BAM file
        sort_cmd = ['samtools', 'sort', '-@', str(threads), temp_bam, '-o', output_sorted_bam]
        
        # Sort BAM file (quiet mode)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(sort_cmd, stderr=devnull, check=True)
        
        # Index BAM file
        index_cmd = ['samtools', 'index', output_sorted_bam]
        with open(os.devnull, 'w') as devnull:
            subprocess.run(index_cmd, stderr=devnull, check=True)
        
        # Alignment completed
            
    finally:
        # Clean up temporary file
        if os.path.exists(temp_bam):
            os.remove(temp_bam)


# =============================================================================
# DEMULTIPLEXING FUNCTIONS
# =============================================================================

def load_config_yaml(config_file):
    """Load configuration from YAML file."""
    if not THIRD_PARTY_AVAILABLE:
        raise ImportError("PyYAML is required for YAML parsing")
    
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def load_sample_sheet_csv(sample_sheet_file):
    """Load sample sheet from CSV file."""
    sample_map = {}
    with open(sample_sheet_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_map[(row["StartBarcode"], row["EndBarcode"])] = row["Sample"]
    return sample_map


def extract_soft_clipped_sequence(read):
    """
    Extract soft-clipped sequences from a BAM read.
    
    Args:
        read: pysam AlignedSegment object
        
    Returns:
        tuple: (left_soft_clipped, right_soft_clipped, full_sequence)
    """
    if not read.cigartuples:
        return "", "", read.query_sequence
    
    left_soft = ""
    right_soft = ""
    
    # Check for soft clipping at the beginning
    if read.cigartuples[0][0] == 4:  # 4 = BAM_CSOFT_CLIP
        left_soft = read.query_sequence[:read.cigartuples[0][1]]
    
    # Check for soft clipping at the end
    if read.cigartuples[-1][0] == 4:  # 4 = BAM_CSOFT_CLIP
        right_soft = read.query_sequence[-read.cigartuples[-1][1]:]
    
    return left_soft, right_soft, read.query_sequence


def create_fastq_from_bam_read(read):
    """
    Create a FASTQ record from a BAM read.
    
    Args:
        read: pysam AlignedSegment object
        
    Returns:
        Bio.SeqIO.SeqRecord: FASTQ record
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # Get quality scores (convert from BAM format to FASTQ format)
    if read.query_qualities:
        quality_scores = "".join(chr(q + 33) for q in read.query_qualities)
    else:
        quality_scores = "I" * len(read.query_sequence)  # Default quality
    
    return SeqRecord(
        Seq(read.query_sequence),
        id=read.query_name,
        description="",
        letter_annotations={"phred_quality": [ord(q) - 33 for q in quality_scores]}
    )


def reorganize_bam_by_samples(bam_file, sample_assignments, output_dir, threads=DEFAULT_THREADS, verbose=True):
    """
    Reorganize a BAM file into separate sample BAM files based on read assignments.
    
    Args:
        bam_file (str): Path to input BAM file
        sample_assignments (dict): Dictionary mapping read names to sample names
        output_dir (str): Path to output directory
        threads (int): Number of threads to use
        verbose (bool): Print progress messages
        
    Returns:
        dict: Paths to output BAM files and statistics
    """
    if not THIRD_PARTY_AVAILABLE:
        raise ImportError("pysam is required for BAM processing")
    
    if verbose:
        logger.log_message(f"Reorganizing BAM file: {bam_file}")
    
    # Get unique sample names
    unique_samples = set(sample_assignments.values())
    
    # Create output BAM files
    sample_bam_files = {}
    sample_writers = {}
    
    for sample_name in unique_samples:
        output_bam = os.path.join(output_dir, f"{sample_name}.bam")
        sample_bam_files[sample_name] = output_bam
        
        # Create BAM writer for each sample
        with pysam.AlignmentFile(bam_file, "rb") as input_bam:
            header = input_bam.header.copy()
            sample_writers[sample_name] = pysam.AlignmentFile(output_bam, "wb", header=header)
    
    # Process reads and write to appropriate sample files
    read_counts = {sample: 0 for sample in unique_samples}
    total_reads = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as input_bam:
        for read in input_bam:
            if read.is_unmapped:
                continue
                
            total_reads += 1
            read_name = read.query_name
            
            # Get sample assignment for this read
            sample_name = sample_assignments.get(read_name, "unassigned")
            
            # Write to appropriate sample BAM file
            if sample_name in sample_writers:
                sample_writers[sample_name].write(read)
                read_counts[sample_name] += 1
    
    # Close all BAM writers
    for writer in sample_writers.values():
        writer.close()
    
    # Sort and index all output BAM files
    if verbose:
        logger.log_message("Sorting and indexing BAM files...")
    
    for sample_name, bam_path in sample_bam_files.items():
        if read_counts[sample_name] > 0:  # Only process files with reads
            # Sort BAM file
            sorted_bam = bam_path.replace('.bam', '_sorted.bam')
            sort_cmd = ['samtools', 'sort', '-@', str(threads), bam_path, '-o', sorted_bam]
            
            with open(os.devnull, 'w') as devnull:
                subprocess.run(sort_cmd, stderr=devnull, check=True)
            
            # Replace unsorted with sorted
            os.replace(sorted_bam, bam_path)
            
            # Index BAM file
            index_cmd = ['samtools', 'index', bam_path]
            with open(os.devnull, 'w') as devnull:
                subprocess.run(index_cmd, stderr=devnull, check=True)
    
    if verbose:
        logger.log_message(f"Reorganized {total_reads:,} reads into {len(unique_samples)} sample BAM files")
        for sample_name, count in read_counts.items():
            if count > 0:
                logger.log_message(f"  {sample_name}: {count:,} reads")
    
    return {
        'sample_bam_files': sample_bam_files,
        'read_counts': read_counts,
        'total_reads': total_reads
    }


def calculate_levenshtein_distance(seq1, seq2):
    """
    Calculate the Levenshtein distance (edit distance) between two sequences.
    This handles substitutions, insertions, and deletions.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        int: Levenshtein distance (minimum number of edits needed)
    """
    if len(seq1) == 0:
        return len(seq2)
    if len(seq2) == 0:
        return len(seq1)
    
    # Create a matrix to store distances
    matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    
    # Initialize first row and column
    for i in range(len(seq1) + 1):
        matrix[i][0] = i
    for j in range(len(seq2) + 1):
        matrix[0][j] = j
    
    # Fill the matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if seq1[i-1] == seq2[j-1]:
                cost = 0
            else:
                cost = 1
            
            matrix[i][j] = min(
                matrix[i-1][j] + 1,      # deletion
                matrix[i][j-1] + 1,      # insertion
                matrix[i-1][j-1] + cost  # substitution
            )
    
    return matrix[len(seq1)][len(seq2)]


def find_best_barcode_match(extracted_barcode, barcodes, max_mismatches):
    """
    Find the best matching barcode using Levenshtein distance.
    
    Args:
        extracted_barcode (str): The barcode sequence extracted from the read
        barcodes (dict): Dictionary mapping barcode names to sequences
        max_mismatches (int): Maximum Levenshtein distance allowed
        
    Returns:
        tuple: (best_match_name, levenshtein_distance) or (None, None) if no good match
    """
    best_match = None
    best_distance = max_mismatches + 1
    
    for barcode_name, barcode_seq in barcodes.items():
        distance = calculate_levenshtein_distance(extracted_barcode, barcode_seq)
        if distance <= max_mismatches and distance < best_distance:
            best_match = barcode_name
            best_distance = distance
    
    if best_match is not None:
        return best_match, best_distance
    else:
        return None, None


def checkMotif(sequence, side):
    """
    Extract barcode using motif-based logic.
    
    Args:
        sequence (str): The soft-clipped sequence to search
        side (str): "left" or "right" to determine extraction logic
        
    Returns:
        str or None: Extracted barcode sequence or None if no motif found
    """
    if side == "left":
        motif_index = sequence.find("AAGGTTAA")
        if motif_index != -1:
            remaining_seq = sequence[motif_index + 8:]
            if len(remaining_seq) >= 20:
                return remaining_seq[:20]
            else:
                return remaining_seq + 'N' * (20 - len(remaining_seq))
        else:
            return None
    elif side == "right":
        motif_index = sequence.rfind("AAGGTTAA")
        if motif_index != -1:
            remaining_seq = sequence[:motif_index]
            if len(remaining_seq) >= 20:
                return remaining_seq[-20:]
            else:
                return 'N' * (20 - len(remaining_seq)) + remaining_seq
        else:
            return None


def assign_read_to_sample_from_barcodes(start_barcode, end_barcode, barcodes, sample_map, mismatch_tolerance=DEFAULT_MISMATCH_THRESHOLD):
    """
    Assign a read to a sample based on extracted barcode sequences with mismatch tolerance.
    
    Args:
        start_barcode (str): Start barcode sequence (20bp from left soft-clipped region)
        end_barcode (str): End barcode sequence (20bp from right soft-clipped region)
        barcodes (dict): Dictionary mapping barcode names to sequences
        sample_map (dict): Dictionary mapping (start_barcode, end_barcode) to sample names
        mismatch_tolerance (int): Maximum number of mismatches allowed
        
    Returns:
        tuple: (sample_name, assignment_reason, barcode_info)
    """
    # Find best matches with mismatch tolerance
    start_name, start_mismatches = find_best_barcode_match(start_barcode, barcodes, mismatch_tolerance)
    end_name, end_mismatches = find_best_barcode_match(end_barcode, barcodes, mismatch_tolerance)
    
    if not start_name or not end_name:
        return "unassigned", "barcode_not_found", {
            'start_barcode': start_barcode,
            'end_barcode': end_barcode,
            'start_name': start_name or 'unknown',
            'end_name': end_name or 'unknown',
            'start_mismatches': start_mismatches if start_name else None,
            'end_mismatches': end_mismatches if end_name else None
        }
    
    # Check if this barcode combination exists in sample map
    sample_name = sample_map.get((start_name, end_name), "unassigned")
    
    if sample_name == "unassigned":
        return "unassigned", "sample_not_found", {
            'start_barcode': start_barcode,
            'end_barcode': end_barcode,
            'start_name': start_name,
            'end_name': end_name,
            'start_mismatches': start_mismatches,
            'end_mismatches': end_mismatches
        }
    
    return sample_name, "success", {
        'start_barcode': start_barcode,
        'end_barcode': end_barcode,
        'start_name': start_name,
        'end_name': end_name,
        'start_mismatches': start_mismatches,
        'end_mismatches': end_mismatches
    }


def extract_barcodes_from_bam(bam_file, config_file, sample_sheet_file, output_dir):
    """
    Extract barcodes from soft-clipped regions of mapped reads and reorganize into sample BAM files.
    This optimized version works directly with BAM files instead of converting to FASTQ.
    
    Args:
        bam_file (str): Path to mapped BAM file
        config_file (str): Path to config YAML file
        sample_sheet_file (str): Path to sample sheet CSV file
        output_dir (str): Path to output directory
        
    Returns:
        dict: Demultiplexing results
    """
    if not THIRD_PARTY_AVAILABLE:
        raise ImportError("BioPython and pysam are required for BAM processing")
    
    # Load configuration and sample sheet
    config = load_config_yaml(config_file)
    sample_map = load_sample_sheet_csv(sample_sheet_file)
    barcodes = config["barcodes"]
    
    # Create output directories
    demux_dir = os.path.join(output_dir, DEMULTIPLEXED_DIR)
    os.makedirs(demux_dir, exist_ok=True)
    
    # Get all unique sample names from sample sheet
    unique_samples = set(sample_map.values())
    
    # Create CSV results file
    csv_file = open(os.path.join(demux_dir, MOTIF_CSV_NAME), 'w')
    csv_writer = csv.writer(csv_file)
    
    # Write CSV header
    header = [
        'read_name', 'sequence_length', 'assigned_sample', 'assignment_reason',
        'start_barcode_seq', 'end_barcode_seq', 'start_barcode_name', 'end_barcode_name',
        'start_mismatches', 'end_mismatches', 'left_soft_clipped', 'right_soft_clipped'
    ]
    csv_writer.writerow(header)
    
    # Process reads from BAM file and collect assignments
    total_reads = 0
    assigned_reads = 0
    reads_without_softclip = 0
    
    # Track assignment reasons
    assignment_stats = {
        'success': 0,
        'sample_not_found': 0,
        'barcode_not_found': 0,
        'motif_success': 0,
        'motif_not_found': 0
    }
    
    # Track barcode detection statistics
    barcode_stats = {
        'both_barcodes': 0,
        'only_5p_barcode': 0,
        'only_3p_barcode': 0,
        'no_barcodes': 0,
        'no_softclip': 0
    }
    
    # Store sample assignments for BAM reorganization
    sample_assignments = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        for read in bamfile:
            if read.is_unmapped:
                continue
                
            total_reads += 1
            motif_forward = None
            motif_reverse = None
            # Extract soft-clipped sequences
            left_soft, right_soft, full_sequence = extract_soft_clipped_sequence(read)
            
            # Check if read has any soft-clipped regions
            has_softclip = len(left_soft) > 0 or len(right_soft) > 0
            if not has_softclip:
                reads_without_softclip += 1
                barcode_stats['no_softclip'] += 1
                # Assign to unassigned
                sample_assignments[read.query_name] = "unassigned"
                continue
            
            # Extract 20bp from each soft-clipped region (end of left, beginning of right)
            forward_barcode = left_soft[-20:] if len(left_soft) >= 20 else left_soft
            reverse_barcode = right_soft[:20] if len(right_soft) >= 20 else right_soft

            # Try to find matches using Levenshtein distance for initial extraction
            forward_match, forward_distance = find_best_barcode_match(forward_barcode, barcodes, max_mismatches=DEFAULT_MISMATCH_THRESHOLD)
            reverse_match, reverse_distance = find_best_barcode_match(reverse_barcode, barcodes, max_mismatches=DEFAULT_MISMATCH_THRESHOLD)

            # If either barcode doesn't match well enough, try motif-based extraction
            if not forward_match or not reverse_match:
                # Try motif extraction for barcodes that didn't match
                if not forward_match:
                    motif_forward = checkMotif(left_soft, "left")
                    if motif_forward:
                        # Check if motif-extracted forward barcode matches with Levenshtein distance
                        motif_forward_match, motif_forward_distance = find_best_barcode_match(motif_forward, barcodes, max_mismatches=DEFAULT_MISMATCH_THRESHOLD)
                        if True:
                            forward_barcode = motif_forward
                            forward_match = motif_forward_match
                            forward_distance = motif_forward_distance
                if not reverse_match:
                    motif_reverse = checkMotif(right_soft, "right")
                    if motif_reverse:
                        # Check if motif-extracted reverse barcode matches with Levenshtein distance
                        motif_reverse_match, motif_reverse_distance = find_best_barcode_match(motif_reverse, barcodes, max_mismatches=DEFAULT_MISMATCH_THRESHOLD)
                        if True:
                            reverse_barcode = motif_reverse
                            reverse_match = motif_reverse_match
                            reverse_distance = motif_reverse_distance
            # Track barcode detection patterns
            has_5p = (forward_match is not None) and (forward_match != '')
            has_3p = (reverse_match is not None) and (reverse_match != '')
            
            if has_5p and has_3p:
                barcode_stats['both_barcodes'] += 1
            elif has_5p and not has_3p:
                barcode_stats['only_5p_barcode'] += 1
            elif not has_5p and has_3p:
                barcode_stats['only_3p_barcode'] += 1
            else:
                barcode_stats['no_barcodes'] += 1
            
            # Now assign to sample using the final barcodes
            sample_name, assignment_reason, barcode_info = assign_read_to_sample_from_barcodes(
                forward_barcode, reverse_barcode, barcodes, sample_map
            )
            
            # Store assignment for BAM reorganization
            sample_assignments[read.query_name] = sample_name
            
            # Update assignment reason if motif extraction was used
            original_forward = left_soft[-20:] if len(left_soft) >= 20 else left_soft
            original_reverse = right_soft[:20] if len(right_soft) >= 20 else right_soft
            
            if (forward_barcode != original_forward) or (reverse_barcode != original_reverse):
                if sample_name != "unassigned":
                    assignment_reason = f"motif_{assignment_reason}"
            
            # Track assignment statistics
            if assignment_reason == "success":
                assignment_stats['success'] += 1
                assigned_reads += 1
            elif assignment_reason.startswith("motif_"):
                assignment_stats['motif_success'] += 1
                assigned_reads += 1
            elif assignment_reason == "sample_not_found":
                assignment_stats['sample_not_found'] += 1
            elif assignment_reason == "barcode_not_found":
                assignment_stats['barcode_not_found'] += 1
            
            # Write CSV row - use the barcodes that were actually used for assignment
            if barcode_info:
                # Use the actual sequences that were used for assignment (forward_barcode and reverse_barcode)
                row = [
                    read.query_name,
                    len(full_sequence),
                    sample_name,
                    assignment_reason,
                    motif_forward,  # The actual sequence that was used for assignment
                    motif_reverse,  # The actual sequence that was used for assignment
                    barcode_info.get('start_name', ''),
                    barcode_info.get('end_name', ''),
                    barcode_info.get('start_mismatches', ''),
                    barcode_info.get('end_mismatches', ''),
                    left_soft,
                    right_soft
                ]
            else:
                row = [read.query_name, len(full_sequence), sample_name, assignment_reason, '', '', '', '', '', '', left_soft, right_soft]
            
            csv_writer.writerow(row)
    
    # Close CSV file
    csv_file.close()
    
    # Reorganize BAM file into sample-specific BAM files
    print(f"\n🔄 Reorganizing BAM files by sample assignments...")
    bam_reorg_results = reorganize_bam_by_samples(bam_file, sample_assignments, demux_dir, verbose=True)
    
    # Print demultiplexing statistics (always show)
    print(f"\n" + "="*80)
    print(f"📊 DEMULTIPLEXING STATISTICS")
    print(f"="*80)
    print(f"   Total reads processed: {total_reads:,}")
    print(f"   Reads without soft-clip: {reads_without_softclip:,} ({reads_without_softclip/total_reads*100:.1f}%)")
    print(f"   Reads with soft-clip: {total_reads - reads_without_softclip:,} ({(total_reads - reads_without_softclip)/total_reads*100:.1f}%)")
    print(f"   Successfully assigned: {assigned_reads:,} ({assigned_reads/(total_reads - reads_without_softclip)*100:.1f}% of soft-clipped reads)")
    print(f"   Unassigned reads: {(total_reads - reads_without_softclip) - assigned_reads:,} ({((total_reads - reads_without_softclip) - assigned_reads)/(total_reads - reads_without_softclip)*100:.1f}% of soft-clipped reads)")
    print(f"   Unique samples found: {len(unique_samples)}")
    print(f"   Output directory: {demux_dir}")
    print()
    
    # Show barcode detection breakdown
    softclip_reads = total_reads - reads_without_softclip
    print(f"🔍 BARCODE DETECTION BREAKDOWN (reads with soft-clip only):")
    print(f"-" * 80)
    print(f"   Both 5' and 3' barcodes found:     {barcode_stats['both_barcodes']:6,} reads ({barcode_stats['both_barcodes']/softclip_reads*100:5.1f}%)")
    print(f"   Only 5' barcode found:            {barcode_stats['only_5p_barcode']:6,} reads ({barcode_stats['only_5p_barcode']/softclip_reads*100:5.1f}%)")
    print(f"   Only 3' barcode found:            {barcode_stats['only_3p_barcode']:6,} reads ({barcode_stats['only_3p_barcode']/softclip_reads*100:5.1f}%)")
    print(f"   No barcodes found:                {barcode_stats['no_barcodes']:6,} reads ({barcode_stats['no_barcodes']/softclip_reads*100:5.1f}%)")
    print(f"-" * 80)
    print(f"   No soft-clip regions:             {barcode_stats['no_softclip']:6,} reads ({barcode_stats['no_softclip']/total_reads*100:5.1f}%)")
    print()
    
    # Show assignment reason statistics
    print(f"📋 ASSIGNMENT REASON BREAKDOWN:")
    print(f"-" * 80)
    print(f"   Success (soft-clip):   {assignment_stats['success']:6,} reads ({assignment_stats['success']/softclip_reads*100:5.1f}%)")
    print(f"   Success (motif):      {assignment_stats['motif_success']:6,} reads ({assignment_stats['motif_success']/softclip_reads*100:5.1f}%)")
    print(f"   Sample not found:      {assignment_stats['sample_not_found']:6,} reads ({assignment_stats['sample_not_found']/softclip_reads*100:5.1f}%)")
    print(f"   Barcode not found:    {assignment_stats['barcode_not_found']:6,} reads ({assignment_stats['barcode_not_found']/softclip_reads*100:5.1f}%)")
    print(f"   Motif not found:      {assignment_stats['motif_not_found']:6,} reads ({assignment_stats['motif_not_found']/softclip_reads*100:5.1f}%)")
    print()
    
    # Show sample assignments
    print(f"📋 SAMPLE ASSIGNMENTS:")
    print(f"-" * 60)
    
    # Get read counts from BAM reorganization results
    sample_counts = bam_reorg_results['read_counts']
    
    # Sort samples by read count (descending)
    sorted_samples = sorted(sample_counts.items(), key=lambda x: x[1], reverse=True)
    
    for sample_name, count in sorted_samples:
        percentage = (count / total_reads) * 100 if total_reads > 0 else 0
        if sample_name == "unassigned":
            print(f"   {sample_name:25} {count:6,} reads ({percentage:5.1f}%)")
        else:
            print(f"   {sample_name:25} {count:6,} reads ({percentage:5.1f}%)")
    
    print(f"="*60)
    print()
    
    return {
        'total_reads': total_reads,
        'assigned_reads': assigned_reads,
        'unassigned_reads': (total_reads - reads_without_softclip) - assigned_reads,
        'reads_without_softclip': reads_without_softclip,
        'barcode_stats': barcode_stats,
        'assignment_stats': assignment_stats,
        'demultiplexed_dir': demux_dir,
        'unique_samples': list(unique_samples),
        'sample_bam_files': bam_reorg_results['sample_bam_files'],
        'sample_assignments': sample_assignments
    }


# =============================================================================
# MAIN PIPELINE FUNCTIONS
# =============================================================================

def map_pipeline(reference, fastq, output_dir, min_length=DEFAULT_MIN_LENGTH, min_quality=DEFAULT_MIN_QUALITY, 
                threads=DEFAULT_THREADS, config_file=None, sample_sheet=None, verbose=False):
    """
    Main mapping pipeline - maps FASTQ files with optional demultiplexing.
    
    Args:
        reference (str): Path to reference FASTA file
        fastq (str): Path to input FASTQ file
        output_dir (str): Path to output directory
        min_length (int): Minimum read length
        min_quality (int): Minimum mean quality score
        threads (int): Number of threads to use
        config_file (str): Optional path to config YAML file for demultiplexing
        sample_sheet (str): Optional path to sample sheet CSV file for demultiplexing
        verbose (bool): Print progress messages
        
    Returns:
        dict: Pipeline results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    results = {
        'reference': reference,
        'input_fastq': fastq,
        'output_dir': output_dir,
        'filtered_fastq': None,
        'bam_file': None,
        'demultiplexed_files': []
    }
    
    # Step 1: Filter FASTQ
    filtered_fastq = os.path.join(output_dir, FILTERED_FASTQ_NAME)
    filter_fastq(fastq, filtered_fastq, min_length, min_quality, verbose)
    results['filtered_fastq'] = filtered_fastq
    
    # Step 2: Map first, then demultiplex if config and sample sheet provided
    if config_file and sample_sheet:
        if verbose:
            logger.log_message("Mapping reads to reference first...")
        
        # Map the filtered FASTQ to reference
        output_bam = os.path.join(output_dir, "mapped.bam")
        align_and_sort_reads(reference, filtered_fastq, output_bam, threads, verbose)
        results['bam_file'] = output_bam
        
        if verbose:
            logger.log_message("Extracting barcodes from soft-clipped regions...")
        
        # Extract barcodes from soft-clipped regions and reorganize BAM files directly
        demux_results = extract_barcodes_from_bam(output_bam, config_file, sample_sheet, output_dir)
        
        # Get the sample BAM files that were created during reorganization
        sample_bam_files = demux_results.get('sample_bam_files', {})
        
        if verbose:
            logger.log_message(f"Created {len(sample_bam_files)} sample BAM files (already sorted and indexed)")
        
        # Add sample BAM files to results
        for sample_name, bam_path in sample_bam_files.items():
            results['demultiplexed_files'].append(bam_path)
        
        results['demux_results'] = demux_results
    else:
        # Step 3: Direct mapping
        output_bam = os.path.join(output_dir, "mapped.bam")
        align_and_sort_reads(reference, filtered_fastq, output_bam, threads, verbose)
        results['bam_file'] = output_bam
    
    if verbose:
        logger.log_message("Pipeline completed successfully!")
    
    return results


def main():
    """Main entry point for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description="BELT Mapping Pipeline")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH, help="Minimum read length")
    parser.add_argument("--min-quality", type=int, default=DEFAULT_MIN_QUALITY, help="Minimum quality score")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads")
    parser.add_argument("--config", help="Config YAML file for demultiplexing")
    parser.add_argument("--sample-sheet", help="Sample sheet CSV file for demultiplexing")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    try:
        results = map_pipeline(
            reference=args.reference,
            fastq=args.fastq,
            output_dir=args.output_dir,
            min_length=args.min_length,
            min_quality=args.min_quality,
            threads=args.threads,
            config_file=args.config,
            sample_sheet=args.sample_sheet,
            verbose=args.verbose
        )
        
        print("✅ Mapping pipeline completed successfully!")
        print(f"Results: {results}")
        
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
