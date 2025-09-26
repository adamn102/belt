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
DEFAULT_MISMATCH_THRESHOLD = 2

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


def compile_barcode_regex(config):
    """
    Compile regex patterns for barcode detection.
    
    Args:
        config (dict): Configuration dictionary with motifs
        
    Returns:
        tuple: (start_regex, end_regex) compiled regex objects
    """
    motifs = config["motifs"]
    
    # Start barcode pattern: 7bp prefix + 20bp barcode + 10bp suffix
    start_pattern = rf"{motifs['5p_prefix']}([ACGT]{{20}}){motifs['5p_suffix']}"
    start_regex = re.compile(start_pattern)
    
    # End barcode pattern: 10bp prefix + 20bp barcode + 7bp suffix  
    end_pattern = rf"{motifs['3p_prefix']}([ACGT]{{20}}){motifs['3p_suffix']}"
    end_regex = re.compile(end_pattern)
    
    return start_regex, end_regex


def calculate_mismatches(seq1, seq2):
    """
    Calculate the number of mismatches between two sequences.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        int: Number of mismatches
    """
    if len(seq1) != len(seq2):
        return max(len(seq1), len(seq2))  # Return max length as penalty for length mismatch
    
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches


def find_best_barcode_match(extracted_barcode, barcodes, max_mismatches):
    """
    Find the best matching barcode with mismatch tolerance.
    
    Args:
        extracted_barcode (str): The 20bp barcode sequence extracted from the read
        barcodes (dict): Dictionary mapping barcode names to sequences
        max_mismatches (int): Maximum number of mismatches allowed
        
    Returns:
        tuple: (best_match_name, mismatch_count) or (None, None) if no good match
    """
    best_match = None
    best_mismatches = max_mismatches + 1
    
    for barcode_name, barcode_seq in barcodes.items():
        mismatches = calculate_mismatches(extracted_barcode, barcode_seq)
        if mismatches <= max_mismatches and mismatches < best_mismatches:
            best_match = barcode_name
            best_mismatches = mismatches
    
    if best_match is not None:
        return best_match, best_mismatches
    else:
        return None, None


def assign_read_to_sample(sequence, start_regex, end_regex, barcodes, sample_map, mismatch_tolerance=DEFAULT_MISMATCH_THRESHOLD):
    """
    Assign a read to a sample based on barcode detection with mismatch tolerance.
    
    Args:
        sequence (str): DNA sequence to analyze
        start_regex: Compiled regex for start barcode detection
        end_regex: Compiled regex for end barcode detection
        barcodes (dict): Dictionary mapping barcode names to sequences
        sample_map (dict): Dictionary mapping (start_barcode, end_barcode) to sample names
        mismatch_tolerance (int): Maximum number of mismatches allowed in barcode sequences
        
    Returns:
        tuple: (sample_name, assignment_reason, barcode_info)
    """
    # Find start and end barcodes (only check first/last ~100 bases for efficiency)
    start_match = start_regex.search(sequence[:100])
    end_match = end_regex.search(sequence[-100:])
    
    if not start_match or not end_match:
        return "unassigned", "no_barcode_match", None
    
    start_barcode_seq = start_match.group(1)
    end_barcode_seq = end_match.group(1)
    
    # Match to known barcodes with mismatch tolerance
    start_name, start_mismatches = find_best_barcode_match(start_barcode_seq, barcodes, mismatch_tolerance)
    end_name, end_mismatches = find_best_barcode_match(end_barcode_seq, barcodes, mismatch_tolerance)
    
    if not start_name or not end_name:
        return "unassigned", "unknown_barcode", {
            'start_barcode': start_barcode_seq,
            'end_barcode': end_barcode_seq,
            'start_name': start_name,
            'end_name': end_name,
            'start_mismatches': start_mismatches,
            'end_mismatches': end_mismatches
        }
    
    # Look up sample
    sample_name = sample_map.get((start_name, end_name), None)
    if sample_name:
        return sample_name, "complete_barcode_match", {
            'start_barcode': start_barcode_seq,
            'end_barcode': end_barcode_seq,
            'start_name': start_name,
            'end_name': end_name,
            'start_mismatches': start_mismatches,
            'end_mismatches': end_mismatches
        }
    else:
        return "unassigned", "no_sample_match", {
            'start_barcode': start_barcode_seq,
            'end_barcode': end_barcode_seq,
            'start_name': start_name,
            'end_name': end_name,
            'start_mismatches': start_mismatches,
            'end_mismatches': end_mismatches
        }


def demultiplex_fastq_barcode_based(input_fastq, output_dir, config_file, sample_sheet_file):
    """
    Demultiplex FASTQ file based on barcode detection using regex patterns.
    
    Args:
        input_fastq (str): Path to input FASTQ file
        output_dir (str): Path to output directory
        config_file (str): Path to config YAML file
        sample_sheet_file (str): Path to sample sheet CSV file
        
    Returns:
        dict: Results summary
    """
    if not THIRD_PARTY_AVAILABLE:
        raise ImportError("BioPython is required for FASTQ processing")
    
    # Load configuration and sample sheet
    config = load_config_yaml(config_file)
    sample_map = load_sample_sheet_csv(sample_sheet_file)
    
    # Compile regex patterns
    start_regex, end_regex = compile_barcode_regex(config)
    barcodes = config["barcodes"]
    
    # Create output directories
    demux_dir = os.path.join(output_dir, DEMULTIPLEXED_DIR)
    os.makedirs(demux_dir, exist_ok=True)
    
    # Get all unique sample names from sample sheet
    unique_samples = set(sample_map.values())
    
    # Create output files
    sample_files = {}
    for sample_name in unique_samples:
        sample_files[sample_name] = open(os.path.join(demux_dir, f"{sample_name}.fastq"), 'w')
    
    # Create unassigned file
    unassigned_file = open(os.path.join(demux_dir, f"{UNASSIGNED_NAME}.fastq"), 'w')
    
    # Create CSV results file
    csv_file = open(os.path.join(demux_dir, MOTIF_CSV_NAME), 'w')
    csv_writer = csv.writer(csv_file)
    
    # Write CSV header
    header = [
        'read_name', 'sequence_length', 'assigned_sample', 'assignment_reason',
        'start_barcode_seq', 'end_barcode_seq', 'start_barcode_name', 'end_barcode_name',
        'start_mismatches', 'end_mismatches'
    ]
    csv_writer.writerow(header)
    
    # Process reads
    total_reads = 0
    assigned_reads = 0
    sample_counts = {}  # Track reads per sample
    
    for record in SeqIO.parse(input_fastq, "fastq"):
        total_reads += 1
        sequence = str(record.seq)
        
        # Assign read to sample
        sample_name, assignment_reason, barcode_info = assign_read_to_sample(
            sequence, start_regex, end_regex, barcodes, sample_map, DEFAULT_MISMATCH_THRESHOLD
        )
        
        # Write to appropriate file
        if sample_name != "unassigned":
            if sample_name in sample_files:
                SeqIO.write(record, sample_files[sample_name], "fastq")
                assigned_reads += 1
                # Track sample counts
                sample_counts[sample_name] = sample_counts.get(sample_name, 0) + 1
        else:
            SeqIO.write(record, unassigned_file, "fastq")
            sample_counts["unassigned"] = sample_counts.get("unassigned", 0) + 1
        
        # Write CSV row
        if barcode_info:
            row = [
                record.id, len(sequence), sample_name, assignment_reason,
                barcode_info.get('start_barcode', ''),
                barcode_info.get('end_barcode', ''),
                barcode_info.get('start_name', ''),
                barcode_info.get('end_name', ''),
                barcode_info.get('start_mismatches', ''),
                barcode_info.get('end_mismatches', '')
            ]
        else:
            row = [record.id, len(sequence), sample_name, assignment_reason, '', '', '', '', '', '']
        
        csv_writer.writerow(row)
    
    # Close files
    for file_handle in sample_files.values():
        file_handle.close()
    unassigned_file.close()
    csv_file.close()
    
    # Print demultiplexing statistics (always show)
    print(f"\n" + "="*60)
    print(f"📊 DEMULTIPLEXING STATISTICS")
    print(f"="*60)
    print(f"   Total reads processed: {total_reads:,}")
    print(f"   Successfully assigned: {assigned_reads:,} ({assigned_reads/total_reads*100:.1f}%)")
    print(f"   Unassigned reads: {total_reads - assigned_reads:,} ({(total_reads - assigned_reads)/total_reads*100:.1f}%)")
    print(f"   Unique samples found: {len(unique_samples)}")
    print(f"   Output directory: {demux_dir}")
    print()
    
    # Show sample assignments
    print(f"📋 SAMPLE ASSIGNMENTS:")
    print(f"-" * 60)
    
    # Sort samples by read count (descending)
    sorted_samples = sorted(sample_counts.items(), key=lambda x: x[1], reverse=True)
    
    for sample_name, count in sorted_samples:
        percentage = (count / total_reads) * 100
        if sample_name == "unassigned":
            print(f"   {sample_name:25} {count:6,} reads ({percentage:5.1f}%)")
        else:
            print(f"   {sample_name:25} {count:6,} reads ({percentage:5.1f}%)")
    
    print(f"="*60)
    print()
    
    return {
        'total_reads': total_reads,
        'assigned_reads': assigned_reads,
        'unassigned_reads': total_reads - assigned_reads,
        'demultiplexed_dir': demux_dir,
        'unique_samples': list(unique_samples)
    }


# =============================================================================
# MAIN PIPELINE FUNCTIONS
# =============================================================================

def map_pipeline(reference, fastq, output_dir, min_length=DEFAULT_MIN_LENGTH, min_quality=DEFAULT_MIN_QUALITY, 
                threads=DEFAULT_THREADS, config_file=None, sample_sheet=None, verbose=True):
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
    
    # Step 2: Demultiplex if config and sample sheet provided
    if config_file and sample_sheet:
        if verbose:
            logger.log_message(f"Demultiplexing with config: {config_file} and sample sheet: {sample_sheet}")
        
        try:
            demux_results = demultiplex_fastq_barcode_based(filtered_fastq, output_dir, config_file, sample_sheet)
            
            if verbose:
                logger.log_message(f"Demultiplexing completed: {demux_results['assigned_reads']} assigned, {demux_results['unassigned_reads']} unassigned")
                logger.log_message(f"Demultiplexed directory: {demux_results['demultiplexed_dir']}")
        except Exception as e:
            logger.log_message(f"Error during demultiplexing: {e}")
            raise
        
        # Wait for demultiplexing to complete fully before mapping
        if verbose:
            logger.log_message("Demultiplexing completed, starting mapping phase...")
        
        # Map each demultiplexed file
        demux_dir = demux_results['demultiplexed_dir']
        demux_files = [f for f in os.listdir(demux_dir) if f.endswith('.fastq')]
        
        if verbose:
            logger.log_message(f"Found {len(demux_files)} demultiplexed files to map")
        
        for sample_file in demux_files:
            sample_name = sample_file.replace('.fastq', '')
            sample_fastq = os.path.join(demux_dir, sample_file)
            sample_bam = os.path.join(output_dir, f"{sample_name}.bam")
            
        if verbose:
            logger.log_message(f"Mapping {sample_name}...")
        
        align_and_sort_reads(reference, sample_fastq, sample_bam, threads, verbose)
        results['demultiplexed_files'].append(sample_bam)
        
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
