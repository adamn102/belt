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
DEFAULT_MISMATCH_THRESHOLD = 6

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


# Note: compile_barcode_regex function removed - no longer needed for BAM-based approach


# Note: assign_read_to_sample function removed - replaced by BAM-based approach


# Note: demultiplex_fastq_barcode_based function removed - replaced by BAM-based approach


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


def assign_read_to_sample_from_barcodes(start_barcode, end_barcode, barcodes, sample_map):
    """
    Assign a read to a sample based on extracted barcode sequences.
    
    Args:
        start_barcode (str): Start barcode sequence (20bp from left soft-clipped region)
        end_barcode (str): End barcode sequence (20bp from right soft-clipped region)
        barcodes (dict): Dictionary mapping barcode names to sequences
        sample_map (dict): Dictionary mapping (start_barcode, end_barcode) to sample names
        
    Returns:
        tuple: (sample_name, assignment_reason, barcode_info)
    """
    # Match to known barcodes
    start_name = next((k for k, v in barcodes.items() if v == start_barcode), None)
    end_name = next((k for k, v in barcodes.items() if v == end_barcode), None)
    
    if not start_name or not end_name:
        return "unassigned", "barcode_not_found", {
            'start_barcode': start_barcode,
            'end_barcode': end_barcode,
            'start_name': start_name or 'unknown',
            'end_name': end_name or 'unknown'
        }
    
    # Check if this barcode combination exists in sample map
    sample_name = sample_map.get((start_name, end_name), "unassigned")
    
    if sample_name == "unassigned":
        return "unassigned", "sample_not_found", {
            'start_barcode': start_barcode,
            'end_barcode': end_barcode,
            'start_name': start_name,
            'end_name': end_name
        }
    
    return sample_name, "success", {
        'start_barcode': start_barcode,
        'end_barcode': end_barcode,
        'start_name': start_name,
        'end_name': end_name
    }


def extract_barcodes_from_bam(bam_file, config_file, sample_sheet_file, output_dir):
    """
    Extract barcodes from soft-clipped regions of mapped reads.
    
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
        'left_soft_clipped', 'right_soft_clipped'
    ]
    csv_writer.writerow(header)
    
    # Process reads from BAM file
    total_reads = 0
    assigned_reads = 0
    
    # Track assignment reasons
    assignment_stats = {
        'success': 0,
        'sample_not_found': 0,
        'barcode_not_found': 0
    }
    
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        for read in bamfile:
            if read.is_unmapped:
                continue
                
            total_reads += 1
            
            # Extract soft-clipped sequences
            left_soft, right_soft, full_sequence = extract_soft_clipped_sequence(read)
            
            # Extract 20bp from each soft-clipped region
            start_barcode = left_soft[-20:] if len(left_soft) >= 20 else left_soft
            end_barcode = right_soft[:20] if len(right_soft) >= 20 else right_soft
            
            # Assign read to sample based on extracted barcodes
            sample_name, assignment_reason, barcode_info = assign_read_to_sample_from_barcodes(
                start_barcode, end_barcode, barcodes, sample_map
            )
            
            # Track assignment statistics
            if assignment_reason == "success":
                assignment_stats['success'] += 1
            elif assignment_reason == "sample_not_found":
                assignment_stats['sample_not_found'] += 1
            elif assignment_reason == "barcode_not_found":
                assignment_stats['barcode_not_found'] += 1
            
            # Create FASTQ record from BAM read
            fastq_record = create_fastq_from_bam_read(read)
            
            # Write to appropriate file
            if sample_name != "unassigned":
                if sample_name in sample_files:
                    SeqIO.write(fastq_record, sample_files[sample_name], "fastq")
                    assigned_reads += 1
            else:
                SeqIO.write(fastq_record, unassigned_file, "fastq")
            
            # Write CSV row
            if barcode_info:
                row = [
                    read.query_name,
                    len(full_sequence),
                    sample_name,
                    assignment_reason,
                    barcode_info.get('start_barcode', ''),
                    barcode_info.get('end_barcode', ''),
                    barcode_info.get('start_name', ''),
                    barcode_info.get('end_name', ''),
                    left_soft,
                    right_soft
                ]
            else:
                row = [read.query_name, len(full_sequence), sample_name, assignment_reason, '', '', '', '', left_soft, right_soft]
            
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
    
    # Show assignment reason statistics
    print(f"🔍 ASSIGNMENT REASON BREAKDOWN:")
    print(f"-" * 60)
    print(f"   Success:              {assignment_stats['success']:6,} reads ({assignment_stats['success']/total_reads*100:5.1f}%)")
    print(f"   Sample not found:      {assignment_stats['sample_not_found']:6,} reads ({assignment_stats['sample_not_found']/total_reads*100:5.1f}%)")
    print(f"   Barcode not found:    {assignment_stats['barcode_not_found']:6,} reads ({assignment_stats['barcode_not_found']/total_reads*100:5.1f}%)")
    print()
    
    # Show sample assignments
    print(f"📋 SAMPLE ASSIGNMENTS:")
    print(f"-" * 60)
    
    # Count reads per sample
    sample_counts = {}
    for sample_name in unique_samples:
        sample_file = os.path.join(demux_dir, f"{sample_name}.fastq")
        if os.path.exists(sample_file):
            with open(sample_file, 'r') as f:
                # Count FASTQ records (4 lines per record)
                line_count = sum(1 for line in f)
                sample_counts[sample_name] = line_count // 4
        else:
            sample_counts[sample_name] = 0
    
    # Add unassigned count
    unassigned_file_path = os.path.join(demux_dir, f"{UNASSIGNED_NAME}.fastq")
    if os.path.exists(unassigned_file_path):
        with open(unassigned_file_path, 'r') as f:
            line_count = sum(1 for line in f)
            sample_counts["unassigned"] = line_count // 4
    else:
        sample_counts["unassigned"] = 0
    
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
        'unassigned_reads': total_reads - assigned_reads,
        'demultiplexed_dir': demux_dir,
        'unique_samples': list(unique_samples)
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
        
        # Extract barcodes from soft-clipped regions and demultiplex
        demux_results = extract_barcodes_from_bam(output_bam, config_file, sample_sheet, output_dir)
        
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
