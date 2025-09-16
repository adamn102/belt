#!/usr/bin/env python3
"""
BELT Mapping Pipeline

This module provides functionality for mapping FASTQ reads to reference sequences
using minimap2 with filtering, sorting, and indexing.
"""

import os
import subprocess
from pathlib import Path
from utils import log


def filter_fastq(input_fastq, output_fastq, min_length=1000, min_quality=7):
    """
    Filter FASTQ file based on read length and quality
    
    Args:
        input_fastq (str): Path to input FASTQ file
        output_fastq (str): Path to output filtered FASTQ file
        min_length (int): Minimum read length
        min_quality (int): Minimum average quality score
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger = log()
    logger.log_message(f"Filtering FASTQ file: {input_fastq}")
    
    try:
        from Bio import SeqIO
        
        total_reads = 0
        filtered_reads = 0
        
        with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
            for record in SeqIO.parse(infile, 'fastq'):
                total_reads += 1
                
                # Check read length
                if len(record.seq) < min_length:
                    continue
                
                # Check average quality
                if record.letter_annotations['phred_quality']:
                    avg_quality = sum(record.letter_annotations['phred_quality']) / len(record.letter_annotations['phred_quality'])
                    if avg_quality < min_quality:
                        continue
                
                # Write filtered read
                SeqIO.write(record, outfile, 'fastq')
                filtered_reads += 1
        
        logger.log_message(f"Filtered {filtered_reads} reads from {total_reads} total reads")
        return True
        
    except ImportError:
        logger.log_message("BioPython not available, skipping FASTQ filtering")
        # Copy file without filtering
        import shutil
        shutil.copy2(input_fastq, output_fastq)
        return True
        
    except Exception as e:
        logger.log_message(f"Error filtering FASTQ: {e}")
        return False


def run_minimap2_to_sorted(reference, fastq, output_sorted_bam, threads=4):
    """
    Run minimap2 alignment and pipe directly to samtools sort
    
    Args:
        reference (str): Path to reference genome
        fastq (str): Path to input FASTQ file
        output_sorted_bam (str): Path to output sorted BAM file
        threads (int): Number of threads to use
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger = log()
    logger.log_message(f"Running minimap2 alignment and sorting: {fastq} -> {output_sorted_bam}")
    
    try:
        # Build minimap2 command optimized for N-base mapping
        minimap2_cmd = [
            "minimap2",
            "-ax", "map-ont",  # Oxford Nanopore preset
            "-t", str(threads),
            "--score-N=0",  # Don't penalize N bases (important for barcodes)
            "--end-bonus=0",  # Don't penalize end alignments
            "--min-occ-floor=0",  # Don't filter based on occurrence
            "--secondary=no",  # Don't output secondary alignments
            reference, fastq
        ]
        
        samtools_cmd = [
            "samtools", "sort",
            "-o", output_sorted_bam,
            "-"
        ]
        
        # Run minimap2 and pipe to samtools sort
        minimap2_process = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        samtools_process = subprocess.Popen(samtools_cmd, stdin=minimap2_process.stdout, stderr=subprocess.PIPE)
        
        # Close minimap2 stdout to allow proper pipe handling
        minimap2_process.stdout.close()
        
        # Wait for processes to complete
        minimap2_stderr = minimap2_process.communicate()[1]
        samtools_stderr = samtools_process.communicate()[1]
        
        # Check return codes
        if minimap2_process.returncode != 0:
            logger.log_message(f"minimap2 failed: {minimap2_stderr.decode()}")
            return False
        
        if samtools_process.returncode != 0:
            logger.log_message(f"samtools sort failed: {samtools_stderr.decode()}")
            return False
        
        logger.log_message("minimap2 alignment and sorting completed successfully")
        return True
        
    except Exception as e:
        logger.log_message(f"Error running minimap2 and sorting: {e}")
        return False


def index_bam(bam_file):
    """
    Index BAM file using samtools
    
    Args:
        bam_file (str): Path to BAM file to index
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger = log()
    logger.log_message(f"Indexing BAM file: {bam_file}")
    
    try:
        # Run samtools index
        cmd = ["samtools", "index", bam_file]
        
        result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            logger.log_message(f"samtools index failed: {result.stderr}")
            return False
        
        logger.log_message(f"BAM file indexed successfully")
        return True
        
    except Exception as e:
        logger.log_message(f"Error indexing BAM: {e}")
        return False


def load_barcodes_from_csv(csv_file):
    """
    Load barcodes from CSV file
    
    Expected CSV format:
    barcode_5p,barcode_3p,sample_name
    ATCGATCGATCGATC,TTTTTTTTTTTTTT,sample1
    GCATGCATGCATGC,CCCCCCCCCCCCCC,sample2
    
    Returns:
        dict: Dictionary with barcode information
    """
    import csv
    
    barcodes = {
        '5p': [],
        '3p': [],
        'samples': []
    }
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'barcode_5p' in row and row['barcode_5p'].strip():
                    barcodes['5p'].append(row['barcode_5p'].strip())
                if 'barcode_3p' in row and row['barcode_3p'].strip():
                    barcodes['3p'].append(row['barcode_3p'].strip())
                if 'sample_name' in row and row['sample_name'].strip():
                    barcodes['samples'].append(row['sample_name'].strip())
        
        return barcodes
    except Exception as e:
        logger = log()
        logger.log_message(f"Error loading barcodes from CSV: {e}")
        return None


def demultiplex_reads_fuzzy(bam_file, output_dir, barcodes_5p=None, barcodes_3p=None, 
                           sample_names=None, search_window_5p=100, search_window_3p=100, mismatch_threshold=2):
    """
    Demultiplex reads using fuzzy matching within specified windows
    
    Args:
        bam_file (str): Path to input BAM file
        output_dir (str): Path to output directory
        barcodes_5p (list): List of 5' barcode sequences
        barcodes_3p (list): List of 3' barcode sequences
        search_window_5p (int): Number of bases to search from start for 5' barcodes
        search_window_3p (int): Number of bases to search from end for 3' barcodes
        mismatch_threshold (int): Maximum mismatches allowed for barcode matching
    
    Returns:
        dict: Demultiplexing statistics and output files
    """
    logger = log()
    logger.log_message(f"Demultiplexing reads from {bam_file} using fuzzy matching")
    
    try:
        import pysam
        
        # Create demultiplexed output directory
        demux_dir = os.path.join(output_dir, "demultiplexed")
        os.makedirs(demux_dir, exist_ok=True)
        
        # Initialize statistics
        stats = {
            'total_reads': 0,
            'demultiplexed_reads': 0,
            'unassigned_reads': 0,
            'barcode_combinations': {},
            'output_files': {}
        }
        
        # Create output BAM files for each barcode combination
        bam_files = {}
        if barcodes_5p and barcodes_3p:
            # Both 5' and 3' barcodes
            for i, (b5p, b3p) in enumerate(zip(barcodes_5p, barcodes_3p)):
                # Use sample name from CSV if available, otherwise use index
                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                bam_path = os.path.join(demux_dir, f"{sample_name}.bam")
                bam_files[sample_name] = pysam.AlignmentFile(bam_path, "wb", template=pysam.AlignmentFile(bam_file, "rb"))
                stats['barcode_combinations'][sample_name] = 0
                stats['output_files'][sample_name] = bam_path
        elif barcodes_5p:
            # Only 5' barcodes
            for i, b5p in enumerate(barcodes_5p):
                # Use sample name from CSV if available, otherwise use index
                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                bam_path = os.path.join(demux_dir, f"{sample_name}.bam")
                bam_files[sample_name] = pysam.AlignmentFile(bam_path, "wb", template=pysam.AlignmentFile(bam_file, "rb"))
                stats['barcode_combinations'][sample_name] = 0
                stats['output_files'][sample_name] = bam_path
        elif barcodes_3p:
            # Only 3' barcodes
            for i, b3p in enumerate(barcodes_3p):
                # Use sample name from CSV if available, otherwise use index
                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                bam_path = os.path.join(demux_dir, f"{sample_name}.bam")
                bam_files[sample_name] = pysam.AlignmentFile(bam_path, "wb", template=pysam.AlignmentFile(bam_file, "rb"))
                stats['barcode_combinations'][sample_name] = 0
                stats['output_files'][sample_name] = bam_path
        
        # Unassigned reads file
        unassigned_bam = os.path.join(demux_dir, "unassigned.bam")
        unassigned_file = pysam.AlignmentFile(unassigned_bam, "wb", template=pysam.AlignmentFile(bam_file, "rb"))
        stats['output_files']['unassigned'] = unassigned_bam
        
        def calculate_mismatches(seq1, seq2):
            """Calculate number of mismatches between two sequences"""
            if len(seq1) != len(seq2):
                return float('inf')
            return sum(1 for a, b in zip(seq1, seq2) if a != b)
        
        def find_best_barcode_match(sequence, barcodes, search_window, is_3p=False):
            """Find the best matching barcode with minimum mismatches within search window"""
            best_match = None
            min_mismatches = float('inf')
            best_position = -1
            
            # Determine search region
            if is_3p:
                # Search from the end
                search_start = max(0, len(sequence) - search_window)
                search_sequence = sequence[search_start:]
                search_offset = search_start
            else:
                # Search from the beginning
                search_sequence = sequence[:search_window]
                search_offset = 0
            
            for barcode in barcodes:
                # Check all possible positions in the search region
                for i in range(len(search_sequence) - len(barcode) + 1):
                    subsequence = search_sequence[i:i+len(barcode)]
                    mismatches = calculate_mismatches(subsequence, barcode)
                    
                    if mismatches < min_mismatches and mismatches <= mismatch_threshold:
                        min_mismatches = mismatches
                        best_match = barcode
                        best_position = search_offset + i
            
            return best_match if min_mismatches <= mismatch_threshold else None
        
        # Process reads
        with pysam.AlignmentFile(bam_file, "rb") as infile:
            for read in infile:
                stats['total_reads'] += 1
                sequence = read.query_sequence
                
                if not sequence:
                    unassigned_file.write(read)
                    stats['unassigned_reads'] += 1
                    continue
                
                assigned = False
                
                if barcodes_5p and barcodes_3p:
                    # Both 5' and 3' barcodes
                    barcode_5p_match = find_best_barcode_match(sequence, barcodes_5p, search_window_5p, is_3p=False)
                    barcode_3p_match = find_best_barcode_match(sequence, barcodes_3p, search_window_3p, is_3p=True)
                    
                    if barcode_5p_match and barcode_3p_match:
                        # Find the index of the matching combination
                        for i, (b5p, b3p) in enumerate(zip(barcodes_5p, barcodes_3p)):
                            if b5p == barcode_5p_match and b3p == barcode_3p_match:
                                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                                if sample_name in bam_files:
                                    bam_files[sample_name].write(read)
                                    stats['barcode_combinations'][sample_name] += 1
                                    stats['demultiplexed_reads'] += 1
                                    assigned = True
                                    break
                
                elif barcodes_5p:
                    # Only 5' barcodes
                    barcode_match = find_best_barcode_match(sequence, barcodes_5p, search_window_5p, is_3p=False)
                    if barcode_match:
                        for i, b5p in enumerate(barcodes_5p):
                            if b5p == barcode_match:
                                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                                if sample_name in bam_files:
                                    bam_files[sample_name].write(read)
                                    stats['barcode_combinations'][sample_name] += 1
                                    stats['demultiplexed_reads'] += 1
                                    assigned = True
                                    break
                
                elif barcodes_3p:
                    # Only 3' barcodes
                    barcode_match = find_best_barcode_match(sequence, barcodes_3p, search_window_3p, is_3p=True)
                    if barcode_match:
                        for i, b3p in enumerate(barcodes_3p):
                            if b3p == barcode_match:
                                sample_name = sample_names[i] if sample_names and i < len(sample_names) else f"sample_{i+1}"
                                if sample_name in bam_files:
                                    bam_files[sample_name].write(read)
                                    stats['barcode_combinations'][sample_name] += 1
                                    stats['demultiplexed_reads'] += 1
                                    assigned = True
                                    break
                
                if not assigned:
                    unassigned_file.write(read)
                    stats['unassigned_reads'] += 1
        
        # Close all BAM files
        for bam_file_obj in bam_files.values():
            bam_file_obj.close()
        unassigned_file.close()
        
        # Index all demultiplexed BAM files
        for bam_path in stats['output_files'].values():
            if bam_path != unassigned_bam:  # Don't index unassigned file
                index_bam(bam_path)
        
        # Log statistics
        logger.log_message(f"Demultiplexing completed:")
        logger.log_message(f"  Total reads: {stats['total_reads']}")
        logger.log_message(f"  Demultiplexed reads: {stats['demultiplexed_reads']}")
        logger.log_message(f"  Unassigned reads: {stats['unassigned_reads']}")
        logger.log_message(f"  Demultiplexing rate: {stats['demultiplexed_reads']/stats['total_reads']:.2%}")
        
        for combo, count in stats['barcode_combinations'].items():
            logger.log_message(f"  {combo}: {count} reads")
        
        return stats
        
    except Exception as e:
        logger.log_message(f"Error during demultiplexing: {e}")
        return None


def map_pipeline(reference, fastq, output_dir, min_length=1000, min_quality=7, threads=4,
                demultiplex_csv=None, search_window_5p=100, search_window_3p=100, mismatch_threshold=2):
    """
    Complete mapping pipeline: filter -> map+sort -> index -> demultiplex (optional)
    
    Args:
        reference (str): Path to reference genome
        fastq (str): Path to input FASTQ file
        output_dir (str): Path to output directory
        min_length (int): Minimum read length for filtering
        min_quality (int): Minimum average quality score for filtering
        threads (int): Number of threads to use
        demultiplex_csv (str): Path to CSV file with barcode information (optional)
        search_window_5p (int): Number of bases to search from start for 5' barcodes
        search_window_3p (int): Number of bases to search from end for 3' barcodes
        mismatch_threshold (int): Maximum mismatches allowed for barcode matching
    
    Returns:
        dict: Dictionary with output file paths, or None if failed
    """
    logger = log()
    logger.log_message(f"Starting mapping pipeline for {fastq}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate output file paths
    sample_name = Path(fastq).stem
    filtered_fastq = os.path.join(output_dir, f"{sample_name}_filtered.fastq")
    bam_file = os.path.join(output_dir, f"{sample_name}.bam")
    temp_bam = os.path.join(output_dir, f"{sample_name}_temp.bam")
    
    # Step 1: Filter FASTQ
    logger.log_message("Step 1: Filtering FASTQ file")
    if not filter_fastq(fastq, filtered_fastq, min_length, min_quality):
        logger.log_message("FASTQ filtering failed")
        return None
    
    # Step 2: Run minimap2 and sort in one step
    logger.log_message("Step 2: Running minimap2 alignment and sorting")
    if not run_minimap2_to_sorted(reference, filtered_fastq, temp_bam, threads):
        logger.log_message("minimap2 alignment and sorting failed")
        return None
    
    # Step 3: Move temp file to final location
    logger.log_message("Step 3: Finalizing BAM file")
    os.rename(temp_bam, bam_file)
    
    # Step 4: Index BAM
    logger.log_message("Step 4: Indexing BAM file")
    if not index_bam(bam_file):
        logger.log_message("BAM indexing failed")
        return None
    
    # Step 5: Demultiplex (if CSV provided)
    demux_stats = None
    if demultiplex_csv:
        logger.log_message("Step 5: Loading barcodes and demultiplexing reads")
        barcodes = load_barcodes_from_csv(demultiplex_csv)
        if barcodes:
            demux_stats = demultiplex_reads_fuzzy(
                bam_file, output_dir, 
                barcodes_5p=barcodes['5p'] if barcodes['5p'] else None,
                barcodes_3p=barcodes['3p'] if barcodes['3p'] else None,
                sample_names=barcodes['samples'] if barcodes['samples'] else None,
                search_window_5p=search_window_5p,
                search_window_3p=search_window_3p,
                mismatch_threshold=mismatch_threshold
            )
            if not demux_stats:
                logger.log_message("Demultiplexing failed")
                return None
        else:
            logger.log_message("Failed to load barcodes from CSV")
            return None
    
    # Return output file paths
    output_files = {
        'filtered_fastq': filtered_fastq,
        'bam_file': bam_file,
        'bam_index': bam_file + '.bai'
    }
    
    if demux_stats:
        output_files['demultiplexed_files'] = demux_stats['output_files']
        output_files['demultiplexing_stats'] = demux_stats
    
    logger.log_message(f"Mapping pipeline completed successfully")
    logger.log_message(f"Output files:")
    for key, path in output_files.items():
        if key not in ['demultiplexed_files', 'demultiplexing_stats']:
            logger.log_message(f"  {key}: {path}")
    
    if demux_stats:
        logger.log_message(f"  Demultiplexed files: {len(demux_stats['output_files'])} files in demultiplexed/ directory")
    
    return output_files


def main():
    """Main function for standalone usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="BELT Mapping Pipeline")
    parser.add_argument("--reference", required=True, help="Reference genome fasta file")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--min-length", type=int, default=1000, help="Minimum read length for filtering")
    parser.add_argument("--min-quality", type=int, default=7, help="Minimum average quality score for filtering")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    
    args = parser.parse_args()
    
    result = map_pipeline(
        reference=args.reference,
        fastq=args.fastq,
        output_dir=args.output_dir,
        min_length=args.min_length,
        min_quality=args.min_quality,
        threads=args.threads
    )
    
    if result:
        print("Mapping pipeline completed successfully!")
        for key, path in result.items():
            print(f"{key}: {path}")
    else:
        print("Mapping pipeline failed!")
        exit(1)


if __name__ == "__main__":
    main()
