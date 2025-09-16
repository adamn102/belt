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
        # Build minimap2 command - simple and basic
        minimap2_cmd = [
            "minimap2",
            "-ax", "map-ont",  # Oxford Nanopore preset
            "-t", str(threads),
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


def map_pipeline(reference, fastq, output_dir, min_length=1000, min_quality=7, threads=4):
    """
    Complete mapping pipeline: filter -> map+sort -> index
    
    Args:
        reference (str): Path to reference genome
        fastq (str): Path to input FASTQ file
        output_dir (str): Path to output directory
        min_length (int): Minimum read length for filtering
        min_quality (int): Minimum average quality score for filtering
        threads (int): Number of threads to use
    
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
    
    # Return output file paths
    output_files = {
        'filtered_fastq': filtered_fastq,
        'bam_file': bam_file,
        'bam_index': bam_file + '.bai'
    }
    
    logger.log_message(f"Mapping pipeline completed successfully")
    logger.log_message(f"Output files:")
    for key, path in output_files.items():
        logger.log_message(f"  {key}: {path}")
    
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