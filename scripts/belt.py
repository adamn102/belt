import click
import sys
import os

# Add the current directory to the path so we can import local modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Also add the parent directory in case we're running from the top level
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from extract import extract_umis
from whitelist import umi_whitelist
from umi_map import umi_extraction_map
from umi_group import long_read_umi_grouping_fast
from consensus import process_grouped_bam
from map import map_pipeline
from utils import *

import yaml


@click.group()
def main():
    """BELT - BalazsLab ERRSeq Long-Read Toolkit"""
    pass


@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file(s) or directory - can use glob patterns")
@click.option('--output', required=True, help="Output directory")
@click.option('--min_read_length', default=3900, type=int, help="Minimum expected read length")
@click.option('--motif', required=False, default=["NNNNNNNNNNNN"], help="Barcode motif sequence")
def extract(reference, alignment, min_read_length, motif, output):
    """Extract reads based on motifs and UMIs from aligned BAM files"""
    import glob
    from pathlib import Path
    
    # Parse motif if it's a string representation of a list
    if isinstance(motif, str):
        import ast
        try:
            motif = ast.literal_eval(motif)
        except (ValueError, SyntaxError):
            # If parsing fails, treat as single motif
            motif = [motif]
    elif not isinstance(motif, list):
        motif = [motif]
    
    # Expand glob patterns and find BAM files
    expanded_bam_files = []
    
    # Check if alignment is a directory
    if os.path.isdir(alignment):
        # Find all BAM files in the directory
        bam_pattern = os.path.join(alignment, "*.bam")
        matched_files = glob.glob(bam_pattern)
        if not matched_files:
            print(f"Warning: No BAM files found in directory: {alignment}")
        else:
            expanded_bam_files.extend(matched_files)
    else:
        # Check if the pattern contains glob characters
        if any(char in alignment for char in ['*', '?', '[', ']']):
            # Expand glob pattern
            matched_files = glob.glob(alignment)
            if not matched_files:
                print(f"Warning: No files found matching pattern: {alignment}")
            else:
                # Filter to only BAM files
                bam_files = [f for f in matched_files if f.endswith('.bam')]
                if not bam_files:
                    print(f"Warning: No BAM files found matching pattern: {alignment}")
                else:
                    expanded_bam_files.extend(bam_files)
        else:
            # Single file - check if it exists and is a BAM file
            if os.path.exists(alignment):
                if alignment.endswith('.bam'):
                    expanded_bam_files.append(alignment)
                else:
                    print(f"Warning: File {alignment} is not a BAM file")
            else:
                print(f"Warning: File not found: {alignment}")
    
    if not expanded_bam_files:
        print("❌ No valid BAM files found!")
        exit(1)
    
    print(f"Found {len(expanded_bam_files)} BAM file(s) to process:")
    for f in expanded_bam_files:
        print(f"  - {f}")
    
    # Process each BAM file
    for bam_file in expanded_bam_files:
        print(f"Processing BAM file: {bam_file}")
        # Generate output path for this BAM file - output CSV directly to output directory
        sample_name = Path(bam_file).stem
        bam_output = os.path.join(output, f"{sample_name}_extracted.csv")
        
        extract_umis(reference, bam_file, min_read_length, motif, bam_output)


@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--fastq', required=True, multiple=True, help="Input FASTQ file(s) - can be specified multiple times")
@click.option('--output-dir', required=True, help="Output directory")
@click.option('--min-length', default=1000, type=int, help="Minimum read length for filtering")
@click.option('--min-quality', default=7, type=int, help="Minimum average quality score for filtering")
@click.option('--threads', default=4, type=int, help="Number of threads to use")
@click.option('--config', help="Config YAML file for barcode demultiplexing (optional)")
@click.option('--sample-sheet', help="Sample sheet CSV file for demultiplexing (optional)")
def map(reference, fastq, output_dir, min_length, min_quality, threads, config, sample_sheet):
    """Map FASTQ reads to reference using minimap2 with filtering, sorting, indexing, and optional demultiplexing"""
    import glob
    from pathlib import Path
    
    # Expand glob patterns in FASTQ file arguments
    expanded_fastq_files = []
    for fastq_pattern in fastq:
        # Check if the pattern contains glob characters
        if any(char in fastq_pattern for char in ['*', '?', '[', ']']):
            # Expand glob pattern
            matched_files = glob.glob(fastq_pattern)
            if not matched_files:
                print(f"Warning: No files found matching pattern: {fastq_pattern}")
            else:
                expanded_fastq_files.extend(matched_files)
        else:
            # No glob pattern, use as-is
            if os.path.exists(fastq_pattern):
                expanded_fastq_files.append(fastq_pattern)
            else:
                print(f"Warning: File not found: {fastq_pattern}")
    
    if not expanded_fastq_files:
        print("❌ No valid FASTQ files found!")
        exit(1)
    
    print(f"Found {len(expanded_fastq_files)} FASTQ file(s) to process:")
    for f in expanded_fastq_files:
        print(f"  - {f}")
    
    # Process each FASTQ file
    all_results = []
    
    for fastq_file in expanded_fastq_files:
        print(f"Processing FASTQ file: {fastq_file}")
        result = map_pipeline(
            reference=reference,
            fastq=fastq_file,
            output_dir=output_dir,
            min_length=min_length,
            min_quality=min_quality,
            threads=threads,
            config_file=config,
            sample_sheet=sample_sheet
        )
        all_results.append(result)
    
    # Check if all results were successful
    successful_results = [r for r in all_results if r is not None]
    
    if successful_results:
        print("✅ Mapping pipeline completed successfully!")
        print(f"Processed {len(successful_results)} out of {len(all_results)} FASTQ files")
        
        # Show statistics for each successful result
        for i, result in enumerate(successful_results):
            if result:
                print(f"\n📊 Results for FASTQ file {i+1}:")
                
                # Show filtering statistics if available
                if 'filtering_stats' in result:
                    filter_stats = result['filtering_stats']
                    print(f"  Filtering Statistics:")
                    print(f"    Original reads: {filter_stats['total_reads']}")
                    print(f"    Filtered reads: {filter_stats['filtered_reads']}")
                    print(f"    Length filtered: {filter_stats['length_filtered']}")
                    print(f"    Quality filtered: {filter_stats['quality_filtered']}")
                    print(f"    Filter rate: {filter_stats['filter_rate']:.2%}")
                
                # Show demultiplexing statistics if available
                if 'demultiplexing_stats' in result:
                    stats = result['demultiplexing_stats']
                    # Calculate corrected statistics
                    original_total = result['filtering_stats']['total_reads'] if 'filtering_stats' in result else stats['total_reads']
                    filtered_total = result['filtering_stats']['filtered_reads'] if 'filtering_stats' in result else stats['total_reads']
                    
                    print(f"  Demultiplexing Statistics:")
                    print(f"    Original reads: {original_total}")
                    print(f"    Filtered reads: {filtered_total}")
                    print(f"    Demultiplexed reads: {stats['demultiplexed_reads']}")
                    print(f"    Unassigned reads: {stats['unassigned_reads']}")
                    print(f"    Demultiplexing rate (of filtered): {stats['demultiplexed_reads']/stats['total_reads']:.2%}")
                    print(f"    Overall success rate: {stats['demultiplexed_reads']/original_total:.2%}")
                    print(f"  Demultiplexed files:")
                    # Use the new sample_counts format
                    sample_counts = stats.get('sample_counts', {})
                    for sample, count in sample_counts.items():
                        print(f"    {sample}: {count} reads")
                    print(f"    unassigned: {stats['unassigned_reads']} reads")
    else:
        print("❌ Mapping pipeline failed for all FASTQ files!")
        exit(1)


@main.command()
@click.option('--whitelist', required=True, help="Path to the whitelist file")
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--min_read_length', default=2000, type=int, help="Minimum expected read length")
@click.option('--motif', required=False, default=["NNNNNNNNNNNN"], help="Barcode motif sequence(s) - can be a single motif or list of motifs")
@click.option('--output_dir', required=True, help="Output directory")
@click.option('--max_distance', required=False, default=1, type=int, help="Maximum hamming distance for whitelist mapping")
def umi_map(reference, alignment, min_read_length, motif, whitelist, output_dir, max_distance):
    """Map UMIs to a whitelist and generate paired UMI output from BAM files"""
    # Parse motif if it's a string representation of a list
    if isinstance(motif, str):
        import ast
        try:
            motif = ast.literal_eval(motif)
        except (ValueError, SyntaxError):
            # If parsing fails, treat as single motif
            motif = [motif]
    elif not isinstance(motif, list):
        motif = [motif]
    
    umi_extraction_map(reference, alignment, min_read_length, motif, whitelist, output_dir, max_distance)


@main.command(name='umi-group')
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--corrected_umis', required=True, help="Path to the corrected UMIs file")
@click.option('--output_dir', required=True, help="Output directory")
@click.option('--grouping', required=True, help="UMI column to group by (e.g., UMI1, UMI2, 5pUMI, 3pUMI)")
@click.option('--threshold', default=10, type=int, help="Threshold for grouping")


def umi_group(reference, alignment, corrected_umis, output_dir, grouping, threshold):
    """Group long-read UMIs based on corrected UMI data from BAM files"""
    long_read_umi_grouping_fast(reference, alignment, corrected_umis, output_dir, grouping, threshold)


@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--min_read_length', default=250, type=int, help="Minimum expected read length")
@click.option('--motif', required=False, default=["NNNNNNNNNNNN"], help="Barcode motif sequence(s) - can be a single motif or list of motifs")
@click.option('--cluster_method', default="default", help="Clustering method for whitelist generation")
@click.option('--output_dir', required=True, help="Output directory")
def whitelist(reference, alignment, min_read_length, motif, cluster_method, output_dir):
    """Generate a whitelist of UMIs"""
    # Parse motif if it's a string representation of a list
    if isinstance(motif, str):
        import ast
        try:
            motif = ast.literal_eval(motif)
        except (ValueError, SyntaxError):
            # If parsing fails, treat as single motif
            motif = [motif]
    elif not isinstance(motif, list):
        motif = [motif]
    
    umi_whitelist(reference, alignment, min_read_length, motif, cluster_method, output_dir)


@main.command()
@click.option('--grouped-bam', required=True, help="Input grouped BAM file with UMI tags")
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--output-dir', required=True, help="Output directory")
@click.option('--sample', required=True, help="Sample name")
@click.option('--tag', default="ZU", help="UMI tag name")
@click.option('--threads', default=2, type=int, help="Total threads")
@click.option('--procs', default=1, type=int, help="Parallel UMI processes")
@click.option('--umi-threads', default=1, type=int, help="Threads per UMI for consensus calling")
@click.option('--model', default="r1041_e82_400bps_sup_variant_v5.0.0", help="Medaka model")
@click.option('--manifest', help="Output manifest file")
@click.option('--log-file', help="Log file path")
def consensus(grouped_bam, reference, output_dir, sample, tag, threads, procs, umi_threads, model, manifest, log_file):
    """Generate consensus sequences from grouped BAM file"""
    produced = process_grouped_bam(
        grouped_bam=grouped_bam,
        reference=reference,
        output_dir=output_dir,
        sample=sample,
        tag=tag,
        threads=threads,
        procs=procs,
        umi_threads=umi_threads,
        model=model,
        manifest=manifest,
        log_file=log_file
    )
    for path in produced:
        print(path)

## TODO Compare long read data to short read whitelist - return the % overlap with whitelist

## TODO Compare long read data to long read whitelist - return the concordence rate


@main.command()
@click.option('--config', required=True, help="Path to the YAML configuration file")
def run_pipeline(config):
    """Run the entire pipeline using a configuration file"""
    config_data = load_config(config)
    extract_umis(
        config_data['reference'],
        config_data['alignment'],
        config_data['min_read_length'],
        config_data['motif'],
        config_data['output']
    )
    umi_extraction_map(
        config_data['reference'],
        config_data['alignment'],
        config_data['min_read_length'],
        config_data['motif'],
        config_data['whitelist'],
        config_data['output_dir']
    )
    long_read_umi_grouping_fast(
        config_data['reference'],
        config_data['alignment'],
        config_data['corrected_umis'],
        config_data['output_dir'],
        config_data['grouping'],
        config_data['threshold']
    )


def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f) or {}


if __name__ == "__main__":
    main()