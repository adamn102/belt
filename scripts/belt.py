import click

import sys


import os

# Add the current directory to the path so we can import local modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from extract import extract_UMIs
    from whitelist import umi_whitelist
    from umi_map import umi_extraction_map
    from umi_group import long_read_umi_grouping_fast
    from utils import *
except ImportError as e:
    # If imports fail, create mock functions for CLI testing
    print(f"Warning: Could not import some modules: {e}")
    print("Creating mock functions for CLI testing...")
    
    def extract_UMIs(*args, **kwargs):
        print("Mock extract_UMIs called with:", args, kwargs)
    
    def umi_whitelist(*args, **kwargs):
        print("Mock umi_whitelist called with:", args, kwargs)
    
    def umi_extraction_map(*args, **kwargs):
        print("Mock umi_extraction_map called with:", args, kwargs)
    
    def long_read_umi_grouping_fast(*args, **kwargs):
        print("Mock long_read_umi_grouping_fast called with:", args, kwargs)

import yaml


@click.group()
def main():
    """BELT - BalazsLab ERRSeq Long-Read Toolkit"""
    pass


@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--output', required=True, help="Output BAM file")
@click.option('--min_read_length', default=3900, type=int, help="Minimum expected read length")
@click.option('--motif', required=False, default="NNNNNNNNNNNN", help="Barcode motif sequence")
def extract(reference, alignment, min_read_length, motif, output):
    """Extract reads based on motifs and UMIs"""
    extract_UMIs(reference, alignment, min_read_length, motif, output)


@main.command()
@click.option('--whitelist', required=True, help="Path to the whitelist file")
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--min_read_length', default=2000, type=int, help="Minimum expected read length")
@click.option('--motif', required=True, help="Motif for UMI extraction")
@click.option('--output_dir', required=True, help="Output directory")
@click.option('--analysis_out', required=False, default="analysis", help="Analysis output directory")
def umi_map(reference, alignment, min_read_length, motif, whitelist, output_dir, analysis_out):
    """Map UMIs to a whitelist and generate paired UMI output"""
    umi_extraction_map(reference, alignment, min_read_length, motif, whitelist, output_dir, analysis_out)


@main.command(name='umi-group')
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--corrected_umis', required=True, help="Path to the corrected UMIs file")
@click.option('--output_dir', required=True, help="Output directory")
@click.option('--grouping', required=True, type=click.Choice(['5pUMI', '3pUMI']), help="Grouping method (5pUMI or 3pUMI)")
@click.option('--threshold', default=10, type=int, help="Threshold for grouping")
def umi_group(reference, alignment, corrected_umis, output_dir, grouping, threshold):
    """Group long-read UMIs based on corrected UMI data"""
    long_read_umi_grouping_fast(reference, alignment, corrected_umis, output_dir, grouping, threshold)


@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment', required=True, help="Input BAM file")
@click.option('--min_read_length', default=250, type=int, help="Minimum expected read length")
@click.option('--motif', required=True, help="Motif for UMI extraction")
@click.option('--cluster_method', default="default", help="Clustering method for whitelist generation")
@click.option('--output_dir', required=True, help="Output directory")
def whitelist(reference, alignment, min_read_length, motif, cluster_method, output_dir):
    """Generate a whitelist of UMIs"""
    umi_whitelist(reference, alignment, min_read_length, motif, cluster_method, output_dir)

## TODO Compare long read data to short read whitelist - return the % overlap with whitelist

## TODO Compare long read data to long read whitelist - return the concordence rate


@main.command()
@click.option('--config', required=True, help="Path to the YAML configuration file")
def run_pipeline(config):
    """Run the entire pipeline using a configuration file"""
    config_data = load_config(config)
    extract_UMIs(
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