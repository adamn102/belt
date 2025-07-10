import click
from extract import extract_UMIs
from whitelist import umi_whitelist
from umi_map import umi_extraction_map
import yaml


@click.group()
def main():
    """BELT - BalazsLab ERRSeq Long-Read Toolkit"""
    pass

@main.command()
@click.option('--reference', required=True, help="Reference genome fasta file")
@click.option('--alignment',required=True, help="Input bam file")
@click.option('--output', required=True, help="Output bam file")
@click.option('--min_read_length', default=3900, type=int, help="Minimum expected read length")
@click.option('--motif', required=False, default = "NNNNNNNNNNNN", help="barcode motif sequence")
def extract(reference, alignment, min_read_length ,motif, output):
    """Extract reads based on motifs and UMIs"""
    extract_UMIs(reference, alignment, min_read_length, motif, output)

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f) or {}
if __name__ == "__main__": 
    main()