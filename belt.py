import click
from extract import extract_reads

@click.group()
def main():
    """BEBT - BalazsLab ERRSeq Bioinformatics Toolkit"""
    pass

@main.command()
@click.argument('input_fastq')
@click.option('--output', required=True, help="Output FASTQ file")
@click.option('--motif5', required=True, help="5' motif sequence")
@click.option('--motif3', required=True, help="3' motif sequence")
@click.option('--fuzzy', default=2, help="Number of fuzzy mismatches allowed")
@click.option('--window', default=50, help="Motif search window size")
def extract(input_fastq, output, motif5, motif3, fuzzy, window):
    """Extract reads based on motifs and UMIs"""
    extract_reads(input_fastq, output, motif5, motif3, fuzzy, window)

# Additional commands will follow this pattern


from group_umis import group

@main.command()
@click.argument('input_bam')
@click.option('--output_dir', required=True, help="Directory to store UMI-specific BAM files")
def group(input_bam, output_dir):
    """Group reads by 5' UMI sequence and output UMI-specific BAM files."""
    group(input_bam, output_dir)

if __name__ == "__main__":
    main()