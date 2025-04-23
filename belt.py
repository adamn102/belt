import click
from extract import extract_reads
import yaml

@click.group()
def main():
    """BELT - BalazsLab ERRSeq Long-Read Toolkit"""
    pass

@main.command()
@click.argument('input_fastq')
@click.option('--output', required=False, help="Output FASTQ file")
@click.option('--min_read_length', default=None, type=int, help="Minimum expected read length")
@click.option('--motif5', required=False, help="5' motif sequence")
@click.option('--motif3', required=False, help="3' motif sequence")
@click.option('--fuzzy', default=None, type=int, help="Number of fuzzy mismatches allowed")
@click.option('--window', default=None, type=int, help="Motif search window size")
@click.option('--config', type=click.Path(exists=True), help="YAML config file")
def extract(input_fastq, output, min_read_length, motif5, motif3, fuzzy, window, config):
    """Extract reads based on motifs and UMIs"""

    # Load config if provided
    cfg = {}
    if config:
        cfg = load_config(config)

    # Merge CLI overrides
    def get_param(name, cli_val):
        return cli_val if cli_val is not None else cfg.get(name)

    final_config = {
        'input_fastq': input_fastq,
        'output': get_param('output', output),
        'min_read_length': get_param('min_read_length', min_read_length) or 3900,
        'motif5': get_param('motif5', motif5),
        'motif3': get_param('motif3', motif3),
        'fuzzy': get_param('fuzzy', fuzzy) or 2,
        'window': get_param('window', window) or 50,
    }
    # Sanity check
    required_keys = ['output', 'motif5', 'motif3']
    for key in required_keys:
        if final_config[key] is None:
            raise click.UsageError(f"Missing required parameter: '{key}' (not set in CLI or config)")

    # Call extraction logic
    print("Running extract_reads with:")
    print(f"Input: {final_config['input_fastq']}")
    print(f"Output: {final_config['output']}")
    print(f"min_read_length: {final_config['min_read_length']}")
    print(f"motif5: {final_config['motif5']}")
    print(f"motif3: {final_config['motif3']}")
    print(f"fuzzy: {final_config['fuzzy']}")
    print(f"window: {final_config['window']}")
    
    extract_reads(
    final_config['input_fastq'],
    final_config['min_read_length'],
    final_config['output'],
    final_config['motif5'],
    final_config['motif3'],
    final_config['fuzzy'],
    final_config['window']
)



# Additional commands will follow this pattern


from group_umis import group

@main.command()
@click.argument('input_bam')
@click.option('--output_dir', required=True, help="Directory to store UMI-specific BAM files")
def group(input_bam, output_dir):
    """Group reads by 5' UMI sequence and output UMI-specific BAM files."""
    group(input_bam, output_dir)


def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f) or {}
    
    
if __name__ == "__main__":
    main()