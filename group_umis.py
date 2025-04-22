import click
import pysam
from collections import defaultdict
from utils import log

@click.command()
@click.argument('input_bam')
@click.option('--output_dir', required=True, help="Directory to store UMI-specific BAM files")
def cli(input_bam, output_dir):
    """Extract UMIs from reads based on 5' and 3' motifs."""
    group(input_bam, output_dir)


def group(input_bam, output_dir):
    """Group reads by 5' UMI sequence and output UMI-specific BAM files."""
    log(f"Reading input BAM file: {input_bam}")
    umi_groups = defaultdict(list)

    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as bam_file:
        for read in bam_file:
            # Parse the 5' UMI from the read name
            if ":UMI5-" in read.query_name:
                umi_5p = read.query_name.split(":UMI5-")[1].split("-UMI3")[0]
                umi_groups[umi_5p].append(read)

        log(f"Found {len(umi_groups)} unique 5' UMIs. Writing grouped BAM files...")

        # Write UMI-specific BAM files
        for umi, reads in umi_groups.items():
            output_bam_path = f"{output_dir}/{umi}.bam"
            with pysam.AlignmentFile(output_bam_path, "wb", header=bam_file.header) as out_bam:
                for read in reads:
                    out_bam.write(read)

    log(f"UMI-specific BAM files written to {output_dir}")