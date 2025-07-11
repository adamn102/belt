import pysam
import pandas as pd
from belt.scripts.utils import log
from collections import Counter
import os

def long_read_umi_grouping(reference, alignedReads, correctedUMIs, outputDir, grouping="5p", threshold=5):
    """
    Generate grouped FASTQ files from long reads using corrected UMIs.

    :param reference: Path to the reference genome.
    :param alignedReads: Path to the BAM file with aligned reads.
    :param correctedUMIs: Path to the CSV file with corrected UMIs.
    :param outputDir: Directory where the output files will be saved.
    """
    os.makedirs(outputDir, exist_ok=True)

    sampleName = os.path.basename(alignedReads).split(".")[0]
    logger = log()
    logger.log_message("Starting long read UMI grouping...")

    # Read in corrected UMIs CSV
    df = pd.read_csv(correctedUMIs)

    # Clean up 5pUMI
    df["5pUMI"] = df["5pUMI"].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    # Group by 5pUMI
    grouped = df.groupby("5pUMI").agg({
        "Read Name": ["count", lambda x: list(x)]
    })

    grouped.columns = ["ReadCount", "ReadNames"]

    # Filter groups by threshold
    filtered = grouped[grouped["ReadCount"] >= threshold]

    filtered_csv = os.path.join(outputDir, f"{sampleName}_filtered_umi_groups.csv")
    filtered.reset_index().to_csv(filtered_csv, index=False)
    logger.log_message(f"Filtered UMI groups saved to {filtered_csv}")

    # Create output dir per sample
    sample_output_dir = os.path.join(outputDir, sampleName)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Read name -> UMI lookup
    read_to_umis = df.set_index("Read Name")[["5pUMI", "3pUMI"]].to_dict(orient="index")

    # Iterate over UMI groups
    for umi, row in filtered.iterrows():
        out_fastq = os.path.join(sample_output_dir, f"{umi}_reads.fastq")

        extract_reads_to_fastq(
            read_names=row["ReadNames"],
            bam_path=alignedReads,
            output_fastq_path=out_fastq,
            read_umi_lookup=read_to_umis
        )

    logger.log_message(f"Grouped FASTQ files written to {sample_output_dir}")

def extract_reads_to_fastq(
        read_names,
        bam_path,
        output_fastq_path,
        read_umi_lookup
):
    """
    Extract reads from BAM and write FASTQ with UMI annotations.

    Parameters:
        read_names (Iterable[str]): Read names to extract.
        bam_path (str): Path to BAM file.
        output_fastq_path (str): Output FASTQ file path.
        read_umi_lookup (dict): Read name -> UMI info.
    """
    read_names_set = set(read_names)
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    n_written = 0

    with open(output_fastq_path, "w") as fastq_out:
        for read in bamfile.fetch(until_eof=True):
            if read.query_name in read_names_set and not read.is_secondary and not read.is_supplementary:
                seq = read.query_sequence
                qual = read.qual

                if seq and qual:
                    umi_info = read_umi_lookup.get(read.query_name, {})
                    umi_5p = umi_info.get("5pUMI", "NA")
                    umi_3p = umi_info.get("3pUMI", "NA")

                    # FASTQ format
                    header = f"{read.query_name}|5pUMI:{umi_5p}|3pUMI:{umi_3p}"
                    fastq_out.write(f"@{header}\n{seq}\n+\n{qual}\n")
                    n_written += 1

    bamfile.close()

    if n_written == 0:
        print(f"Warning: No reads written for {output_fastq_path}")
