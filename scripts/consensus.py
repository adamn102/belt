import pysam
import pandas as pd
from belt.scripts.utils import log
from collections import Counter
import os

def long_read_consensus(reference, alignedReads, correctedUMIs, outputDir, grouping="5p", threshold=5):
    """
    Generate consensus sequences from long reads using corrected UMIs.

    :param reference: Path to the reference genome.
    :param alignedReads: Path to the BAM file with aligned reads.
    :param correctedUMIs: Path to the CSV file with corrected UMIs.
    :param outputDir: Directory where the output files will be saved.
    """
    # Create output directory if not exists
    os.makedirs(outputDir, exist_ok=True)

    sampleName = os.path.basename(alignedReads).split(".")[0]
    logger = log()
    logger.log_message("Starting long read consensus generation...")

    # Read in correctedUMIs CSV file
    df = pd.read_csv(correctedUMIs)
    print(df.head())

    # Replace empty or 'NoMatch' 5pUMI with NaN to exclude them
    df["5pUMI"] = df["5pUMI"].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    # Group by 5pUMI
    grouped = df.groupby("5pUMI").agg({
        "Read Name": ["count", lambda x: list(x)]
    })

    # Flatten multi-index columns
    grouped.columns = ["ReadCount", "ReadNames"]

    # Filter groups with >= threshold reads
    filtered = grouped[grouped["ReadCount"] >= threshold]

    # Save filtered groups
    filtered_csv = os.path.join(outputDir, f"{sampleName}_filtered_umi_groups.csv")
    filtered.reset_index().to_csv(filtered_csv, index=False)
    logger.log_message(f"Filtered UMI groups saved to {filtered_csv}")

    # Prepare path for combined consensus FASTA
    consensus_fasta_path = os.path.join(outputDir, f"{sampleName}_all_consensus.fasta")

    # Create directory for sample if it doesn't exist outputDir/sampleName
    sample_output_dir = os.path.join(outputDir, f"{sampleName}")
    os.makedirs(sample_output_dir, exist_ok=True)

    # Build a lookup dict: read name -> (5pUMI, 3pUMI)
    read_to_umis = df.set_index("Read Name")[["5pUMI", "3pUMI"]].to_dict(orient="index")

    with open(consensus_fasta_path, "w") as consensus_fasta:
        # Iterate through UMI groups and write consensus sequences
        for umi, row in filtered.iterrows():
            out_fasta = os.path.join(sample_output_dir, f"{umi}_reads.fasta")

            consensus_sequence = extract_reads_and_call_consensus(
                read_names=row["ReadNames"],
                bam_path=alignedReads,
                output_fasta_path=out_fasta,
                read_umi_lookup=read_to_umis
            )

            if consensus_sequence:
                # Write consensus to combined FASTA
                consensus_fasta.write(f">{sampleName}|{umi}\n{consensus_sequence}\n")

    logger.log_message(f"All consensus sequences written to {consensus_fasta_path}")


def extract_reads_and_call_consensus(
        read_names,
        bam_path,
        output_fasta_path,
        read_umi_lookup
):
    """
    Extract reads by read name from BAM, write to FASTA with UMI annotations, and return consensus sequence.

    Parameters:
        read_names (Iterable[str]): Read names to extract.
        bam_path (str): Path to BAM file.
        output_fasta_path (str): Path to write FASTA of individual reads.
        read_umi_lookup (dict): Mapping from read name -> dict with '5pUMI' and '3pUMI'.

    Returns:
        str: consensus sequence, or None if no sequences were found.
    """
    read_names_set = set(read_names)
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    sequences = []
    with open(output_fasta_path, "w") as fasta_out:
        for read in bamfile.fetch(until_eof=True):
            if read.query_name in read_names_set:
                seq = read.query_sequence
                if seq:
                    sequences.append(seq)

                    # Look up UMIs
                    umi_info = read_umi_lookup.get(read.query_name, {})
                    umi_5p = umi_info.get("5pUMI", "NA")
                    umi_3p = umi_info.get("3pUMI", "NA")

                    # Write FASTA header with UMIs
                    fasta_out.write(f">{read.query_name}|5pUMI:{umi_5p}|3pUMI:{umi_3p}\n{seq}\n")

    bamfile.close()

    # Skip consensus if no sequences
    if not sequences:
        print(f"No sequences found for {output_fasta_path}. Skipping consensus.")
        return None

    # Trim all sequences to the shortest length
    min_len = min(len(s) for s in sequences)
    trimmed_sequences = [s[:min_len] for s in sequences]

    # Build consensus: simple majority per column
    consensus = ""
    for i in range(min_len):
        column_bases = [s[i] for s in trimmed_sequences]
        base_counts = Counter(column_bases)
        most_common_base, _ = base_counts.most_common(1)[0]
        consensus += most_common_base

    return consensus
