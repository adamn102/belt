import pysam
import pandas as pd
import sys
import os
from collections import Counter

# Add the current directory to the path so we can import local modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from utils import log
except ImportError:
    # Create a simple mock log function if utils is not available
    class MockLog:
        def log_message(self, message):
            print(f"[LOG] {message}")
    
    log = MockLog


def long_read_umi_grouping(reference, alignedReads, correctedUMIs, outputDir, grouping="5pUMI", threshold=5):
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
    # df = pd.read_csv(correctedUMIs)
    df = pd.read_csv(
        correctedUMIs,
        dtype={"Read Name": "string", "5pUMI": "string", "3pUMI": "string"},
        low_memory=False,
        keep_default_na=False  # so "NA" stays "NA" instead of NaN, optional
    )

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


import os, subprocess, pysam, pandas as pd

def long_read_umi_grouping_fast(
        reference, alignedReads, correctedUMIs, outputDir,
        grouping="5pUMI", threshold=5, threads=8
):
    os.makedirs(outputDir, exist_ok=True)
    sample = os.path.basename(alignedReads).split(".")[0]
    logger = log()
    logger.log_message("UMI grouping via tag+sort...")

    # 1) Load just what we need
    df = pd.read_csv(
        correctedUMIs, usecols=["Read Name", "5pUMI", "3pUMI"],
        dtype={"Read Name": "string", "5pUMI": "string", "3pUMI": "string"},
        low_memory=False, keep_default_na=False
    )
    if grouping not in ("5pUMI", "3pUMI"):
        raise ValueError("grouping must be '5pUMI' or '3pUMI'")
    gcol = grouping
    df[gcol] = df[gcol].replace(["", "NoMatch", "NoBarcode"], pd.NA)

    counts = df[gcol].value_counts(dropna=True)
    keepers = counts[counts >= threshold]
    filt_csv = os.path.join(outputDir, f"{sample}_filtered_umi_groups.csv")
    keepers.rename("ReadCount").rename_axis(gcol).reset_index().to_csv(filt_csv, index=False)
    if keepers.empty:
        logger.log_message("No UMIs pass threshold; done.");
        return

    # map read name -> UMI for kept groups
    name_to_umi = dict(df.loc[df[gcol].isin(keepers.index), ["Read Name", gcol]].itertuples(index=False, name=None))
    del df

    # 2) Tag BAM with ZU:<string> = UMI (only for kept reads)
    tagged_bam = os.path.join(outputDir, f"{sample}.tagged.bam")
    with pysam.AlignmentFile(alignedReads, "rb") as bam_in, \
            pysam.AlignmentFile(tagged_bam, "wb", template=bam_in) as bam_out:
        for r in bam_in.fetch(until_eof=True):
            umi = name_to_umi.get(r.query_name)
            if umi is None:
                continue
            r.set_tag("ZU", str(umi), value_type='Z')
            bam_out.write(r)
    # index optional, but sort will read sequentially
    # 3) Sort by tag
    grouped_bam = os.path.join(outputDir, f"{sample}.grouped.bam")
    subprocess.run(
        ["samtools", "sort", "-@",
         str(threads), "-t", "ZU", "-o", grouped_bam, tagged_bam],
        check=True
    )
    # optional: os.remove(tagged_bam)
    # remove tagged bam
    os.remove(tagged_bam)

    #
    # # 4) Stream grouped.bam; write one UMI at a time
    # out_dir = os.path.join(outputDir, sample)
    # os.makedirs(out_dir, exist_ok=True)
    # curr_umi, fh = None, None
    # written = 0
    # with pysam.AlignmentFile(grouped_bam, "rb") as bam:
    #     for r in bam.fetch(until_eof=True):
    #         umi = r.get_tag("ZU")
    #         if umi != curr_umi:
    #             if fh: fh.close()
    #             curr_umi = umi
    #             path = _umi_subdir(out_dir, umi)  # sharded subfolders
    #             fh = open(path, "w", buffering=1024 * 1024)
    #         rn = r.query_name
    #         seq = r.query_sequence or ""
    #         qual = r.qual if r.qual is not None else ("I" * len(seq))
    #         fh.write(f"@{rn}\n{seq}\n+\n{qual}\n")
    #         written += 1
    # if fh: fh.close()
    # logger.log_message(f"Done. Wrote {written} reads across {len(keepers)} UMIs.")
