import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from utils import log
from utils import calculator



def group_and_extract(input_bam, reference,motif5, motif3 ,output_dir):
    logger = log()
    calc = calculator()

    #Find the index of the motif5 and motif3 in the refernece sequence
    motif5_index = None
    motif3_index = None
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    logger.log_message(f"Reference sequence loaded from {reference}")
    logger.log_message(f"Reference sequence length: {len(ref_seq)}")

    motif5_index = ref_seq.find(motif5) + len(motif5)+ 1
    if motif5_index == -1:
        logger.log_message(f"Motif5 '{motif5}' not found in reference sequence.")
        return
    motif3_index = ref_seq.find(motif3) - len(motif3) - 1
    if motif3_index == -1:
        logger.log_message(f"Motif3 '{motif3}' not found in reference sequence.")
        return
    logger.log_message(f"Motif5 found at index {motif5_index}, Motif3 found at index {motif3_index} in reference sequence.")

    # Open the BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    logger.log_message(f"BAM file opened: {input_bam}")
    #Extract 5' UMI and 3' UMI from the reads
    raw_umi5_list = []
    raw_umi3_list = []

    # Prepare to write to a new BAM file with updated query names
    output_bam_path = f"{output_dir}/grouped_reads.bam"
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bamfile)

    region_start = motif5_index  # 1-based reference coordinate
    region_end = motif3_index + 12  # exclusive

    for read in bamfile.fetch(start=region_start - 1, end=region_end):
        # Map reference positions to query positions
        ref_positions = read.get_reference_positions(full_length=True)
        # Find query indices corresponding to the region
        query_start = None
        query_end = None
        for i, ref_pos in enumerate(ref_positions):
            if ref_pos == region_start - 1 and query_start is None:
                query_start = i
            if ref_pos == region_end - 1:
                query_end = i + 1  # end is exclusive
                break
        if query_start is not None and query_end is not None:
            seq = read.query_sequence[query_start:query_end]
            raw_umi5_list.append(seq)
            continue
    
    print(len(raw_umi5_list))
   
   
    bamfile.close()
    output_bam.close()





    logger.log_message(f"Written updated BAM file with new query names to {output_bam_path}")
    logger.log_message(f"Extracted {len(raw_umi5_list)} 5' UMIs and {len(raw_umi3_list)} 3' UMIs from the reads.")





if __name__ == "__main__":
    input_bam = "/Users/nomur/GitHub/test_data_belt/test_data/NL43-Test-Aligned.bam"
    reference = "/Users/nomur/GitHub/test_data_belt/test_data/NL43-Test-Aligned.fasta"  # Example reference sequence
    motif5 = "ACCGTCGGAT"
    motif3 = "TGTTGGCCTC"
    output_dir = "/Users/nomur/GitHub/belt/output_files"

    group_and_extract(input_bam, reference, motif5, motif3, output_dir)
