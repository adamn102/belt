import random
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv

def random_dna(length):
    return ''.join(random.choices("ACGT", k=length))

def random_qual(length, min_q=30, max_q=40):
    return [random.randint(min_q, max_q) for _ in range(length)]

def generate_reference_with_Ns(read_length, umi_length, umi5_pos, umi3_pos):
    ref = list(random_dna(read_length))
    ref[umi5_pos:umi5_pos+umi_length] = ['N'] * umi_length
    ref[umi3_pos:umi3_pos+umi_length] = ['N'] * umi_length
    return ''.join(ref)

def generate_read_with_umis(ref_seq, umi_length, umi5_pos, umi3_pos, ins_rate=0.01, del_rate=0.01):
    read = list(ref_seq)
    qual = random_qual(len(ref_seq))

    umi5 = random_dna(umi_length)
    umi3 = random_dna(umi_length)

    read[umi5_pos:umi5_pos+umi_length] = list(umi5)
    read[umi3_pos:umi3_pos+umi_length] = list(umi3)

    qual[umi5_pos:umi5_pos+umi_length] = [40] * umi_length
    qual[umi3_pos:umi3_pos+umi_length] = [40] * umi_length

    # Apply indels around UMI positions
    read, qual = apply_indels(read, qual, umi5_pos, window=5, ins_rate=ins_rate, del_rate=del_rate)
    read, qual = apply_indels(read, qual, umi3_pos, window=5, ins_rate=ins_rate, del_rate=del_rate)

    sequence = ''.join(read)
    quality_str = ''.join([chr(q + 33) for q in qual])

    return sequence, quality_str, umi5, umi3


def apply_indels(seq, qual, pos, window=5, ins_rate=0.01, del_rate=0.01):
    new_seq = []
    new_qual = []

    for i in range(len(seq)):
        in_window = abs(i - pos) <= window

        # Random deletion
        if in_window and random.random() < del_rate:
            continue  # Skip this base (deleted)

        # Random insertion
        if in_window and random.random() < ins_rate:
            inserted_base = random.choice("ACGT")
            inserted_qual = random.randint(30, 40)
            new_seq.append(inserted_base)
            new_qual.append(inserted_qual)

        # Add original base
        new_seq.append(seq[i])
        new_qual.append(qual[i])

    return new_seq, new_qual


def main(args):
    # Generate reference with Ns
    ref_seq = generate_reference_with_Ns(args.read_length, args.umi_length, args.umi5_pos, args.umi3_pos)

    # Write reference FASTA
    with open(args.reference_output, 'w') as ref_fh:
        ref_record = SeqRecord(Seq(ref_seq), id="synthetic_ref", description="Reference with Ns at UMI positions")
        SeqIO.write(ref_record, ref_fh, "fasta")

    # Prepare for CSV
    umi_records = []

    # Generate reads with UMIs
    with open(args.output, 'w') as fastq_fh:
        for i in range(args.num_reads):
            seq, qual, umi5, umi3 = generate_read_with_umis(
                ref_seq,
                args.umi_length,
                args.umi5_pos,
                args.umi3_pos
            )
            read_id = f"read_{i}_UMI5:{umi5}_UMI3:{umi3}"
            record = SeqRecord(
                Seq(seq),
                id=read_id,
                description=""
            )
            record.letter_annotations["phred_quality"] = [ord(c) - 33 for c in qual]
            SeqIO.write(record, fastq_fh, "fastq")

            umi_records.append({'read_id': read_id, 'umi5': umi5, 'umi3': umi3})

    # Write UMI CSV
    with open(args.umi_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['read_id', 'umi5', 'umi3'])
        writer.writeheader()
        writer.writerows(umi_records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic FASTQ reads with 2 UMIs and a masked reference")
    parser.add_argument("--read_length", type=int, default=100, help="Total read length")
    parser.add_argument("--umi_length", type=int, default=8, help="UMI length")
    parser.add_argument("--umi5_pos", type=int, default=0, help="5' UMI position")
    parser.add_argument("--umi3_pos", type=int, default=92, help="3' UMI position")
    parser.add_argument("--num_reads", type=int, default=10, help="Number of reads to generate")
    parser.add_argument("--output", type=str, default="output.fastq", help="Output FASTQ file")
    parser.add_argument("--reference_output", type=str, default="reference_with_Ns.fasta", help="Output FASTA reference with Ns")
    parser.add_argument("--umi_csv", type=str, default="umi_table.csv", help="CSV file to write UMI sequences")
    parser.add_argument("--ins_rate", type=float, default=0.01, help="Insertion rate around UMIs")
    parser.add_argument("--del_rate", type=float, default=0.01, help="Deletion rate around UMIs")


    
    args = parser.parse_args()
    main(args)
