#!/usr/bin/env python3
"""
Test script to fix the extract logic
"""

import sys
import os
import pysam
from Bio import SeqIO

# Add the scripts directory to the path
sys.path.insert(0, '../scripts')

def test_extract_fixed():
    """Test the extract logic with fixed approach"""
    
    print("🔧 TESTING FIXED EXTRACT LOGIC")
    print("=" * 50)
    
    # Parameters
    reference = "data/0-reference_sequence/B145-XTMMF9-12.fa"
    bam_file = "data/4-sorted_reads/B145-XTMMF9-12.bam"
    barcode_motif = "NNNNNNNNNNNNNNM"
    min_read_length = 1000
    
    print(f"Reference: {reference}")
    print(f"BAM file: {bam_file}")
    print(f"Barcode motif: {barcode_motif}")
    
    # Load reference
    with open(reference, "r") as ref_file:
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
    
    # Find UMI motifs
    motif5_index = ref_seq.find(barcode_motif)
    motif3_index = ref_seq.rfind(barcode_motif)
    
    print(f"\nUMI motifs found:")
    print(f"  5' motif at position: {motif5_index}")
    print(f"  3' motif at position: {motif3_index}")
    
    # Calculate UMI regions
    barcode_length = len(barcode_motif)
    umi5_start = motif5_index + 1
    umi5_end = motif5_index + barcode_length
    umi3_start = motif3_index + 1
    umi3_end = motif3_index + barcode_length
    
    print(f"\nUMI regions:")
    print(f"  5' UMI region: {umi5_start}-{umi5_end}")
    print(f"  3' UMI region: {umi3_start}-{umi3_end}")
    
    # Test 1: Original logic (only reads spanning both regions)
    print(f"\n🧪 Test 1: Original logic (reads spanning both regions)")
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    reads_spanning_both = 0
    for read in bamfile.fetch(start=umi5_start - 1, end=umi3_end):
        reads_spanning_both += 1
        if reads_spanning_both <= 3:
            print(f"  Read {reads_spanning_both}: {read.query_name} (start: {read.reference_start}, end: {read.reference_end})")
    bamfile.close()
    print(f"  Total reads spanning both regions: {reads_spanning_both}")
    
    # Test 2: Fixed logic (reads in 5' region)
    print(f"\n🧪 Test 2: Fixed logic (reads in 5' region)")
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    reads_in_5p = 0
    umis_5p = []
    for read in bamfile.fetch(start=umi5_start - 10, end=umi5_end + 10):
        reads_in_5p += 1
        if reads_in_5p <= 3:
            print(f"  Read {reads_in_5p}: {read.query_name} (start: {read.reference_start}, end: {read.reference_end})")
            
            if read.query_sequence and len(read.query_sequence) >= min_read_length:
                ref_positions = read.get_reference_positions(full_length=True)
                
                # Find 5' UMI
                query_start_5p = None
                query_end_5p = None
                for i, ref_pos in enumerate(ref_positions):
                    if ref_pos == umi5_start - 1 and query_start_5p is None:
                        query_start_5p = i
                    if ref_pos == umi5_end - 1:
                        query_end_5p = i + 1
                        break
                
                if query_start_5p is not None and query_end_5p is not None:
                    seq5p = read.query_sequence[query_start_5p:query_end_5p]
                    print(f"    5' UMI: {seq5p} (length: {len(seq5p)})")
                    if len(seq5p) == barcode_length:
                        umis_5p.append(seq5p)
                        print(f"    ✅ 5' UMI extracted!")
                    else:
                        print(f"    ❌ 5' UMI length mismatch")
                else:
                    print(f"    ❌ Could not find 5' UMI region")
    bamfile.close()
    print(f"  Total reads in 5' region: {reads_in_5p}")
    print(f"  UMIs extracted from 5' region: {len(umis_5p)}")
    
    # Test 3: Fixed logic (reads in 3' region)
    print(f"\n🧪 Test 3: Fixed logic (reads in 3' region)")
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    reads_in_3p = 0
    umis_3p = []
    for read in bamfile.fetch(start=umi3_start - 10, end=umi3_end + 10):
        reads_in_3p += 1
        if reads_in_3p <= 3:
            print(f"  Read {reads_in_3p}: {read.query_name} (start: {read.reference_start}, end: {read.reference_end})")
            
            if read.query_sequence and len(read.query_sequence) >= min_read_length:
                ref_positions = read.get_reference_positions(full_length=True)
                
                # Find 3' UMI
                query_start_3p = None
                query_end_3p = None
                for i, ref_pos in enumerate(ref_positions):
                    if ref_pos == umi3_start - 1 and query_start_3p is None:
                        query_start_3p = i
                    if ref_pos == umi3_end - 1:
                        query_end_3p = i + 1
                        break
                
                if query_start_3p is not None and query_end_3p is not None:
                    seq3p = read.query_sequence[query_start_3p:query_end_3p]
                    print(f"    3' UMI: {seq3p} (length: {len(seq3p)})")
                    if len(seq3p) == barcode_length:
                        umis_3p.append(seq3p)
                        print(f"    ✅ 3' UMI extracted!")
                    else:
                        print(f"    ❌ 3' UMI length mismatch")
                else:
                    print(f"    ❌ Could not find 3' UMI region")
    bamfile.close()
    print(f"  Total reads in 3' region: {reads_in_3p}")
    print(f"  UMIs extracted from 3' region: {len(umis_3p)}")
    
    print(f"\n📊 SUMMARY:")
    print(f"  Reads spanning both regions: {reads_spanning_both}")
    print(f"  Reads in 5' region: {reads_in_5p}")
    print(f"  Reads in 3' region: {reads_in_3p}")
    print(f"  UMIs from 5' region: {len(umis_5p)}")
    print(f"  UMIs from 3' region: {len(umis_3p)}")

if __name__ == "__main__":
    test_extract_fixed()
