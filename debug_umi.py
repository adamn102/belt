#!/usr/bin/env python3
"""
Debug script for UMI extraction
"""

import pysam
from Bio import SeqIO
import sys
import os

def debug_umi_extraction():
    """Debug UMI extraction step by step"""
    
    # Parameters
    reference = "data/0-reference_sequence/B145-XTMMF9-12.fa"
    bam_file = "data/4-sorted_reads/B145-XTMMF9-12.bam"
    barcode_motif = "NNNNNNNNNNNNN"
    min_read_length = 1000
    
    print("🔍 DEBUGGING UMI EXTRACTION")
    print("=" * 50)
    
    # Step 1: Check if files exist
    print("Step 1: Checking input files...")
    if not os.path.exists(reference):
        print(f"❌ Reference file not found: {reference}")
        return
    if not os.path.exists(bam_file):
        print(f"❌ BAM file not found: {bam_file}")
        return
    print("✅ Input files exist")
    
    # Step 2: Load reference sequence
    print("\nStep 2: Loading reference sequence...")
    try:
        with open(reference, "r") as ref_file:
            ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
        print(f"✅ Reference sequence loaded, length: {len(ref_seq)}")
        print(f"   Reference name: {SeqIO.read(reference, 'fasta').id}")
    except Exception as e:
        print(f"❌ Error loading reference: {e}")
        return
    
    # Step 3: Find UMI motifs in reference
    print(f"\nStep 3: Finding UMI motifs '{barcode_motif}' in reference...")
    motif5_index = ref_seq.find(barcode_motif)
    motif3_index = ref_seq.rfind(barcode_motif)
    
    print(f"   5' motif found at position: {motif5_index}")
    print(f"   3' motif found at position: {motif3_index}")
    
    if motif5_index == -1:
        print("❌ First barcode region not found in reference sequence")
        return
    if motif3_index == -1:
        print("❌ Second barcode region not found in reference sequence")
        return
    if motif5_index == motif3_index:
        print("❌ Only one barcode region found (need two)")
        return
    
    print("✅ Both UMI motifs found in reference")
    
    # Step 4: Check BAM file
    print(f"\nStep 4: Checking BAM file...")
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        read_count = 0
        for read in bamfile.fetch():
            read_count += 1
            if read_count >= 5:  # Just check first few reads
                break
        bamfile.close()
        print(f"✅ BAM file is valid, contains reads")
    except Exception as e:
        print(f"❌ Error reading BAM file: {e}")
        return
    
    # Step 5: Check read alignments around UMI regions
    print(f"\nStep 5: Checking read alignments around UMI regions...")
    barcode_length = len(barcode_motif)
    umi5_start = motif5_index + 1
    umi5_end = motif5_index + barcode_length
    umi3_start = motif3_index + 1
    umi3_end = motif3_index + barcode_length
    
    print(f"   5' UMI region: {umi5_start}-{umi5_end}")
    print(f"   3' UMI region: {umi3_start}-{umi3_end}")
    
    # Step 6: Extract UMIs from first few reads
    print(f"\nStep 6: Extracting UMIs from first 5 reads...")
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    reads_processed = 0
    reads_with_5p_umi = 0
    reads_with_3p_umi = 0
    
    for read in bamfile.fetch(start=umi5_start - 1, end=umi3_end):
        reads_processed += 1
        if reads_processed > 5:
            break
            
        print(f"\n   Read {reads_processed}: {read.query_name}")
        print(f"     Alignment start: {read.reference_start}")
        print(f"     Alignment end: {read.reference_end}")
        print(f"     Read length: {len(read.query_sequence) if read.query_sequence else 'None'}")
        
        if read.query_sequence is None or len(read.query_sequence) < min_read_length:
            print(f"     ❌ Skipped: No sequence or too short")
            continue
            
        # Get reference positions
        ref_positions = read.get_reference_positions(full_length=True)
        print(f"     Reference positions array length: {len(ref_positions)}")
        
        # Check 5' UMI
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
            print(f"     5' UMI: {seq5p} (length: {len(seq5p)})")
            if len(seq5p) == barcode_length:
                reads_with_5p_umi += 1
                print(f"     ✅ 5' UMI extracted successfully")
            else:
                print(f"     ❌ 5' UMI length mismatch")
        else:
            print(f"     ❌ Could not find 5' UMI region in read")
            
        # Check 3' UMI
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
            print(f"     3' UMI: {seq3p} (length: {len(seq3p)})")
            if len(seq3p) == barcode_length:
                reads_with_3p_umi += 1
                print(f"     ✅ 3' UMI extracted successfully")
            else:
                print(f"     ❌ 3' UMI length mismatch")
        else:
            print(f"     ❌ Could not find 3' UMI region in read")
    
    bamfile.close()
    
    print(f"\n📊 SUMMARY:")
    print(f"   Reads processed: {reads_processed}")
    print(f"   Reads with 5' UMI: {reads_with_5p_umi}")
    print(f"   Reads with 3' UMI: {reads_with_3p_umi}")

if __name__ == "__main__":
    debug_umi_extraction()
