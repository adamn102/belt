#!/usr/bin/env python3
"""
Test script to demonstrate UMI extraction with synthetic data
"""

import pysam
from Bio import SeqIO
import tempfile
import os

def create_test_data():
    """Create test data with actual UMI sequences"""
    
    # Create a simple reference with UMI motifs
    ref_seq = "ATCG" * 100 + "NNNNNNNNNNNNNNM" + "GCTA" * 100 + "NNNNNNNNNNNNNNM" + "ATCG" * 100
    
    # Create reference file
    with open("test_ref.fa", "w") as f:
        f.write(">test_ref\n")
        f.write(ref_seq + "\n")
    
    # Create a test read with actual UMI sequences
    read_seq = "ATCG" * 100 + "ATCGATCGATCGAT" + "GCTA" * 100 + "GCTAGCTAGCTAG" + "ATCG" * 100
    
    # Create BAM file with test read
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': len(ref_seq), 'SN': 'test_ref'}]
    }
    
    with pysam.AlignmentFile("test.bam", "wb", header=header) as bam:
        a = pysam.AlignedSegment()
        a.query_name = "test_read"
        a.query_sequence = read_seq
        a.reference_id = 0
        a.reference_start = 0
        a.cigarstring = f"{len(read_seq)}M"
        a.mapping_quality = 60
        bam.write(a)
    
    # Index the BAM file
    pysam.index("test.bam")
    
    return "test_ref.fa", "test.bam"

def test_umi_extraction():
    """Test UMI extraction with synthetic data"""
    
    print("🧪 TESTING UMI EXTRACTION WITH SYNTHETIC DATA")
    print("=" * 60)
    
    # Create test data
    ref_file, bam_file = create_test_data()
    
    # Parameters
    barcode_motif = "NNNNNNNNNNNNNNM"
    min_read_length = 100
    
    print(f"Created test reference: {ref_file}")
    print(f"Created test BAM: {bam_file}")
    print(f"Barcode motif: {barcode_motif}")
    
    # Load reference
    with open(ref_file, "r") as f:
        ref_seq = str(SeqIO.read(f, "fasta").seq)
    
    # Find UMI motifs
    motif5_index = ref_seq.find(barcode_motif)
    motif3_index = ref_seq.rfind(barcode_motif)
    
    print(f"\nUMI motifs found:")
    print(f"  5' motif at position: {motif5_index}")
    print(f"  3' motif at position: {motif3_index}")
    
    # Extract UMIs
    barcode_length = len(barcode_motif)
    umi5_start = motif5_index + 1
    umi5_end = motif5_index + barcode_length
    umi3_start = motif3_index + 1
    umi3_end = motif3_index + barcode_length
    
    print(f"\nUMI regions:")
    print(f"  5' UMI region: {umi5_start}-{umi5_end}")
    print(f"  3' UMI region: {umi3_start}-{umi3_end}")
    
    # Process reads
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    raw_UMI_dict = {}
    
    for read in bamfile.fetch():
        print(f"\nProcessing read: {read.query_name}")
        print(f"  Read length: {len(read.query_sequence)}")
        
        if read.query_sequence is None or len(read.query_sequence) < min_read_length:
            print("  ❌ Skipped: No sequence or too short")
            continue
        
        # Get reference positions
        ref_positions = read.get_reference_positions(full_length=True)
        
        # Extract 5' UMI
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
            print(f"  5' UMI: {seq5p} (length: {len(seq5p)})")
            if len(seq5p) == barcode_length:
                print(f"  ✅ 5' UMI extracted successfully")
                raw_UMI_dict[read.query_name] = [seq5p, ""]
            else:
                print(f"  ❌ 5' UMI length mismatch")
        else:
            print(f"  ❌ Could not find 5' UMI region")
        
        # Extract 3' UMI
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
            print(f"  3' UMI: {seq3p} (length: {len(seq3p)})")
            if len(seq3p) == barcode_length:
                print(f"  ✅ 3' UMI extracted successfully")
                if read.query_name in raw_UMI_dict:
                    raw_UMI_dict[read.query_name] = [raw_UMI_dict[read.query_name][0], seq3p]
                else:
                    raw_UMI_dict[read.query_name] = ["", seq3p]
            else:
                print(f"  ❌ 3' UMI length mismatch")
        else:
            print(f"  ❌ Could not find 3' UMI region")
    
    bamfile.close()
    
    print(f"\n📊 RESULTS:")
    print(f"  UMIs extracted: {len(raw_UMI_dict)}")
    for read_name, (umi5, umi3) in raw_UMI_dict.items():
        print(f"  {read_name}: 5'={umi5}, 3'={umi3}")
    
    # Clean up
    os.remove(ref_file)
    os.remove(bam_file)
    os.remove(bam_file + ".bai")
    
    print(f"\n✅ Test completed successfully!")

if __name__ == "__main__":
    test_umi_extraction()
