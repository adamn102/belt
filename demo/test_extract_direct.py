#!/usr/bin/env python3
"""
Test script to call extract.py function directly
"""

import sys
import os

# Add the scripts directory to the path
sys.path.insert(0, '../scripts')

# Import the extract function
from extract import extract_UMIs

def test_extract_direct():
    """Test the extract_UMIs function directly"""
    
    print("🧪 TESTING EXTRACT.PY FUNCTION DIRECTLY")
    print("=" * 50)
    
    # Parameters
    reference = "data/0-reference_sequence/B145-XTMMF9-12.fa"
    bam_file = "data/4-sorted_reads/B145-XTMMF9-12.bam"
    barcode_motif = "NNNNNNNNNNNNNNM"
    min_read_length = 1000
    output_dir = "data/test_extract_output"
    
    print(f"Reference: {reference}")
    print(f"BAM file: {bam_file}")
    print(f"Barcode motif: {barcode_motif}")
    print(f"Min read length: {min_read_length}")
    print(f"Output directory: {output_dir}")
    
    # Check if files exist
    if not os.path.exists(reference):
        print(f"❌ Reference file not found: {reference}")
        return
    if not os.path.exists(bam_file):
        print(f"❌ BAM file not found: {bam_file}")
        return
    
    print("\n✅ Input files exist, calling extract_UMIs...")
    
    try:
        # Call the extract function
        extract_UMIs(reference, bam_file, min_read_length, barcode_motif, output_dir)
        
        print("\n✅ Extract function completed successfully!")
        
        # Check output
        output_file = f"{output_dir}/B145-XTMMF9-12_rawUMIs_read_names.csv"
        if os.path.exists(output_file):
            print(f"\n📄 Output file created: {output_file}")
            print("Content:")
            with open(output_file, 'r') as f:
                content = f.read()
                print(content)
        else:
            print(f"\n❌ Output file not found: {output_file}")
            
    except Exception as e:
        print(f"\n❌ Error calling extract function: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_extract_direct()
