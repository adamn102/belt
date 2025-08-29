#!/usr/bin/env python3
"""
Demo script to show real extract functionality working with test data
"""

import sys
import os
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

from extract import extract_UMIs, findUMI

def main():
    """Demo the extract functionality with real test data"""
    
    # Paths to real test data
    reference_file = "scripts/test-data/test-umi_map/reference_sequence/B145.fa"
    bam_file = "scripts/test-data/test-umi_map/sorted_reads/B145-XTMMF9-12.bam"
    output_dir = "test_output"
    barcode_motif = "NNNNNNNNNNNNN"  # 13 N's as seen in reference
    min_read_length = 2000
    
    print("🧪 BELT Extract Function Demo")
    print("=" * 50)
    print(f"Reference file: {reference_file}")
    print(f"BAM file: {bam_file}")
    print(f"Barcode motif: {barcode_motif}")
    print(f"Min read length: {min_read_length}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Test 1: Check if files exist
    print("✅ Checking input files...")
    if not os.path.exists(reference_file):
        print(f"❌ Reference file not found: {reference_file}")
        return
    if not os.path.exists(bam_file):
        print(f"❌ BAM file not found: {bam_file}")
        return
    print("✅ All input files found!")
    print()
    
    # Test 2: Run findUMI function
    print("🔍 Testing findUMI function...")
    try:
        result = findUMI(reference_file, bam_file, barcode_motif, min_read_length)
        raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = result
        
        print(f"✅ findUMI completed successfully!")
        print(f"   5' motif index: {motif5_index}")
        print(f"   3' motif index: {motif3_index}")
        print(f"   5' UMIs found: {len(raw_umi5_list)}")
        print(f"   3' UMIs found: {len(raw_umi3_list)}")
        print(f"   Total reads with UMIs: {len(raw_UMI_dict)}")
        
        # Show some example UMIs
        if raw_UMI_dict:
            print("\n📋 Example UMIs:")
            count = 0
            for read_name, umis in raw_UMI_dict.items():
                if count < 5:  # Show first 5
                    print(f"   {read_name}: 5'={umis[0]}, 3'={umis[1]}")
                    count += 1
                else:
                    break
        print()
        
    except Exception as e:
        print(f"❌ findUMI failed: {e}")
        return
    
    # Test 3: Run extract_UMIs function
    print("📝 Testing extract_UMIs function...")
    try:
        extract_UMIs(reference_file, bam_file, min_read_length, barcode_motif, output_dir)
        
        # Check output file
        expected_output_file = os.path.join(output_dir, "B145-XTMMF9-12_rawUMIs_read_names.csv")
        if os.path.exists(expected_output_file):
            print(f"✅ Output file created: {expected_output_file}")
            
            # Read and show some stats
            import pandas as pd
            df = pd.read_csv(expected_output_file)
            print(f"   Total rows: {len(df)}")
            print(f"   Columns: {list(df.columns)}")
            
            # Show some example data
            print("\n📋 Example output data:")
            print(df.head().to_string(index=False))
            
        else:
            print(f"❌ Output file not found: {expected_output_file}")
            
    except Exception as e:
        print(f"❌ extract_UMIs failed: {e}")
        return
    
    print("\n🎉 Demo completed successfully!")

if __name__ == "__main__":
    main()
