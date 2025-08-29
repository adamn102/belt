#!/usr/bin/env python3
"""
Demo script to show real umi-group functionality working with test data
"""

import sys
import os
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

from umi_group import long_read_umi_grouping, long_read_umi_grouping_fast

def main():
    """Demo the umi-group functionality with real test data"""
    
    # Paths to real test data
    reference_file = "scripts/test-data/test-umi_map/reference_sequence/B145.fa"
    bam_file = "scripts/test-data/test-umi_map/sorted_reads/B145-XTMMF9-12.bam"
    corrected_umis_file = "scripts/test-data/5-umi_maps/B145-XTMMF9-12_Paired_UMIs_read_names.csv"
    output_dir = "test_umi_group_output"
    grouping = "5pUMI"
    threshold = 5
    threads = 4
    
    print("🧪 BELT UMI-Group Function Demo")
    print("=" * 50)
    print(f"Reference file: {reference_file}")
    print(f"BAM file: {bam_file}")
    print(f"Corrected UMIs file: {corrected_umis_file}")
    print(f"Grouping: {grouping}")
    print(f"Threshold: {threshold}")
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
    if not os.path.exists(corrected_umis_file):
        print(f"❌ Corrected UMIs file not found: {corrected_umis_file}")
        return
    print("✅ All input files found!")
    print()
    
    # Test 2: Analyze the corrected UMIs data
    print("📊 Analyzing corrected UMIs data...")
    import pandas as pd
    
    df = pd.read_csv(corrected_umis_file)
    print(f"   Total reads: {len(df)}")
    print(f"   Valid 5' UMIs: {df['5pUMI'].notna().sum()}")
    print(f"   Valid 3' UMIs: {df['3pUMI'].notna().sum()}")
    
    # Clean up 5pUMI (same as in the function)
    df["5pUMI"] = df["5pUMI"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
    
    # Group by 5pUMI
    grouped = df.groupby("5pUMI").agg({
        "Read Name": ["count", lambda x: list(x)]
    })
    grouped.columns = ["ReadCount", "ReadNames"]
    
    # Filter groups by threshold
    filtered = grouped[grouped["ReadCount"] >= threshold]
    
    print(f"   Total UMI groups: {len(grouped)}")
    print(f"   Groups meeting threshold ({threshold}): {len(filtered)}")
    
    if len(filtered) > 0:
        print(f"   Largest group size: {filtered['ReadCount'].max()}")
        print(f"   Average group size: {filtered['ReadCount'].mean():.1f}")
        
        # Show some example groups
        print("\n📋 Example UMI groups:")
        for i, (umi, row) in enumerate(filtered.head(5).iterrows()):
            print(f"   {umi}: {row['ReadCount']} reads")
            if i >= 4:  # Show only first 5
                break
    print()
    
    # Test 3: Run the main umi-grouping function
    print("🔍 Testing long_read_umi_grouping function...")
    try:
        long_read_umi_grouping(
            reference_file,
            bam_file,
            corrected_umis_file,
            output_dir,
            grouping,
            threshold
        )
        
        # Check output files
        sample_name = "B145-XTMMF9-12"
        expected_filtered_csv = os.path.join(output_dir, f"{sample_name}_filtered_umi_groups.csv")
        expected_sample_dir = os.path.join(output_dir, sample_name)
        
        if os.path.exists(expected_filtered_csv):
            print(f"✅ Filtered UMI groups CSV created: {expected_filtered_csv}")
            
            # Read and show stats
            filtered_df = pd.read_csv(expected_filtered_csv)
            print(f"   Groups in output: {len(filtered_df)}")
            print(f"   Total reads in groups: {filtered_df['ReadCount'].sum()}")
            
            # Show some example groups
            print("\n📋 Output UMI groups:")
            for _, row in filtered_df.head(5).iterrows():
                print(f"   {row['5pUMI']}: {row['ReadCount']} reads")
        else:
            print(f"❌ Filtered CSV not found: {expected_filtered_csv}")
        
        if os.path.exists(expected_sample_dir):
            print(f"✅ Sample output directory created: {expected_sample_dir}")
            
            # Count FASTQ files
            fastq_files = list(Path(expected_sample_dir).glob("*.fastq"))
            print(f"   FASTQ files created: {len(fastq_files)}")
            
            # Show some example FASTQ files
            if fastq_files:
                print("\n📋 Example FASTQ files:")
                for fastq_file in fastq_files[:3]:
                    file_size = fastq_file.stat().st_size
                    print(f"   {fastq_file.name}: {file_size} bytes")
                    
                    # Show first few lines of first file
                    if fastq_file == fastq_files[0]:
                        with open(fastq_file, 'r') as f:
                            lines = f.readlines()[:8]  # First 2 reads (8 lines)
                            print("   First 2 reads:")
                            for line in lines:
                                print(f"     {line.strip()}")
        else:
            print(f"❌ Sample directory not found: {expected_sample_dir}")
            
    except Exception as e:
        print(f"❌ long_read_umi_grouping failed: {e}")
        return
    
    # Test 4: Try the fast function (if samtools is available)
    print("\n⚡ Testing long_read_umi_grouping_fast function...")
    import subprocess
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True)
        samtools_available = True
        print("✅ samtools is available")
    except (subprocess.CalledProcessError, FileNotFoundError):
        samtools_available = False
        print("⚠️  samtools not available - skipping fast function")
    
    if samtools_available:
        try:
            fast_output_dir = f"{output_dir}_fast"
            long_read_umi_grouping_fast(
                reference_file,
                bam_file,
                corrected_umis_file,
                fast_output_dir,
                grouping,
                threshold,
                threads
            )
            
            # Check fast output
            fast_filtered_csv = os.path.join(fast_output_dir, f"{sample_name}_filtered_umi_groups.csv")
            if os.path.exists(fast_filtered_csv):
                print(f"✅ Fast function completed: {fast_filtered_csv}")
                
                # Compare with regular output
                fast_df = pd.read_csv(fast_filtered_csv)
                print(f"   Groups in fast output: {len(fast_df)}")
                print(f"   Total reads in fast output: {fast_df['ReadCount'].sum()}")
                
                # Should match regular output
                if os.path.exists(expected_filtered_csv):
                    regular_df = pd.read_csv(expected_filtered_csv)
                    if len(fast_df) == len(regular_df):
                        print("✅ Fast and regular outputs have same number of groups")
                    else:
                        print(f"⚠️  Fast output has {len(fast_df)} groups, regular has {len(regular_df)}")
            else:
                print(f"❌ Fast output not found: {fast_filtered_csv}")
                
        except Exception as e:
            print(f"❌ long_read_umi_grouping_fast failed: {e}")
    
    print("\n🎉 Demo completed successfully!")

if __name__ == "__main__":
    main()
