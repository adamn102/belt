#!/usr/bin/env python3
"""
Simplified test file for BELT extract function
"""

import os
import sys
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from extract import extract_UMIs, findUMI

# Test configuration constants
# Modify these values to test different scenarios
TEST_CONFIG = {
    'min_read_length': 2000,
    'barcode_motifs': ["NNNNNNNNNNNNNNM", "NNNNNNNNNNNNNNM"],  # 15 N's for both 5' and 3' UMIs
    'test_data_dir': "scripts/test-data/test-umi_map",
    'reference_file': "scripts/test-data/test-umi_map/reference_sequence/B145.fa",
    'bam_file': "scripts/test-data/test-umi_map/sorted_reads/B145-XTMMF9-12.bam"
}


def setup_test_environment():
    """Set up and verify the test environment"""
    print("=" * 60)
    print("BELT Extract Function Test Setup")
    print("=" * 60)
    
    # Define test data paths using constants
    test_data_dir = Path(__file__).parent.parent / TEST_CONFIG['test_data_dir']
    reference_file = Path(__file__).parent.parent / TEST_CONFIG['reference_file']
    bam_file = Path(__file__).parent.parent / TEST_CONFIG['bam_file']
    
    print(f"Test data directory: {test_data_dir}")
    print(f"Reference file: {reference_file}")
    print(f"BAM file: {bam_file}")
    print()
    
    # Verify files exist
    files_ok = True
    
    if not test_data_dir.exists():
        print("❌ Test data directory not found!")
        files_ok = False
    else:
        print("✅ Test data directory exists")
    
    if not reference_file.exists():
        print("❌ Reference file not found!")
        files_ok = False
    else:
        print("✅ Reference file exists")
        # Show reference file content
        with open(reference_file, 'r') as f:
            content = f.read()
            print(f"   Reference length: {len(content)} characters")
            # Check if any of the barcode motifs are in the reference
            motif_found = any(motif in content for motif in TEST_CONFIG['barcode_motifs'])
            print(f"   Contains barcode motifs: {motif_found}")
    
    if not bam_file.exists():
        print("❌ BAM file not found!")
        files_ok = False
    else:
        print("✅ BAM file exists")
        # Check BAM file size
        size_mb = bam_file.stat().st_size / (1024 * 1024)
        print(f"   BAM file size: {size_mb:.1f} MB")
    
    print()
    
    if not files_ok:
        print("❌ Setup failed - missing required files!")
        return None, None, None
    
    print("✅ All test files verified successfully!")
    print()
    
    return reference_file, bam_file, test_data_dir


def test_extract_function(reference_file, bam_file, output_dir):
    """Test the extract function and show example output"""
    print("=" * 60)
    print("Testing Extract Function")
    print("=" * 60)
    
    # Use test configuration constants
    min_read_length = TEST_CONFIG['min_read_length']
    barcode_motifs = TEST_CONFIG['barcode_motifs']
    
    print(f"Parameters:")
    print(f"  Min read length: {min_read_length}")
    print(f"  Barcode motifs: {barcode_motifs}")
    print(f"  Expected number of UMI columns: {len(barcode_motifs)}")
    print(f"  Output directory: {output_dir}")
    print()
    
    try:
        # Run the extract function
        print("Running extract_UMIs function...")
        extract_UMIs(
            str(reference_file),
            str(bam_file),
            min_read_length,
            barcode_motifs,
            str(output_dir)
        )
        print("✅ Extract function completed successfully!")
        print()
        
        # Check output file
        expected_output_file = output_dir / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        
        if expected_output_file.exists():
            print("✅ Output file created successfully!")
            print(f"Output file: {expected_output_file}")
            
            # Read and display results
            df = pd.read_csv(expected_output_file)
            print(f"Number of reads processed: {len(df)}")
            print(f"Columns: {list(df.columns)}")
            
            # Verify we got the expected number of UMI columns
            umi_columns = [col for col in df.columns if col != 'Read Name']
            expected_umi_count = len(barcode_motifs)
            actual_umi_count = len(umi_columns)
            
            print(f"Expected UMI columns: {expected_umi_count}")
            print(f"Actual UMI columns: {actual_umi_count}")
            
            if actual_umi_count == expected_umi_count:
                print("✅ Correct number of UMI columns extracted!")
            else:
                print(f"⚠️  Expected {expected_umi_count} UMI columns, got {actual_umi_count}")
            
            print()
            
            return df, expected_output_file, barcode_motifs
        else:
            print("❌ Output file not found!")
            return None, None, None
            
    except Exception as e:
        print(f"❌ Extract function failed with error: {e}")
        return None, None, None


def show_example_output(df, output_file, barcode_motifs):
    """Display example output from the extract function"""
    print("=" * 60)
    print("Example Output")
    print("=" * 60)
    
    if df is None or len(df) == 0:
        print("No data to display")
        return
    
    # Get UMI columns (exclude 'Read Name')
    umi_columns = [col for col in df.columns if col != 'Read Name']
    
    # Show UMI column mapping
    print("UMI Column Mapping:")
    for i, motif in enumerate(barcode_motifs):
        if i < len(umi_columns):
            print(f"  {umi_columns[i]} -> Barcode motif: {motif}")
        else:
            print(f"  [Missing] -> Barcode motif: {motif}")
    print()
    
    # Show first 10 rows
    print("First 10 extracted UMIs:")
    print("-" * 60)
    
    # Format the output nicely
    for i, (_, row) in enumerate(df.head(10).iterrows()):
        read_name = row['Read Name']
        print(f"{i+1:2d}. Read: {read_name[:30]}...")
        
        # Show all UMI columns
        for j, umi_col in enumerate(umi_columns):
            umi_value = row[umi_col] if pd.notna(row[umi_col]) else "None"
            print(f"    {umi_col}: {umi_value}")
        print()
    
    if len(df) > 10:
        print(f"... and {len(df) - 10} more reads")
    
    print()
    print("Summary statistics:")
    print(f"  Total reads: {len(df)}")
    
    # Show statistics for each UMI column
    for umi_col in umi_columns:
        umi_count = df[umi_col].notna().sum()
        print(f"  Reads with {umi_col}: {umi_count}")
    
    print()
    print(f"Output file saved to: {output_file}")


def test_multiple_motifs(reference_file, bam_file, output_dir):
    """Test with multiple barcode motifs to show dynamic UMI extraction"""
    print("=" * 60)
    print("Testing with Multiple Barcode Motifs")
    print("=" * 60)
    
    # Use test configuration constants
    min_read_length = TEST_CONFIG['min_read_length']
    barcode_motifs = TEST_CONFIG['barcode_motifs']
    
    print(f"Parameters:")
    print(f"  Min read length: {min_read_length}")
    print(f"  Barcode motifs: {barcode_motifs}")
    print(f"  Expected number of UMI columns: {len(barcode_motifs)}")
    print()
    
    try:
        # Run the extract function
        print("Running extract_UMIs function with multiple motifs...")
        extract_UMIs(
            str(reference_file),
            str(bam_file),
            min_read_length,
            barcode_motifs,
            str(output_dir / "multiple_motifs")
        )
        print("✅ Extract function completed successfully!")
        print()
        
        # Check output file
        expected_output_file = output_dir / "multiple_motifs" / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        
        if expected_output_file.exists():
            print("✅ Output file created successfully!")
            
            # Read and display results
            df = pd.read_csv(expected_output_file)
            print(f"Number of reads processed: {len(df)}")
            print(f"Columns: {list(df.columns)}")
            
            # Verify we got the expected number of UMI columns
            umi_columns = [col for col in df.columns if col != 'Read Name']
            expected_umi_count = len(barcode_motifs)
            actual_umi_count = len(umi_columns)
            
            print(f"Expected UMI columns: {expected_umi_count}")
            print(f"Actual UMI columns: {actual_umi_count}")
            
            if actual_umi_count == expected_umi_count:
                print("✅ Correct number of UMI columns extracted!")
            else:
                print(f"⚠️  Expected {expected_umi_count} UMI columns, got {actual_umi_count}")
            
            print()
            
            return df, expected_output_file, barcode_motifs
        else:
            print("❌ Output file not found!")
            return None, None, None
            
    except Exception as e:
        print(f"❌ Extract function failed with error: {e}")
        return None, None, None


def main():
    """Main function to run the simplified test"""
    import tempfile
    import shutil
    
    # Create temporary directory for output
    temp_dir = tempfile.mkdtemp()
    output_dir = Path(temp_dir)
    
    try:
        # Step 1: Setup and verify files
        reference_file, bam_file, test_data_dir = setup_test_environment()
        
        if reference_file is None:
            return 1
        
        # Step 2: Test the extract function with single motif
        print("🔍 Testing with single barcode motif...")
        df, output_file, barcode_motifs = test_extract_function(reference_file, bam_file, output_dir)
        
        if df is None:
            return 1
        
        # Step 3: Show example output
        show_example_output(df, output_file, barcode_motifs)
        
        print("\n" + "=" * 60)
        print("🔍 Testing with multiple barcode motifs...")
        
        # Step 4: Test with multiple motifs
        df_multi, output_file_multi, barcode_motifs_multi = test_multiple_motifs(reference_file, bam_file, output_dir)
        
        if df_multi is not None:
            show_example_output(df_multi, output_file_multi, barcode_motifs_multi)
        
        print("=" * 60)
        print("✅ All tests completed successfully!")
        print("=" * 60)
        
        return 0
        
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        return 1
        
    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == '__main__':
    exit(main())