#!/usr/bin/env python3
"""
Real tests for BELT extract function using actual test data
"""

import unittest
import tempfile
import os
import sys
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from extract import extract_UMIs, findUMI

# Global variable to control UMI output display
SHOW_UMI_SAMPLES = True


class TestExtractReal(unittest.TestCase):
    """Test extract functionality with real test data"""
    
    def setUp(self):
        """Set up test fixtures with real test data"""
        # Paths to real test data
        self.test_data_dir = Path(__file__).parent.parent / "scripts" / "test-data" / "test-umi_map"
        self.reference_file = self.test_data_dir / "reference_sequence" / "B145.fa"
        self.bam_file = self.test_data_dir / "sorted_reads" / "B145-XTMMF9-12.bam"
        
        # Create temporary directory for outputs
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir)
        
        # Test parameters
        self.barcode_motif = "NNNNNNNNNNNNN"  # 13 N's as seen in reference
        self.min_read_length = 2000
        
        # Verify test data exists
        self.assertTrue(self.reference_file.exists(), f"Reference file not found: {self.reference_file}")
        self.assertTrue(self.bam_file.exists(), f"BAM file not found: {self.bam_file}")
        
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_reference_file_contains_barcode_motif(self):
        """Test that the reference file contains the expected barcode motif"""
        with open(self.reference_file, 'r') as f:
            content = f.read()
        
        # Check that the barcode motif appears in the reference
        self.assertIn(self.barcode_motif, content, 
                     f"Barcode motif '{self.barcode_motif}' not found in reference file")
        
        # Count occurrences (should be at least 2 for 5' and 3' UMIs)
        motif_count = content.count(self.barcode_motif)
        self.assertGreaterEqual(motif_count, 2, 
                               f"Expected at least 2 barcode motifs, found {motif_count}")
    
    def test_bam_file_is_valid(self):
        """Test that the BAM file is valid and contains reads"""
        import pysam
        
        with pysam.AlignmentFile(str(self.bam_file), "rb") as bam:
            read_count = 0
            for read in bam.fetch():
                read_count += 1
                if read_count >= 10:  # Just check first 10 reads
                    break
            
            self.assertGreater(read_count, 0, "BAM file contains no reads")
    
    def test_findUMI_function(self):
        """Test the findUMI function with real data"""
        # Call findUMI function
        result = findUMI(
            str(self.reference_file),
            str(self.bam_file),
            self.barcode_motif,
            self.min_read_length
        )
        
        # Unpack results
        raw_umi5_list, raw_umi3_list, motif5_index, motif3_index, raw_UMI_dict = result
        
        # Basic validation
        self.assertIsNotNone(raw_umi5_list, "5' UMI list should not be None")
        self.assertIsNotNone(raw_umi3_list, "3' UMI list should not be None")
        self.assertIsNotNone(motif5_index, "5' motif index should not be None")
        self.assertIsNotNone(motif3_index, "3' motif index should not be None")
        self.assertIsNotNone(raw_UMI_dict, "UMI dictionary should not be None")
        
        # Check motif indices are valid
        self.assertGreaterEqual(motif5_index, 0, "5' motif index should be >= 0")
        self.assertGreaterEqual(motif3_index, 0, "3' motif index should be >= 0")
        self.assertNotEqual(motif5_index, motif3_index, "5' and 3' motif indices should be different")
        
        # Check UMI lists
        if len(raw_umi5_list) > 0:
            self.assertEqual(len(raw_umi5_list[0]), len(self.barcode_motif), 
                           "5' UMI length should match barcode motif length")
        
        if len(raw_umi3_list) > 0:
            self.assertEqual(len(raw_umi3_list[0]), len(self.barcode_motif), 
                           "3' UMI length should match barcode motif length")
        
        # Check UMI dictionary
        if raw_UMI_dict:
            for read_name, umis in raw_UMI_dict.items():
                self.assertIsInstance(read_name, str, "Read name should be string")
                self.assertIsInstance(umis, list, "UMIs should be list")
                self.assertEqual(len(umis), 2, "Should have 2 UMIs (5' and 3')")
                
                # Check UMI sequences are valid
                for umi in umis:
                    if umi is not None:
                        self.assertIsInstance(umi, str, "UMI should be string")
                        self.assertIn(len(umi), [0, len(self.barcode_motif)], 
                                    f"UMI length should be 0 or {len(self.barcode_motif)}")
    
    def test_extract_UMIs_function(self):
        """Test the extract_UMIs function with real data"""
        # Call extract_UMIs function
        extract_UMIs(
            str(self.reference_file),
            str(self.bam_file),
            self.min_read_length,
            self.barcode_motif,
            str(self.output_dir)
        )
        
        # Check that output file was created
        expected_output_file = self.output_dir / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        self.assertTrue(expected_output_file.exists(), 
                       f"Expected output file not created: {expected_output_file}")
        
        # Read and validate the output CSV
        df = pd.read_csv(expected_output_file)
        
        # Check column structure
        expected_columns = ['Read Name', '5pUMI', '3pUMI']
        self.assertEqual(list(df.columns), expected_columns, 
                        f"Expected columns {expected_columns}, got {list(df.columns)}")
        
        # Check that we have some data
        self.assertGreater(len(df), 0, "Output CSV should contain data")
        
        # Show sample of extracted UMIs if enabled
        if SHOW_UMI_SAMPLES:
            # Find all UMI columns (exclude 'Read Name')
            umi_columns = [col for col in df.columns if col != 'Read Name']
            num_umis = len(umi_columns)
            
            print(f"\nExtracted {len(df)} reads with {num_umis} UMIs:")
            
            # Calculate column widths dynamically
            read_name_width = 30
            umi_width = 15
            total_width = read_name_width + (num_umis * umi_width) + (num_umis - 1) * 2  # +2 for spaces between columns
            
            print("=" * total_width)
            
            # Create header row
            header_parts = [f"{'Read Name':<{read_name_width}}"]
            for i, col in enumerate(umi_columns, 1):
                header_parts.append(f"UMI{i}")
            print(" ".join(header_parts))
            
            print("-" * total_width)
            
            # Show first 10 rows
            for i, (_, row) in enumerate(df.head(10).iterrows()):
                read_name = row['Read Name'][:read_name_width-3] + "..." if len(row['Read Name']) > read_name_width else row['Read Name']
                
                # Format row data
                row_parts = [f"{read_name:<{read_name_width}}"]
                for col in umi_columns:
                    umi_value = row[col] if pd.notna(row[col]) and row[col] != "" else "None"
                    row_parts.append(f"{umi_value:<{umi_width}}")
                
                print(" ".join(row_parts))
            
            if len(df) > 10:
                print(f"... and {len(df) - 10} more reads")
                
            # Show UMI column mapping
            print(f"\nUMI Column Mapping:")
            for i, col in enumerate(umi_columns, 1):
                print(f"  UMI{i} -> {col}")
        
        # Validate UMI sequences
        for _, row in df.iterrows():
            read_name = row['Read Name']
            umi_5p = row['5pUMI']
            umi_3p = row['3pUMI']
            
            # Read name should not be empty
            self.assertIsInstance(read_name, str, "Read name should be string")
            self.assertGreater(len(read_name), 0, "Read name should not be empty")
            
            # At least one UMI should be present
            has_umi = (pd.notna(umi_5p) and umi_5p != "") or (pd.notna(umi_3p) and umi_3p != "")
            self.assertTrue(has_umi, f"At least one UMI should be present for read {read_name}")
            
            # If UMIs are present, they should have correct length
            if pd.notna(umi_5p) and umi_5p != "":
                self.assertEqual(len(umi_5p), len(self.barcode_motif), 
                               f"5' UMI length should be {len(self.barcode_motif)}")
            
            if pd.notna(umi_3p) and umi_3p != "":
                self.assertEqual(len(umi_3p), len(self.barcode_motif), 
                               f"3' UMI length should be {len(self.barcode_motif)}")
    
    def test_extract_UMIs_with_different_min_read_length(self):
        """Test extract_UMIs with different minimum read length"""
        # Test with a higher minimum read length
        high_min_length = 5000
        
        extract_UMIs(
            str(self.reference_file),
            str(self.bam_file),
            high_min_length,
            self.barcode_motif,
            str(self.output_dir / "high_min_length")
        )
        
        # Test with a lower minimum read length
        low_min_length = 1000
        
        extract_UMIs(
            str(self.reference_file),
            str(self.bam_file),
            low_min_length,
            self.barcode_motif,
            str(self.output_dir / "low_min_length")
        )
        
        # Check that both output files were created
        high_output = self.output_dir / "high_min_length" / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        low_output = self.output_dir / "low_min_length" / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        
        self.assertTrue(high_output.exists(), "High min length output should exist")
        self.assertTrue(low_output.exists(), "Low min length output should exist")
        
        # Compare results - lower min length should have more reads
        high_df = pd.read_csv(high_output)
        low_df = pd.read_csv(low_output)
        
        self.assertGreaterEqual(len(low_df), len(high_df), 
                               "Lower min read length should produce more or equal reads")
    
    def test_extract_UMIs_with_different_barcode_motif(self):
        """Test extract_UMIs with a different barcode motif"""
        # Use a shorter motif for testing
        short_motif = "NNNNNNNNNNNN"  # 12 N's
        
        extract_UMIs(
            str(self.reference_file),
            str(self.bam_file),
            self.min_read_length,
            short_motif,
            str(self.output_dir / "short_motif")
        )
        
        output_file = self.output_dir / "short_motif" / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        
        if output_file.exists():
            df = pd.read_csv(output_file)
            
            # If we have data, check UMI lengths
            if len(df) > 0:
                for _, row in df.iterrows():
                    umi_5p = row['5pUMI']
                    umi_3p = row['3pUMI']
                    
                    if pd.notna(umi_5p) and umi_5p != "":
                        self.assertEqual(len(umi_5p), len(short_motif), 
                                       f"5' UMI length should be {len(short_motif)}")
                    
                    if pd.notna(umi_3p) and umi_3p != "":
                        self.assertEqual(len(umi_3p), len(short_motif), 
                                       f"3' UMI length should be {len(short_motif)}")
    
    def test_error_handling_invalid_files(self):
        """Test error handling with invalid input files"""
        # Test with non-existent reference file
        with self.assertRaises(Exception):
            findUMI(
                "nonexistent_reference.fa",
                str(self.bam_file),
                self.barcode_motif,
                self.min_read_length
            )
        
        # Test with non-existent BAM file
        with self.assertRaises(Exception):
            findUMI(
                str(self.reference_file),
                "nonexistent_bam.bam",
                self.barcode_motif,
                self.min_read_length
            )
    
    def test_output_file_format(self):
        """Test that the output CSV file has the correct format"""
        extract_UMIs(
            str(self.reference_file),
            str(self.bam_file),
            self.min_read_length,
            self.barcode_motif,
            str(self.output_dir)
        )
        
        output_file = self.output_dir / "B145-XTMMF9-12_rawUMIs_read_names.csv"
        
        # Check file exists and is readable
        self.assertTrue(output_file.exists())
        self.assertTrue(output_file.is_file())
        
        # Read and validate CSV format
        df = pd.read_csv(output_file)
        
        # Check basic structure
        self.assertGreater(len(df), 0, "CSV should contain data")
        self.assertEqual(len(df.columns), 3, "Should have 3 columns")
        
        # Check data types
        self.assertTrue(df['Read Name'].dtype == 'object', "Read Name should be string/object")
        self.assertTrue(df['5pUMI'].dtype == 'object', "5pUMI should be string/object")
        self.assertTrue(df['3pUMI'].dtype == 'object', "3pUMI should be string/object")
        
        # Check for any completely empty rows
        empty_rows = df[df['Read Name'].isna() & df['5pUMI'].isna() & df['3pUMI'].isna()]
        self.assertEqual(len(empty_rows), 0, "Should not have completely empty rows")


def simple_umi_preview(save_csv_dir=None):
    """Simple function to just show UMI extraction preview without test framework"""
    import tempfile
    import shutil
    
    # Set up paths
    test_data_dir = Path(__file__).parent.parent / "scripts" / "test-data" / "test-umi_map"
    reference_file = test_data_dir / "reference_sequence" / "LP249-Gldx-CircAmp.fasta"
    bam_file = test_data_dir / "sorted_reads" / "LP249-Gldx-CircAmp-sorted.bam"
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Run extraction
        extract_UMIs(
            str(reference_file),
            str(bam_file),
            1000,  # min_read_length
            "NNNNNNNNNNNN",  # barcode_motif
            temp_dir
        )
        
        # Read and display results
        output_file = Path(temp_dir) / "LP249-Gldx-CircAmp-sorted_rawUMIs_read_names.csv"
        if output_file.exists():
            df = pd.read_csv(output_file)
            
            # Find all UMI columns (exclude 'Read Name')
            umi_columns = [col for col in df.columns if col != 'Read Name']
            num_umis = len(umi_columns)
            
            print(f"Extracted {len(df)} reads with {num_umis} UMIs:")
            
            # Calculate column widths dynamically
            read_name_width = 30
            umi_width = 15
            total_width = read_name_width + (num_umis * umi_width) + (num_umis - 1) * 2  # +2 for spaces between columns
            
            print("=" * total_width)
            
            # Create header row
            header_parts = [f"{'Read Name':<{read_name_width}}"]
            for i, col in enumerate(umi_columns, 1):
                header_parts.append(f"UMI{i}")
            print(" ".join(header_parts))
            
            print("-" * total_width)
            
            # Show first 10 rows
            for i, (_, row) in enumerate(df.head(10).iterrows()):
                read_name = row['Read Name'][:read_name_width-3] + "..." if len(row['Read Name']) > read_name_width else row['Read Name']
                
                # Format row data
                row_parts = [f"{read_name:<{read_name_width}}"]
                for col in umi_columns:
                    umi_value = row[col] if pd.notna(row[col]) and row[col] != "" else "None"
                    row_parts.append(f"{umi_value:<{umi_width}}")
                
                print(" ".join(row_parts))
            
            if len(df) > 10:
                print(f"... and {len(df) - 10} more reads")
                
            # Show UMI column mapping
            print(f"\nUMI Column Mapping:")
            for i, col in enumerate(umi_columns, 1):
                print(f"  UMI{i} -> {col}")
            
            # Save CSV if requested
            if save_csv_dir:
                save_dir = Path(save_csv_dir)
                save_dir.mkdir(parents=True, exist_ok=True)
                saved_file = save_dir / output_file.name
                shutil.copy2(output_file, saved_file)
                print(f"\nCSV file saved to: {saved_file}")
                
        else:
            print("No output file created")
            
    finally:
        # Clean up
        shutil.rmtree(temp_dir, ignore_errors=True)


def main():
    """Main function to run the extract tests with detailed output"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run BELT extract function tests with real data')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Verbose output showing individual test results')
    parser.add_argument('--failfast', '-f', action='store_true',
                       help='Stop on first failure')
    parser.add_argument('--buffer', '-b', action='store_true',
                       help='Buffer stdout and stderr during test runs')
    parser.add_argument('--pattern', '-p', type=str,
                       help='Pattern to match test method names')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress UMI sample output')
    parser.add_argument('--simple', '-s', action='store_true',
                       help='Simple mode: only show UMI preview, no test output')
    parser.add_argument('--save-csv', type=str, metavar='DIR',
                       help='Save the CSV file to the specified directory')
    
    args = parser.parse_args()
    
    # Handle simple mode
    if args.simple:
        simple_umi_preview(save_csv_dir=args.save_csv)
        return 0
    
    # Set global variable for UMI sample display
    global SHOW_UMI_SAMPLES
    SHOW_UMI_SAMPLES = not args.quiet
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestExtractReal)
    
    # Filter tests by pattern if provided
    if args.pattern:
        suite = unittest.TestSuite()
        for test in loader.loadTestsFromTestCase(TestExtractReal):
            if args.pattern.lower() in test._testMethodName.lower():
                suite.addTest(test)
    
    # Create test runner with options
    verbosity = 2 if args.verbose else 1
    runner = unittest.TextTestRunner(
        verbosity=verbosity,
        failfast=args.failfast,
        buffer=args.buffer,
        stream=open(os.devnull, 'w') if not args.verbose else None
    )
    
    print("=" * 70)
    print("BELT Extract Function Tests with Real Data")
    print("=" * 70)
    print(f"Running {suite.countTestCases()} test(s)...")
    print()
    
    # Run tests
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped) if hasattr(result, 'skipped') else 0}")
    
    if result.failures:
        print("\nFAILURES:")
        for test, traceback in result.failures:
            error_msg = traceback.split('AssertionError: ')[-1].split('\n')[0]
            print(f"  - {test}: {error_msg}")
    
    if result.errors:
        print("\nERRORS:")
        for test, traceback in result.errors:
            error_msg = traceback.split('\n')[-2]
            print(f"  - {test}: {error_msg}")
    
    # Exit with appropriate code
    exit_code = 0 if result.wasSuccessful() else 1
    print(f"\nExit code: {exit_code}")
    return exit_code


if __name__ == '__main__':
    exit(main())
