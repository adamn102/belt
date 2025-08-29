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


if __name__ == '__main__':
    unittest.main()
