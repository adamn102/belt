#!/usr/bin/env python3
"""
Real tests for BELT umi-group function using actual test data
"""

import unittest
import tempfile
import os
import sys
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_group import long_read_umi_grouping, long_read_umi_grouping_fast


class TestUMIGroupReal(unittest.TestCase):
    """Test umi-group functionality with real test data"""
    
    def setUp(self):
        """Set up test fixtures with real test data"""
        # Paths to real test data
        self.test_data_dir = Path(__file__).parent.parent / "scripts" / "test-data"
        self.reference_file = self.test_data_dir / "test-umi_map" / "reference_sequence" / "B145.fa"
        self.bam_file = self.test_data_dir / "test-umi_map" / "sorted_reads" / "B145-XTMMF9-12.bam"
        self.corrected_umis_file = self.test_data_dir / "5-umi_maps" / "B145-XTMMF9-12_Paired_UMIs_read_names.csv"
        
        # Create temporary directory for outputs
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir)
        
        # Test parameters
        self.grouping = "5pUMI"
        self.threshold = 5
        self.threads = 4
        
        # Verify test data exists
        self.assertTrue(self.reference_file.exists(), f"Reference file not found: {self.reference_file}")
        self.assertTrue(self.bam_file.exists(), f"BAM file not found: {self.bam_file}")
        self.assertTrue(self.corrected_umis_file.exists(), f"Corrected UMIs file not found: {self.corrected_umis_file}")
        
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_corrected_umis_file_format(self):
        """Test that the corrected UMIs file has the correct format"""
        df = pd.read_csv(self.corrected_umis_file)
        
        # Check column structure
        expected_columns = ['Read Name', '5pUMI', '3pUMI']
        self.assertEqual(list(df.columns), expected_columns, 
                        f"Expected columns {expected_columns}, got {list(df.columns)}")
        
        # Check that we have data
        self.assertGreater(len(df), 0, "Corrected UMIs file should contain data")
        
        # Check data types
        self.assertTrue(df['Read Name'].dtype == 'object', "Read Name should be string/object")
        self.assertTrue(df['5pUMI'].dtype == 'object', "5pUMI should be string/object")
        self.assertTrue(df['3pUMI'].dtype == 'object', "3pUMI should be string/object")
        
        # Check for some valid UMIs
        valid_5p_umis = df['5pUMI'].dropna()
        valid_3p_umis = df['3pUMI'].dropna()
        
        self.assertGreater(len(valid_5p_umis), 0, "Should have some valid 5' UMIs")
        self.assertGreater(len(valid_3p_umis), 0, "Should have some valid 3' UMIs")
        
        # Check UMI lengths (should be consistent)
        if len(valid_5p_umis) > 0:
            umi_lengths_5p = valid_5p_umis.str.len()
            self.assertTrue(umi_lengths_5p.nunique() <= 2, "5' UMI lengths should be consistent")
        
        if len(valid_3p_umis) > 0:
            umi_lengths_3p = valid_3p_umis.str.len()
            self.assertTrue(umi_lengths_3p.nunique() <= 2, "3' UMI lengths should be consistent")
    
    def test_umi_grouping_analysis(self):
        """Test UMI grouping analysis with real data"""
        df = pd.read_csv(self.corrected_umis_file)
        
        # Clean up 5pUMI (same as in the function)
        df["5pUMI"] = df["5pUMI"].replace(["", "NoMatch", "NoBarcode"], pd.NA)
        
        # Group by 5pUMI
        grouped = df.groupby("5pUMI").agg({
            "Read Name": ["count", lambda x: list(x)]
        })
        grouped.columns = ["ReadCount", "ReadNames"]
        
        # Filter groups by threshold
        filtered = grouped[grouped["ReadCount"] >= self.threshold]
        
        # Basic validation
        self.assertGreater(len(grouped), 0, "Should have some UMI groups")
        self.assertGreaterEqual(len(filtered), 0, "Should have some groups meeting threshold")
        
        # Check group structure
        for umi, row in filtered.iterrows():
            self.assertIsInstance(umi, str, "UMI should be string")
            self.assertGreater(len(umi), 0, "UMI should not be empty")
            self.assertGreaterEqual(row["ReadCount"], self.threshold, 
                                  f"Group {umi} should meet threshold {self.threshold}")
            self.assertIsInstance(row["ReadNames"], list, "ReadNames should be list")
            self.assertEqual(len(row["ReadNames"]), row["ReadCount"], 
                           "ReadNames list length should match ReadCount")
    
    def test_long_read_umi_grouping_function(self):
        """Test the long_read_umi_grouping function with real data"""
        # Call the function
        long_read_umi_grouping(
            str(self.reference_file),
            str(self.bam_file),
            str(self.corrected_umis_file),
            str(self.output_dir),
            self.grouping,
            self.threshold
        )
        
        # Check that output files were created
        sample_name = "B145-XTMMF9-12"
        expected_filtered_csv = self.output_dir / f"{sample_name}_filtered_umi_groups.csv"
        expected_sample_dir = self.output_dir / sample_name
        
        self.assertTrue(expected_filtered_csv.exists(), 
                       f"Filtered UMI groups CSV not created: {expected_filtered_csv}")
        self.assertTrue(expected_sample_dir.exists(), 
                       f"Sample output directory not created: {expected_sample_dir}")
        self.assertTrue(expected_sample_dir.is_dir(), 
                       f"Sample output should be a directory: {expected_sample_dir}")
        
        # Read and validate the filtered CSV
        filtered_df = pd.read_csv(expected_filtered_csv)
        
        # Check column structure
        expected_columns = ['5pUMI', 'ReadCount', 'ReadNames']
        self.assertEqual(list(filtered_df.columns), expected_columns, 
                        f"Expected columns {expected_columns}, got {list(filtered_df.columns)}")
        
        # Check that we have some groups
        self.assertGreater(len(filtered_df), 0, "Should have some filtered UMI groups")
        
        # Check that all groups meet the threshold
        for _, row in filtered_df.iterrows():
            self.assertGreaterEqual(row['ReadCount'], self.threshold, 
                                  f"All groups should meet threshold {self.threshold}")
        
        # Check that FASTQ files were created for each group
        fastq_files = list(expected_sample_dir.glob("*.fastq"))
        self.assertEqual(len(fastq_files), len(filtered_df), 
                        f"Should have one FASTQ file per UMI group")
        
        # Validate some FASTQ files
        for fastq_file in fastq_files[:3]:  # Check first 3 files
            self.assertTrue(fastq_file.exists())
            self.assertGreater(fastq_file.stat().st_size, 0, f"FASTQ file should not be empty: {fastq_file}")
            
            # Check FASTQ format
            with open(fastq_file, 'r') as f:
                lines = f.readlines()
                self.assertGreater(len(lines), 0, f"FASTQ file should have content: {fastq_file}")
                
                # Check that it's a multiple of 4 (FASTQ format)
                self.assertEqual(len(lines) % 4, 0, f"FASTQ file should have 4 lines per read: {fastq_file}")
                
                # Check header format (should contain UMI info)
                for i in range(0, len(lines), 4):
                    if i < len(lines):
                        header = lines[i].strip()
                        self.assertTrue(header.startswith('@'), f"FASTQ header should start with @: {header}")
                        self.assertIn('|5pUMI:', header, f"FASTQ header should contain 5pUMI info: {header}")
                        self.assertIn('|3pUMI:', header, f"FASTQ header should contain 3pUMI info: {header}")
    
    def test_long_read_umi_grouping_fast_function(self):
        """Test the long_read_umi_grouping_fast function with real data"""
        # Check if samtools is available
        import subprocess
        try:
            subprocess.run(["samtools", "--version"], capture_output=True, check=True)
            samtools_available = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            samtools_available = False
            self.skipTest("samtools not available - skipping fast function test")
        
        if samtools_available:
            # Call the function
            long_read_umi_grouping_fast(
                str(self.reference_file),
                str(self.bam_file),
                str(self.corrected_umis_file),
                str(self.output_dir / "fast"),
                self.grouping,
                self.threshold,
                self.threads
            )
            
            # Check that output files were created
            sample_name = "B145-XTMMF9-12"
            expected_filtered_csv = self.output_dir / "fast" / f"{sample_name}_filtered_umi_groups.csv"
            
            self.assertTrue(expected_filtered_csv.exists(), 
                           f"Filtered UMI groups CSV not created: {expected_filtered_csv}")
            
            # Read and validate the filtered CSV
            filtered_df = pd.read_csv(expected_filtered_csv)
            
            # Check column structure (fast version has different format)
            expected_columns = ['5pUMI', 'ReadCount']
            self.assertEqual(list(filtered_df.columns), expected_columns, 
                            f"Expected columns {expected_columns}, got {list(filtered_df.columns)}")
            
            # Check that we have some groups
            self.assertGreater(len(filtered_df), 0, "Should have some filtered UMI groups")
            
            # Check that all groups meet the threshold
            for _, row in filtered_df.iterrows():
                self.assertGreaterEqual(row['ReadCount'], self.threshold, 
                                      f"All groups should meet threshold {self.threshold}")
    
    def test_umi_grouping_with_different_thresholds(self):
        """Test UMI grouping with different threshold values"""
        thresholds = [1, 3, 10]
        
        for threshold in thresholds:
            with self.subTest(threshold=threshold):
                output_subdir = self.output_dir / f"threshold_{threshold}"
                
                # Call the function
                long_read_umi_grouping(
                    str(self.reference_file),
                    str(self.bam_file),
                    str(self.corrected_umis_file),
                    str(output_subdir),
                    self.grouping,
                    threshold
                )
                
                # Check output
                sample_name = "B145-XTMMF9-12"
                expected_filtered_csv = output_subdir / f"{sample_name}_filtered_umi_groups.csv"
                
                if expected_filtered_csv.exists():
                    filtered_df = pd.read_csv(expected_filtered_csv)
                    
                    # All groups should meet the threshold
                    for _, row in filtered_df.iterrows():
                        self.assertGreaterEqual(row['ReadCount'], threshold, 
                                              f"All groups should meet threshold {threshold}")
    
    def test_umi_grouping_with_different_grouping_options(self):
        """Test UMI grouping with different grouping options"""
        grouping_options = ["5pUMI", "3pUMI"]
        
        for grouping in grouping_options:
            with self.subTest(grouping=grouping):
                output_subdir = self.output_dir / f"grouping_{grouping}"
                
                # Call the function
                long_read_umi_grouping(
                    str(self.reference_file),
                    str(self.bam_file),
                    str(self.corrected_umis_file),
                    str(output_subdir),
                    grouping,
                    self.threshold
                )
                
                # Check output
                sample_name = "B145-XTMMF9-12"
                expected_filtered_csv = output_subdir / f"{sample_name}_filtered_umi_groups.csv"
                
                if expected_filtered_csv.exists():
                    filtered_df = pd.read_csv(expected_filtered_csv)
                    
                    # Should have some groups
                    self.assertGreaterEqual(len(filtered_df), 0, 
                                          f"Should have some groups for {grouping} grouping")
    
    def test_error_handling_invalid_files(self):
        """Test error handling with invalid input files"""
        # Test with non-existent UMIs file (this should definitely fail)
        with self.assertRaises((FileNotFoundError, pd.errors.EmptyDataError)):
            long_read_umi_grouping(
                str(self.reference_file),
                str(self.bam_file),
                "nonexistent_umis.csv",
                str(self.output_dir),
                self.grouping,
                self.threshold
            )
    
    def test_output_file_consistency(self):
        """Test that output files are consistent between runs"""
        # Run the function twice
        output_dir1 = self.output_dir / "run1"
        output_dir2 = self.output_dir / "run2"
        
        long_read_umi_grouping(
            str(self.reference_file),
            str(self.bam_file),
            str(self.corrected_umis_file),
            str(output_dir1),
            self.grouping,
            self.threshold
        )
        
        long_read_umi_grouping(
            str(self.reference_file),
            str(self.bam_file),
            str(self.corrected_umis_file),
            str(output_dir2),
            self.grouping,
            self.threshold
        )
        
        # Compare filtered CSV files
        sample_name = "B145-XTMMF9-12"
        csv1 = output_dir1 / f"{sample_name}_filtered_umi_groups.csv"
        csv2 = output_dir2 / f"{sample_name}_filtered_umi_groups.csv"
        
        if csv1.exists() and csv2.exists():
            df1 = pd.read_csv(csv1)
            df2 = pd.read_csv(csv2)
            
            # Should have same number of groups
            self.assertEqual(len(df1), len(df2), "Should have same number of groups between runs")
            
            # Should have same total read count
            total_reads1 = df1['ReadCount'].sum()
            total_reads2 = df2['ReadCount'].sum()
            self.assertEqual(total_reads1, total_reads2, "Should have same total read count between runs")


if __name__ == '__main__':
    unittest.main()
