#!/usr/bin/env python3
"""
Comprehensive test suite for UMI2 majority filtering in long_read_umi_grouping_fast
Tests the new umi2_majority_threshold parameter and read-level filtering functionality.
"""
import os
import sys
import pytest
import pandas as pd
import pysam
import json
from pathlib import Path
import tempfile

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_group import long_read_umi_grouping_fast


def create_test_bam(bam_path, read_names):
    """Create a minimal test BAM file with the specified read names"""
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }
    
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for read_name in read_names:
            a = pysam.AlignedSegment()
            a.query_name = read_name
            a.query_sequence = "ACGTACGTACGT"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigar = [(0, 12)]  # 12M
            a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIII")
            outf.write(a)


def create_test_umi_csv(csv_path, data):
    """
    Create a test UMI CSV file.
    data should be a list of dicts with keys: 'Read Name', 'UMI1', 'UMI2', etc.
    """
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)


class TestUMI2MajorityFilter:
    """Test suite for UMI2 majority filtering functionality"""
    
    def test_default_50_percent_threshold(self, tmp_path):
        """
        Test with default 50% threshold.
        UMI1-A has 10 reads: 6 with UMI2-X, 4 with UMI2-Y
        Should keep UMI1-A and only the 6 reads with UMI2-X
        """
        # Create test data
        umi_data = []
        read_names = []
        
        # UMI1-A group: 6 reads with UMI2-X (majority), 4 with UMI2-Y
        for i in range(6):
            read_name = f"read_A_X_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        
        for i in range(4):
            read_name = f"read_A_Y_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "YYYY"})
        
        # UMI1-B group: 3 reads with UMI2-Z (all same)
        for i in range(3):
            read_name = f"read_B_Z_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "BBBB", "UMI2": "ZZZZ"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run function with default threshold (0.5 = 50%)
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=3,  # minimum reads per group
            threads=1,
            umi2_majority_threshold=0.5  # default
        )
        
        # Verify results
        grouped_bam = output_dir / "test.grouped.bam"
        assert grouped_bam.exists(), "Grouped BAM should be created"
        
        # Check BAM contents
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        # Should have 6 reads from UMI1-A (only UMI2-X) + 3 reads from UMI1-B (all UMI2-Z) = 9 reads
        assert len(reads_in_bam) == 9, f"Expected 9 reads, got {len(reads_in_bam)}"
        
        # Verify only UMI2-X reads from UMI1-A are kept
        umi_a_reads = [r for r in read_names_in_bam if r.startswith("read_A")]
        assert all("_X_" in r for r in umi_a_reads), "Only UMI2-X reads should be kept from UMI1-A"
        assert len(umi_a_reads) == 6, f"Should have 6 reads from UMI1-A, got {len(umi_a_reads)}"
        
        # Verify UMI1-B reads are all kept
        umi_b_reads = [r for r in read_names_in_bam if r.startswith("read_B")]
        assert len(umi_b_reads) == 3, f"Should have 3 reads from UMI1-B, got {len(umi_b_reads)}"
        
        # Check log file
        assert os.path.exists(log_path), "Log file should be created"
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["umi2_majority_filter_applied"] == True
        assert log_data["umi2_majority_threshold"] == 0.5
        assert log_data["reads_before_umi2_filter"] == 13  # 10 + 3
        assert log_data["reads_after_umi2_filter"] == 9   # 6 + 3
        assert log_data["reads_filtered_by_umi2"] == 4     # 4 UMI2-Y reads removed
        assert log_data["groups_passing_umi2_majority"] == 2  # Both UMI1-A and UMI1-B pass
        
        print("✓ Test passed: Default 50% threshold works correctly")
    
    def test_70_percent_strict_threshold(self, tmp_path):
        """
        Test with stricter 70% threshold.
        UMI1-A has 10 reads: 6 with UMI2-X (60%), 4 with UMI2-Y
        Should be excluded because 60% < 70%
        UMI1-B has 10 reads: 8 with UMI2-Z (80%), 2 with UMI2-W
        Should pass and keep only the 8 UMI2-Z reads
        """
        umi_data = []
        read_names = []
        
        # UMI1-A: 60% majority (should be excluded at 70% threshold)
        for i in range(6):
            read_name = f"read_A_X_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        for i in range(4):
            read_name = f"read_A_Y_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "YYYY"})
        
        # UMI1-B: 80% majority (should pass at 70% threshold)
        for i in range(8):
            read_name = f"read_B_Z_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "BBBB", "UMI2": "ZZZZ"})
        for i in range(2):
            read_name = f"read_B_W_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "BBBB", "UMI2": "WWWW"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run with 70% threshold
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.7  # 70%
        )
        
        # Verify results
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        # Should only have 8 reads from UMI1-B (UMI2-Z only)
        assert len(reads_in_bam) == 8, f"Expected 8 reads, got {len(reads_in_bam)}"
        assert all("read_B_Z_" in r for r in read_names_in_bam), "Only UMI1-B UMI2-Z reads should be kept"
        
        # Check log
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["umi2_majority_threshold"] == 0.7
        assert log_data["groups_passing_umi2_majority"] == 1  # Only UMI1-B passes
        assert log_data["reads_after_umi2_filter"] == 8  # 8 reads from UMI1-B
        
        print("✓ Test passed: 70% strict threshold correctly filters groups")
    
    def test_30_percent_loose_threshold(self, tmp_path):
        """
        Test with loose 30% threshold.
        UMI1-A has 10 reads: 3 with UMI2-X (30%), 7 with other UMIs
        Should just barely pass and keep only the 3 UMI2-X reads
        """
        umi_data = []
        read_names = []
        
        # UMI1-A: exactly 30% majority (should pass at 30% threshold)
        for i in range(3):
            read_name = f"read_A_X_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        for i in range(7):
            read_name = f"read_A_other_{i}"
            read_names.append(read_name)
            # Different UMI2 for each
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": f"UMI{i}"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run with 30% threshold
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.3  # 30%
        )
        
        # Verify results
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        # Should have only 3 reads with UMI2-X
        assert len(reads_in_bam) == 3, f"Expected 3 reads, got {len(reads_in_bam)}"
        assert all("read_A_X_" in r for r in read_names_in_bam), "Only UMI2-X reads should be kept"
        
        # Check log
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["reads_after_umi2_filter"] == 3
        assert log_data["reads_filtered_by_umi2"] == 7  # 7 non-majority reads filtered
        
        print("✓ Test passed: 30% loose threshold works correctly")
    
    def test_no_umi2_column(self, tmp_path):
        """
        Test when UMI2 column doesn't exist - should work normally without filtering
        """
        umi_data = []
        read_names = []
        
        # Only UMI1, no UMI2
        for i in range(10):
            read_name = f"read_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run function
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.5
        )
        
        # Verify results - all reads should be kept
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
        
        assert len(reads_in_bam) == 10, f"Expected 10 reads, got {len(reads_in_bam)}"
        
        # Check log - UMI2 filter should not be applied
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["umi2_majority_filter_applied"] == False
        assert log_data["reads_before_umi2_filter"] == 0
        assert log_data["reads_after_umi2_filter"] == 0
        
        print("✓ Test passed: Works correctly without UMI2 column")
    
    def test_grouping_by_umi2(self, tmp_path):
        """
        Test when grouping by UMI2 instead of UMI1 - UMI2 filter should not apply
        """
        umi_data = []
        read_names = []
        
        # Group by UMI2
        for i in range(10):
            read_name = f"read_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": f"UMI{i}", "UMI2": "XXXX"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run function grouping by UMI2
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI2",  # Group by UMI2, not UMI1
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.5
        )
        
        # Verify UMI2 filter was not applied
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["umi2_majority_filter_applied"] == False
        assert log_data["grouping"] == "UMI2"
        
        print("✓ Test passed: UMI2 filter not applied when grouping by UMI2")
    
    def test_exact_50_percent_split(self, tmp_path):
        """
        Test edge case: exactly 50% split with default 50% threshold
        UMI1-A has 10 reads: 5 with UMI2-X, 5 with UMI2-Y
        With threshold >= 0.5, the first 5 should be considered majority
        """
        umi_data = []
        read_names = []
        
        # Exactly 50/50 split
        for i in range(5):
            read_name = f"read_X_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        for i in range(5):
            read_name = f"read_Y_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "YYYY"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run with default 50% threshold
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.5
        )
        
        # Verify results
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        # At exactly 50%, should pass and keep the 5 reads with the first majority UMI2
        assert len(reads_in_bam) == 5, f"Expected 5 reads, got {len(reads_in_bam)}"
        
        # Check log
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["groups_passing_umi2_majority"] == 1
        assert log_data["reads_after_umi2_filter"] == 5
        
        print("✓ Test passed: Exact 50% split handled correctly")
    
    def test_multiple_groups_different_outcomes(self, tmp_path):
        """
        Test multiple UMI1 groups with different UMI2 distributions
        """
        umi_data = []
        read_names = []
        
        # UMI1-A: 80% majority (pass)
        for i in range(8):
            read_name = f"read_A_pass_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        for i in range(2):
            read_name = f"read_A_filtered_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "YYYY"})
        
        # UMI1-B: 30% majority (fail at 50% threshold)
        for i in range(3):
            read_name = f"read_B_maj_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "BBBB", "UMI2": "ZZZZ"})
        for i in range(7):
            read_name = f"read_B_other_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "BBBB", "UMI2": "WWWW"})
        
        # UMI1-C: 100% unanimous (pass)
        for i in range(6):
            read_name = f"read_C_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "CCCC", "UMI2": "QQQQ"})
        
        # UMI1-D: below threshold count (only 2 reads, threshold is 3)
        for i in range(2):
            read_name = f"read_D_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "DDDD", "UMI2": "PPPP"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run with 50% threshold and minimum 3 reads
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=3,
            threads=1,
            umi2_majority_threshold=0.5
        )
        
        # Verify results
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        # Should have: 8 from A (XXXX only), 0 from B (failed majority), 6 from C, 0 from D (below threshold)
        expected_count = 8 + 6
        assert len(reads_in_bam) == expected_count, f"Expected {expected_count} reads, got {len(reads_in_bam)}"
        
        # Verify correct reads are kept
        assert sum(1 for r in read_names_in_bam if "read_A_pass_" in r) == 8
        assert sum(1 for r in read_names_in_bam if "read_A_filtered_" in r) == 0
        assert sum(1 for r in read_names_in_bam if "read_B_" in r) == 0
        assert sum(1 for r in read_names_in_bam if "read_C_" in r) == 6
        assert sum(1 for r in read_names_in_bam if "read_D_" in r) == 0
        
        # Check log
        with open(log_path, 'r') as f:
            log_data = json.load(f)
        
        assert log_data["total_unique_umis"] == 4  # A, B, C, D
        assert log_data["groups_passing_threshold"] == 3  # A, B, C (D below threshold)
        assert log_data["groups_passing_umi2_majority"] == 2  # A, C (B failed majority)
        assert log_data["reads_before_umi2_filter"] == 10 + 10 + 6  # A + B + C
        assert log_data["reads_after_umi2_filter"] == 8 + 6  # A (filtered) + C
        assert log_data["reads_filtered_by_umi2"] == 2 + 10  # 2 from A, all 10 from B
        
        print("✓ Test passed: Multiple groups with different outcomes handled correctly")
    
    def test_umi2_with_na_values(self, tmp_path):
        """
        Test handling of NA/NoMatch/NoBarcode values in UMI2
        """
        umi_data = []
        read_names = []
        
        # UMI1-A: 6 valid UMI2-X, 2 NoMatch, 2 NoBarcode
        for i in range(6):
            read_name = f"read_valid_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "XXXX"})
        for i in range(2):
            read_name = f"read_nomatch_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "NoMatch"})
        for i in range(2):
            read_name = f"read_nobarcode_{i}"
            read_names.append(read_name)
            umi_data.append({"Read Name": read_name, "UMI1": "AAAA", "UMI2": "NoBarcode"})
        
        # Create test files
        bam_path = tmp_path / "test.bam"
        csv_path = tmp_path / "corrected_umis.csv"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        
        create_test_bam(bam_path, read_names)
        create_test_umi_csv(csv_path, umi_data)
        
        # Run function
        log_path = long_read_umi_grouping_fast(
            alignedReads=str(bam_path),
            correctedUMIs=str(csv_path),
            outputDir=str(output_dir),
            grouping="UMI1",
            threshold=5,
            threads=1,
            umi2_majority_threshold=0.5
        )
        
        # Verify - should keep only the 6 valid UMI2-X reads
        grouped_bam = output_dir / "test.grouped.bam"
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            reads_in_bam = list(bam.fetch(until_eof=True))
            read_names_in_bam = [r.query_name for r in reads_in_bam]
        
        assert len(reads_in_bam) == 6, f"Expected 6 reads, got {len(reads_in_bam)}"
        assert all("read_valid_" in r for r in read_names_in_bam)
        
        print("✓ Test passed: NA/NoMatch/NoBarcode values handled correctly")


def test_comprehensive_scenarios(tmp_path):
    """Run all test scenarios"""
    test_suite = TestUMI2MajorityFilter()
    
    print("\n" + "="*80)
    print("COMPREHENSIVE UMI2 MAJORITY FILTER TEST SUITE")
    print("="*80 + "\n")
    
    test_suite.test_default_50_percent_threshold(tmp_path / "test1")
    test_suite.test_70_percent_strict_threshold(tmp_path / "test2")
    test_suite.test_30_percent_loose_threshold(tmp_path / "test3")
    test_suite.test_no_umi2_column(tmp_path / "test4")
    test_suite.test_grouping_by_umi2(tmp_path / "test5")
    test_suite.test_exact_50_percent_split(tmp_path / "test6")
    test_suite.test_multiple_groups_different_outcomes(tmp_path / "test7")
    test_suite.test_umi2_with_na_values(tmp_path / "test8")
    
    print("\n" + "="*80)
    print("ALL TESTS PASSED!")
    print("="*80 + "\n")


if __name__ == "__main__":
    # Can be run directly for quick testing
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        test_comprehensive_scenarios(Path(tmpdir))



