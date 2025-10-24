#!/usr/bin/env python3
"""
Pytest test file for BELT extract function
"""
import os
import sys
import pytest
import pandas as pd
from pathlib import Path
import tempfile
import shutil

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_map import umi_extraction_map  # type: ignore

@pytest.fixture
def test_config_five_motifs():
    """Fixture for test configuration with five motifs"""
    return {
        'min_read_length': 2000,
        'barcode_motifs': ["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"],
        'test_data_dir': "test-data/test_umi_map",
        'reference_file': "test-data/test_umi_map/extract_umi_5.fa",
        'bam_file': "test-data/test_umi_map/extract_umi_5-large.bam",
        'whitelist': "test-data/test_umi_map/extract_umi_5_whitelist.csv",
        'output_dir': "test-data/test_umi_map/output",

    }

@pytest.fixture
def test_paths_five_motifs(test_config_five_motifs):
    """Fixture to build absolute paths for five motif test"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_five_motifs['reference_file'],
        'bam_file': base_dir / test_config_five_motifs['bam_file'],
        'whitelist_file': base_dir / test_config_five_motifs['whitelist'],
        'output_dir': base_dir / test_config_five_motifs['output_dir'],
        'base_dir': base_dir
    }

def test_umi_extraction_map_five_motifs(test_config_five_motifs, test_paths_five_motifs):
    """Test the umi_extraction_map function with five motifs"""
    reference_file = test_paths_five_motifs['reference_file']
    bam_file = test_paths_five_motifs['bam_file']
    whitelist_file = test_paths_five_motifs['whitelist_file']
    output_dir = test_paths_five_motifs['output_dir']
    
    # Verify test files exist
    assert reference_file.exists(), f"Reference file not found: {reference_file}"
    assert bam_file.exists(), f"BAM file not found: {bam_file}"

    
    # Run the umi_extraction_map function (convert Path objects to strings)
    umi_extraction_map(
        reference=str(reference_file),
        alignedReads=str(bam_file), 
        minReadLength=test_config_five_motifs['min_read_length'], 
        barcodeMotifs=test_config_five_motifs['barcode_motifs'], 
        whitelist=str(whitelist_file), 
        outputDir=str(output_dir)
    )

    print(f"\n✓ Test passed: UMI extraction completed")
    print(f"Reference: {reference_file}")
    print(f"BAM: {bam_file}")
    print(f"Whitelist: {whitelist_file}")
    print(f"Output: {output_dir}")


# class TestExtract(unittest.TestCase):
#     def test_umi_map(self):
#         # Test extracting UMIs from a short read BAM file
#         base_dir = Path(__file__).parent.parent / "test-data" / "test_umi_map"
#         reference = str(base_dir / "LP250.fa")
#         alignment = str(base_dir / "LP250-long.bam")
#         minReadLength = 2000
#         motif = "NHHHHHHHHHHHHHM"
#         whitelist = str(base_dir / "LP250-short_correctedUMIs.csv")
#         outputDir = str(base_dir)
#         analysisOut = str(base_dir)

#         try:
#             umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, outputDir, analysisOut)
#         except Exception as e:
#             self.fail(f"extract_umis raised an exception: {e}")

#     def test_umi_map_sample2(self):
#         # Test extracting UMIs from a short read BAM file
#         base_dir = Path(__file__).parent.parent / "test-data" / "test_umi_map"
#         reference = str(base_dir / "LP250.fa")
#         alignment = str(base_dir / "LP250v1.bam")
#         minReadLength = 2000
#         motif = "NHHHHHHHHHHHHHM"
#         whitelist = str(base_dir / "LP250-short_correctedUMIs.csv")
#         outputDir = str(base_dir)
#         analysisOut = str(base_dir)

#         try:
#             umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, outputDir, analysisOut)
#         except Exception as e:
#             self.fail(f"extract_umis raised an exception: {e}")
