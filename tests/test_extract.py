#!/usr/bin/env python3
"""
Pytest test file for BELT extract function
"""
import os
import sys
import pytest
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from extract import extract_umis, findUMI  # type: ignore

"""
Test the findUMI function with a single motif
"""

@pytest.fixture
def test_config_one_motif():
    """Fixture for test configuration with single motif"""
    return {
        'min_read_length': 2000,
        'barcode_motifs': ["NNNNNNNNNNNN"],  # 12N
        'test_data_dir': "test-data/test_extract",
        'reference_file': "test-data/test_extract/extract_umi_1.fa",
        'bam_file': "test-data/test_extract/extract_umi_1.bam"
    }

@pytest.fixture
def test_paths(test_config_one_motif):
    """Fixture to build absolute paths from test file location"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_one_motif['reference_file'],
        'bam_file': base_dir / test_config_one_motif['bam_file'],
        'base_dir': base_dir
    }

def test_findUMI(test_config_one_motif, test_paths):
    """Test the findUMI function extracts UMIs correctly"""
    reference_file = test_paths['reference_file']
    bam_file = test_paths['bam_file']
    
    # Verify test files exist
    assert reference_file.exists(), f"Reference file not found: {reference_file}"
    assert bam_file.exists(), f"BAM file not found: {bam_file}"

    # Run the findUMI function
    motif_indices, raw_UMI_dict = findUMI(
        str(reference_file), 
        str(bam_file), 
        test_config_one_motif['barcode_motifs'], 
        test_config_one_motif['min_read_length']
    )

    print(motif_indices)
    print(raw_UMI_dict)

    # Verify results
    # Make sure only 1 UMI is found per read (no duplicate read names)
    assert len(raw_UMI_dict) == len(set(raw_UMI_dict.keys())), \
        "Duplicate read names found in UMI dictionary"
    
    # Verify we found at least some UMIs
    assert len(raw_UMI_dict) > 0, "No UMIs were found"
    
    # Verify each read has the correct number of UMIs (one per motif)
    expected_umi_count = len(test_config_one_motif['barcode_motifs'])
    for read_name, umis in raw_UMI_dict.items():
        assert len(umis) == expected_umi_count, \
            f"Read {read_name} has {len(umis)} UMIs, expected {expected_umi_count}"
    
    # Verify the number of motif positions found matches the expected count
    assert len(motif_indices) == expected_umi_count, \
        f"Found {len(motif_indices)} motif positions, expected {expected_umi_count}"

"""
Test the findUMI function with five motifs
"""
@pytest.fixture
def test_config_five_motifs():
    """Fixture for test configuration with five motifs"""
    return {
        'min_read_length': 2000,
        'barcode_motifs': ["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"],
        'test_data_dir': "test-data/test_extract",
        'reference_file': "test-data/test_extract/extract_umi_5.fa",
        'bam_file': "test-data/test_extract/extract_umi_5.bam"
    }

@pytest.fixture
def test_paths_five_motifs(test_config_five_motifs):
    """Fixture to build absolute paths for five motif test"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_five_motifs['reference_file'],
        'bam_file': base_dir / test_config_five_motifs['bam_file'],
        'base_dir': base_dir
    }

def test_findUMI_five_motifs(test_config_five_motifs, test_paths_five_motifs):  
    """Test the findUMI function extracts UMIs correctly with five motifs"""
    reference_file = test_paths_five_motifs['reference_file']
    bam_file = test_paths_five_motifs['bam_file']
    
    # Verify test files exist
    assert reference_file.exists(), f"Reference file not found: {reference_file}"
    assert bam_file.exists(), f"BAM file not found: {bam_file}"
    
    # Run the findUMI function
    motif_indices, raw_UMI_dict = findUMI(
        str(reference_file), 
        str(bam_file), 
        test_config_five_motifs['barcode_motifs'], 
        test_config_five_motifs['min_read_length']
    )

    print(motif_indices)
    print(raw_UMI_dict)
    
    # Verify results
    # Make sure only 1 UMI is found per read (no duplicate read names)
    assert len(raw_UMI_dict) == len(set(raw_UMI_dict.keys())), \
        "Duplicate read names found in UMI dictionary"
    
    # Verify we found at least some UMIs
    assert len(raw_UMI_dict) > 0, "No UMIs were found"
    
    # Verify each read has the correct number of UMIs (one per motif)
    expected_umi_count = len(test_config_five_motifs['barcode_motifs'])
    for read_name, umis in raw_UMI_dict.items():
        assert len(umis) == expected_umi_count, \
            f"Read {read_name} has {len(umis)} UMIs, expected {expected_umi_count}"
    
    # Verify the number of motif positions found matches the expected count
    assert len(motif_indices) == expected_umi_count, \
        f"Found {len(motif_indices)} motif positions, expected {expected_umi_count}"

### TODO: Add tests for reads where only 4/5 UMIs are found
### TODO: Add tests for reads where no UMIs are found

"""
Test the extract_umis function with a larger dataset
"""
@pytest.fixture
def test_config_extract_large():
    """Fixture for test configuration with larger dataset"""
    return {
        'min_read_length': 2000,
        'barcode_motifs': ["NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN", "NNNNNNNNNNNN"],
        'reference_file': "test-data/test_extract/extract_umi_5-large.fasta",
        'bam_file': "test-data/test_extract/extract_umi_5-large.bam",
        'output_file': "test-data/test_extract/test_output_large.csv"
    }

@pytest.fixture
def test_paths_extract_large(test_config_extract_large):
    """Fixture to build absolute paths for extract_umis test"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_extract_large['reference_file'],
        'bam_file': base_dir / test_config_extract_large['bam_file'],
        'output_file': base_dir / test_config_extract_large['output_file'],
        'base_dir': base_dir
    }

def test_extract_umis_large(test_config_extract_large, test_paths_extract_large):
    """Test the extract_umis function on a larger dataset and validate CSV output"""
    reference_file = test_paths_extract_large['reference_file']
    bam_file = test_paths_extract_large['bam_file']
    output_file = test_paths_extract_large['output_file']
    
    # Verify test files exist
    assert reference_file.exists(), f"Reference file not found: {reference_file}"
    assert bam_file.exists(), f"BAM file not found: {bam_file}"
    
    # Clean up output file if it exists from previous run
    if output_file.exists():
        output_file.unlink()
    
    # Run the extract_umis function
    extract_umis(
        str(reference_file),
        str(bam_file),
        test_config_extract_large['min_read_length'],
        test_config_extract_large['barcode_motifs'],
        str(output_file)
    )
    
    # Verify output file was created
    assert output_file.exists(), f"Output CSV file was not created: {output_file}"
    
    # Read and validate the CSV file
    df = pd.read_csv(output_file)
    
    print(f"\n{'='*60}")
    print(f"CSV Output Summary:")
    print(f"  - Total reads with UMIs: {len(df)}")
    print(f"  - Columns: {list(df.columns)}")
    print(f"  - First 3 rows:")
    print(df.head(3).to_string(index=False))
    print(f"{'='*60}\n")
    
    # Validate CSV structure
    expected_columns = ['Read Name'] + [f'UMI{i+1}' for i in range(len(test_config_extract_large['barcode_motifs']))]
    assert list(df.columns) == expected_columns, \
        f"CSV columns mismatch. Expected {expected_columns}, got {list(df.columns)}"
    
    # Verify we have some reads
    assert len(df) > 0, "No reads found in output CSV"
    
    # Verify each UMI column contains valid DNA sequences
    umi_columns = [col for col in df.columns if col.startswith('UMI')]
    for col in umi_columns:
        # Check that UMIs are non-null and are strings
        assert df[col].notna().all(), f"Column {col} contains null values"
        
        # Check UMI length (should be 12 for NNNNNNNNNNNN motif)
        umi_lengths = df[col].str.len()
        assert (umi_lengths == 12).all(), \
            f"Column {col} has UMIs with incorrect length. Expected 12, got lengths: {umi_lengths.unique()}"
        
        # Check UMIs contain only valid DNA bases
        valid_bases = set('ACGTN')
        for umi in df[col]:
            assert all(base in valid_bases for base in umi), \
                f"Column {col} contains invalid DNA bases in UMI: {umi}"
    
    # Verify Read Name column is unique (no duplicates)
    assert df['Read Name'].is_unique, "Duplicate read names found in CSV"
    
    # Clean up test output file
    if output_file.exists():
        output_file.unlink()
    
    print(f"✓ Test passed: {len(df)} reads processed successfully")
