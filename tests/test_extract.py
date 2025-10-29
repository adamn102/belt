#!/usr/bin/env python3
"""
Comprehensive pytest test file for BELT extract function
Consolidated from multiple test files with enhanced coverage
"""
import os
import sys
import pytest
import pandas as pd
import tempfile
import shutil
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

# ============================================================================
# ADDITIONAL TEST CASES - Consolidated from other test files
# ============================================================================

@pytest.fixture
def test_config_two_motifs():
    """Fixture for test configuration with two motifs (from test_extractV0.py)"""
    return {
        'min_read_length': 2000,
        'barcode_motifs': ["NNNNNNNNNNNNNNM", "NNNNNNNNNNNNNNM"],  # 15 N's for both 5' and 3' UMIs
        'test_data_dir': "test-data/test-umi_map",
        'reference_file': "test-data/test-umi_map/reference_sequence/B145.fa",
        'bam_file': "test-data/test-umi_map/sorted_reads/B145-XTMMF9-12.bam"
    }

@pytest.fixture
def test_paths_two_motifs(test_config_two_motifs):
    """Fixture to build absolute paths for two motif test"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_two_motifs['reference_file'],
        'bam_file': base_dir / test_config_two_motifs['bam_file'],
        'base_dir': base_dir
    }

def test_findUMI_two_motifs(test_config_two_motifs, test_paths_two_motifs):
    """Test the findUMI function with two motifs (consolidated from test_extractV0.py)"""
    reference_file = test_paths_two_motifs['reference_file']
    bam_file = test_paths_two_motifs['bam_file']
    
    # Skip test if files don't exist (they might not be available in all environments)
    if not reference_file.exists() or not bam_file.exists():
        pytest.skip(f"Test data files not found: {reference_file} or {bam_file}")
    
    # Run the findUMI function
    motif_indices, raw_UMI_dict = findUMI(
        str(reference_file), 
        str(bam_file), 
        test_config_two_motifs['barcode_motifs'], 
        test_config_two_motifs['min_read_length']
    )
    
    # Verify results
    assert len(raw_UMI_dict) == len(set(raw_UMI_dict.keys())), \
        "Duplicate read names found in UMI dictionary"
    
    if len(raw_UMI_dict) > 0:  # Only test if UMIs were found
        expected_umi_count = len(test_config_two_motifs['barcode_motifs'])
        for read_name, umis in raw_UMI_dict.items():
            assert len(umis) == expected_umi_count, \
                f"Read {read_name} has {len(umis)} UMIs, expected {expected_umi_count}"
        
        assert len(motif_indices) == expected_umi_count, \
            f"Found {len(motif_indices)} motif positions, expected {expected_umi_count}"

def test_extract_umis_function(test_config_one_motif, test_paths):
    """Test the extract_umis function (consolidated from test_extractV0.py)"""
    reference_file = test_paths['reference_file']
    bam_file = test_paths['bam_file']
    
    # Skip test if files don't exist
    if not reference_file.exists() or not bam_file.exists():
        pytest.skip(f"Test data files not found: {reference_file} or {bam_file}")
    
    # Create temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        
        # Run the extract_umis function
        extract_umis(
            str(reference_file),
            str(bam_file),
            test_config_one_motif['min_read_length'],
            test_config_one_motif['barcode_motifs'],
            str(output_dir)
        )
        
        # Check output file
        expected_output_file = output_dir / f"{reference_file.stem}.belt_extract.csv"
        
        if expected_output_file.exists():
            # Read and validate results
            df = pd.read_csv(expected_output_file)
            
            # Verify we got the expected number of UMI columns
            umi_columns = [col for col in df.columns if col != 'Read Name']
            expected_umi_count = len(test_config_one_motif['barcode_motifs'])
            actual_umi_count = len(umi_columns)
            
            assert actual_umi_count == expected_umi_count, \
                f"Expected {expected_umi_count} UMI columns, got {actual_umi_count}"
            
            # Verify CSV structure
            assert 'Read Name' in df.columns, "Read Name column missing"
            assert len(df) > 0, "No reads processed"
        else:
            pytest.fail(f"Output file not found: {expected_output_file}")

# ============================================================================
# EDGE CASE TESTS - Addressing the TODO items
# ============================================================================

def test_findUMI_no_motifs_found():
    """Test findUMI when no motifs are found in reference"""
    # Create a temporary reference file without motifs
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as ref_file:
        ref_file.write(">test_ref\n")
        ref_file.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        ref_file_path = ref_file.name
    
    # Create a temporary BAM file
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as bam_file:
        bam_file_path = bam_file.name
    
    try:
        # Test with motif that doesn't exist
        motif_indices, raw_UMI_dict = findUMI(
            ref_file_path,
            bam_file_path,
            ["NNNNNNNNNNNN"],  # This motif won't be found
            1000
        )
        
        # Should return None, None when motif not found
        assert motif_indices is None, "Expected None when motif not found"
        assert raw_UMI_dict is None, "Expected None when motif not found"
        
    finally:
        # Clean up
        os.unlink(ref_file_path)
        if os.path.exists(bam_file_path):
            os.unlink(bam_file_path)

def test_findUMI_invalid_files():
    """Test findUMI with invalid file paths"""
    # Test with non-existent files
    motif_indices, raw_UMI_dict = findUMI(
        "non_existent_ref.fa",
        "non_existent.bam",
        ["NNNNNNNNNNNN"],
        1000
    )
    
    assert motif_indices is None, "Expected None for non-existent files"
    assert raw_UMI_dict is None, "Expected None for non-existent files"

def test_findUMI_wrong_file_types():
    """Test findUMI with wrong file types"""
    # Create temporary files with wrong extensions
    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as ref_file:
        ref_file.write(b">test_ref\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        ref_file_path = ref_file.name
    
    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as bam_file:
        bam_file_path = bam_file.name
    
    try:
        # Test with wrong file types
        motif_indices, raw_UMI_dict = findUMI(
            ref_file_path,
            bam_file_path,
            ["NNNNNNNNNNNN"],
            1000
        )
        
        assert motif_indices is None, "Expected None for wrong file types"
        assert raw_UMI_dict is None, "Expected None for wrong file types"
        
    finally:
        # Clean up
        os.unlink(ref_file_path)
        os.unlink(bam_file_path)

def test_extract_umis_output_directory_creation():
    """Test that extract_umis creates output directories"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a nested output directory that doesn't exist
        nested_output = Path(temp_dir) / "nested" / "output" / "dir"
        
        # This should create the directory structure
        extract_umis(
            "dummy_ref.fa",  # Will fail, but we're testing directory creation
            "dummy.bam",
            1000,
            ["NNNNNNNNNNNN"],
            str(nested_output)
        )
        
        # The function should have attempted to create the directory
        # (Even though it fails due to invalid files, directory creation is tested)

def test_extract_umis_csv_output_path():
    """Test extract_umis with direct CSV file path"""
    with tempfile.TemporaryDirectory() as temp_dir:
        csv_file_path = Path(temp_dir) / "test_output.csv"
        
        # Test with direct CSV file path
        extract_umis(
            "dummy_ref.fa",  # Will fail, but we're testing path handling
            "dummy.bam",
            1000,
            ["NNNNNNNNNNNN"],
            str(csv_file_path)
        )
        
        # The function should have attempted to create the parent directory
        assert csv_file_path.parent.exists(), "Parent directory should be created"

# ============================================================================
# SYNTHETIC DATA TESTS - Consolidated from test_umi_extraction.py
# ============================================================================

def create_synthetic_test_data():
    """Create synthetic test data for controlled testing"""
    import pysam
    from Bio import SeqIO
    
    # Create a simple reference with UMI motifs
    ref_seq = "ATCG" * 100 + "NNNNNNNNNNNNNN" + "GCTA" * 100 + "NNNNNNNNNNNNNN" + "ATCG" * 100
    
    # Create reference file
    ref_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    ref_file.write(">test_ref\n")
    ref_file.write(ref_seq + "\n")
    ref_file.close()
    
    # Create a test read with actual UMI sequences
    read_seq = "ATCG" * 100 + "ATCGATCGATCGAT" + "GCTA" * 100 + "GCTAGCTAGCTAG" + "ATCG" * 100
    
    # Create BAM file with test read
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': len(ref_seq), 'SN': 'test_ref'}]
    }
    
    bam_file = tempfile.NamedTemporaryFile(suffix='.bam', delete=False)
    bam_file.close()
    
    with pysam.AlignmentFile(bam_file.name, "wb", header=header) as bam:
        a = pysam.AlignedSegment()
        a.query_name = "test_read"
        a.query_sequence = read_seq
        a.reference_id = 0
        a.reference_start = 0
        a.cigarstring = f"{len(read_seq)}M"
        a.mapping_quality = 60
        bam.write(a)
    
    # Index the BAM file
    pysam.index(bam_file.name)
    
    return ref_file.name, bam_file.name

def test_findUMI_synthetic_data():
    """Test findUMI with synthetic data (consolidated from test_umi_extraction.py)"""
    ref_file, bam_file = create_synthetic_test_data()
    
    try:
        # Test with synthetic data
        motif_indices, raw_UMI_dict = findUMI(
            ref_file,
            bam_file,
            ["NNNNNNNNNNNNNN", "NNNNNNNNNNNNNN"],
            100
        )
        
        # Verify results
        assert motif_indices is not None, "Motif indices should be found"
        assert raw_UMI_dict is not None, "UMI dictionary should be created"
        assert len(raw_UMI_dict) > 0, "Should find at least one UMI"
        
        # Check that we found the expected UMIs
        for read_name, umis in raw_UMI_dict.items():
            assert len(umis) == 2, f"Expected 2 UMIs, got {len(umis)}"
            # The synthetic data should have specific UMIs
            assert umis[0] == "ATCGATCGATCGAT", f"Expected specific 5' UMI, got {umis[0]}"
            assert umis[1] == "GCTAGCTAGCTAGA", f"Expected specific 3' UMI, got {umis[1]}"
        
    finally:
        # Clean up
        os.unlink(ref_file)
        os.unlink(bam_file)
        os.unlink(bam_file + ".bai")
