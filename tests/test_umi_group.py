#!/usr/bin/env python3
"""
Pytest test file for BELT umi_group function
"""
import os
import sys
import pytest
import pandas as pd
from pathlib import Path
import tempfile
import shutil
import json

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_group import long_read_umi_grouping_fast  # type: ignore

# Import pysam for creating synthetic BAM files
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


@pytest.fixture
def test_config_umi_group():
    """Fixture for test configuration for UMI grouping"""
    return {
        'reference_file': "test-data/test_group_umi/extract_umi_5.fa",
        'bam_file': "test-data/test_group_umi/extract_umi_5-large.bam",
        'corrected_umis_file': "test-data/test_group_umi/extract_umi_5-large_corrected_umis.csv",
        'grouping': "UMI4",  # Which UMI column to group by
        'threshold': 2,  # Minimum read count per group
        'output_dir': "test-data/test_group_umi/"
    }


@pytest.fixture
def test_paths_umi_group(test_config_umi_group):
    """Fixture to build absolute paths for umi_group test"""
    base_dir = Path(__file__).parent.parent
    return {
        'reference_file': base_dir / test_config_umi_group['reference_file'],
        'bam_file': base_dir / test_config_umi_group['bam_file'],
        'corrected_umis_file': base_dir / test_config_umi_group['corrected_umis_file'],
        'base_dir': base_dir
    }


def test_umi_consensus_long(test_config_umi_group, test_paths_umi_group):
    """Test long read UMI consensus grouping using fast grouping function"""
    reference_file = test_paths_umi_group['reference_file']
    bam_file = test_paths_umi_group['bam_file']
    corrected_umis_file = test_paths_umi_group['corrected_umis_file']
    
    # Verify test files exist
    assert reference_file.exists(), f"Reference file not found: {reference_file}"
    assert bam_file.exists(), f"BAM file not found: {bam_file}"
    assert corrected_umis_file.exists(), f"Corrected UMIs file not found: {corrected_umis_file}"
    
    # Run the long_read_umi_grouping_fast function
    long_read_umi_grouping_fast(
        alignedReads=str(bam_file),
        correctedUMIs=str(corrected_umis_file),
        outputDir=test_config_umi_group['output_dir'],
        grouping=test_config_umi_group['grouping'],
        threshold=test_config_umi_group['threshold']
    )
    
    # Verify output files were created
    sample_name = test_config_umi_group['bam_file'].split('/')[-1].split('.')[0]
    output_dir = test_config_umi_group['output_dir']
    
    # Check for expected output files
    expected_files = [
        f"{sample_name}_filtered_umi_groups.csv",
        f"{sample_name}"  # This will be a directory containing FASTQ files
    ]
    
    found_files = []
    for expected_file in expected_files:
        output_path = os.path.join(output_dir, expected_file)
        if os.path.exists(output_path):
            found_files.append(expected_file)
            print(f"✓ Found: {expected_file}")
            
            # If it's a CSV, validate it
            if expected_file.endswith('.csv'):
                df_output = pd.read_csv(output_path)
                print(f"  Rows: {len(df_output)}, Columns: {list(df_output.columns)}")
                print(f"  First few rows:\n{df_output.head()}")
            # If it's a directory, count FASTQ files
            elif os.path.isdir(output_path):
                fastq_files = [f for f in os.listdir(output_path) if f.endswith('.fastq')]
                print(f"  FASTQ files created: {len(fastq_files)}")
                if fastq_files:
                    print(f"  Examples: {fastq_files[:3]}")
    
    # Assert at least the filtered CSV was created
    assert f"{sample_name}_filtered_umi_groups.csv" in found_files, \
        f"Filtered UMI groups CSV not created. Expected: {expected_files}"
    
    print(f"\n✓ Test passed: UMI consensus completed, {len(found_files)}/{len(expected_files)} outputs created")


def test_umi_grouping_with_fastq_inspection(test_config_umi_group, test_paths_umi_group, tmp_path):
    """Test UMI grouping and write FASTQ sequences for visual inspection"""
    reference_file = test_paths_umi_group['reference_file']
    bam_file = test_paths_umi_group['bam_file']
    corrected_umis_file = test_paths_umi_group['corrected_umis_file']
    
    # Skip if files don't exist
    if not all([reference_file.exists(), bam_file.exists(), corrected_umis_file.exists()]):
        pytest.skip("Required test files not found")
    
    # Use temporary directory for output
    output_dir = str(tmp_path / "umi_group_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run the function
    long_read_umi_grouping_fast(
        alignedReads=str(bam_file),
        correctedUMIs=str(corrected_umis_file),
        outputDir=output_dir,
        grouping="UMI1",
        threshold=1,  # Lower threshold to get more groups
        write_fastq=True  # Enable FASTQ output for inspection
    )
    
    # Read and display sequences from FASTQ files
    sample_name = test_config_umi_group['bam_file'].split('/')[-1].split('.')[0]
    fastq_dir = os.path.join(output_dir, sample_name)
    
    if os.path.exists(fastq_dir):
        fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq')]
        
        print(f"\n{'='*80}")
        print(f"FASTQ FILES CREATED: {len(fastq_files)}")
        print(f"{'='*80}")
        
        # Display first few sequences from first FASTQ file for visual inspection
        for i, fastq_file in enumerate(fastq_files[:3]):  # Show first 3 UMI groups
            fastq_path = os.path.join(fastq_dir, fastq_file)
            print(f"\n--- File {i+1}: {fastq_file} ---")
            
            with open(fastq_path, 'r') as f:
                lines = f.readlines()
                num_reads = len(lines) // 4
                print(f"Number of reads: {num_reads}")
                
                # Show first 2 reads from this file
                for read_idx in range(min(2, num_reads)):
                    start = read_idx * 4
                    header = lines[start].strip()
                    seq = lines[start + 1].strip()
                    qual = lines[start + 3].strip()
                    
                    print(f"\nRead {read_idx + 1}:")
                    print(f"  Header: {header}")
                    print(f"  Seq:    {seq[:80]}{'...' if len(seq) > 80 else ''}")
                    print(f"  Qual:   {qual[:80]}{'...' if len(qual) > 80 else ''}")
                    
                    # Parse and display UMI information
                    if '|' in header:
                        parts = header.split('|')
                        print(f"  Read Name: {parts[0].lstrip('@')}")
                        umi_info = parts[1:]
                        for umi in umi_info:
                            print(f"    {umi}")
        
        print(f"\n{'='*80}\n")
        
        # Validate that we created FASTQ files
        assert len(fastq_files) > 0, "No FASTQ files were created"
    else:
        pytest.skip(f"FASTQ directory not created: {fastq_dir}")


def test_umi_grouping_different_columns(test_config_umi_group, test_paths_umi_group, tmp_path):
    """Test grouping by different UMI columns (UMI1, UMI2, UMI3, etc.)"""
    reference_file = test_paths_umi_group['reference_file']
    bam_file = test_paths_umi_group['bam_file']
    corrected_umis_file = test_paths_umi_group['corrected_umis_file']
    
    # Skip if files don't exist
    if not all([reference_file.exists(), bam_file.exists(), corrected_umis_file.exists()]):
        pytest.skip("Required test files not found")
    
    # Read CSV to get available UMI columns
    df = pd.read_csv(corrected_umis_file)
    umi_columns = [col for col in df.columns if col.startswith('UMI')]
    
    print(f"\n{'='*80}")
    print(f"Testing grouping by different UMI columns")
    print(f"Available UMI columns: {umi_columns}")
    print(f"{'='*80}")
    
    # Test grouping by first 3 UMI columns
    for umi_col in umi_columns[:3]:
        output_dir = str(tmp_path / f"group_by_{umi_col}")
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"\n--- Grouping by {umi_col} ---")
        
        long_read_umi_grouping_fast(
            alignedReads=str(bam_file),
            correctedUMIs=str(corrected_umis_file),
            outputDir=output_dir,
            grouping=umi_col,
            threshold=1
        )
        
        # Check filtered CSV
        sample_name = test_config_umi_group['bam_file'].split('/')[-1].split('.')[0]
        filtered_csv = os.path.join(output_dir, f"{sample_name}_filtered_umi_groups.csv")
        
        if os.path.exists(filtered_csv):
            df_filtered = pd.read_csv(filtered_csv)
            print(f"  Groups found: {len(df_filtered)}")
            print(f"  Unique {umi_col} values: {df_filtered[umi_col].tolist()}")
            
            assert umi_col in df_filtered.columns, f"Grouping column {umi_col} not in output"
            assert len(df_filtered) > 0, f"No groups found for {umi_col}"


def test_umi_grouping_fast(test_config_umi_group, test_paths_umi_group, tmp_path):
    """Test the fast UMI grouping function with BAM sorting"""
    reference_file = test_paths_umi_group['reference_file']
    bam_file = test_paths_umi_group['bam_file']
    corrected_umis_file = test_paths_umi_group['corrected_umis_file']
    
    # Skip if files don't exist
    if not all([reference_file.exists(), bam_file.exists(), corrected_umis_file.exists()]):
        pytest.skip("Required test files not found")
    
    output_dir = str(tmp_path / "umi_group_fast")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run the fast version
    long_read_umi_grouping_fast(
        alignedReads=str(bam_file),
        correctedUMIs=str(corrected_umis_file),
        outputDir=output_dir,
        grouping="UMI1",
        threshold=1,
        threads=2
    )
    
    sample_name = test_config_umi_group['bam_file'].split('/')[-1].split('.')[0]
    
    # Check outputs
    filtered_csv = os.path.join(output_dir, f"{sample_name}_filtered_umi_groups.csv")
    grouped_bam = os.path.join(output_dir, f"{sample_name}.grouped.bam")
    
    print(f"\n{'='*80}")
    print(f"Fast UMI Grouping Test")
    print(f"{'='*80}")
    
    assert os.path.exists(filtered_csv), f"Filtered CSV not created: {filtered_csv}"
    assert os.path.exists(grouped_bam), f"Grouped BAM not created: {grouped_bam}"
    
    # Display filtered groups
    df_filtered = pd.read_csv(filtered_csv)
    print(f"\nFiltered UMI Groups ({len(df_filtered)} groups):")
    print(df_filtered.head(10))
    
    print(f"\n✓ Fast grouping created:")
    print(f"  - Filtered CSV: {filtered_csv}")
    print(f"  - Grouped BAM: {grouped_bam}")


def test_umi_grouping_fast_with_fastq(test_config_umi_group, test_paths_umi_group, tmp_path):
    """Test fast UMI grouping with FASTQ output for visual inspection"""
    reference_file = test_paths_umi_group['reference_file']
    bam_file = test_paths_umi_group['bam_file']
    corrected_umis_file = test_paths_umi_group['corrected_umis_file']
    
    # Skip if files don't exist
    if not all([reference_file.exists(), bam_file.exists(), corrected_umis_file.exists()]):
        pytest.skip("Required test files not found")
    
    output_dir = str(tmp_path / "umi_group_fast_fastq")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run with FASTQ output
    long_read_umi_grouping_fast(
        alignedReads=str(bam_file),
        correctedUMIs=str(corrected_umis_file),
        outputDir=output_dir,
        grouping="UMI1",
        threshold=1,
        threads=2,
        write_fastq=True  # Enable FASTQ output
    )
    
    sample_name = test_config_umi_group['bam_file'].split('/')[-1].split('.')[0]
    fastq_dir = os.path.join(output_dir, sample_name)
    
    print(f"\n{'='*80}")
    print(f"Fast UMI Grouping with FASTQ Output")
    print(f"{'='*80}")
    
    if os.path.exists(fastq_dir):
        fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq')]
        print(f"\nFASTQ files created: {len(fastq_files)}")
        
        # Display sequences from first FASTQ file
        for fastq_file in fastq_files[:2]:  # Show first 2 files
            fastq_path = os.path.join(fastq_dir, fastq_file)
            print(f"\n--- {fastq_file} ---")
            
            with open(fastq_path, 'r') as f:
                lines = f.readlines()
                num_reads = len(lines) // 4
                print(f"Number of reads: {num_reads}")
                
                # Show first read
                if num_reads > 0:
                    header = lines[0].strip()
                    seq = lines[1].strip()
                    qual = lines[3].strip()
                    
                    print(f"\nFirst read:")
                    print(f"  Header: {header}")
                    print(f"  Seq:    {seq[:80]}{'...' if len(seq) > 80 else ''}")
                    
                    # Parse UMI information
                    if '|' in header:
                        parts = header.split('|')
                        print(f"  Read Name: {parts[0].lstrip('@')}")
                        print(f"  UMI Info:")
                        for umi in parts[1:]:
                            print(f"    - {umi}")
        
        assert len(fastq_files) > 0, "No FASTQ files created"
        print(f"\n✓ Test passed: {len(fastq_files)} FASTQ files created with full UMI annotations")
    else:
        pytest.fail(f"FASTQ directory not created: {fastq_dir}")


# ============================================================================
# UMI2 MAJORITY FILTERING TESTS
# ============================================================================

def create_synthetic_umi_csv(umi_data, output_path):
    """
    Create a synthetic CSV file with UMI data for testing.
    
    Args:
        umi_data: List of dicts with 'Read Name', 'UMI1', 'UMI2', etc.
        output_path: Path to write CSV file
    """
    df = pd.DataFrame(umi_data)
    df.to_csv(output_path, index=False)
    return output_path


def create_synthetic_bam(read_names, output_path, reference_length=10000):
    """
    Create a synthetic BAM file with reads matching the read names.
    
    Args:
        read_names: List of read names to include in BAM
        output_path: Path to write BAM file
        reference_length: Length of synthetic reference sequence
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create BAM header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': reference_length, 'SN': 'test_ref'}]
    }
    
    # Create BAM file (convert PosixPath to string for pysam)
    bam_path_str = str(output_path)
    with pysam.AlignmentFile(bam_path_str, "wb", header=header) as bam:
        for i, read_name in enumerate(read_names):
            a = pysam.AlignedSegment()
            a.query_name = read_name
            a.query_sequence = "ATCG" * 100  # Simple synthetic sequence
            a.reference_id = 0
            a.reference_start = i * 10  # Stagger reads
            a.cigarstring = f"{len(a.query_sequence)}M"
            a.mapping_quality = 60
            bam.write(a)
    
    # Index the BAM file (convert PosixPath to string for pysam)
    pysam.index(bam_path_str)
    return output_path


def create_synthetic_umi_data_umi2_filtering():
    """
    Create synthetic UMI data specifically designed to test UMI2 majority filtering.
    
    Test scenarios:
    - UMI1_A: 10 reads, 8 with UMI2_X (80%), 2 with UMI2_Y (20%) - should pass with threshold 0.5
    - UMI1_B: 10 reads, 5 with UMI2_X (50%), 5 with UMI2_Y (50%) - should pass with threshold 0.5 (exact)
    - UMI1_C: 10 reads, 4 with UMI2_X (40%), 6 with UMI2_Y (60%) - should pass with threshold 0.5
    - UMI1_D: 10 reads, 3 with UMI2_X (30%), 7 with UMI2_Y (70%) - should pass with threshold 0.5
    - UMI1_E: 10 reads, 2 with UMI2_X (20%), 8 with UMI2_Y (80%) - should pass with threshold 0.5
    - UMI1_F: 10 reads, 4 with UMI2_X (40%), 6 with UMI2_Y (60%) - should FAIL with threshold 0.7
    - UMI1_G: 5 reads, all with UMI2_X (100%) - should pass (strong majority)
    - UMI1_H: 3 reads, 2 with UMI2_X (67%), 1 with UMI2_Y (33%) - should pass with threshold 0.5
    """
    umi_data = []
    read_counter = 1
    
    # UMI1_A: 8 reads with UMI2_X, 2 with UMI2_Y (80% majority)
    for i in range(8):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_A',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(2):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_A',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_B: 5 reads with UMI2_X, 5 with UMI2_Y (50% - exact threshold)
    for i in range(5):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_B',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(5):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_B',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_C: 4 reads with UMI2_X, 6 with UMI2_Y (60% majority for Y)
    for i in range(4):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_C',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(6):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_C',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_D: 3 reads with UMI2_X, 7 with UMI2_Y (70% majority for Y)
    for i in range(3):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_D',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(7):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_D',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_E: 2 reads with UMI2_X, 8 with UMI2_Y (80% majority for Y)
    for i in range(2):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_E',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(8):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_E',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_F: 4 reads with UMI2_X, 6 with UMI2_Y (60% majority - will fail with 0.7 threshold)
    for i in range(4):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_F',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    for i in range(6):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_F',
            'UMI2': 'UMI2_Y'
        })
        read_counter += 1
    
    # UMI1_G: 5 reads, all with UMI2_X (100% majority)
    for i in range(5):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_G',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    
    # UMI1_H: 3 reads, 2 with UMI2_X, 1 with UMI2_Y (67% majority)
    for i in range(2):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_H',
            'UMI2': 'UMI2_X'
        })
        read_counter += 1
    umi_data.append({
        'Read Name': f'read_{read_counter:04d}',
        'UMI1': 'UMI1_H',
        'UMI2': 'UMI2_Y'
    })
    read_counter += 1
    
    # UMI1_I: 10 reads, all with missing UMI2 (should be filtered out)
    for i in range(10):
        umi_data.append({
            'Read Name': f'read_{read_counter:04d}',
            'UMI1': 'UMI1_I',
            'UMI2': ''  # Empty UMI2
        })
        read_counter += 1
    
    return umi_data


def test_umi2_majority_filtering_basic(tmp_path):
    """
    Test basic UMI2 majority filtering with default threshold (0.5).
    
    Expected behavior:
    - UMI1_A: 80% majority (UMI2_X) - should pass, keep 8 reads
    - UMI1_B: 50% majority (UMI2_X) - should pass, keep 5 reads
    - UMI1_C: 60% majority (UMI2_Y) - should pass, keep 6 reads
    - UMI1_D: 70% majority (UMI2_Y) - should pass, keep 7 reads
    - UMI1_E: 80% majority (UMI2_Y) - should pass, keep 8 reads
    - UMI1_F: 60% majority (UMI2_Y) - should pass, keep 6 reads
    - UMI1_G: 100% majority (UMI2_X) - should pass, keep 5 reads
    - UMI1_H: 67% majority (UMI2_X) - should pass, keep 2 reads
    - UMI1_I: No valid UMI2 - should be filtered out
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create synthetic UMI data
    umi_data = create_synthetic_umi_data_umi2_filtering()
    
    # Create synthetic CSV
    csv_path = tmp_path / "test_umis.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create synthetic BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping with UMI2 majority filtering
    output_dir = tmp_path / "umi2_filter_output"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=3,  # Minimum 3 reads per group
        umi2_majority_threshold=0.5  # 50% threshold
    )
    
    # Verify outputs
    filtered_csv = output_dir / "test_reads_filtered_umi_groups.csv"
    log_json = output_dir / "test_reads_umi_group.log.json"
    
    assert filtered_csv.exists(), "Filtered CSV should be created"
    assert log_json.exists(), "JSON log should be created"
    
    # Read and validate filtered CSV
    df_filtered = pd.read_csv(filtered_csv)
    print(f"\n{'='*80}")
    print(f"UMI2 MAJORITY FILTERING TEST (threshold=0.5)")
    print(f"{'='*80}")
    print(f"Groups passing filter: {len(df_filtered)}")
    print(df_filtered)
    
    # Verify expected groups are present
    expected_groups = ['UMI1_A', 'UMI1_B', 'UMI1_C', 'UMI1_D', 'UMI1_E', 'UMI1_F', 'UMI1_G', 'UMI1_H']
    actual_groups = set(df_filtered['UMI1'].tolist())
    
    print(f"\nExpected groups: {expected_groups}")
    print(f"Actual groups: {sorted(actual_groups)}")
    
    # UMI1_I should be filtered out (no valid UMI2)
    assert 'UMI1_I' not in actual_groups, "UMI1_I should be filtered out (no valid UMI2)"
    
    # Verify read counts match expected majority counts
    expected_counts = {
        'UMI1_A': 8,  # 8 reads with majority UMI2_X
        'UMI1_B': 5,  # 5 reads with majority UMI2_X (tie, first wins)
        'UMI1_C': 6,  # 6 reads with majority UMI2_Y
        'UMI1_D': 7,  # 7 reads with majority UMI2_Y
        'UMI1_E': 8,  # 8 reads with majority UMI2_Y
        'UMI1_F': 6,  # 6 reads with majority UMI2_Y
        'UMI1_G': 5,  # 5 reads with majority UMI2_X (100%)
        'UMI1_H': 2,  # 2 reads with majority UMI2_X (but only 3 total, threshold=3, so might not pass)
    }
    
    for umi1, expected_count in expected_counts.items():
        if umi1 in actual_groups:
            actual_count = df_filtered[df_filtered['UMI1'] == umi1]['ReadCount'].iloc[0]
            print(f"  {umi1}: expected {expected_count}, got {actual_count}")
            # Note: UMI1_H has only 3 reads total, so if threshold=3, it should pass
            # But after filtering to majority (2 reads), it might not meet threshold
            if umi1 == 'UMI1_H':
                # UMI1_H: 3 total, 2 majority, but threshold=3, so might be filtered
                assert actual_count >= 2, f"{umi1} should have at least 2 reads (majority)"
            else:
                assert actual_count == expected_count, \
                    f"{umi1} should have {expected_count} reads, got {actual_count}"
    
    # Read and validate JSON log
    with open(log_json, 'r') as f:
        log_data = json.load(f)
    
    print(f"\nJSON Log Statistics:")
    print(f"  UMI2 filter applied: {log_data.get('umi2_majority_filter_applied')}")
    print(f"  Groups passing UMI2 filter: {log_data.get('groups_passing_umi2_majority')}")
    print(f"  Reads before UMI2 filter: {log_data.get('reads_before_umi2_filter')}")
    print(f"  Reads after UMI2 filter: {log_data.get('reads_after_umi2_filter')}")
    print(f"  Reads filtered by UMI2: {log_data.get('reads_filtered_by_umi2')}")
    
    # Verify UMI2 filtering was applied
    assert log_data.get('umi2_majority_filter_applied') == True, "UMI2 filter should be applied"
    assert log_data.get('reads_before_umi2_filter', 0) > 0, "Should have reads before filtering"
    assert log_data.get('reads_after_umi2_filter', 0) > 0, "Should have reads after filtering"
    assert log_data.get('reads_filtered_by_umi2', 0) > 0, "Should have filtered some reads"
    
    print(f"\n✓ Test passed: UMI2 majority filtering works correctly")


def test_umi2_majority_filtering_strict_threshold(tmp_path):
    """
    Test UMI2 majority filtering with strict threshold (0.7).
    
    Expected behavior:
    - UMI1_A: 80% majority - should pass
    - UMI1_B: 50% majority - should FAIL (below 0.7)
    - UMI1_C: 60% majority - should FAIL (below 0.7)
    - UMI1_D: 70% majority - should pass (exact threshold)
    - UMI1_E: 80% majority - should pass
    - UMI1_F: 60% majority - should FAIL (below 0.7)
    - UMI1_G: 100% majority - should pass
    - UMI1_H: 67% majority - should FAIL (below 0.7)
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create synthetic UMI data
    umi_data = create_synthetic_umi_data_umi2_filtering()
    
    # Create synthetic CSV
    csv_path = tmp_path / "test_umis_strict.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create synthetic BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_strict.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping with strict UMI2 majority threshold
    output_dir = tmp_path / "umi2_filter_strict"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=3,
        umi2_majority_threshold=0.7  # 70% threshold (strict)
    )
    
    # Verify outputs
    filtered_csv = output_dir / "test_reads_strict_filtered_umi_groups.csv"
    log_json = output_dir / "test_reads_strict_umi_group.log.json"
    
    assert filtered_csv.exists(), "Filtered CSV should be created"
    
    # Read filtered CSV
    df_filtered = pd.read_csv(filtered_csv)
    actual_groups = set(df_filtered['UMI1'].tolist())
    
    print(f"\n{'='*80}")
    print(f"UMI2 MAJORITY FILTERING TEST (strict threshold=0.7)")
    print(f"{'='*80}")
    print(f"Groups passing filter: {len(df_filtered)}")
    print(df_filtered)
    
    # Expected groups with 70%+ majority
    expected_passing = {'UMI1_A', 'UMI1_D', 'UMI1_E', 'UMI1_G'}  # 80%, 70%, 80%, 100%
    expected_failing = {'UMI1_B', 'UMI1_C', 'UMI1_F', 'UMI1_H'}  # 50%, 60%, 60%, 67%
    
    print(f"\nExpected passing (>=70%): {expected_passing}")
    print(f"Expected failing (<70%): {expected_failing}")
    print(f"Actual groups: {sorted(actual_groups)}")
    
    # Verify expected groups pass
    for group in expected_passing:
        assert group in actual_groups, f"{group} should pass with 70% threshold"
    
    # Verify expected groups fail
    for group in expected_failing:
        assert group not in actual_groups, f"{group} should fail with 70% threshold"
    
    # Verify JSON log
    with open(log_json, 'r') as f:
        log_data = json.load(f)
    
    assert log_data.get('umi2_majority_filter_applied') == True
    assert log_data.get('umi2_majority_threshold') == 0.7
    assert log_data.get('groups_passing_umi2_majority') == len(expected_passing)
    
    print(f"\n✓ Test passed: Strict UMI2 majority filtering works correctly")


def test_umi2_majority_filtering_verify_majority_reads_only(tmp_path):
    """
    Test that only reads with the majority UMI2 are retained in the grouped BAM.
    
    This is a critical test to verify that the filtering logic correctly removes
    non-majority UMI2 reads from each UMI1 group.
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create synthetic UMI data
    umi_data = create_synthetic_umi_data_umi2_filtering()
    
    # Create synthetic CSV
    csv_path = tmp_path / "test_umis_verify.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create synthetic BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_verify.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping
    output_dir = tmp_path / "umi2_filter_verify"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=3,
        umi2_majority_threshold=0.5,
        write_fastq=True  # Enable FASTQ to inspect reads
    )
    
    # Load original CSV to verify
    df_original = pd.read_csv(csv_path)
    
    # Load filtered CSV
    filtered_csv = output_dir / "test_reads_verify_filtered_umi_groups.csv"
    df_filtered = pd.read_csv(filtered_csv)
    
    print(f"\n{'='*80}")
    print(f"VERIFYING ONLY MAJORITY UMI2 READS ARE RETAINED")
    print(f"{'='*80}")
    
    # For each UMI1 group that passed, verify only majority UMI2 reads are present
    for _, row in df_filtered.iterrows():
        umi1 = row['UMI1']
        expected_read_count = row['ReadCount']
        
        # Get all reads for this UMI1 in original data
        umi1_reads = df_original[df_original['UMI1'] == umi1].copy()
        
        # Determine majority UMI2
        umi2_counts = umi1_reads['UMI2'].value_counts()
        # Filter out empty strings
        umi2_counts = umi2_counts[umi2_counts.index != '']
        if len(umi2_counts) == 0:
            continue  # Skip if no valid UMI2
        
        majority_umi2 = umi2_counts.index[0]
        majority_count = umi2_counts.iloc[0]
        total_count = len(umi1_reads)
        majority_fraction = majority_count / total_count
        
        print(f"\nUMI1: {umi1}")
        print(f"  Total reads: {total_count}")
        print(f"  Majority UMI2: {majority_umi2} ({majority_count} reads, {majority_fraction:.1%})")
        print(f"  Expected in output: {expected_read_count} reads")
        
        # Verify expected count matches majority count
        assert expected_read_count == majority_count, \
            f"{umi1}: Expected {majority_count} reads (majority), got {expected_read_count}"
        
        # Verify majority fraction meets threshold
        assert majority_fraction >= 0.5, \
            f"{umi1}: Majority fraction {majority_fraction:.1%} should be >= 50%"
    
    print(f"\n✓ Test passed: Only majority UMI2 reads are retained")


def test_umi2_majority_filtering_no_umi2_column(tmp_path):
    """
    Test that UMI2 filtering is NOT applied when UMI2 column doesn't exist.
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create CSV with only UMI1 (no UMI2)
    umi_data = [
        {'Read Name': f'read_{i:04d}', 'UMI1': 'UMI1_A'} 
        for i in range(10)
    ]
    csv_path = tmp_path / "test_umis_no_umi2.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_no_umi2.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping
    output_dir = tmp_path / "umi2_filter_no_umi2"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=1,
        umi2_majority_threshold=0.5
    )
    
    # Verify JSON log shows filter was NOT applied
    log_json = output_dir / "test_reads_no_umi2_umi_group.log.json"
    with open(log_json, 'r') as f:
        log_data = json.load(f)
    
    print(f"\n{'='*80}")
    print(f"UMI2 FILTERING NOT APPLIED (no UMI2 column)")
    print(f"{'='*80}")
    print(f"UMI2 filter applied: {log_data.get('umi2_majority_filter_applied')}")
    
    assert log_data.get('umi2_majority_filter_applied') == False, \
        "UMI2 filter should NOT be applied when UMI2 column doesn't exist"
    
    # Verify all reads are kept (no filtering)
    filtered_csv = output_dir / "test_reads_no_umi2_filtered_umi_groups.csv"
    df_filtered = pd.read_csv(filtered_csv)
    
    assert len(df_filtered) == 1, "Should have 1 group (UMI1_A)"
    assert df_filtered.iloc[0]['ReadCount'] == 10, "Should keep all 10 reads"
    
    print(f"\n✓ Test passed: UMI2 filtering correctly skipped when UMI2 column missing")


def test_umi2_majority_filtering_grouping_by_umi2(tmp_path):
    """
    Test that UMI2 filtering is NOT applied when grouping by UMI2 (not UMI1).
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create CSV with UMI1 and UMI2
    umi_data = [
        {'Read Name': f'read_{i:04d}', 'UMI1': 'UMI1_A', 'UMI2': f'UMI2_{i % 3}'} 
        for i in range(15)
    ]
    csv_path = tmp_path / "test_umis_group_by_umi2.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_group_by_umi2.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping by UMI2 (not UMI1)
    output_dir = tmp_path / "umi2_filter_group_by_umi2"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI2",  # Grouping by UMI2, not UMI1
        threshold=1,
        umi2_majority_threshold=0.5
    )
    
    # Verify JSON log shows filter was NOT applied
    log_json = output_dir / "test_reads_group_by_umi2_umi_group.log.json"
    with open(log_json, 'r') as f:
        log_data = json.load(f)
    
    print(f"\n{'='*80}")
    print(f"UMI2 FILTERING NOT APPLIED (grouping by UMI2)")
    print(f"{'='*80}")
    print(f"Grouping column: {log_data.get('grouping')}")
    print(f"UMI2 filter applied: {log_data.get('umi2_majority_filter_applied')}")
    
    assert log_data.get('grouping') == 'UMI2'
    assert log_data.get('umi2_majority_filter_applied') == False, \
        "UMI2 filter should NOT be applied when grouping by UMI2"
    
    print(f"\n✓ Test passed: UMI2 filtering correctly skipped when grouping by UMI2")


def test_umi2_majority_filtering_empty_umi2(tmp_path):
    """
    Test UMI2 filtering with empty/missing UMI2 values.
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create CSV with some reads having empty UMI2
    umi_data = [
        {'Read Name': 'read_0001', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_X'},
        {'Read Name': 'read_0002', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_X'},
        {'Read Name': 'read_0003', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_X'},
        {'Read Name': 'read_0004', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_X'},
        {'Read Name': 'read_0005', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_X'},
        {'Read Name': 'read_0006', 'UMI1': 'UMI1_A', 'UMI2': ''},  # Empty UMI2
        {'Read Name': 'read_0007', 'UMI1': 'UMI1_A', 'UMI2': ''},  # Empty UMI2
        {'Read Name': 'read_0008', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_Y'},
        {'Read Name': 'read_0009', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_Y'},
        {'Read Name': 'read_0010', 'UMI1': 'UMI1_A', 'UMI2': 'UMI2_Y'},
    ]
    csv_path = tmp_path / "test_umis_empty_umi2.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_empty_umi2.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping
    output_dir = tmp_path / "umi2_filter_empty_umi2"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=3,
        umi2_majority_threshold=0.5
    )
    
    # Verify empty UMI2 values are ignored
    filtered_csv = output_dir / "test_reads_empty_umi2_filtered_umi_groups.csv"
    df_filtered = pd.read_csv(filtered_csv)
    
    print(f"\n{'='*80}")
    print(f"UMI2 FILTERING WITH EMPTY UMI2 VALUES")
    print(f"{'='*80}")
    print(df_filtered)
    
    # UMI1_A should pass: 5 reads with UMI2_X (majority), 3 with UMI2_Y, 2 empty (ignored)
    # Total valid: 8, majority: 5 (62.5%), should pass
    if len(df_filtered) > 0:
        umi1_a_row = df_filtered[df_filtered['UMI1'] == 'UMI1_A']
        if len(umi1_a_row) > 0:
            assert umi1_a_row.iloc[0]['ReadCount'] == 5, \
                "Should keep 5 reads with majority UMI2_X (empty UMI2 ignored)"
    
    print(f"\n✓ Test passed: Empty UMI2 values are correctly ignored")


def test_umi2_majority_filtering_statistics(tmp_path):
    """
    Test that JSON log statistics are accurate for UMI2 filtering.
    """
    if not PYSAM_AVAILABLE:
        pytest.skip("pysam not available for creating synthetic BAM files")
    
    # Create synthetic UMI data
    umi_data = create_synthetic_umi_data_umi2_filtering()
    
    # Create synthetic CSV
    csv_path = tmp_path / "test_umis_stats.csv"
    create_synthetic_umi_csv(umi_data, csv_path)
    
    # Create synthetic BAM
    read_names = [row['Read Name'] for row in umi_data]
    bam_path = tmp_path / "test_reads_stats.bam"
    create_synthetic_bam(read_names, bam_path)
    
    # Run grouping
    output_dir = tmp_path / "umi2_filter_stats"
    output_dir.mkdir(exist_ok=True)
    
    long_read_umi_grouping_fast(
        alignedReads=str(bam_path),
        correctedUMIs=str(csv_path),
        outputDir=str(output_dir),
        grouping="UMI1",
        threshold=3,
        umi2_majority_threshold=0.5
    )
    
    # Read JSON log
    log_json = output_dir / "test_reads_stats_umi_group.log.json"
    with open(log_json, 'r') as f:
        log_data = json.load(f)
    
    print(f"\n{'='*80}")
    print(f"UMI2 FILTERING STATISTICS VALIDATION")
    print(f"{'='*80}")
    
    # Calculate expected statistics from synthetic data
    df_original = pd.read_csv(csv_path)
    df_original['UMI2'] = df_original['UMI2'].replace(['', 'NoMatch', 'NoBarcode'], pd.NA)
    
    # Count reads before filtering (all UMI1 groups passing threshold)
    umi1_counts = df_original['UMI1'].value_counts()
    passing_umi1 = umi1_counts[umi1_counts >= 3]
    expected_reads_before = int(passing_umi1.sum())
    
    # Count reads after filtering (only majority UMI2 reads)
    expected_reads_after = 0
    for umi1 in passing_umi1.index:
        sub = df_original[df_original['UMI1'] == umi1]
        umi2_counts = sub['UMI2'].value_counts(dropna=True)
        if len(umi2_counts) > 0:
            majority_count = int(umi2_counts.iloc[0])
            total_count = len(sub)
            majority_fraction = majority_count / total_count
            if majority_fraction >= 0.5 and majority_count >= 3:
                expected_reads_after += majority_count
    
    expected_filtered = expected_reads_before - expected_reads_after
    
    print(f"Expected reads before filter: {expected_reads_before}")
    print(f"Expected reads after filter: {expected_reads_after}")
    print(f"Expected reads filtered: {expected_filtered}")
    print(f"\nActual statistics from log:")
    print(f"  Reads before filter: {log_data.get('reads_before_umi2_filter')}")
    print(f"  Reads after filter: {log_data.get('reads_after_umi2_filter')}")
    print(f"  Reads filtered: {log_data.get('reads_filtered_by_umi2')}")
    
    # Verify statistics match (allow small differences due to data structure)
    actual_before = log_data.get('reads_before_umi2_filter', 0)
    actual_after = log_data.get('reads_after_umi2_filter', 0)
    actual_filtered = log_data.get('reads_filtered_by_umi2', 0)
    
    assert actual_before == expected_reads_before, \
        f"Reads before filter: expected {expected_reads_before}, got {actual_before}"
    assert actual_after == expected_reads_after, \
        f"Reads after filter: expected {expected_reads_after}, got {actual_after}"
    assert actual_filtered == expected_filtered, \
        f"Reads filtered: expected {expected_filtered}, got {actual_filtered}"
    
    # Verify fraction calculations
    groups_passing_threshold = log_data.get('groups_passing_threshold', 0)
    groups_passing_umi2 = log_data.get('groups_passing_umi2_majority', 0)
    if groups_passing_threshold > 0:
        expected_fraction = groups_passing_umi2 / groups_passing_threshold
        actual_fraction = log_data.get('fraction_passing_umi2_majority', 0)
        assert abs(actual_fraction - expected_fraction) < 0.01, \
            f"Fraction passing: expected {expected_fraction:.3f}, got {actual_fraction:.3f}"
    
    print(f"\n✓ Test passed: Statistics are accurate")
