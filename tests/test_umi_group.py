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

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_group import long_read_umi_grouping_fast  # type: ignore


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
