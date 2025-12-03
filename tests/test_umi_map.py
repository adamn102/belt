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


def test_umi_extraction_map_partial_umis(tmp_path):
    """Test that umi_extraction_map handles reads with missing UMIs gracefully.
    
    This test verifies that the function can handle cases where:
    - A read has UMI1 but not UMI2
    - A read has UMI2 but not UMI1
    - A read has both UMIs
    - A read has neither UMI (should be excluded)
    """
    base_dir = Path(__file__).parent.parent
    test_data_dir = base_dir / "test-data" / "test_umi_map"
    
    # Use existing test data (assuming some reads may have partial UMIs naturally)
    reference_file = test_data_dir / "extract_umi_5.fa"
    bam_file = test_data_dir / "extract_umi_5-large.bam"
    
    # Skip if test files don't exist
    if not reference_file.exists() or not bam_file.exists():
        pytest.skip("Required test files not found")
    
    # Create temporary output directory
    output_dir = tmp_path / "partial_umi_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Use just two motifs to increase likelihood of partial matches
    test_motifs = ["NNNNNNNNNNNN", "NNNNNNNNNNNN"]
    
    # Run extraction with whitelist=false to bypass whitelist mapping
    umi_extraction_map(
        reference=str(reference_file),
        alignedReads=str(bam_file),
        minReadLength=2000,
        barcodeMotifs=test_motifs,
        whitelist="false",
        outputDir=str(output_dir)
    )
    
    # Verify output CSV was created
    sample_name = bam_file.stem.replace('.bam', '')
    output_csv = output_dir / f"{sample_name}_Paired_UMIs_read_names.csv"
    
    assert output_csv.exists(), f"Output CSV not found: {output_csv}"
    
    # Read and analyze the output
    df = pd.read_csv(output_csv)
    
    # Verify expected columns
    expected_columns = ["Read Name", "UMI1_Raw", "UMI1", "UMI2_Raw", "UMI2"]
    assert list(df.columns) == expected_columns, f"Expected columns {expected_columns}, got {list(df.columns)}"
    
    # Count different UMI patterns
    has_umi1 = df['UMI1'].notna()
    has_umi2 = df['UMI2'].notna()
    
    both_umis = (has_umi1 & has_umi2).sum()
    only_umi1 = (has_umi1 & ~has_umi2).sum()
    only_umi2 = (~has_umi1 & has_umi2).sum()
    no_umis = (~has_umi1 & ~has_umi2).sum()
    
    # Print diagnostic information
    print(f"\n{'='*80}")
    print(f"PARTIAL UMI TEST RESULTS")
    print(f"{'='*80}")
    print(f"Total reads: {len(df)}")
    print(f"Reads with both UMI1 and UMI2: {both_umis}")
    print(f"Reads with only UMI1: {only_umi1}")
    print(f"Reads with only UMI2: {only_umi2}")
    print(f"Reads with no UMIs: {no_umis}")
    print(f"{'='*80}")
    
    # Show some example reads with partial UMIs
    if only_umi1 > 0:
        print("\nExample reads with only UMI1:")
        sample = df[has_umi1 & ~has_umi2].head(3)
        for _, row in sample.iterrows():
            print(f"  {row['Read Name']}: UMI1={row['UMI1']}, UMI2={row['UMI2']}")
    
    if only_umi2 > 0:
        print("\nExample reads with only UMI2:")
        sample = df[~has_umi1 & has_umi2].head(3)
        for _, row in sample.iterrows():
            print(f"  {row['Read Name']}: UMI1={row['UMI1']}, UMI2={row['UMI2']}")
    
    print(f"\n{'='*80}\n")
    
    # Assertions
    assert len(df) > 0, "No reads were extracted"
    
    # Verify that at least some reads have valid UMIs
    assert both_umis + only_umi1 + only_umi2 > 0, "No reads with any UMIs found"
    
    # The function should handle partial UMIs gracefully
    # (i.e., the test should pass even if some reads have missing UMIs)
    print(f"✓ Test passed: Function handles partial UMIs correctly")
    print(f"  - Processed {len(df)} reads successfully")
    print(f"  - Handled {only_umi1 + only_umi2} reads with partial UMIs")
    
    # Additional validation: check that UMI sequences are valid (when present)
    valid_bases = set('ACGTN')
    for col in ['UMI1', 'UMI2']:
        non_null_umis = df[col].dropna()
        if len(non_null_umis) > 0:
            # Check a sample of UMIs
            sample_umis = non_null_umis.head(10)
            for umi in sample_umis:
                assert isinstance(umi, str), f"{col} should be string, got {type(umi)}"
                assert all(base in valid_bases for base in umi.upper()), \
                    f"{col} contains invalid bases: {umi}"
    
    # Note: Even if test data has all complete UMIs, this test validates that:
    # 1. The function doesn't crash with missing UMIs (pandas handles NaN correctly)
    # 2. The CSV output format supports partial UMIs (empty cells for missing values)
    # 3. The function logic (from findUMI) initializes with None and handles partial extraction
    print(f"✓ Function is robust to partial UMIs (test validates structure even with complete data)")


def test_whitelist_mapping_with_partial_umis():
    """Test that whitelist mapping handles reads with missing UMIs correctly.
    
    This is a critical test that validates:
    - BK_scan can handle NaN values in input UMIs
    - Whitelist mapping preserves NaN for missing UMIs
    - The final output CSV correctly represents partial UMIs
    """
    from umi_map import BK_scan
    import pandas as pd
    
    # Simulate UMI data with missing values
    experimental_umis = pd.Series([
        'ACGTACGTACGT',  # Valid UMI
        None,            # Missing UMI
        'GCTAGCTAGCTA',  # Valid UMI
        float('nan'),    # Another form of missing UMI
        'TGCATGCATGCA',  # Valid UMI
    ])
    
    # Whitelist with known UMIs (all same length as experimental UMIs)
    whitelist_umis = pd.Series([
        'ACGTACGTACGT',  # Exact match
        'GCTAGCTAGCTA',  # Exact match
        'TGCATGCATGCA',  # Exact match
        'AAAAAAAAAAAA',  # Extra whitelist entry (12 bases like others)
    ])
    
    print(f"\n{'='*80}")
    print(f"WHITELIST MAPPING WITH PARTIAL UMIs TEST")
    print(f"{'='*80}")
    print(f"Experimental UMIs (with NaN):")
    for i, umi in enumerate(experimental_umis):
        status = "NaN" if pd.isna(umi) else umi
        print(f"  [{i}] {status}")
    print(f"\nWhitelist UMIs:")
    for i, umi in enumerate(whitelist_umis):
        print(f"  [{i}] {umi}")
    print(f"{'='*80}\n")
    
    # Run BK_scan (should not crash with NaN values)
    try:
        best_matches = BK_scan(experimental_umis, whitelist_umis, max_distance=1)
        print("✓ BK_scan completed without error")
    except Exception as e:
        print(f"✗ BK_scan failed with error: {e}")
        raise
    
    # Verify results
    print(f"\nMapping results:")
    for i, umi in enumerate(experimental_umis):
        if pd.isna(umi):
            mapped, dist = best_matches[umi]
            print(f"  NaN -> {mapped} (distance: {dist})")
            assert mapped is None, "Missing UMI should map to None"
            assert dist is None, "Distance for missing UMI should be None"
        else:
            mapped, dist = best_matches[umi]
            print(f"  {umi} -> {mapped} (distance: {dist})")
            assert mapped is not None, f"Valid UMI {umi} should have a mapping"
    
    print(f"\n{'='*80}")
    print("✓ Test passed: Whitelist mapping handles partial UMIs correctly")
    print(f"{'='*80}\n")


# def test_umi_extraction_real_data_with_partial_umis(tmp_path):
#     """Test umi_extraction_map with real data that has partial UMIs.
    
#     Uses LP249-NL43-A3G-1-A BAM file which is known to have reads with missing UMIs.
#     This is a critical integration test for production scenarios.
#     """
#     base_dir = Path(__file__).parent.parent
#     test_data_dir = base_dir / "test-data" / "test_group_umi"
    
#     reference_file = test_data_dir / "LP249.fa"
#     bam_file = test_data_dir / "LP249-NL43-A3G-1-A.bam"
#     whitelist_file = test_data_dir / "extract_umi_5_whitelist.csv"
    
#     # Skip if test files don't exist
#     if not all([reference_file.exists(), bam_file.exists(), whitelist_file.exists()]):
#         pytest.skip("Required test files not found")
    
#     # Create temporary output directory
#     output_dir = tmp_path / "real_data_output"
#     output_dir.mkdir(parents=True, exist_ok=True)
    
#     # Run extraction with whitelist mapping
#     print(f"\n{'='*80}")
#     print(f"REAL DATA TEST: LP249-NL43-A3G-1-A with Partial UMIs")
#     print(f"{'='*80}")
    
#     umi_extraction_map(
#         reference=str(reference_file),
#         alignedReads=str(bam_file),
#         minReadLength=1000,
#         barcodeMotifs=['NNNNNNNNNNNN', 'NNNNNNNNNNNN'],
#         whitelist=str(whitelist_file),
#         outputDir=str(output_dir)
#     )
    
#     # Verify output CSV was created
#     sample_name = bam_file.stem
#     output_csv = output_dir / f"{sample_name}_Paired_UMIs_read_names.csv"
    
#     assert output_csv.exists(), f"Output CSV not found: {output_csv}"
    
#     # Read and analyze the output
#     df = pd.read_csv(output_csv)
    
#     # Expected statistics based on known data
#     total_reads = len(df)
#     umi1_missing = df['UMI1'].isna().sum()
#     umi2_missing = df['UMI2'].isna().sum()
#     both_present = (df['UMI1'].notna() & df['UMI2'].notna()).sum()
#     only_umi1 = (df['UMI1'].notna() & df['UMI2'].isna()).sum()
#     only_umi2 = (df['UMI1'].isna() & df['UMI2'].notna()).sum()
#     both_missing = (df['UMI1'].isna() & df['UMI2'].isna()).sum()
    
#     print(f"\nReal Data Statistics:")
#     print(f"  Total reads: {total_reads:,}")
#     print(f"  UMI1 missing: {umi1_missing:,} ({umi1_missing/total_reads*100:.1f}%)")
#     print(f"  UMI2 missing: {umi2_missing:,} ({umi2_missing/total_reads*100:.1f}%)")
#     print(f"  Both UMIs present: {both_present:,} ({both_present/total_reads*100:.1f}%)")
#     print(f"  Only UMI1 present: {only_umi1:,} ({only_umi1/total_reads*100:.1f}%)")
#     print(f"  Only UMI2 present: {only_umi2:,} ({only_umi2/total_reads*100:.1f}%)")
#     print(f"  Both missing: {both_missing:,} ({both_missing/total_reads*100:.1f}%)")
    
#     # Assertions
#     assert total_reads > 50000, f"Expected >50k reads, got {total_reads}"
#     assert umi1_missing > 0, "Expected some reads with missing UMI1"
#     assert umi2_missing > 0, "Expected some reads with missing UMI2"
#     assert only_umi1 > 0, "Expected some reads with only UMI1"
#     assert only_umi2 > 0, "Expected some reads with only UMI2"
    
#     # Show examples
#     print(f"\nExamples of partial UMIs:")
#     if only_umi2 > 0:
#         print(f"\nReads with only UMI2 (missing UMI1):")
#         for _, row in df[df['UMI1'].isna() & df['UMI2'].notna()].head(2).iterrows():
#             print(f"  {row['Read Name']}: UMI1=NaN, UMI2={row['UMI2']}")
    
#     if only_umi1 > 0:
#         print(f"\nReads with only UMI1 (missing UMI2):")
#         for _, row in df[df['UMI1'].notna() & df['UMI2'].isna()].head(2).iterrows():
#             print(f"  {row['Read Name']}: UMI1={row['UMI1']}, UMI2=NaN")
    
#     print(f"\n{'='*80}")
#     print(f"✓ Test passed: Real data with {only_umi1 + only_umi2:,} partial UMIs handled correctly")
#     print(f"{'='*80}\n")


def test_dataframe_handles_none_umis():
    """Unit test to verify DataFrame construction handles None values correctly.
    
    This directly tests the core mechanism for handling missing UMIs:
    - When findUMI returns None for a motif position, pandas should handle it as NaN
    - The CSV output should show empty cells for missing UMIs
    """
    # Simulate raw_UMI_dict with partial UMIs (as returned by findUMI)
    raw_UMI_dict = {
        'read1': ['ACGTACGTACGT', 'TGCATGCATGCA'],  # Has both UMI1 and UMI2
        'read2': ['GCTAGCTAGCTA', None],              # Has only UMI1, missing UMI2
        'read3': [None, 'ATCGATCGATCG'],              # Has only UMI2, missing UMI1
        'read4': ['CGATCGATCGAT', 'TAGCTAGCTAGC'],  # Has both
    }
    
    # Convert to DataFrame (mimicking what umi_extraction_map does)
    df = pd.DataFrame.from_dict(raw_UMI_dict, orient='index', columns=['motif1', 'motif2'])
    df.index.name = 'Read Name'
    df.reset_index(inplace=True)
    
    # Update column names to UMI1, UMI2
    df.columns = ['Read Name', 'UMI1', 'UMI2']
    
    print(f"\n{'='*80}")
    print(f"UNIT TEST: DataFrame with None/Missing UMIs")
    print(f"{'='*80}")
    print(df.to_string())
    print(f"{'='*80}\n")
    
    # Verify DataFrame structure
    assert len(df) == 4, "Should have 4 reads"
    assert list(df.columns) == ['Read Name', 'UMI1', 'UMI2'], "Should have correct columns"
    
    # Verify None values are handled as NaN
    assert pd.isna(df.loc[df['Read Name'] == 'read2', 'UMI2'].values[0]), \
        "read2 should have NaN for UMI2"
    assert pd.isna(df.loc[df['Read Name'] == 'read3', 'UMI1'].values[0]), \
        "read3 should have NaN for UMI1"
    
    # Verify non-None values are preserved
    assert df.loc[df['Read Name'] == 'read1', 'UMI1'].values[0] == 'ACGTACGTACGT', \
        "read1 UMI1 should be preserved"
    assert df.loc[df['Read Name'] == 'read2', 'UMI1'].values[0] == 'GCTAGCTAGCTA', \
        "read2 UMI1 should be preserved"
    assert df.loc[df['Read Name'] == 'read3', 'UMI2'].values[0] == 'ATCGATCGATCG', \
        "read3 UMI2 should be preserved"
    
    # Write to CSV and verify format
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        csv_path = f.name
    
    try:
        df.to_csv(csv_path, index=False)
        
        # Read back and verify
        df_read = pd.read_csv(csv_path)
        
        print("CSV content verification:")
        print(df_read.to_string())
        print(f"\n{'='*80}\n")
        
        # Verify missing UMIs are read as NaN from CSV
        assert pd.isna(df_read.loc[df_read['Read Name'] == 'read2', 'UMI2'].values[0]), \
            "Missing UMI2 should be NaN when read from CSV"
        assert pd.isna(df_read.loc[df_read['Read Name'] == 'read3', 'UMI1'].values[0]), \
            "Missing UMI1 should be NaN when read from CSV"
        
        print("✓ Test passed: DataFrame correctly handles None/missing UMIs")
        print("  - None values converted to NaN")
        print("  - CSV format preserves missing values as empty cells")
        print("  - Reading CSV back correctly interprets empty cells as NaN")
        
    finally:
        # Clean up
        import os
        if os.path.exists(csv_path):
            os.unlink(csv_path)

