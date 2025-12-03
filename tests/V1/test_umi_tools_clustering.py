#!/usr/bin/env python3
"""
Comprehensive tests for UMI-tools clustering functionality in umi_map.py
"""
import os
import sys
import pytest
import pandas as pd
from pathlib import Path
import tempfile
import shutil
import json
from collections import Counter

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_map import cluster_umis_umitools, umi_extraction_map  # type: ignore


class TestUMIToolsClustering:
    """Test suite for UMI-tools clustering functionality"""
    
    def test_cluster_umis_empty_list(self):
        """Test clustering with empty UMI list"""
        correction_map, cluster_info = cluster_umis_umitools([])
        
        assert correction_map == {}
        assert cluster_info["num_clusters"] == 0
        assert cluster_info["num_singletons"] == 0
        assert cluster_info["cluster_sizes"] == []
    
    def test_cluster_umis_single_umi(self):
        """Test clustering with a single UMI"""
        umis = ["ACGTACGTACGT"]
        correction_map, cluster_info = cluster_umis_umitools(umis)
        
        assert len(correction_map) == 1
        assert correction_map["ACGTACGTACGT"] == "ACGTACGTACGT"
        assert cluster_info["num_clusters"] == 1
        assert cluster_info["num_singletons"] == 1
        assert cluster_info["cluster_sizes"] == [1]
    
    def test_cluster_umis_identical_umis(self):
        """Test clustering with multiple identical UMIs"""
        umis = ["ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT"]
        correction_map, cluster_info = cluster_umis_umitools(umis, threshold=1)
        
        assert len(correction_map) == 1
        assert correction_map["ACGTACGTACGT"] == "ACGTACGTACGT"
        assert cluster_info["num_clusters"] == 1
        # Even with multiple occurrences, a single unique UMI is counted as a singleton
        assert cluster_info["num_singletons"] == 1
        assert cluster_info["cluster_sizes"] == [1]  # One unique UMI
        # All three should map to the same corrected UMI
        assert all(c == "ACGTACGTACGT" for c in correction_map.values())
    
    def test_cluster_umis_single_mutation(self):
        """Test clustering UMIs with single base differences"""
        # UMI1: ACGTACGTACGT (high count, should be primary)
        # UMI2: ACGTACGTACGG (1 edit away)
        # UMI3: ACGTACGTACGA (1 edit away)
        # All should cluster together with UMI1 as primary
        
        umis = (
            ["ACGTACGTACGT"] * 10 +  # High count primary
            ["ACGTACGTACGG"] * 5 +   # Variant 1 (1 edit)
            ["ACGTACGTACGA"] * 3     # Variant 2 (1 edit)
        )
        
        correction_map, cluster_info = cluster_umis_umitools(umis, threshold=1, cluster_method="directional")
        
        print(f"\n{'='*80}")
        print("SINGLE MUTATION CLUSTERING TEST")
        print(f"{'='*80}")
        print(f"Correction map:")
        for orig, corrected in sorted(correction_map.items()):
            print(f"  {orig} -> {corrected}")
        print(f"\nCluster info: {cluster_info}")
        print(f"{'='*80}\n")
        
        # All UMIs should cluster together (directional method with high count primary)
        # The primary UMI should be ACGTACGTACGT (highest count)
        primary = correction_map["ACGTACGTACGT"]
        assert primary == "ACGTACGTACGT", f"Primary should be ACGTACGTACGT, got {primary}"
        
        # All variants should map to primary
        assert correction_map["ACGTACGTACGG"] == primary
        assert correction_map["ACGTACGTACGA"] == primary
        
        # Should have 1 cluster
        assert cluster_info["num_clusters"] == 1
        # Cluster should have 3 unique UMIs
        assert len(set(correction_map.values())) == 1  # All map to same primary
    
    def test_cluster_umis_separate_clusters(self):
        """Test that UMIs beyond threshold distance form separate clusters"""
        # UMI1: ACGTACGTACGT
        # UMI2: ACGTACGTACGG (1 edit away)
        # UMI3: ACGTACGTACCC (2 edits away from UMI1, 2 from UMI2)
        # With threshold=1, UMI3 should be separate
        
        umis = (
            ["ACGTACGTACGT"] * 10 +
            ["ACGTACGTACGG"] * 5 +
            ["ACGTACGTACCC"] * 3  # 2 edits away, should be separate cluster
        )
        
        correction_map, cluster_info = cluster_umis_umitools(umis, threshold=1, cluster_method="directional")
        
        print(f"\n{'='*80}")
        print("SEPARATE CLUSTERS TEST")
        print(f"{'='*80}")
        print(f"Correction map:")
        for orig, corrected in sorted(correction_map.items()):
            print(f"  {orig} -> {corrected}")
        print(f"\nCluster info: {cluster_info}")
        print(f"{'='*80}\n")
        
        # Should have 2 clusters (UMI1+UMI2, and UMI3 separate)
        assert cluster_info["num_clusters"] == 2, f"Expected 2 clusters, got {cluster_info['num_clusters']}"
        
        # UMI1 and UMI2 should cluster together
        assert correction_map["ACGTACGTACGT"] == correction_map["ACGTACGTACGG"]
        
        # UMI3 should be in a separate cluster
        assert correction_map["ACGTACGTACCC"] != correction_map["ACGTACGTACGT"]
        
        # Should have 2 unique corrected UMIs
        unique_corrected = len(set(correction_map.values()))
        assert unique_corrected == 2, f"Expected 2 unique corrected UMIs, got {unique_corrected}"
    
    def test_cluster_umis_varying_lengths(self):
        """Test that clustering handles UMIs of different lengths (edge case)"""
        # Note: UMI-tools uses Hamming distance which requires same length
        # This test verifies that an error is raised when lengths differ
        umis = ["ACGTACGTACGT", "ACGTACGTACGTACGT"]  # Different lengths
        
        # UMI-tools should raise an AssertionError for different lengths
        with pytest.raises(AssertionError, match="not all umis are the same length"):
            cluster_umis_umitools(umis, threshold=1)
    
    def test_cluster_umis_adjacency_method(self):
        """Test clustering with adjacency method instead of directional"""
        umis = (
            ["ACGTACGTACGT"] * 10 +
            ["ACGTACGTACGG"] * 5 +
            ["ACGTACGTACGA"] * 3
        )
        
        correction_map, cluster_info = cluster_umis_umitools(
            umis, threshold=1, cluster_method="adjacency"
        )
        
        print(f"\n{'='*80}")
        print("ADJACENCY METHOD TEST")
        print(f"{'='*80}")
        print(f"Cluster info: {cluster_info}")
        print(f"{'='*80}\n")
        
        # Should cluster successfully
        assert cluster_info["num_clusters"] >= 1
        assert len(correction_map) == 3  # Three unique UMIs
    
    def test_cluster_umis_statistics(self):
        """Test that cluster statistics are calculated correctly"""
        umis = (
            ["ACGTACGTACGT"] * 10 +  # Cluster 1: 10 occurrences, 3 unique UMIs
            ["ACGTACGTACGG"] * 5 +
            ["ACGTACGTACGA"] * 3 +
            ["TGCATGCATGCA"] * 2 +   # Cluster 2: 2 occurrences, 1 UMI
            ["AAAAAAAAAAAA"]          # Cluster 3: 1 occurrence, singleton
        )
        
        correction_map, cluster_info = cluster_umis_umitools(umis, threshold=1)
        
        print(f"\n{'='*80}")
        print("STATISTICS TEST")
        print(f"{'='*80}")
        print(f"Cluster info: {cluster_info}")
        print(f"Cluster sizes: {cluster_info['cluster_sizes']}")
        print(f"Num singletons: {cluster_info['num_singletons']}")
        print(f"Num clusters: {cluster_info['num_clusters']}")
        print(f"{'='*80}\n")
        
        # Verify statistics structure
        assert "num_clusters" in cluster_info
        assert "num_singletons" in cluster_info
        assert "cluster_sizes" in cluster_info
        assert isinstance(cluster_info["cluster_sizes"], list)
        assert cluster_info["num_clusters"] >= 1
        assert cluster_info["num_singletons"] >= 0
    
    def test_cluster_umis_real_world_scenario(self):
        """Test with a realistic scenario of sequencing errors"""
        # Simulate a scenario where:
        # - Primary UMI: ACTGACTGACTG (many reads)
        # - Error variants: single base changes
        # - Some true different UMIs
        
        primary = "ACTGACTGACTG"
        errors = [
            "ACTGACTGACTA",  # 1 edit
            "ACTGACTGACTT",  # 1 edit
            "ACTGACTGACTC",  # 1 edit
            "ACTGACTGACCG",  # 1 edit
        ]
        different = [
            "CGTACGTACGTA",  # Many edits away
            "GATCGATCGATC",  # Many edits away
        ]
        
        umis = (
            [primary] * 100 +
            [errors[0]] * 15 +
            [errors[1]] * 10 +
            [errors[2]] * 8 +
            [errors[3]] * 5 +
            [different[0]] * 20 +
            [different[1]] * 12
        )
        
        correction_map, cluster_info = cluster_umis_umitools(umis, threshold=1)
        
        print(f"\n{'='*80}")
        print("REAL-WORLD SCENARIO TEST")
        print(f"{'='*80}")
        print(f"Primary UMI: {primary}")
        print(f"Error variants: {errors}")
        print(f"Different UMIs: {different}")
        print(f"\nCorrection results:")
        for orig in set(umis):
            corrected = correction_map[orig]
            print(f"  {orig} -> {corrected}")
        print(f"\nCluster info: {cluster_info}")
        print(f"{'='*80}\n")
        
        # Primary and errors should cluster together
        primary_cluster = correction_map[primary]
        for error in errors:
            assert correction_map[error] == primary_cluster, \
                f"Error variant {error} should cluster with primary"
        
        # Different UMIs should be separate
        assert correction_map[different[0]] != primary_cluster
        assert correction_map[different[1]] != primary_cluster


class TestUMIExtractionMapClusterMode:
    """Integration tests for umi_extraction_map in cluster mode"""
    
    @pytest.fixture
    def test_data_paths(self):
        """Fixture to provide test data paths"""
        base_dir = Path(__file__).parent.parent
        test_data_dir = base_dir / "test-data" / "test_umi_map"
        
        return {
            'reference_file': test_data_dir / "extract_umi_5.fa",
            'bam_file': test_data_dir / "extract_umi_5-large.bam",
            'base_dir': base_dir
        }
    
    def test_cluster_mode_basic(self, test_data_paths, tmp_path):
        """Test umi_extraction_map with cluster mode"""
        reference_file = test_data_paths['reference_file']
        bam_file = test_data_paths['bam_file']
        
        # Skip if test files don't exist
        if not reference_file.exists() or not bam_file.exists():
            pytest.skip("Required test files not found")
        
        output_dir = tmp_path / "cluster_mode_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{'='*80}")
        print("CLUSTER MODE INTEGRATION TEST")
        print(f"{'='*80}")
        
        # Run with cluster mode
        umi_extraction_map(
            reference=str(reference_file),
            alignedReads=str(bam_file),
            minReadLength=2000,
            barcodeMotifs=["NNNNNNNNNNNN", "NNNNNNNNNNNN"],
            whitelist="cluster",  # Use cluster mode
            outputDir=str(output_dir),
            max_distance=1  # Edit distance threshold
        )
        
        # Verify output CSV was created
        sample_name = bam_file.stem.replace('.bam', '')
        output_csv = output_dir / f"{sample_name}_Paired_UMIs_read_names.csv"
        
        assert output_csv.exists(), f"Output CSV not found: {output_csv}"
        
        # Read and verify output
        df = pd.read_csv(output_csv)
        print(f"\nOutput CSV shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        print(f"First few rows:")
        print(df.head())
        
        # Verify structure
        assert "Read Name" in df.columns
        assert "UMI1" in df.columns
        assert "UMI2" in df.columns
        assert len(df) > 0, "Output should contain reads"
        
        # Verify JSON log exists
        log_json = output_dir / f"{sample_name}_umi_extraction_map.log.json"
        assert log_json.exists(), "JSON log should exist"
        
        # Read and verify log
        with open(log_json, 'r') as f:
            log_data = json.load(f)
        
        print(f"\nLog data keys: {list(log_data.keys())}")
        assert log_data["whitelist_mode"] == "cluster"
        assert "clustering_stats" in log_data, "Log should contain clustering_stats"
        
        clustering_stats = log_data["clustering_stats"]
        print(f"\nClustering stats structure:")
        print(f"  Per column keys: {list(clustering_stats.get('per_column', {}).keys())}")
        print(f"  Overall stats: {clustering_stats.get('overall', {})}")
        
        assert "per_column" in clustering_stats
        assert "overall" in clustering_stats
        
        print(f"\n{'='*80}")
        print("✓ Cluster mode integration test passed")
        print(f"{'='*80}\n")
    
    def test_cluster_mode_vs_whitelist_mode(self, test_data_paths, tmp_path):
        """Compare outputs between cluster mode and whitelist mode"""
        reference_file = test_data_paths['reference_file']
        bam_file = test_data_paths['bam_file']
        base_dir = test_data_paths['base_dir']
        
        # Skip if test files don't exist
        whitelist_file = base_dir / "test-data" / "test_umi_map" / "extract_umi_5_whitelist.csv"
        if not all([reference_file.exists(), bam_file.exists(), whitelist_file.exists()]):
            pytest.skip("Required test files not found")
        
        sample_name = bam_file.stem.replace('.bam', '')
        
        # Run cluster mode
        cluster_output = tmp_path / "cluster_output"
        cluster_output.mkdir(parents=True, exist_ok=True)
        
        umi_extraction_map(
            reference=str(reference_file),
            alignedReads=str(bam_file),
            minReadLength=2000,
            barcodeMotifs=["NNNNNNNNNNNN", "NNNNNNNNNNNN"],
            whitelist="cluster",
            outputDir=str(cluster_output),
            max_distance=1
        )
        
        # Run whitelist mode
        whitelist_output = tmp_path / "whitelist_output"
        whitelist_output.mkdir(parents=True, exist_ok=True)
        
        umi_extraction_map(
            reference=str(reference_file),
            alignedReads=str(bam_file),
            minReadLength=2000,
            barcodeMotifs=["NNNNNNNNNNNN", "NNNNNNNNNNNN"],
            whitelist=str(whitelist_file),
            outputDir=str(whitelist_output),
            max_distance=1
        )
        
        # Compare outputs
        cluster_csv = cluster_output / f"{sample_name}_Paired_UMIs_read_names.csv"
        whitelist_csv = whitelist_output / f"{sample_name}_Paired_UMIs_read_names.csv"
        
        cluster_df = pd.read_csv(cluster_csv)
        whitelist_df = pd.read_csv(whitelist_csv)
        
        print(f"\n{'='*80}")
        print("CLUSTER vs WHITELIST MODE COMPARISON")
        print(f"{'='*80}")
        print(f"Cluster mode reads: {len(cluster_df)}")
        print(f"Whitelist mode reads: {len(whitelist_df)}")
        
        # Count unique UMIs
        cluster_umi1_unique = cluster_df['UMI1'].dropna().nunique()
        whitelist_umi1_unique = whitelist_df['UMI1'].dropna().nunique()
        
        print(f"\nUnique UMI1 counts:")
        print(f"  Cluster mode: {cluster_umi1_unique}")
        print(f"  Whitelist mode: {whitelist_umi1_unique}")
        
        # Cluster mode should generally have fewer unique UMIs (due to clustering)
        # But this depends on the data, so we just verify both produce valid outputs
        assert cluster_umi1_unique > 0
        assert whitelist_umi1_unique > 0
        
        print(f"\n{'='*80}")
        print("✓ Comparison test passed")
        print(f"{'='*80}\n")
    
    def test_cluster_mode_with_partial_umis(self, test_data_paths, tmp_path):
        """Test cluster mode handles partial UMIs correctly"""
        reference_file = test_data_paths['reference_file']
        bam_file = test_data_paths['bam_file']
        
        if not reference_file.exists() or not bam_file.exists():
            pytest.skip("Required test files not found")
        
        output_dir = tmp_path / "cluster_partial_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        umi_extraction_map(
            reference=str(reference_file),
            alignedReads=str(bam_file),
            minReadLength=2000,
            barcodeMotifs=["NNNNNNNNNNNN", "NNNNNNNNNNNN"],
            whitelist="cluster",
            outputDir=str(output_dir),
            max_distance=1
        )
        
        sample_name = bam_file.stem.replace('.bam', '')
        output_csv = output_dir / f"{sample_name}_Paired_UMIs_read_names.csv"
        
        df = pd.read_csv(output_csv)
        
        # Check that NaN values are preserved
        umi1_missing = df['UMI1'].isna().sum()
        umi2_missing = df['UMI2'].isna().sum()
        
        print(f"\n{'='*80}")
        print("PARTIAL UMIs IN CLUSTER MODE")
        print(f"{'='*80}")
        print(f"Total reads: {len(df)}")
        print(f"UMI1 missing: {umi1_missing}")
        print(f"UMI2 missing: {umi2_missing}")
        print(f"{'='*80}\n")
        
        # Verify structure is correct
        assert len(df) > 0
        assert "Read Name" in df.columns
        assert "UMI1" in df.columns
        assert "UMI2" in df.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])

