#!/usr/bin/env python3
"""
Comprehensive test for max_distance parameter in UMI mapping.
Tests that different max_distance values result in improved mapping to the whitelist.
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

from umi_map import BK_scan, load_whitelist

class TestMaxDistanceUMIMapping:
    """Test class for max_distance parameter testing"""
    
    @pytest.fixture
    def test_data_paths(self):
        """Fixture for test data paths"""
        base_dir = Path(__file__).parent.parent
        test_data_dir = base_dir / "test-data" / "test_max_distance"
        return {
            'whitelist_file': test_data_dir / "max_distance_whitelist.csv",
            'test_data_file': test_data_dir / "max_distance_test_data.csv",
            'test_data_dir': test_data_dir
        }
    
    @pytest.fixture
    def whitelist_data(self, test_data_paths):
        """Load whitelist data"""
        return load_whitelist(str(test_data_paths['whitelist_file']))
    
    @pytest.fixture
    def test_data(self, test_data_paths):
        """Load test data"""
        return pd.read_csv(test_data_paths['test_data_file'])
    
    def test_max_distance_0_exact_matches_only(self, whitelist_data, test_data):
        """Test that max_distance=0 only matches exact sequences"""
        max_distance = 0
        
        # Test each UMI column
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name == "exact_match_1":
                    # Exact match should map to itself
                    assert match == umi, f"Exact match {umi} should map to itself, got {match}"
                    assert distance == 0, f"Exact match should have distance 0, got {distance}"
                else:
                    # Non-exact matches should not map (return original)
                    assert match == umi, f"Non-exact match {umi} should not map with distance=0, got {match}"
                    assert distance == "NoMatch", f"Non-exact match should have 'NoMatch', got {distance}"
    
    def test_max_distance_1_single_mutations(self, whitelist_data, test_data):
        """Test that max_distance=1 matches single mutations"""
        max_distance = 1
        
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name == "exact_match_1":
                    assert match == umi and distance == 0
                elif read_name == "distance_1_1":
                    # Should map to base UMI with distance 1
                    base_umi = whitelist_data[col].iloc[0]
                    assert match == base_umi, f"Distance 1 UMI {umi} should map to {base_umi}, got {match}"
                    assert distance == 1, f"Distance 1 UMI should have distance 1, got {distance}"
                else:
                    # Higher distance UMIs should not map
                    assert match == umi, f"Higher distance UMI {umi} should not map with distance=1"
                    assert distance == "NoMatch", f"Higher distance UMI should have 'NoMatch'"
    
    def test_max_distance_2_double_mutations(self, whitelist_data, test_data):
        """Test that max_distance=2 matches double mutations"""
        max_distance = 2
        
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name in ["exact_match_1", "distance_1_1", "distance_2_1"]:
                    # Should map to base UMI
                    base_umi = whitelist_data[col].iloc[0]
                    assert match == base_umi, f"UMI {umi} should map to {base_umi}, got {match}"
                    if read_name == "exact_match_1":
                        assert distance == 0
                    elif read_name == "distance_1_1":
                        assert distance == 1
                    elif read_name == "distance_2_1":
                        assert distance == 2
                else:
                    # Higher distance UMIs should not map
                    assert match == umi, f"Higher distance UMI {umi} should not map with distance=2"
                    assert distance == "NoMatch"
    
    def test_max_distance_3_triple_mutations(self, whitelist_data, test_data):
        """Test that max_distance=3 matches triple mutations"""
        max_distance = 3
        
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name in ["exact_match_1", "distance_1_1", "distance_2_1", "distance_3_1"]:
                    # Should map to base UMI
                    base_umi = whitelist_data[col].iloc[0]
                    assert match == base_umi, f"UMI {umi} should map to {base_umi}, got {match}"
                    if read_name == "exact_match_1":
                        assert distance == 0
                    elif read_name == "distance_1_1":
                        assert distance == 1
                    elif read_name == "distance_2_1":
                        assert distance == 2
                    elif read_name == "distance_3_1":
                        assert distance == 3
                else:
                    # Distance 4 UMIs should not map
                    assert match == umi, f"Distance 4 UMI {umi} should not map with distance=3"
                    assert distance == "NoMatch"
    
    def test_mapping_improvement_with_higher_distance(self, whitelist_data, test_data):
        """Test that higher max_distance values improve mapping accuracy"""
        
        # Test a specific UMI that should benefit from higher distance
        test_umi = "GTAGCCTTAACT"  # distance_3_1 UMI1
        col = 'UMI1'
        base_umi = whitelist_data[col].iloc[0]
        
        # Test with different max_distance values
        for max_distance in [0, 1, 2, 3]:
            results = BK_scan([test_umi], whitelist_data[col], max_distance)
            match, distance = results[test_umi]
            
            if max_distance < 3:
                # Should not map with lower distances
                assert match == test_umi, f"Should not map with distance={max_distance}"
                assert distance == "NoMatch", f"Should have 'NoMatch' with distance={max_distance}"
            else:
                # Should map with distance=3
                assert match == base_umi, f"Should map to {base_umi} with distance=3, got {match}"
                assert distance == 3, f"Should have distance 3, got {distance}"
    
    def test_mixed_distance_cases(self, whitelist_data, test_data):
        """Test mixed distance cases"""
        max_distance = 3
        
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name in ["mixed_1", "mixed_2"]:
                    # Mixed cases should map to base UMI
                    base_umi = whitelist_data[col].iloc[0]
                    assert match == base_umi, f"Mixed case UMI {umi} should map to {base_umi}"
                    assert isinstance(distance, int), f"Distance should be integer, got {distance}"
                    assert 0 <= distance <= 3, f"Distance should be 0-3, got {distance}"
    
    def test_no_match_cases(self, whitelist_data, test_data):
        """Test cases that should not match even with high distance"""
        max_distance = 3
        
        for col in ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']:
            results = BK_scan(test_data[col], whitelist_data[col], max_distance)
            
            for i, umi in enumerate(test_data[col]):
                read_name = test_data.iloc[i]['Read Name']
                match, distance = results[umi]
                
                if read_name == "distance_4_1":
                    # Distance 4 should not map even with max_distance=3
                    assert match == umi, f"Distance 4 UMI {umi} should not map"
                    assert distance == "NoMatch", f"Distance 4 UMI should have 'NoMatch'"

if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v"])
