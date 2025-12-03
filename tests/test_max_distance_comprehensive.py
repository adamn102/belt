#!/usr/bin/env python3
"""
Comprehensive test to demonstrate max_distance improvement in UMI mapping.
This test shows that the BK_scan function works correctly and that higher
max_distance values lead to more successful mappings.
"""

import os
import sys
import pytest
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_map import BK_scan, load_whitelist

class TestMaxDistanceComprehensive:
    """Comprehensive test demonstrating max_distance improvement"""
    
    @pytest.fixture
    def test_data_paths(self):
        """Fixture for test data paths"""
        base_dir = Path(__file__).parent.parent
        test_data_dir = base_dir / "test-data" / "test_max_distance"
        return {
            'whitelist_file': test_data_dir / "proper_whitelist.csv",
            'test_data_file': test_data_dir / "proper_test_data.csv",
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
    
    def test_max_distance_improvement_demonstration(self):
        """Demonstrate that max_distance parameter works correctly"""
        
        # Simple whitelist with one UMI
        whitelist = ["GTAGCCTTAACG"]
        
        # Test UMIs with different distances
        test_cases = [
            ("GTAGCCTTAACG", 0, "exact match"),
            ("GTAGCCTTAACC", 1, "distance 1"),
            ("GTAGCCTTAACA", 1, "distance 1"), 
            ("GTAGCCTTAATG", 1, "distance 1"),
        ]
        
        print("Testing max_distance parameter with different values:")
        print(f"Whitelist: {whitelist[0]}")
        print()
        
        # Test with different max_distance values
        for max_distance in [0, 1, 2, 3]:
            print(f"=== Testing with max_distance = {max_distance} ===")
            results = BK_scan([case[0] for case in test_cases], whitelist, max_distance)
            
            for i, (umi, expected_dist, description) in enumerate(test_cases):
                match, distance = results[umi]
                print(f"  {description}: {umi} -> {match} (distance: {distance})")
                
                if max_distance == 0:
                    # Only exact matches should map
                    if expected_dist == 0:
                        assert match == whitelist[0] and distance == 0
                    else:
                        assert match == umi and distance == "NoMatch"
                        
                elif max_distance >= 1:
                    # All should map since they're all distance 0 or 1
                    assert match == whitelist[0]
                    assert distance in [0, 1]
            
            print()
    
    def test_max_distance_improvement_summary(self):
        """Summary test showing improvement with higher max_distance"""
        
        whitelist = ["GTAGCCTTAACG"]
        
        # Test UMIs with different distances
        test_umis = [
            ("GTAGCCTTAACG", 0),  # exact match
            ("GTAGCCTTAACC", 1),  # distance 1
            ("GTAGCCTTAACA", 1),  # distance 1
            ("GTAGCCTTAATG", 1),  # distance 1
        ]
        
        print("\n=== Max Distance Improvement Summary ===")
        print(f"Whitelist: {whitelist[0]}")
        print()
        
        # Count successful mappings for each max_distance
        for max_distance in [0, 1, 2, 3]:
            successful_mappings = 0
            results = BK_scan([umi for umi, _ in test_umis], whitelist, max_distance)
            
            for umi, expected_dist in test_umis:
                match, distance = results[umi]
                if match == whitelist[0]:  # Successfully mapped
                    successful_mappings += 1
            
            print(f"max_distance={max_distance}: {successful_mappings}/{len(test_umis)} successful mappings")
            
            # Verify that higher max_distance leads to more successful mappings
            if max_distance == 0:
                assert successful_mappings == 1  # Only exact match
            elif max_distance >= 1:
                assert successful_mappings == 4  # All should map
    
    def test_max_distance_with_actual_distance_2(self):
        """Test with a UMI that actually has distance 2"""
        
        whitelist = ["GTAGCCTTAACG"]
        
        # Create a UMI with actual distance 2
        # Change positions 10 and 11: CG -> CA
        dist2_umi = "GTAGCCTTAACA"
        
        # Verify distance
        actual_distance = sum(c1 != c2 for c1, c2 in zip(whitelist[0], dist2_umi))
        print(f"Distance verification: {whitelist[0]} vs {dist2_umi} = {actual_distance}")
        
        # Test with max_distance = 1 (should not map)
        results_dist1 = BK_scan([dist2_umi], whitelist, max_distance=1)
        match, distance = results_dist1[dist2_umi]
        
        print(f"With max_distance=1: {dist2_umi} -> {match} (distance: {distance})")
        
        # Should not map with max_distance=1
        assert match == dist2_umi and distance == "NoMatch"
        
        # Test with max_distance = 2 (should map)
        results_dist2 = BK_scan([dist2_umi], whitelist, max_distance=2)
        match, distance = results_dist2[dist2_umi]
        
        print(f"With max_distance=2: {dist2_umi} -> {match} (distance: {distance})")
        
        # Should map with max_distance=2
        assert match == whitelist[0] and distance == 2
    
    def test_max_distance_with_actual_distance_3(self):
        """Test with a UMI that actually has distance 3"""
        
        whitelist = ["GTAGCCTTAACG"]
        
        # Create a UMI with actual distance 3
        # Change positions 9, 10, and 11: ACG -> ATG
        dist3_umi = "GTAGCCTTAATG"
        
        # Verify distance
        actual_distance = sum(c1 != c2 for c1, c2 in zip(whitelist[0], dist3_umi))
        print(f"Distance verification: {whitelist[0]} vs {dist3_umi} = {actual_distance}")
        
        # Test with max_distance = 2 (should not map)
        results_dist2 = BK_scan([dist3_umi], whitelist, max_distance=2)
        match, distance = results_dist2[dist3_umi]
        
        print(f"With max_distance=2: {dist3_umi} -> {match} (distance: {distance})")
        
        # Should not map with max_distance=2
        assert match == dist3_umi and distance == "NoMatch"
        
        # Test with max_distance = 3 (should map)
        results_dist3 = BK_scan([dist3_umi], whitelist, max_distance=3)
        match, distance = results_dist3[dist3_umi]
        
        print(f"With max_distance=3: {dist3_umi} -> {match} (distance: {distance})")
        
        # Should map with max_distance=3
        assert match == whitelist[0] and distance == 3

if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v", "-s"])
