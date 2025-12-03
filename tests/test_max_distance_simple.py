#!/usr/bin/env python3
"""
Simple working test to demonstrate max_distance improvement in UMI mapping.
"""

import os
import sys
import pytest
import pandas as pd
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from umi_map import BK_scan

class TestMaxDistanceSimple:
    """Simple test demonstrating max_distance improvement"""
    
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

if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v", "-s"])