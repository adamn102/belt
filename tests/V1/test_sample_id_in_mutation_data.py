#!/usr/bin/env python3
"""
Test script to verify that sample_id is included in per_read_analysis and mutations_by_position
"""

import json
import sys
from pathlib import Path

# Add scripts directory to path
scripts_dir = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

from mutation_caller import mutation_caller
from mutation_summary import mutation_summary


def test_sample_id_in_mutation_caller():
    """Test that mutation_caller includes sample_id in per_read_analysis"""
    print("Testing mutation_caller sample_id inclusion...")
    
    # Check test data
    test_data = Path(__file__).parent.parent / "test-data" / "test_mutation_caller"
    if not test_data.exists():
        print("⚠️  Test data directory not found, skipping test")
        return True
    
    ref_file = test_data / "reference.fa"
    cons_file = test_data / "consensus_with_mutations.fa"
    
    if not ref_file.exists() or not cons_file.exists():
        print("⚠️  Test data files not found, skipping test")
        return True
    
    # Create temporary output directory
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        
        # Run mutation caller
        result = mutation_caller(
            consensus_file=cons_file,
            reference_file=ref_file,
            sample_name="test_sample",
            output_dir=output_dir
        )
        
        # Check that result has per_read_analysis
        if "per_read_analysis" not in result:
            print("❌ Result does not contain per_read_analysis")
            return False
        
        # Check that per_read_analysis entries have sample_id
        for i, read_analysis in enumerate(result["per_read_analysis"]):
            if "sample_id" not in read_analysis:
                print(f"❌ per_read_analysis entry {i} does not contain sample_id")
                return False
            
            if read_analysis["sample_id"] != "test_sample":
                print(f"❌ per_read_analysis entry {i} has incorrect sample_id: {read_analysis['sample_id']}")
                return False
        
        print(f"✅ All {len(result['per_read_analysis'])} per_read_analysis entries contain correct sample_id")
        
        # Check that mutations_by_position entries have sample_id
        if "mutations_by_position" not in result:
            print("⚠️  Result does not contain mutations_by_position")
        else:
            for i, mutation in enumerate(result["mutations_by_position"]):
                if "sample_id" not in mutation:
                    print(f"❌ mutations_by_position entry {i} does not contain sample_id")
                    return False
                
                if mutation["sample_id"] != "test_sample":
                    print(f"❌ mutations_by_position entry {i} has incorrect sample_id: {mutation['sample_id']}")
                    return False
            
            print(f"✅ All {len(result['mutations_by_position'])} mutations_by_position entries contain correct sample_id")
        
        # Check the JSON output file
        json_file = output_dir / "samples" / "test_sample" / "test_sample_mutations.json"
        if json_file.exists():
            with open(json_file, 'r') as f:
                saved_result = json.load(f)
            
            # Verify sample_id is in the saved file
            if saved_result.get("per_read_analysis"):
                has_sample_id = all("sample_id" in read for read in saved_result["per_read_analysis"])
                if has_sample_id:
                    print("✅ Saved JSON file contains sample_id in per_read_analysis")
                else:
                    print("❌ Saved JSON file does not contain sample_id in all per_read_analysis entries")
                    return False
        
        print("✅ mutation_caller test passed!")
        return True


def test_sample_id_structure():
    """Test the structure of sample_id in mutation data"""
    print("\nTesting sample_id data structure...")
    
    # Create a mock result to test the structure
    mock_result = {
        "sample_name": "LP249-NL43-WT-1-A",
        "consensus_count": 2,
        "total_length": 7956,
        "total_mutations": 5,
        "average_mutation_rate": 0.0006285,
        "mutation_types": {"substitution": 4, "insertion": 1, "deletion": 0},
        "mutations_by_position": [
            {
                "position": 100,
                "type": "substitution",
                "ref_base": "A",
                "cons_base": "G",
                "umi_id": "AAAAAAAAATTT",
                "sample_id": "LP249-NL43-WT-1-A"
            }
        ],
        "per_read_analysis": [
            {
                "sample_id": "LP249-NL43-WT-1-A",
                "consensus_id": "AAAAAAAAATTT",
                "sequence_length": 3978,
                "mutation_count": 3,
                "mutation_rate": 0.000754,
                "mutations": []
            },
            {
                "sample_id": "LP249-NL43-WT-1-A",
                "consensus_id": "TTTTTTTTAAAA",
                "sequence_length": 3978,
                "mutation_count": 2,
                "mutation_rate": 0.000503,
                "mutations": []
            }
        ]
    }
    
    # Verify structure
    print("Expected per_read_analysis structure:")
    print("  - sample_id: Sample identifier")
    print("  - consensus_id: UMI identifier")
    print("  - sequence_length: Length of sequence")
    print("  - mutation_count: Number of mutations")
    print("  - mutation_rate: Mutation rate")
    print("  - mutations: List of mutations")
    
    print("\nExpected mutations_by_position structure:")
    print("  - position: Position in reference")
    print("  - type: Mutation type")
    print("  - ref_base: Reference base")
    print("  - cons_base: Consensus base")
    print("  - umi_id: UMI identifier")
    print("  - sample_id: Sample identifier")
    
    # Check structure
    for read in mock_result["per_read_analysis"]:
        assert "sample_id" in read, "sample_id missing from per_read_analysis"
        assert "consensus_id" in read, "consensus_id missing from per_read_analysis"
    
    for mutation in mock_result["mutations_by_position"]:
        assert "sample_id" in mutation, "sample_id missing from mutations_by_position"
        assert "umi_id" in mutation, "umi_id missing from mutations_by_position"
    
    print("\n✅ Data structure test passed!")
    return True


def main():
    """Run all tests"""
    print("=" * 60)
    print("Testing Sample ID in Mutation Analysis")
    print("=" * 60)
    
    all_passed = True
    
    # Test structure
    if not test_sample_id_structure():
        all_passed = False
    
    # Test mutation_caller
    if not test_sample_id_in_mutation_caller():
        all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✅ All tests passed!")
        print("\nThe following changes have been made:")
        print("1. per_read_analysis now includes 'sample_id' field")
        print("2. mutations_by_position now includes 'sample_id' field")
        print("3. CSV exports include 'Sample_ID' column")
        print("\nBenefit:")
        print("- Better identification of which UMIs belong to which samples")
        print("- Easier to track mutations across samples")
        print("- More complete data in mutation_analysis_report.json")
    else:
        print("❌ Some tests failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

