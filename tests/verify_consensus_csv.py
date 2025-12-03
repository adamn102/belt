#!/usr/bin/env python3
"""
Simple verification script to test consensus sequence data CSV generation.
"""

import sys
import os
import tempfile
import json
from pathlib import Path
import pandas as pd

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_summary import generate_consensus_sequence_data

def test_consensus_csv_generation():
    """Test consensus sequence data CSV generation"""
    print("=" * 70)
    print("TESTING CONSENSUS SEQUENCE DATA CSV GENERATION")
    print("=" * 70)
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp())
    
    try:
        # Create mock sample results
        sample_results = [
            {
                'sample_name': 'Sample_A',
                'per_read_analysis': [
                    {
                        'consensus_id': 'UMI_001',
                        'sequence_length': 1000,
                        'mutation_count': 5,
                        'mutation_rate': 0.005,
                        'mutations': [
                            {'type': 'substitution'},
                            {'type': 'substitution'},
                            {'type': 'insertion'},
                            {'type': 'deletion'},
                            {'type': 'deletion'}
                        ]
                    },
                    {
                        'consensus_id': 'UMI_002',
                        'sequence_length': 950,
                        'mutation_count': 2,
                        'mutation_rate': 0.0021,
                        'mutations': [
                            {'type': 'substitution'},
                            {'type': 'insertion'}
                        ]
                    }
                ]
            },
            {
                'sample_name': 'Sample_B',
                'per_read_analysis': [
                    {
                        'consensus_id': 'UMI_003',
                        'sequence_length': 1200,
                        'mutation_count': 8,
                        'mutation_rate': 0.0067,
                        'mutations': [
                            {'type': 'substitution'},
                            {'type': 'substitution'},
                            {'type': 'substitution'},
                            {'type': 'insertion'},
                            {'type': 'insertion'},
                            {'type': 'deletion'},
                            {'type': 'deletion'},
                            {'type': 'deletion'}
                        ]
                    }
                ]
            }
        ]
        
        print(f"Created mock data with {len(sample_results)} samples")
        total_consensus = sum(len(s['per_read_analysis']) for s in sample_results)
        print(f"Total consensus sequences: {total_consensus}")
        
        # Generate CSV
        print("\nGenerating consensus sequence data CSV...")
        generate_consensus_sequence_data(sample_results, temp_dir)
        
        # Check file was created
        csv_file = temp_dir / "consensus_sequence_data.csv"
        if not csv_file.exists():
            print("✗ FAIL: CSV file was not created")
            return False
        
        print(f"✓ CSV file created: {csv_file}")
        
        # Load and display CSV content
        df = pd.read_csv(csv_file)
        print(f"\nCSV contains {len(df)} rows and {len(df.columns)} columns")
        print(f"Columns: {list(df.columns)}")
        
        print("\nCSV Content:")
        print("-" * 70)
        print(df.to_string(index=False))
        
        # Verify data
        print("\nVerification:")
        print("-" * 70)
        
        # Check Sample_A UMI_001
        row1 = df[df['UMI_ID'] == 'UMI_001'].iloc[0]
        print(f"Sample_A UMI_001:")
        print(f"  Total Bases: {row1['Total_Bases']} (expected: 1000)")
        print(f"  Total Mutations: {row1['Total_Mutations']} (expected: 5)")
        print(f"  Substitutions: {row1['Substitutions']} (expected: 2)")
        print(f"  Insertions: {row1['Insertions']} (expected: 1)")
        print(f"  Deletions: {row1['Deletions']} (expected: 2)")
        
        # Check Sample_A UMI_002
        row2 = df[df['UMI_ID'] == 'UMI_002'].iloc[0]
        print(f"\nSample_A UMI_002:")
        print(f"  Total Bases: {row2['Total_Bases']} (expected: 950)")
        print(f"  Total Mutations: {row2['Total_Mutations']} (expected: 2)")
        print(f"  Substitutions: {row2['Substitutions']} (expected: 1)")
        print(f"  Insertions: {row2['Insertions']} (expected: 1)")
        print(f"  Deletions: {row2['Deletions']} (expected: 0)")
        
        # Check Sample_B UMI_003
        row3 = df[df['UMI_ID'] == 'UMI_003'].iloc[0]
        print(f"\nSample_B UMI_003:")
        print(f"  Total Bases: {row3['Total_Bases']} (expected: 1200)")
        print(f"  Total Mutations: {row3['Total_Mutations']} (expected: 8)")
        print(f"  Substitutions: {row3['Substitutions']} (expected: 3)")
        print(f"  Insertions: {row3['Insertions']} (expected: 2)")
        print(f"  Deletions: {row3['Deletions']} (expected: 3)")
        
        # Verify all expected values
        checks = [
            (row1['Total_Bases'], 1000, "Sample_A UMI_001 Total_Bases"),
            (row1['Total_Mutations'], 5, "Sample_A UMI_001 Total_Mutations"),
            (row1['Substitutions'], 2, "Sample_A UMI_001 Substitutions"),
            (row1['Insertions'], 1, "Sample_A UMI_001 Insertions"),
            (row1['Deletions'], 2, "Sample_A UMI_001 Deletions"),
            (row2['Total_Bases'], 950, "Sample_A UMI_002 Total_Bases"),
            (row2['Total_Mutations'], 2, "Sample_A UMI_002 Total_Mutations"),
            (row2['Substitutions'], 1, "Sample_A UMI_002 Substitutions"),
            (row2['Insertions'], 1, "Sample_A UMI_002 Insertions"),
            (row2['Deletions'], 0, "Sample_A UMI_002 Deletions"),
            (row3['Total_Bases'], 1200, "Sample_B UMI_003 Total_Bases"),
            (row3['Total_Mutations'], 8, "Sample_B UMI_003 Total_Mutations"),
            (row3['Substitutions'], 3, "Sample_B UMI_003 Substitutions"),
            (row3['Insertions'], 2, "Sample_B UMI_003 Insertions"),
            (row3['Deletions'], 3, "Sample_B UMI_003 Deletions")
        ]
        
        all_passed = True
        for actual, expected, description in checks:
            if actual == expected:
                print(f"✓ {description}: {actual}")
            else:
                print(f"✗ {description}: {actual} (expected {expected})")
                all_passed = False
        
        if all_passed:
            print("\n✓ ALL TESTS PASSED! Consensus sequence data CSV generation is working correctly.")
            return True
        else:
            print("\n✗ SOME TESTS FAILED!")
            return False
            
    finally:
        # Clean up
        import shutil
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

def main():
    """Run the test"""
    success = test_consensus_csv_generation()
    return 0 if success else 1

if __name__ == '__main__':
    sys.exit(main())
