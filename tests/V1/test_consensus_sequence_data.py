#!/usr/bin/env python3
"""
Test suite for consensus sequence data CSV generation in mutation_summary.py
"""

import unittest
import sys
import os
import tempfile
import json
from pathlib import Path
import pandas as pd

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_summary import generate_consensus_sequence_data


class TestConsensusSequenceDataCSV(unittest.TestCase):
    """Test consensus sequence data CSV generation"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_consensus_sequence_data_generation(self):
        """Test generation of consensus sequence data CSV"""
        # Create mock sample results
        sample_results = [
            {
                'sample_name': 'Sample1',
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
                'sample_name': 'Sample2',
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
        
        # Generate CSV
        generate_consensus_sequence_data(sample_results, self.temp_dir)
        
        # Check file was created
        csv_file = self.temp_dir / "consensus_sequence_data.csv"
        self.assertTrue(csv_file.exists(), "CSV file should be created")
        
        # Load and verify CSV content
        df = pd.read_csv(csv_file)
        
        # Check columns
        expected_columns = [
            'Sample', 'UMI_ID', 'Total_Bases', 'Total_Mutations', 
            'Mutation_Rate', 'Substitutions', 'Insertions', 'Deletions'
        ]
        self.assertEqual(list(df.columns), expected_columns)
        
        # Check number of rows
        self.assertEqual(len(df), 3, "Should have 3 consensus sequences")
        
        # Check Sample1 UMI_001 data
        row1 = df[df['UMI_ID'] == 'UMI_001'].iloc[0]
        self.assertEqual(row1['Sample'], 'Sample1')
        self.assertEqual(row1['Total_Bases'], 1000)
        self.assertEqual(row1['Total_Mutations'], 5)
        self.assertEqual(row1['Substitutions'], 2)
        self.assertEqual(row1['Insertions'], 1)
        self.assertEqual(row1['Deletions'], 2)
        
        # Check Sample1 UMI_002 data
        row2 = df[df['UMI_ID'] == 'UMI_002'].iloc[0]
        self.assertEqual(row2['Sample'], 'Sample1')
        self.assertEqual(row2['Total_Bases'], 950)
        self.assertEqual(row2['Total_Mutations'], 2)
        self.assertEqual(row2['Substitutions'], 1)
        self.assertEqual(row2['Insertions'], 1)
        self.assertEqual(row2['Deletions'], 0)
        
        # Check Sample2 UMI_003 data
        row3 = df[df['UMI_ID'] == 'UMI_003'].iloc[0]
        self.assertEqual(row3['Sample'], 'Sample2')
        self.assertEqual(row3['Total_Bases'], 1200)
        self.assertEqual(row3['Total_Mutations'], 8)
        self.assertEqual(row3['Substitutions'], 3)
        self.assertEqual(row3['Insertions'], 2)
        self.assertEqual(row3['Deletions'], 3)
    
    def test_empty_sample_results(self):
        """Test handling of empty sample results"""
        sample_results = []
        
        generate_consensus_sequence_data(sample_results, self.temp_dir)
        
        csv_file = self.temp_dir / "consensus_sequence_data.csv"
        self.assertTrue(csv_file.exists(), "Empty CSV file should be created")
        
        df = pd.read_csv(csv_file)
        self.assertEqual(len(df), 0, "Empty CSV should have 0 rows")
        
        # Check columns are still present
        expected_columns = [
            'Sample', 'UMI_ID', 'Total_Bases', 'Total_Mutations', 
            'Mutation_Rate', 'Substitutions', 'Insertions', 'Deletions'
        ]
        self.assertEqual(list(df.columns), expected_columns)
    
    def test_sample_with_no_per_read_analysis(self):
        """Test handling of sample with no per_read_analysis"""
        sample_results = [
            {
                'sample_name': 'Sample1',
                'per_read_analysis': []  # Empty per_read_analysis
            },
            {
                'sample_name': 'Sample2',
                'per_read_analysis': [
                    {
                        'consensus_id': 'UMI_001',
                        'sequence_length': 1000,
                        'mutation_count': 3,
                        'mutation_rate': 0.003,
                        'mutations': [
                            {'type': 'substitution'},
                            {'type': 'substitution'},
                            {'type': 'insertion'}
                        ]
                    }
                ]
            }
        ]
        
        generate_consensus_sequence_data(sample_results, self.temp_dir)
        
        csv_file = self.temp_dir / "consensus_sequence_data.csv"
        df = pd.read_csv(csv_file)
        
        # Should only have data for Sample2
        self.assertEqual(len(df), 1, "Should have 1 consensus sequence")
        self.assertEqual(df.iloc[0]['Sample'], 'Sample2')
        self.assertEqual(df.iloc[0]['UMI_ID'], 'UMI_001')
        self.assertEqual(df.iloc[0]['Substitutions'], 2)
        self.assertEqual(df.iloc[0]['Insertions'], 1)
        self.assertEqual(df.iloc[0]['Deletions'], 0)
    
    def test_missing_mutation_data(self):
        """Test handling of missing mutation data"""
        sample_results = [
            {
                'sample_name': 'Sample1',
                'per_read_analysis': [
                    {
                        'consensus_id': 'UMI_001',
                        'sequence_length': 1000,
                        'mutation_count': 0,
                        'mutation_rate': 0.0,
                        'mutations': []  # No mutations
                    },
                    {
                        'consensus_id': 'UMI_002',
                        'sequence_length': 950,
                        'mutation_count': 2,
                        'mutation_rate': 0.0021,
                        # Missing 'mutations' key
                    }
                ]
            }
        ]
        
        generate_consensus_sequence_data(sample_results, self.temp_dir)
        
        csv_file = self.temp_dir / "consensus_sequence_data.csv"
        df = pd.read_csv(csv_file)
        
        # Check UMI_001 (no mutations)
        row1 = df[df['UMI_ID'] == 'UMI_001'].iloc[0]
        self.assertEqual(row1['Total_Mutations'], 0)
        self.assertEqual(row1['Substitutions'], 0)
        self.assertEqual(row1['Insertions'], 0)
        self.assertEqual(row1['Deletions'], 0)
        
        # Check UMI_002 (missing mutations data)
        row2 = df[df['UMI_ID'] == 'UMI_002'].iloc[0]
        self.assertEqual(row2['Total_Mutations'], 2)
        self.assertEqual(row2['Substitutions'], 0)  # Should default to 0
        self.assertEqual(row2['Insertions'], 0)
        self.assertEqual(row2['Deletions'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
