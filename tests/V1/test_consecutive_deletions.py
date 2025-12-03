#!/usr/bin/env python3
"""
Test suite for consecutive deletion grouping in mutation_caller.py

This test verifies that consecutive deletions are properly grouped together
as a single deletion event with size information.
"""

import unittest
import sys
import os
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_caller import load_sequences, align_sequences, find_mutations
from Bio.Seq import Seq


class TestConsecutiveDeletions(unittest.TestCase):
    """Test consecutive deletion handling"""
    
    def test_single_deletion(self):
        """Test that a single deletion is identified correctly"""
        # Reference: ATCG
        # Consensus: ATC (G deleted at position 3)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("ATC")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        print(f"\n--- Single Deletion Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        for d in deletions:
            print(f"  Position {d['ref_position']}: {d['ref_base']} deleted")
        
        # Currently this will find 1 deletion
        self.assertEqual(len(deletions), 1)
        self.assertEqual(deletions[0]["ref_base"], "G")
    
    def test_consecutive_double_deletion_grouped(self):
        """Test that consecutive deletions are grouped as one event"""
        # Reference: ATCG
        # Consensus: AT (CG deleted at positions 2-3)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("AT")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        print(f"\n--- Double Deletion Test (After Fix) ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        for d in deletions:
            print(f"  Position {d['ref_position']}: {d.get('deleted_sequence', d['ref_base'])} deleted (size: {d.get('size', 1)})")
        
        # NEW BEHAVIOR: Should find 1 grouped deletion
        self.assertEqual(len(deletions), 1, 
                        "Should group consecutive deletions into one event")
        self.assertEqual(deletions[0]["deleted_sequence"], "CG")
        self.assertEqual(deletions[0]["size"], 2)
    
    def test_consecutive_triple_deletion_grouped(self):
        """Test that triple consecutive deletions are grouped"""
        # Reference: ATCGATCG
        # Consensus: ATCG (ATCG deleted at positions 4-7)
        ref_seq = Seq("ATCGATCG")
        cons_seq = Seq("ATCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        print(f"\n--- Multiple Deletion Test (After Fix) ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        for i, d in enumerate(deletions):
            print(f"  Deletion {i+1}: Position {d['ref_position']}, Sequence: {d.get('deleted_sequence', 'N/A')}, Size: {d.get('size', 1)}")
        
        # NEW BEHAVIOR: Should find 1 grouped deletion
        self.assertEqual(len(deletions), 1,
                        "Should group consecutive deletions into one event")
        self.assertEqual(deletions[0]["size"], 4, "Should be 4 consecutive deleted bases")
    
    def test_non_consecutive_deletions(self):
        """Test that non-consecutive deletions remain separate"""
        # Reference: ATCGATCG
        # Consensus: TCGTCG (A deleted at pos 0, A deleted at pos 4)
        ref_seq = Seq("ATCGATCG")
        cons_seq = Seq("TCGTCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        print(f"\n--- Non-Consecutive Deletions Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        for i, d in enumerate(deletions):
            print(f"  Deletion {i+1}: Position {d['ref_position']}, Base {d['ref_base']}")
        
        # These deletions are NOT consecutive, so they should remain separate
        self.assertEqual(len(deletions), 2, 
                        "Non-consecutive deletions should remain separate")
    
    def test_large_consecutive_deletion(self):
        """Test a large consecutive deletion (5+ bases)"""
        # Reference: ATCGATCGATCGATCG
        # Consensus: ATCGATCG (ATCGATCG deleted - 8 bases)
        ref_seq = Seq("ATCGATCGATCGATCG")
        cons_seq = Seq("ATCGATCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        print(f"\n--- Large Consecutive Deletion Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        if deletions:
            total_size = sum(d.get('size', 1) for d in deletions)
            print(f"Total deleted bases: {total_size}")
            for d in deletions:
                print(f"  Position: {d['ref_position']}, Size: {d.get('size', 1)}, Sequence: {d.get('deleted_sequence', 'N/A')[:20]}...")
        
        # NEW BEHAVIOR: Should be grouped as 1 deletion of size 8
        self.assertEqual(len(deletions), 1,
                        "Should group all consecutive deletions into one event")
        self.assertEqual(deletions[0]["size"], 8, "Should have size of 8")
        self.assertEqual(deletions[0]["deleted_sequence"], "ATCGATCG")


class TestConsecutiveDeletionsWithRealData(unittest.TestCase):
    """Test with real FASTA files"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_data_dir = Path(__file__).parent.parent / "test-data" / "test_mutation_caller"
    
    def test_with_fasta_files(self):
        """Test consecutive deletions with FASTA files"""
        ref_file = self.test_data_dir / "consecutive_deletions_reference.fa"
        cons_file = self.test_data_dir / "consecutive_deletions_consensus.fa"
        
        if not ref_file.exists() or not cons_file.exists():
            self.skipTest("Test data files not found")
        
        # Load sequences
        ref_seqs = load_sequences(ref_file)
        cons_seqs = load_sequences(cons_file)
        
        ref_seq = list(ref_seqs.values())[0].seq
        
        print(f"\n--- Testing with FASTA files ---")
        print(f"Reference: {ref_seq} (length: {len(ref_seq)})")
        
        for cons_id, cons_record in cons_seqs.items():
            cons_seq = cons_record.seq
            print(f"\n{cons_id}:")
            print(f"  Consensus: {cons_seq} (length: {len(cons_seq)})")
            print(f"  Expected deletion size: {len(ref_seq) - len(cons_seq)}")
            
            ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
            mutations = find_mutations(ref_aligned, cons_aligned)
            
            deletions = [m for m in mutations if m["type"] == "deletion"]
            print(f"  Deletions found: {len(deletions)}")
            
            if deletions:
                print(f"  Deleted bases: {[d['ref_base'] for d in deletions]}")
                print(f"  Positions: {[d['ref_position'] for d in deletions]}")


class TestDesiredDeletionBehavior(unittest.TestCase):
    """Test cases for desired behavior after implementing grouping"""
    
    def test_consecutive_deletions_grouped(self):
        """Consecutive deletions should be grouped as one event"""
        # Reference: ATCG
        # Consensus: AT (CG deleted - size 2)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("AT")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        
        print(f"\n--- Grouped Deletion Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        for d in deletions:
            print(f"  Position: {d['ref_position']}, Size: {d.get('size', 'N/A')}, Sequence: {d.get('deleted_sequence', 'N/A')}")
        
        # NEW BEHAVIOR: Should find 1 grouped deletion
        self.assertEqual(len(deletions), 1, "Should group consecutive deletions")
        self.assertEqual(deletions[0]["size"], 2, "Should report deletion size")
        self.assertEqual(deletions[0]["deleted_sequence"], "CG", 
                        "Should report deleted sequence")
    
    def test_deletion_has_size_field(self):
        """Deletions should have a 'size' field"""
        ref_seq = Seq("ATCGATCG")
        cons_seq = Seq("ATCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        deletions = [m for m in mutations if m["type"] == "deletion"]
        
        print(f"\n--- Size Field Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Deletions found: {len(deletions)}")
        
        # After fix, should have these fields
        for deletion in deletions:
            print(f"  Position: {deletion['ref_position']}, Size: {deletion['size']}, Sequence: {deletion['deleted_sequence']}")
            self.assertIn("size", deletion, "Deletion should have 'size' field")
            self.assertIn("deleted_sequence", deletion, 
                         "Deletion should have 'deleted_sequence' field")
            self.assertGreater(deletion["size"], 0, "Deletion size should be positive")
    
    def test_consecutive_insertions_grouped(self):
        """Consecutive insertions should also be grouped"""
        # Reference: ATCG
        # Consensus: ATGGCG (GG inserted - size 2)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("ATGGCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        insertions = [m for m in mutations if m["type"] == "insertion"]
        
        print(f"\n--- Grouped Insertion Test ---")
        print(f"Reference: {ref_seq}")
        print(f"Consensus: {cons_seq}")
        print(f"Insertions found: {len(insertions)}")
        for ins in insertions:
            print(f"  Position: {ins['ref_position']}, Size: {ins.get('size', 'N/A')}, Sequence: {ins.get('inserted_sequence', 'N/A')}")
        
        # Should find 1 grouped insertion
        self.assertEqual(len(insertions), 1, "Should group consecutive insertions")
        self.assertIn("size", insertions[0], "Insertion should have 'size' field")
        self.assertIn("inserted_sequence", insertions[0], "Insertion should have 'inserted_sequence' field")


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)

