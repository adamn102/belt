#!/usr/bin/env python3
"""
Test suite for mutation_caller.py
Tests to verify correct identification of reference vs consensus nucleotides
and proper handling of N bases (UMI sequences) in the reference
"""

import unittest
import sys
import os
import json
import tempfile
import shutil
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_caller import (
    read_fasta, align_sequences, _find_mutations_skip_ref_N, 
    _callable_sites, run, call_mutations_for_sample
)
from Bio.Seq import Seq


class TestMutationCaller(unittest.TestCase):
    """Test mutation caller functions"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_data_dir = Path(__file__).parent.parent / "test-data" / "test_mutation_caller"
        self.temp_dir = Path(tempfile.mkdtemp())
        
    def tearDown(self):
        """Clean up temporary files"""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_substitution_nucleotides(self):
        """Test that substitution mutations correctly identify ref and cons bases"""
        # Reference: ATCG
        # Consensus: TTCG (A->T substitution at position 0)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("TTCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Should find one substitution
        self.assertEqual(len(mutations), 1)
        mutation = mutations[0]
        
        # Verify mutation type
        self.assertEqual(mutation["type"], "substitution")
        
        # CRITICAL TEST: Reference should be A, consensus should be T
        # If these are swapped, the bug exists
        self.assertEqual(mutation["ref_base"], "A", 
                        f"Expected ref_base='A' but got '{mutation['ref_base']}'. Reference and consensus bases may be swapped!")
        self.assertEqual(mutation["cons_base"], "T",
                        f"Expected cons_base='T' but got '{mutation['cons_base']}'. Reference and consensus bases may be swapped!")
        self.assertEqual(mutation["ref_position"], 0)
    
    def test_multiple_substitutions(self):
        """Test multiple substitutions are correctly identified"""
        # Reference: ATCGATCG
        # Consensus: TTCGTTCG (A->T at pos 0, A->T at pos 4)
        ref_seq = Seq("ATCGATCG")
        cons_seq = Seq("TTCGTTCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Should find two substitutions
        self.assertEqual(len(mutations), 2)
        
        # First mutation: A->T at position 0
        self.assertEqual(mutations[0]["type"], "substitution")
        self.assertEqual(mutations[0]["ref_base"], "A")
        self.assertEqual(mutations[0]["cons_base"], "T")
        self.assertEqual(mutations[0]["ref_position"], 0)
        
        # Second mutation: A->T at position 4
        self.assertEqual(mutations[1]["type"], "substitution")
        self.assertEqual(mutations[1]["ref_base"], "A")
        self.assertEqual(mutations[1]["cons_base"], "T")
        self.assertEqual(mutations[1]["ref_position"], 4)
    
    def test_insertion_identification(self):
        """Test that insertions are correctly identified"""
        # Reference: ATCG
        # Consensus: ATGCG (G inserted after position 2)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("ATGCG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Should find one insertion
        insertions = [m for m in mutations if m["type"] == "insertion"]
        self.assertEqual(len(insertions), 1)
        
        insertion = insertions[0]
        self.assertEqual(insertion["ref_base"], "-")
        self.assertNotEqual(insertion["cons_base"], "-")
    
    def test_deletion_identification(self):
        """Test that deletions are correctly identified"""
        # Reference: ATCG
        # Consensus: ACG (T deleted at position 1)
        ref_seq = Seq("ATCG")
        cons_seq = Seq("ACG")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Should find one deletion
        deletions = [m for m in mutations if m["type"] == "deletion"]
        self.assertEqual(len(deletions), 1)
        
        deletion = deletions[0]
        self.assertNotEqual(deletion["ref_base"], "-")
        self.assertEqual(deletion["cons_base"], "-")
    
    def test_known_mutation_pattern(self):
        """Test with a known mutation pattern to verify correct base assignment"""
        # Reference: AAAAA
        # Consensus: TTTTT (all A->T)
        ref_seq = Seq("AAAAA")
        cons_seq = Seq("TTTTT")
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Should find 5 substitutions
        self.assertEqual(len(mutations), 5)
        
        # All should be A->T (reference A, consensus T)
        for i, mutation in enumerate(mutations):
            self.assertEqual(mutation["type"], "substitution")
            self.assertEqual(mutation["ref_base"], "A", 
                           f"Position {i}: Expected ref_base='A' but got '{mutation['ref_base']}'")
            self.assertEqual(mutation["cons_base"], "T",
                           f"Position {i}: Expected cons_base='T' but got '{mutation['cons_base']}'")
    
    def test_mutation_caller_with_test_data(self):
        """Test mutation_caller function with test data files"""
        if not self.test_data_dir.exists():
            self.skipTest("Test data directory not found")
        
        ref_file = self.test_data_dir / "reference.fa"
        cons_file = self.test_data_dir / "consensus_with_mutations.fa"
        
        if not ref_file.exists() or not cons_file.exists():
            self.skipTest("Test data files not found")
        
        # Run mutation caller
        result = mutation_caller(
            consensus_file=cons_file,
            reference_file=ref_file,
            sample_name="test_sample",
            output_dir=self.temp_dir
        )
        
        # Verify result structure
        self.assertIn("sample_name", result)
        self.assertIn("total_mutations", result)
        self.assertIn("mutation_types", result)
        self.assertIn("per_read_analysis", result)
        
        # Check output files were created
        sample_dir = self.temp_dir / "samples" / "test_sample"
        self.assertTrue(sample_dir.exists())
        self.assertTrue((sample_dir / "test_sample_mutations.json").exists())
        self.assertTrue((sample_dir / "test_sample_per_read.csv").exists())


class TestAlignmentOrder(unittest.TestCase):
    """Test alignment order to ensure sequences are correctly assigned"""
    
    def test_alignment_returns_correct_order(self):
        """Verify that align_sequences returns (reference, consensus) in correct order"""
        # Create clearly different sequences
        ref_seq = Seq("AAAAAAAAAA")  # All A's
        cons_seq = Seq("TTTTTTTTTT")  # All T's
        
        ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
        
        # Reference should be all A's (or gaps)
        # Consensus should be all T's (or gaps)
        ref_bases = set(ref_aligned.replace('-', ''))
        cons_bases = set(cons_aligned.replace('-', ''))
        
        self.assertEqual(ref_bases, {'A'}, "Reference sequence should only contain A's")
        self.assertEqual(cons_bases, {'T'}, "Consensus sequence should only contain T's")


class TestNSkippingBehavior(unittest.TestCase):
    """Test that mutations are not called at N positions (UMI sequences) in reference"""
    
    def test_substitution_at_N_position_is_skipped(self):
        """Test that substitutions at N positions in reference are skipped"""
        # Reference: ANCG (N at position 1)
        # Consensus: ATCG (T at position 1, but ref is N so should be skipped)
        ref_raw = "ANCG"
        cons_seq = "ATCG"  # T instead of N at position 1
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should skip the substitution at position 1 (N->T) because ref is N
        # Check that no mutation is reported at position 1
        n_mutations = [m for m in mutations if m.ref_position == 1]
        self.assertEqual(len(n_mutations), 0, 
                        "Should not call mutations at N positions in reference")
    
    def test_substitution_at_ACGT_position_is_called(self):
        """Test that substitutions at ACGT positions are still called"""
        # Reference: ATCG
        # Consensus: AGCG (T->G at position 1)
        ref_raw = "ATCG"
        cons_seq = "AGCG"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should find one substitution at position 1
        substitutions = [m for m in mutations if m.type == "substitution" and m.ref_position == 1]
        self.assertEqual(len(substitutions), 1, "Should call mutations at ACGT positions")
        self.assertEqual(substitutions[0].ref_base, "T")
        self.assertEqual(substitutions[0].cons_base, "G")
    
    def test_umi_sequence_skipped(self):
        """Test that UMI sequence (NNNNNNNNNNNN) in reference is completely skipped"""
        # Reference: ATCG + NNNNNNNNNNNN + GCTA (UMI in middle)
        # Consensus: ATCG + AAAAAAAAAAAA + GCTA (different sequence in UMI region)
        ref_raw = "ATCG" + "N" * 12 + "GCTA"
        cons_seq = "ATCG" + "A" * 12 + "GCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should find NO mutations because all differences are in the N region
        self.assertEqual(len(mutations), 0, 
                        "Should not call mutations in UMI (N) regions")
    
    def test_deletion_spanning_N_is_skipped(self):
        """Test that deletions spanning N positions are skipped"""
        # Reference: ATCGNGCTA (N at position 4)
        # Consensus: ATCGCTA (deletion of NG, which spans the N)
        ref_raw = "ATCGNGCTA"
        cons_seq = "ATCGCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should skip deletion because it spans an N
        deletions = [m for m in mutations if m.type == "deletion"]
        self.assertEqual(len(deletions), 0, 
                        "Should skip deletions that span N positions")
    
    def test_deletion_not_spanning_N_is_called(self):
        """Test that deletions not spanning N are still called"""
        # Reference: ATCGGCCTA
        # Consensus: ATCGCTA (deletion of GC)
        ref_raw = "ATCGGCCTA"
        cons_seq = "ATCGCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should find deletion
        deletions = [m for m in mutations if m.type == "deletion"]
        self.assertGreater(len(deletions), 0, 
                          "Should call deletions that don't span N")
    
    def test_insertion_at_N_is_skipped(self):
        """Test that insertions anchored at N positions are skipped"""
        # Reference: ATCGNGCTA (N at position 4)
        # Consensus: ATCGNAAGCTA (AA inserted after N)
        ref_raw = "ATCGNGCTA"
        cons_seq = "ATCGNAAGCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw, skip_insertions_at_refN=True)
        
        # Should skip insertion because it's anchored at N
        insertions = [m for m in mutations if m.type == "insertion"]
        self.assertEqual(len(insertions), 0, 
                        "Should skip insertions anchored at N positions")
    
    def test_umi_insertions_skipped(self):
        """Test that insertions in UMI region (NNNNNNNNNNNN) are skipped"""
        # Reference: ATCG + NNNNNNNNNNNN + GCTA (UMI in middle)
        # Consensus: ATCG + ATCGATCGATCG + GCTA (actual bases in UMI region)
        # If alignment creates insertions, they should be skipped
        ref_raw = "ATCG" + "N" * 12 + "GCTA"
        cons_seq = "ATCG" + "ATCG" * 3 + "GCTA"  # 12 bases
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw, skip_insertions_at_refN=True)
        
        # Should have NO insertions because they're in the UMI (N) region
        insertions = [m for m in mutations if m.type == "insertion"]
        self.assertEqual(len(insertions), 0,
                        "Should skip insertions in UMI (N) regions")
        
        # Also verify no mutations at all in the UMI region (positions 4-15)
        umi_region_mutations = [m for m in mutations if 4 <= m.ref_position <= 15]
        self.assertEqual(len(umi_region_mutations), 0,
                        "Should not call any mutations in UMI region (positions 4-15)")
    
    def test_insertion_at_ACGT_is_called(self):
        """Test that insertions at ACGT positions are still called"""
        # Reference: ATCGGCCTA
        # Consensus: ATCGGAACTA (AA inserted after G)
        ref_raw = "ATCGGCCTA"
        cons_seq = "ATCGGAACTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw, skip_insertions_at_refN=True)
        
        # Should find insertion
        insertions = [m for m in mutations if m.type == "insertion"]
        self.assertGreater(len(insertions), 0, 
                          "Should call insertions at ACGT positions")
    
    def test_callable_sites_excludes_N(self):
        """Test that callable sites calculation excludes N positions"""
        # Reference: ATCGNGCTA (N at position 4)
        # Consensus: ATCGAGCTA (matches except N position)
        ref_raw = "ATCGNGCTA"
        cons_seq = "ATCGAGCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        callable_count = _callable_sites(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should count: A(0), T(1), C(2), G(3), skip N(4), G(5), C(6), T(7), A(8)
        # Total: 8 positions (excluding position 4 which is N)
        expected = 8
        self.assertEqual(callable_count, expected, 
                        f"Callable sites should exclude N positions. Expected {expected}, got {callable_count}")
    
    def test_callable_sites_with_umi_sequence(self):
        """Test callable sites with UMI sequence (12 N's)"""
        # Reference: ATCG + NNNNNNNNNNNN + GCTA
        # Consensus: ATCG + AAAAAAAAAAAA + GCTA
        ref_raw = "ATCG" + "N" * 12 + "GCTA"
        cons_seq = "ATCG" + "A" * 12 + "GCTA"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        callable_count = _callable_sites(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should count: A(0), T(1), C(2), G(3), skip 12 N's, G, C, T, A
        # Total: 4 + 4 = 8 positions (excluding the 12 N positions)
        expected = 8
        self.assertEqual(callable_count, expected,
                        f"Callable sites should exclude UMI (N) region. Expected {expected}, got {callable_count}")
    
    def test_end_to_end_with_umi_in_reference(self):
        """End-to-end test with UMI sequence in reference"""
        import tempfile
        
        # Create temporary FASTA files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as ref_file:
            ref_file.write(">ref1\n")
            ref_file.write("ATCG" + "N" * 12 + "GCTA\n")  # UMI in middle
            ref_path = ref_file.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as cons_file:
            cons_file.write(">UMI1\n")
            cons_file.write("ATCG" + "A" * 12 + "GCTA\n")  # Different sequence in UMI region
            cons_path = cons_file.name
        
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = Path(tmpdir) / "mutations.csv"
            
            try:
                # Run mutation caller
                run(
                    reference_fa=ref_path,
                    consensus_fa=cons_path,
                    out_csv=out_csv,
                    per_mutation_tsv=None,
                    ref_contig=None,
                    skip_insertions_at_refN=True
                )
                
                # Read results
                import csv
                with open(out_csv, 'r') as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                
                self.assertEqual(len(rows), 1, "Should have one consensus sequence")
                row = rows[0]
                
                # Should have NO mutations because all differences are in N region
                self.assertEqual(int(row['mutation_count']), 0,
                               "Should not call mutations in UMI (N) regions")
                
                # Callable bases should exclude the 12 N positions
                # Total length is 4 + 12 + 4 = 20, but only 8 are callable (4 + 4)
                callable = int(row['callable_bases'])
                self.assertEqual(callable, 8,
                               f"Callable bases should exclude N positions. Expected 8, got {callable}")
                
            finally:
                # Cleanup
                os.unlink(ref_path)
                os.unlink(cons_path)
    
    def test_mixed_N_and_ACGT_positions(self):
        """Test with mixed N and ACGT positions - only ACGT mutations should be called"""
        # Reference: ANCGNTCGN (N at positions 1, 4, 8)
        # Consensus: ATCGATCGN (substitutions at positions 1, 4, 8)
        ref_raw = "ANCGNTCGN"
        cons_seq = "ATCGATCGN"
        
        # Align sequences
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should have NO mutations because:
        # - Position 1: N->T (skip, ref is N)
        # - Position 2: C->C (match)
        # - Position 3: G->G (match)
        # - Position 4: N->A (skip, ref is N)
        # - Position 5: T->T (match)
        # - Position 6: C->C (match)
        # - Position 7: G->G (match)
        # - Position 8: N->N (match, but ref is N so wouldn't count anyway)
        self.assertEqual(len(mutations), 0,
                       "Should not call mutations at N positions, even if consensus differs")
    
    def test_real_umi_scenario(self):
        """Test realistic scenario: reference with UMI barcode region"""
        # Reference: 50bp + NNNNNNNNNNNN (UMI) + 50bp
        ref_start = "A" * 25 + "C" * 25
        umi = "N" * 12
        ref_end = "G" * 25 + "T" * 25
        ref_raw = ref_start + umi + ref_end
        
        # Consensus: same but with mutations in UMI region (should be ignored)
        # and one real mutation in the callable region
        cons_seq = ref_start + "A" * 12 + ref_end[:-1] + "A"  # T->A mutation at last position
        
        # Align
        aln = align_sequences(ref_raw, cons_seq)
        
        # Find mutations
        mutations = _find_mutations_skip_ref_N(aln.ref_aln, aln.qry_aln, ref_raw)
        
        # Should find 1 mutation (T->A at the end), NOT the differences in UMI region
        substitutions = [m for m in mutations if m.type == "substitution"]
        # The mutation should be at the last position (position 111, 0-based: 25+25+12+25+24 = 111)
        # Actually, let's just check we have exactly 1 mutation and it's not in the UMI region
        self.assertEqual(len(substitutions), 1,
                        "Should find exactly one mutation (not in UMI region)")
        
        # Verify the mutation is NOT in the UMI region (positions 50-61)
        mut_pos = substitutions[0].ref_position
        self.assertNotIn(mut_pos, range(50, 62),
                        f"Mutation at position {mut_pos} should not be in UMI region (50-61)")


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)

