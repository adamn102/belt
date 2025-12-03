#!/usr/bin/env python3
"""
Simple verification script to test that the mutation caller fix works correctly.
This creates inline test sequences and verifies the mutations are correctly identified.
"""

import sys
import os
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_caller import align_sequences, find_mutations
from Bio.Seq import Seq

def test_substitution():
    """Test A->T substitution is correctly identified"""
    print("Test 1: Substitution (Reference A -> Consensus T)")
    print("=" * 60)
    
    # Reference: ATCG
    # Consensus: TTCG (A->T at position 0)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("TTCG")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    print(f"Reference aligned: {ref_aligned}")
    print(f"Consensus aligned: {cons_aligned}")
    
    mutations = find_mutations(ref_aligned, cons_aligned)
    print(f"Mutations found: {len(mutations)}")
    
    if mutations:
        m = mutations[0]
        print(f"  Type: {m['type']}")
        print(f"  Position: {m['ref_position']}")
        print(f"  Reference base: {m['ref_base']}")
        print(f"  Consensus base: {m['cons_base']}")
        
        # Verify correct identification
        if m['ref_base'] == 'A' and m['cons_base'] == 'T':
            print("✓ PASS: Nucleotides correctly identified!")
            return True
        else:
            print("✗ FAIL: Nucleotides are still swapped!")
            return False
    else:
        print("✗ FAIL: No mutations found!")
        return False

def test_multiple_substitutions():
    """Test multiple substitutions"""
    print("\nTest 2: Multiple Substitutions")
    print("=" * 60)
    
    # Reference: AAAAA
    # Consensus: TTTTT (all A->T)
    ref_seq = Seq("AAAAA")
    cons_seq = Seq("TTTTT")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    
    print(f"Mutations found: {len(mutations)}")
    
    all_correct = True
    for i, m in enumerate(mutations):
        if m['ref_base'] != 'A' or m['cons_base'] != 'T':
            print(f"✗ Position {i}: Expected ref='A' cons='T', got ref='{m['ref_base']}' cons='{m['cons_base']}'")
            all_correct = False
    
    if all_correct and len(mutations) == 5:
        print("✓ PASS: All 5 substitutions correctly identified!")
        return True
    else:
        print("✗ FAIL: Some mutations incorrectly identified!")
        return False

def test_insertion():
    """Test insertion identification"""
    print("\nTest 3: Insertion")
    print("=" * 60)
    
    # Reference: ATCG
    # Consensus: ATGCG (G inserted)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("ATGCG")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    
    insertions = [m for m in mutations if m['type'] == 'insertion']
    print(f"Insertions found: {len(insertions)}")
    
    if insertions:
        m = insertions[0]
        print(f"  Reference base: {m['ref_base']}")
        print(f"  Consensus base: {m['cons_base']}")
        
        if m['ref_base'] == '-' and m['cons_base'] != '-':
            print("✓ PASS: Insertion correctly identified!")
            return True
        else:
            print("✗ FAIL: Insertion incorrectly identified!")
            return False
    else:
        print("✗ FAIL: No insertion found!")
        return False

def test_deletion():
    """Test deletion identification"""
    print("\nTest 4: Deletion")
    print("=" * 60)
    
    # Reference: ATCG
    # Consensus: ACG (T deleted)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("ACG")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    
    deletions = [m for m in mutations if m['type'] == 'deletion']
    print(f"Deletions found: {len(deletions)}")
    
    if deletions:
        m = deletions[0]
        print(f"  Reference base: {m['ref_base']}")
        print(f"  Consensus base: {m['cons_base']}")
        
        if m['ref_base'] != '-' and m['cons_base'] == '-':
            print("✓ PASS: Deletion correctly identified!")
            return True
        else:
            print("✗ FAIL: Deletion incorrectly identified!")
            return False
    else:
        print("✗ FAIL: No deletion found!")
        return False

def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("MUTATION CALLER VERIFICATION TESTS")
    print("=" * 60 + "\n")
    
    results = []
    results.append(("Substitution", test_substitution()))
    results.append(("Multiple Substitutions", test_multiple_substitutions()))
    results.append(("Insertion", test_insertion()))
    results.append(("Deletion", test_deletion()))
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    total_passed = sum(1 for _, passed in results if passed)
    print(f"\nTotal: {total_passed}/{len(results)} tests passed")
    
    if total_passed == len(results):
        print("\n✓ All tests passed! The mutation caller fix is working correctly.")
        return 0
    else:
        print("\n✗ Some tests failed. The mutation caller may still have issues.")
        return 1

if __name__ == '__main__':
    sys.exit(main())

