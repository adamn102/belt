#!/usr/bin/env python3
"""
Simple verification script to test consecutive deletion grouping.
"""

import sys
import os
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from mutation_caller import align_sequences, find_mutations
from Bio.Seq import Seq

def test_single_deletion():
    """Test single deletion"""
    print("=" * 70)
    print("TEST 1: Single Deletion")
    print("=" * 70)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("ATC")
    
    print(f"Reference: {ref_seq}")
    print(f"Consensus: {cons_seq}")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    deletions = [m for m in mutations if m["type"] == "deletion"]
    
    print(f"\nDeletions found: {len(deletions)}")
    for d in deletions:
        print(f"  Position: {d['ref_position']}, Size: {d.get('size', 1)}, Sequence: {d.get('deleted_sequence', d['ref_base'])}")
    
    if len(deletions) == 1 and deletions[0].get('size', 1) == 1:
        print("✓ PASS: Single deletion correctly identified\n")
        return True
    else:
        print("✗ FAIL: Single deletion not correctly identified\n")
        return False

def test_double_deletion():
    """Test double consecutive deletion"""
    print("=" * 70)
    print("TEST 2: Double Consecutive Deletion")
    print("=" * 70)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("AT")
    
    print(f"Reference: {ref_seq}")
    print(f"Consensus: {cons_seq}")
    print(f"Expected: 1 deletion of size 2 (CG)")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    deletions = [m for m in mutations if m["type"] == "deletion"]
    
    print(f"\nDeletions found: {len(deletions)}")
    for d in deletions:
        print(f"  Position: {d['ref_position']}, Size: {d.get('size', 1)}, Sequence: {d.get('deleted_sequence', d['ref_base'])}")
    
    if len(deletions) == 1 and deletions[0].get('size') == 2 and deletions[0].get('deleted_sequence') == "CG":
        print("✓ PASS: Double deletion correctly grouped\n")
        return True
    else:
        print("✗ FAIL: Double deletion not correctly grouped\n")
        return False

def test_large_deletion():
    """Test large consecutive deletion"""
    print("=" * 70)
    print("TEST 3: Large Consecutive Deletion (8 bases)")
    print("=" * 70)
    ref_seq = Seq("ATCGATCGATCGATCG")
    cons_seq = Seq("ATCGATCG")
    
    print(f"Reference: {ref_seq} (length: {len(ref_seq)})")
    print(f"Consensus: {cons_seq} (length: {len(cons_seq)})")
    print(f"Expected: 1 deletion of size 8")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    deletions = [m for m in mutations if m["type"] == "deletion"]
    
    print(f"\nDeletions found: {len(deletions)}")
    for d in deletions:
        seq = d.get('deleted_sequence', d['ref_base'])
        print(f"  Position: {d['ref_position']}, Size: {d.get('size', 1)}, Sequence: {seq}")
    
    if len(deletions) == 1 and deletions[0].get('size') == 8:
        print("✓ PASS: Large deletion correctly grouped\n")
        return True
    else:
        print("✗ FAIL: Large deletion not correctly grouped\n")
        return False

def test_insertion_grouping():
    """Test insertion grouping"""
    print("=" * 70)
    print("TEST 4: Consecutive Insertion Grouping")
    print("=" * 70)
    ref_seq = Seq("ATCG")
    cons_seq = Seq("ATGGCG")
    
    print(f"Reference: {ref_seq}")
    print(f"Consensus: {cons_seq}")
    print(f"Expected: 1 insertion of size 2 (GG)")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    insertions = [m for m in mutations if m["type"] == "insertion"]
    
    print(f"\nInsertions found: {len(insertions)}")
    for ins in insertions:
        print(f"  Position: {ins['ref_position']}, Size: {ins.get('size', 1)}, Sequence: {ins.get('inserted_sequence', ins['cons_base'])}")
    
    if len(insertions) >= 1 and any(ins.get('size', 0) >= 2 for ins in insertions):
        print("✓ PASS: Insertion grouping working\n")
        return True
    else:
        print("✗ FAIL: Insertion grouping not working\n")
        return False

def test_non_consecutive_deletions():
    """Test non-consecutive deletions remain separate"""
    print("=" * 70)
    print("TEST 5: Non-Consecutive Deletions (Should Remain Separate)")
    print("=" * 70)
    ref_seq = Seq("ATCGATCG")
    cons_seq = Seq("TCGTCG")
    
    print(f"Reference: {ref_seq}")
    print(f"Consensus: {cons_seq}")
    print(f"Expected: 2 separate deletions")
    
    ref_aligned, cons_aligned, stats = align_sequences(cons_seq, ref_seq)
    mutations = find_mutations(ref_aligned, cons_aligned)
    deletions = [m for m in mutations if m["type"] == "deletion"]
    
    print(f"\nDeletions found: {len(deletions)}")
    for d in deletions:
        print(f"  Position: {d['ref_position']}, Size: {d.get('size', 1)}, Sequence: {d.get('deleted_sequence', d['ref_base'])}")
    
    if len(deletions) == 2:
        print("✓ PASS: Non-consecutive deletions correctly kept separate\n")
        return True
    else:
        print("✗ FAIL: Non-consecutive deletions incorrectly grouped\n")
        return False

def main():
    """Run all tests"""
    print("\n" + "=" * 70)
    print("CONSECUTIVE DELETION GROUPING VERIFICATION")
    print("=" * 70 + "\n")
    
    results = []
    results.append(test_single_deletion())
    results.append(test_double_deletion())
    results.append(test_large_deletion())
    results.append(test_insertion_grouping())
    results.append(test_non_consecutive_deletions())
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    passed = sum(results)
    total = len(results)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("\n✓ All tests passed! Consecutive deletion grouping is working correctly.")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed.")
        return 1

if __name__ == '__main__':
    sys.exit(main())

