# Mutation Caller Fixes - Complete Summary

This document summarizes all fixes applied to the mutation caller module.

## Fix #1: Reference/Consensus Nucleotide Swap (October 28, 2025)

### Problem
Reference and consensus nucleotides were reported incorrectly due to swapped function arguments.

### Root Cause
In `mutation_caller.py` line 161, arguments to `align_sequences()` were passed in wrong order.

### Fix
Changed:
```python
ref_aligned, cons_aligned, alignment_stats = align_sequences(reference_seq, cons_seq)
```
To:
```python
ref_aligned, cons_aligned, alignment_stats = align_sequences(cons_seq, reference_seq)
```

### Documentation
See `MUTATION_CALLER_FIX.md` for detailed information.

---

## Fix #2: Consecutive Deletion Grouping (October 28, 2025)

### Problem
Consecutive deletions and insertions were treated as separate mutation events instead of being grouped together.

### Example
- Reference: `ATCGATCG`
- Consensus: `ATCG` (4 bases deleted)
- **Before**: 4 separate deletion records
- **After**: 1 deletion record with `size=4`

### Fix
Rewrote `find_mutations()` function to:
1. Detect consecutive indels using while loop
2. Group them as single events
3. Add `size` and sequence information fields

### New Data Structure

**Deletions:**
```python
{
    "type": "deletion",
    "ref_position": 5,
    "ref_base": "GATC",           # Multiple bases if consecutive
    "cons_base": "-",
    "size": 4,                    # NEW: Size of deletion
    "deleted_sequence": "GATC"    # NEW: Full deleted sequence
}
```

**Insertions:**
```python
{
    "type": "insertion",
    "ref_position": 5,
    "ref_base": "-",
    "cons_base": "GG",
    "size": 2,                     # NEW: Size of insertion
    "inserted_sequence": "GG"      # NEW: Full inserted sequence
}
```

### Documentation
See `CONSECUTIVE_DELETION_GROUPING.md` for detailed information.

---

## Test Files

### Test Data Location
`/test-data/test_mutation_caller/`

- `reference.fa` - Reference for basic mutation tests
- `consensus_with_mutations.fa` - Various mutation types
- `consecutive_deletions_reference.fa` - Reference for deletion tests
- `consecutive_deletions_consensus.fa` - Consecutive deletion patterns

### Test Suites Location
`/tests/`

1. **`test_mutation_caller.py`** - Tests for Fix #1 (nucleotide swap)
2. **`test_consecutive_deletions.py`** - Tests for Fix #2 (deletion grouping)
3. **`verify_mutation_fix.py`** - Quick verification for Fix #1
4. **`verify_deletion_grouping.py`** - Quick verification for Fix #2

### Running All Tests

```bash
cd /Users/anitido/Documents/GitHub/belt/belt

# Run all mutation caller tests
python3 -m unittest discover -s tests -p "test_mutation*.py" -v
python3 -m unittest discover -s tests -p "test_consecutive*.py" -v

# Run quick verifications
python3 tests/verify_mutation_fix.py
python3 tests/verify_deletion_grouping.py
```

---

## Files Modified

### Core Code
- `scripts/mutation_caller.py`
  - Line 161: Fixed argument order for `align_sequences()`
  - Lines 81-152: Rewrote `find_mutations()` to group consecutive indels

### Test Data Created
- `test-data/test_mutation_caller/reference.fa`
- `test-data/test_mutation_caller/consensus_with_mutations.fa`
- `test-data/test_mutation_caller/consecutive_deletions_reference.fa`
- `test-data/test_mutation_caller/consecutive_deletions_consensus.fa`

### Test Files Created
- `tests/test_mutation_caller.py`
- `tests/test_consecutive_deletions.py`
- `tests/verify_mutation_fix.py`
- `tests/verify_deletion_grouping.py`

### Documentation Created
- `tests/MUTATION_CALLER_FIX.md`
- `tests/CONSECUTIVE_DELETION_GROUPING.md`
- `tests/README_MUTATION_TESTS.md`
- `tests/README_MUTATION_FIXES.md` (this file)

---

## Impact on Pipeline

These fixes affect:

1. **Step 9: mutation_calling** - Individual sample mutation calling
2. **Step 10: mutation_analysis** - Summary analysis across samples

### Expected Changes

1. **More accurate mutation reporting**: Reference and consensus bases now correctly identified
2. **Fewer mutation records**: Consecutive indels grouped as single events
3. **Enhanced data**: New `size` and sequence fields for indels
4. **Better statistics**: More accurate mutation rates and patterns

### Backward Compatibility

The changes are mostly backward compatible:
- Core fields (`type`, `ref_position`, `ref_base`, `cons_base`) still present
- New fields are additions, not replacements
- Single-base indels have `size=1`

However, downstream analyses may need updates to:
- Use new `size` field for statistics
- Handle multi-base `ref_base` and `cons_base` values for grouped indels
- Use `deleted_sequence`/`inserted_sequence` fields

---

## Summary

Both fixes have been implemented, tested, and documented. The mutation caller now:

✅ Correctly identifies reference vs consensus nucleotides  
✅ Groups consecutive deletions and insertions as single events  
✅ Provides size information for all indels  
✅ Maintains backward compatibility with core fields  

All changes are production-ready and integrated into the BELT pipeline.

---

## Date

All fixes completed: October 28, 2025

