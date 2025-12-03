# Mutation Caller Fix Documentation

## Issue Identified

The mutation caller had a bug where reference and consensus nucleotides were swapped in the output.

### Root Cause

In `mutation_caller.py` line 161, the function call to `align_sequences()` had the arguments in the wrong order:

**BEFORE (BUGGY):**
```python
ref_aligned, cons_aligned, alignment_stats = align_sequences(reference_seq, cons_seq)
```

However, the `align_sequences()` function signature expects:
```python
def align_sequences(consensus_seq: Seq, reference_seq: Seq) -> Tuple[str, str, Dict]:
```

This caused the arguments to be swapped internally, resulting in:
- `ref_aligned` actually containing the consensus sequence
- `cons_aligned` actually containing the reference sequence

### The Fix

Changed line 161 to pass arguments in the correct order:

**AFTER (FIXED):**
```python
ref_aligned, cons_aligned, alignment_stats = align_sequences(cons_seq, reference_seq)
```

Now the arguments match the function signature, and:
- `ref_aligned` correctly contains the reference sequence
- `cons_aligned` correctly contains the consensus sequence

### Impact

This fix ensures that mutation calls correctly report:
- `ref_base`: The nucleotide in the reference sequence
- `cons_base`: The nucleotide in the consensus sequence

For example, if the reference has `A` at position 10 and the consensus has `T`:
- **Before fix**: Would report `ref_base='T'` and `cons_base='A'` (SWAPPED)
- **After fix**: Correctly reports `ref_base='A'` and `cons_base='T'`

## Testing

### Test Data Created

Located in `/test-data/test_mutation_caller/`:
- `reference.fa`: Reference sequence (ATCGATCGATCGATCGATCG)
- `consensus_with_mutations.fa`: Consensus sequences with known mutations

### Test Files Created

1. **`tests/test_mutation_caller.py`**: Comprehensive unittest suite
   - Tests substitutions, insertions, and deletions
   - Verifies correct base assignment
   - Tests alignment order

2. **`tests/verify_mutation_fix.py`**: Simple verification script
   - Quick inline tests for common mutation types
   - Easy to run and understand
   - Provides clear pass/fail output

### Running Tests

```bash
# Run comprehensive test suite
cd /Users/anitido/Documents/GitHub/belt/belt
python3 -m unittest tests/test_mutation_caller.py -v

# Run simple verification
python3 tests/verify_mutation_fix.py
```

## Files Modified

- `scripts/mutation_caller.py`: Fixed argument order on line 161

## Files Created

- `test-data/test_mutation_caller/reference.fa`
- `test-data/test_mutation_caller/consensus_with_mutations.fa`
- `tests/test_mutation_caller.py`
- `tests/verify_mutation_fix.py`
- `tests/MUTATION_CALLER_FIX.md` (this file)

## Date

Fixed: October 28, 2025

