# Mutation Caller Tests

## Overview

This directory contains tests for the mutation caller functionality in the BELT pipeline.

## Test Files

### Unit Tests

**`test_mutation_caller.py`**
- Comprehensive unittest suite for mutation calling
- Tests substitutions, insertions, and deletions
- Verifies correct identification of reference vs consensus bases
- Tests alignment order and sequence processing

### Verification Script

**`verify_mutation_fix.py`**
- Simple standalone verification script
- Tests basic mutation types with clear examples
- Provides human-readable output
- Quick way to verify the fix is working

## Test Data

Located in `../test-data/test_mutation_caller/`:

- `reference.fa` - Reference sequence for testing
- `consensus_with_mutations.fa` - Consensus sequences with known mutations

## Running Tests

### Run Full Test Suite

```bash
cd /Users/anitido/Documents/GitHub/belt/belt
python3 -m unittest tests/test_mutation_caller.py -v
```

### Run Quick Verification

```bash
cd /Users/anitido/Documents/GitHub/belt/belt
python3 tests/verify_mutation_fix.py
```

### Run Specific Test

```bash
python3 -m unittest tests/test_mutation_caller.py::TestMutationCaller::test_substitution_nucleotides -v
```

## What the Tests Verify

1. **Correct Base Assignment**: Reference bases are correctly identified as reference, consensus bases as consensus
2. **Substitution Detection**: A→T, C→G, etc. mutations are correctly detected
3. **Insertion Detection**: Bases present in consensus but not in reference  
4. **Deletion Detection**: Bases present in reference but not in consensus
5. **Position Tracking**: Mutations are correctly mapped to reference positions
6. **Alignment Order**: Sequences are aligned in the correct order

## Expected Test Results

All tests should pass after the fix. The key verification is:

- **Reference has `A`, Consensus has `T`**
  - `ref_base` should be `'A'`
  - `cons_base` should be `'T'`
  - If swapped, the bug still exists

## Bug Fix Summary

**Issue**: Reference and consensus nucleotides were swapped in mutation output

**Root Cause**: Function arguments passed in wrong order to `align_sequences()`

**Fix**: Corrected argument order in `mutation_caller.py` line 161

See `MUTATION_CALLER_FIX.md` for detailed documentation.

