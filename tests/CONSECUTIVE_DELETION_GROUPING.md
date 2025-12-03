# Consecutive Deletion Grouping Fix

## Summary

Fixed the mutation caller to group consecutive deletions and insertions as single events with size information, rather than treating each base as a separate mutation.

## Problem

Previously, when multiple consecutive bases were deleted or inserted, the mutation caller would create separate mutation records for each base. For example:

**Before Fix:**
- Reference: `ATCGATCG`
- Consensus: `ATCG` (4 bases deleted)
- Result: 4 separate deletion records

This made it difficult to:
- Understand the true mutation pattern
- Calculate accurate mutation statistics
- Identify larger structural changes

## Solution

Modified the `find_mutations()` function in `mutation_caller.py` to:

1. **Detect consecutive indels**: Use a while loop to collect all consecutive insertions or deletions
2. **Group as single event**: Create one mutation record for the entire indel
3. **Add size information**: Include `size` and `deleted_sequence`/`inserted_sequence` fields

**After Fix:**
- Reference: `ATCGATCG`
- Consensus: `ATCG` (4 bases deleted)
- Result: 1 deletion record with `size=4` and `deleted_sequence="GATC"`

## Changes Made

### Code Changes

**File**: `scripts/mutation_caller.py`

**Function**: `find_mutations(ref_aligned: str, cons_aligned: str) -> List[Dict]`

Key improvements:
- Changed from `for` loop to `while` loop to allow skipping ahead
- Added consecutive detection logic for deletions and insertions
- Enhanced mutation data structure with new fields

### New Mutation Data Structure

**Deletions:**
```python
{
    "type": "deletion",
    "ref_position": 5,           # Starting position in reference
    "ref_base": "GATC",           # Deleted bases (or single base)
    "cons_base": "-",
    "size": 4,                    # Number of consecutive deleted bases
    "deleted_sequence": "GATC"    # Full deleted sequence
}
```

**Insertions:**
```python
{
    "type": "insertion",
    "ref_position": 5,
    "ref_base": "-",
    "cons_base": "GG",             # Inserted bases (or single base)
    "size": 2,                     # Number of consecutive inserted bases
    "inserted_sequence": "GG"      # Full inserted sequence
}
```

**Substitutions** (unchanged):
```python
{
    "type": "substitution",
    "ref_position": 5,
    "ref_base": "A",
    "cons_base": "T"
}
```

## Test Files Created

### Test Data

Located in `/test-data/test_mutation_caller/`:

1. **`consecutive_deletions_reference.fa`**
   - Reference sequence for testing consecutive deletions

2. **`consecutive_deletions_consensus.fa`**
   - Multiple consensus sequences with various deletion patterns
   - Single, double, triple, five-base, and middle deletions

### Test Suites

Located in `/tests/`:

1. **`test_consecutive_deletions.py`**
   - Comprehensive unittest suite
   - Tests single deletions, double deletions, large deletions
   - Tests non-consecutive deletions (should remain separate)
   - Tests with real FASTA files
   - Tests insertion grouping

2. **`verify_deletion_grouping.py`**
   - Simple verification script
   - Clear pass/fail output
   - Tests common scenarios
   - Easy to run and understand

## Running Tests

### Run Comprehensive Test Suite

```bash
cd /Users/anitido/Documents/GitHub/belt/belt
python3 -m unittest tests/test_consecutive_deletions.py -v
```

### Run Quick Verification

```bash
cd /Users/anitido/Documents/GitHub/belt/belt
python3 tests/verify_deletion_grouping.py
```

## Examples

### Example 1: Double Deletion

```python
# Reference: ATCG
# Consensus: AT

# OLD BEHAVIOR:
mutations = [
    {"type": "deletion", "ref_position": 2, "ref_base": "C"},
    {"type": "deletion", "ref_position": 3, "ref_base": "G"}
]

# NEW BEHAVIOR:
mutations = [
    {
        "type": "deletion",
        "ref_position": 2,
        "ref_base": "CG",
        "size": 2,
        "deleted_sequence": "CG"
    }
]
```

### Example 2: Large Deletion

```python
# Reference: ATCGATCGATCGATCG (16 bases)
# Consensus: ATCGATCG (8 bases)

# OLD BEHAVIOR: 8 separate deletion records

# NEW BEHAVIOR:
mutations = [
    {
        "type": "deletion",
        "ref_position": 8,
        "ref_base": "ATCGATCG",
        "size": 8,
        "deleted_sequence": "ATCGATCG"
    }
]
```

### Example 3: Non-Consecutive Deletions (Remain Separate)

```python
# Reference: ATCGATCG
# Consensus: TCGTCG (A deleted at pos 0, A deleted at pos 4)

# BEHAVIOR (correct): 2 separate deletions
mutations = [
    {
        "type": "deletion",
        "ref_position": 0,
        "ref_base": "A",
        "size": 1,
        "deleted_sequence": "A"
    },
    {
        "type": "deletion",
        "ref_position": 4,
        "ref_base": "A",
        "size": 1,
        "deleted_sequence": "A"
    }
]
```

## Backward Compatibility

### Potential Impact

The mutation data structure has been enhanced with new fields. Existing code that processes mutations should still work because:

1. **Core fields preserved**: `type`, `ref_position`, `ref_base`, `cons_base` still exist
2. **New fields added**: `size` and `*_sequence` fields are additions, not replacements
3. **Single mutations**: For single-base indels, `size=1` and sequence fields contain single bases

### Migration Notes

If you have downstream code that processes mutations:

1. **Check for size field**: Use `mutation.get('size', 1)` to handle both old and new formats
2. **Use sequence fields**: Prefer `deleted_sequence`/`inserted_sequence` over `ref_base`/`cons_base` for indels
3. **Update statistics**: Mutation counts may change (fewer mutation records, but same total bases affected)

## Benefits

1. **More accurate representation**: Matches biological reality of indel events
2. **Better statistics**: Easier to calculate true mutation rates
3. **Clearer patterns**: Easier to identify hotspots and structural variants
4. **Size information**: Can now distinguish small vs. large indels
5. **Improved analysis**: Downstream tools can better understand mutation patterns

## Date

Implemented: October 28, 2025

