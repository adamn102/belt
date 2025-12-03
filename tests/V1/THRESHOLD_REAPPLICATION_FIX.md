# Critical Fix: Threshold Re-application After UMI2 Filtering

## Issue Identified

The original implementation applied the minimum read threshold only **before** UMI2 filtering. This created a potential problem:

1. A UMI1 group with 10 reads passes the initial threshold (e.g., `threshold=10`)
2. UMI2 filtering identifies majority UMI2 and removes minority reads
3. Only 6 reads remain (those with the majority UMI2)
4. **Problem**: Group with only 6 reads proceeds to consensus calling despite being below threshold

This could result in low-coverage consensus sequences with poor quality.

## Solution

The threshold is now applied **twice**:

### Before Fix:
```python
if majority_fraction >= umi2_majority_threshold:
    passing_umi1_values.append(umi1_value)
    umi1_to_majority_umi2[umi1_value] = majority_umi2
    reads_after_umi2_filter += majority_count
```

### After Fix:
```python
if majority_fraction >= umi2_majority_threshold and majority_count >= threshold:
    passing_umi1_values.append(umi1_value)
    umi1_to_majority_umi2[umi1_value] = majority_umi2
    reads_after_umi2_filter += majority_count
```

## Code Change

**File**: `belt/scripts/umi_group.py`
**Line**: 123

```python
# Check if majority UMI2 meets threshold AND has enough reads after filtering
if majority_fraction >= umi2_majority_threshold and majority_count >= threshold:
```

## Impact

### Scenario 1: Group Passes Both Checks ✓
- **Input**: 15 total reads, 12 with majority UMI2 (80%)
- **Threshold**: 10 reads
- **Result**: 
  - Majority fraction: 80% >= 50% ✓
  - Majority count: 12 >= 10 ✓
  - **PASS**: Group proceeds with 12 reads

### Scenario 2: Group Fails Second Check ✗
- **Input**: 10 total reads, 6 with majority UMI2 (60%)
- **Threshold**: 10 reads
- **Result**: 
  - Majority fraction: 60% >= 50% ✓
  - Majority count: 6 >= 10 ✗
  - **FAIL**: Group excluded (insufficient coverage after filtering)

### Scenario 3: Group Passes Both with Lower Threshold ✓
- **Input**: 10 total reads, 6 with majority UMI2 (60%)
- **Threshold**: 5 reads
- **Result**: 
  - Majority fraction: 60% >= 50% ✓
  - Majority count: 6 >= 5 ✓
  - **PASS**: Group proceeds with 6 reads

## Test Coverage

Added comprehensive test: `test_threshold_reapplied_after_umi2_filter()`

### Test Case 1: Group Below Threshold After Filtering
```python
# UMI1-A: 10 reads, 6 majority (60%), threshold=10
# Expected: EXCLUDED (6 < 10)
assert "AAAA" not in passing_groups
```

### Test Case 2: Group Above Threshold After Filtering
```python
# UMI1-B: 15 reads, 12 majority (80%), threshold=10
# Expected: PASSED (12 >= 10)
assert "BBBB" in passing_groups
```

### Test Results
```
✓ Group with 10 total reads, 6 majority (60%), threshold=10: EXCLUDED
✓ Group with 15 total reads, 12 majority (80%), threshold=10: PASSED
✓ Threshold correctly applied after UMI2 filtering
✓ Low-coverage groups properly rejected
```

## Benefits

1. **Quality Control**: Prevents low-coverage groups from proceeding
2. **Consistency**: Ensures all groups meet minimum coverage requirements
3. **Data Integrity**: Maintains high-quality consensus sequences
4. **Expected Behavior**: Aligns with user expectations for threshold parameter

## Backward Compatibility

This fix **improves** the filtering behavior without breaking existing functionality:
- Groups that passed before may now be excluded (stricter, better quality)
- No change to API or configuration
- Existing thresholds work as expected
- More stringent quality control

## Example Impact

With `threshold=10` and `umi2_majority_threshold=0.5`:

### Before Fix:
- 100 UMI1 groups pass initial threshold
- After UMI2 filtering, some groups may have < 10 reads
- Those low-coverage groups still proceed ❌

### After Fix:
- 100 UMI1 groups pass initial threshold
- After UMI2 filtering, groups with < 10 majority reads are excluded
- Only groups with >= 10 majority reads proceed ✓

## Recommendations

For users with the config change `threshold: 5` (reduced from 10), this fix is especially important:

```yaml
step6_umi_group:
  grouping: "UMI1"
  threshold: 5  # Lower threshold
  umi2_majority_threshold: 0.5
```

With a lower threshold, UMI2 filtering becomes even more critical for maintaining quality. This fix ensures that even with `threshold=5`, groups must still have at least 5 reads with the majority UMI2 species.

## Files Modified

1. `belt/scripts/umi_group.py` - Line 123 (added threshold check)
2. `belt/tests/test_umi2_filter_logic.py` - Added test case
3. `belt/tests/UMI2_MAJORITY_FILTER_UPDATE.md` - Updated documentation

## Verification

Run tests to verify:
```bash
cd belt
python3 tests/test_umi2_filter_logic.py
```

Expected output:
```
TEST: Threshold Re-applied After UMI2 Filtering
✓ All assertions passed!
✓ Threshold correctly re-applied after UMI2 filtering!
```

## Summary

This critical fix ensures that the minimum read threshold is enforced both before and after UMI2 filtering, maintaining high-quality data throughout the pipeline. It prevents low-coverage UMI groups from proceeding to consensus calling, which is essential for reliable results.



