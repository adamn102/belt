# Partial UMI Handling in BELT

## Overview
This document describes how BELT handles reads with missing or partial UMIs (Unique Molecular Identifiers).

## Problem
When extracting UMIs from aligned reads, not all reads may align perfectly with all motif positions. This can result in:
- Reads with UMI1 but missing UMI2
- Reads with UMI2 but missing UMI1
- Reads with some UMIs but not others (when using multiple motifs)

Previously, this caused crashes when trying to map these partial UMIs to a whitelist, with the error:
```
TypeError: object of type 'float' has no len()
```

## Solution
The fix involves handling `None`/`NaN` values throughout the UMI extraction and mapping pipeline:

### 1. UMI Extraction (`findUMI` in `extract.py`)
- Initializes each read's UMI list with `None` values
- Only populates UMI positions that are successfully extracted
- Returns reads with at least one valid UMI (partial UMIs are kept)

```python
read_umis = [None] * len(barcodeMotifs)  # Initialize with None
# ... extract UMIs ...
if any(umi is not None for umi in read_umis):
    raw_UMI_dict[read.query_name] = read_umis
```

### 2. DataFrame Construction (`umi_extraction_map` in `umi_map.py`)
- Converts `None` values to pandas `NaN` automatically
- Logs statistics about missing UMIs per column
- Preserves missing values in CSV output as empty cells

```python
# Log statistics about missing UMIs
for umi_col in umi_names:
    missing_count = df[umi_col].isna().sum()
    if missing_count > 0:
        logger.log_message(f"{umi_col}: {missing_count}/{total_reads} reads have missing UMI")
```

### 3. Whitelist Mapping (`BK_scan` in `umi_map.py`)
- Preserves original input to track NaN positions
- Filters out `NaN` values before building the BK tree
- Validates all whitelist UMIs have the same length
- Maps `NaN` values to `(None, None)` without attempting hamming distance
- Only processes valid (non-NaN) UMIs through the BK tree
- Maps results back to original input including NaN values

```python
# Keep original input to preserve NaN values and their positions
original_umis = list(experimental_umis)

# Filter out NaN values from whitelist before building tree
wl_unique = [w for w in dict.fromkeys(whitelist) if pd.notna(w)]

# Make sure all strings in whitelist are the same length
if not all(len(umi) == len(wl_unique[0]) for umi in wl_unique):
    raise ValueError("All strings in the whitelist must be the same length.")

# Filter valid (non-NaN) UMIs from experimental data for processing
valid_umis = [umi for umi in original_umis if pd.notna(umi)]

# Process only valid UMIs... then map back to original
for umi in original_umis:
    if pd.isna(umi):
        result[umi] = (None, None)
    else:
        result[umi] = uniq_results[umi]
```

### 4. Output Generation
- CSV files correctly represent missing UMIs as empty cells
- Downstream tools can detect missing UMIs using `pd.isna()`

## Example Output

### Input Reads
```
read1: [ACGTACGTACGT, TGCATGCATGCA]  # Has both UMI1 and UMI2
read2: [GCTAGCTAGCTA, None]           # Has only UMI1, missing UMI2
read3: [None, ATCGATCGATCG]           # Has only UMI2, missing UMI1
```

### CSV Output
```csv
Read Name,UMI1,UMI2
read1,ACGTACGTACGT,TGCATGCATGCA
read2,GCTAGCTAGCTA,
read3,,ATCGATCGATCG
```

### After Whitelist Mapping
```csv
Read Name,UMI1,UMI2
read1,ACGTACGTACGT,TGCATGCATGCA
read2,GCTAGCTAGCTA,
read3,,ATCGATCGATCG
```
Note: Missing UMIs remain as empty cells after mapping.

## Test Coverage
The following tests validate partial UMI handling:

1. **`test_umi_extraction_map_five_motifs`**
   - Integration test with 5 motifs
   - Uses extract_umi_5-large.bam test data
   - Validates basic multi-motif extraction

2. **`test_umi_extraction_map_partial_umis`**
   - Integration test using real BAM data
   - Validates full pipeline with partial UMIs
   - Reports statistics on UMI coverage

3. **`test_whitelist_mapping_with_partial_umis`**
   - Unit test for `BK_scan` with NaN values
   - Verifies NaN handling in whitelist mapping
   - Tests both `None` and `float('nan')` representations

4. **`test_umi_extraction_real_data_with_partial_umis`** ⭐ **CRITICAL TEST**
   - Real production data: LP249-NL43-A3G-1-A (55,342 reads)
   - **18,581 reads (33.6%) have partial UMIs**
   - Tests complete pipeline with whitelist mapping
   - Validates length mismatch handling
   - Demonstrates real-world robustness

5. **`test_dataframe_handles_none_umis`**
   - Unit test for DataFrame construction
   - Validates CSV round-trip with missing UMIs
   - Ensures NaN values are preserved correctly

## Usage Recommendations

### When to Expect Partial UMIs
- Low quality reads with alignment gaps
- Reads that don't span all motif positions
- Truncated reads at the ends of molecules
- Incorrect motif positions in reference

### Downstream Analysis
When using BELT output with partial UMIs:

```python
import pandas as pd

# Read UMI data
df = pd.read_csv("sample_Paired_UMIs_read_names.csv")

# Filter for reads with all UMIs present
complete_umis = df.dropna()

# Or filter for reads with at least one UMI
any_umi = df[df.notna().any(axis=1)]

# Or analyze by specific UMI column
umi1_present = df[df['UMI1'].notna()]
```

### Logging Output
When running `belt umi_map`, you'll see logs like:
```
[BELT] UMI1: 2/10 reads have missing UMI (20.0%)
[BELT] UMI2: All reads have valid UMI
```

This helps you understand the quality and completeness of your UMI extraction.

## Benefits
1. **Robustness**: No crashes when reads have partial UMIs
2. **Transparency**: Users can see which reads are missing which UMIs
3. **Flexibility**: Downstream analysis can choose how to handle partial UMIs
4. **Data preservation**: All extracted data is kept, nothing is silently dropped

## Related Files
- `/belt/scripts/extract.py`: `findUMI()` function
- `/belt/scripts/umi_map.py`: `umi_extraction_map()`, `BK_scan()` functions
- `/belt/tests/test_umi_map.py`: Test suite for partial UMI handling

