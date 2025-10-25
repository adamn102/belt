# UMI Grouping Tests - Visual Inspection Guide

This document describes the comprehensive test suite for the updated UMI grouping functions that now support variable numbers of UMIs (UMI1, UMI2, UMI3, etc.).

## Test Suite Overview

### 1. `test_umi_consensus_long`
**Purpose:** Basic test for the standard `long_read_umi_grouping` function

**What it tests:**
- Reads BAM file and corrected UMIs CSV
- Groups reads by UMI1
- Creates filtered CSV and FASTQ files
- Validates output file structure

**Outputs:**
- `{sample}_filtered_umi_groups.csv` - UMI groups passing threshold
- `{sample}/` directory with FASTQ files per UMI group

---

### 2. `test_umi_grouping_with_fastq_inspection` ⭐
**Purpose:** Visual inspection of FASTQ sequences with full UMI annotations

**What it tests:**
- Generates FASTQ files for each UMI group
- Displays actual sequences from multiple files
- Shows UMI annotations in headers

**Visual Output Example:**
```
================================================================================
FASTQ FILES CREATED: 1
================================================================================

--- File 1: GTAGCCTTAACG_reads.fastq ---
Number of reads: 10

Read 1:
  Header: @1_UMI5_1|UMI1:GTAGCCTTAACG|UMI2:TCCGATGATCGA|UMI3:GGCATTAAGCTT|UMI4:CTGATGCGTTAC|UMI5:AATCGGCTTGAC
  Seq:    TATGGCAGGAGTAGCCTTAACGAGAAGCGGAGACAGAGACGAAGACCTCATCAAGACAGTCAGACTCATCAA...
  Qual:   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...
  Read Name: 1_UMI5_1
    UMI1:GTAGCCTTAACG
    UMI2:TCCGATGATCGA
    UMI3:GGCATTAAGCTT
    UMI4:CTGATGCGTTAC
    UMI5:AATCGGCTTGAC
```

---

### 3. `test_umi_grouping_different_columns`
**Purpose:** Test grouping by different UMI columns

**What it tests:**
- Groups by UMI1, UMI2, UMI3 in separate runs
- Validates that each grouping creates different groups
- Shows unique UMI values for each column

**Visual Output Example:**
```
================================================================================
Testing grouping by different UMI columns
Available UMI columns: ['UMI1', 'UMI2', 'UMI3', 'UMI4', 'UMI5']
================================================================================

--- Grouping by UMI1 ---
  Groups found: 1
  Unique UMI1 values: ['GTAGCCTTAACG']

--- Grouping by UMI2 ---
  Groups found: 1
  Unique UMI2 values: ['TCCGATGATCGA']

--- Grouping by UMI3 ---
  Groups found: 1
  Unique UMI3 values: ['GGCATTAAGCTT']
```

---

### 4. `test_umi_grouping_fast`
**Purpose:** Test the fast BAM-based grouping method

**What it tests:**
- Uses BAM tagging (ZU tag) for efficient grouping
- Creates sorted BAM file grouped by UMI
- Validates filtered CSV output

**Visual Output Example:**
```
================================================================================
Fast UMI Grouping Test
================================================================================

Filtered UMI Groups (1 groups):
          UMI1  ReadCount
0  GTAGCCTTAACG         10

✓ Fast grouping created:
  - Filtered CSV: /tmp/.../extract_umi_5-large_filtered_umi_groups.csv
  - Grouped BAM: /tmp/.../extract_umi_5-large.grouped.bam
```

---

### 5. `test_umi_grouping_fast_with_fastq` ⭐
**Purpose:** Fast grouping WITH FASTQ output for visual inspection

**What it tests:**
- Fast BAM-based grouping with FASTQ generation
- All UMIs included in FASTQ headers
- Efficient processing of large datasets

**Visual Output Example:**
```
================================================================================
Fast UMI Grouping with FASTQ Output
================================================================================

FASTQ files created: 1

--- GTAGCCTTAACG_reads.fastq ---
Number of reads: 10

First read:
  Header: @1_UMI5_1|UMI1:GTAGCCTTAACG|UMI2:TCCGATGATCGA|UMI3:GGCATTAAGCTT|UMI4:CTGATGCGTTAC|UMI5:AATCGGCTTGAC
  Seq:    TATGGCAGGAGTAGCCTTAACGAGAAGCGGAGACAGAGACGAAGACCTCATCAAGACAGTCAGACTCATCAA...
  Read Name: 1_UMI5_1
  UMI Info:
    - UMI1:GTAGCCTTAACG
    - UMI2:TCCGATGATCGA
    - UMI3:GGCATTAAGCTT
    - UMI4:CTGATGCGTTAC
    - UMI5:AATCGGCTTGAC

✓ Test passed: 1 FASTQ files created with full UMI annotations
```

---

## Running the Tests

### Run all tests:
```bash
cd /Users/anitido/Documents/GitHub/belt/belt
pytest tests/test_umi_group.py -v -s
```

### Run specific test for visual inspection:
```bash
# Best for visual sequence inspection
pytest tests/test_umi_group.py::test_umi_grouping_with_fastq_inspection -v -s

# Test fast version with FASTQ output
pytest tests/test_umi_group.py::test_umi_grouping_fast_with_fastq -v -s

# Test different grouping columns
pytest tests/test_umi_group.py::test_umi_grouping_different_columns -v -s
```

### Using the helper script:
```bash
./tests/run_umi_group_tests.sh
```
This will run all tests and save the output to `tests/umi_group_test_output.txt` for easy review.

---

## Understanding the Output

### FASTQ Header Format
```
@{read_name}|UMI1:{seq1}|UMI2:{seq2}|UMI3:{seq3}|...
```

**Example:**
```
@1_UMI5_1|UMI1:GTAGCCTTAACG|UMI2:TCCGATGATCGA|UMI3:GGCATTAAGCTT|UMI4:CTGATGCGTTAC|UMI5:AATCGGCTTGAC
```

**Breakdown:**
- Read Name: `1_UMI5_1`
- UMI1: `GTAGCCTTAACG` (12bp)
- UMI2: `TCCGATGATCGA` (12bp)
- UMI3: `GGCATTAAGCTT` (12bp)
- UMI4: `CTGATGCGTTAC` (12bp)
- UMI5: `AATCGGCTTGAC` (12bp)

---

## Key Features Tested

✅ **Variable UMI Numbers:** Handles any number of UMIs (not just 5pUMI/3pUMI)
✅ **Dynamic Column Detection:** Automatically detects available UMI columns
✅ **Full UMI Annotation:** All UMIs included in FASTQ headers (not just grouping UMI)
✅ **Flexible Grouping:** Can group by any UMI column (UMI1, UMI2, etc.)
✅ **Backwards Compatible:** Still works with legacy 5pUMI/3pUMI format
✅ **Visual Validation:** Tests display actual sequences for manual inspection

---

## Troubleshooting

### If tests are skipped:
The tests require these files to exist:
- `test-data/test_group_umi/extract_umi_5.fa`
- `test-data/test_group_umi/extract_umi_5-large.bam`
- `test-data/test_group_umi/extract_umi_5-large_corrected_umis.csv`

You may need to update the paths in the test config fixture to point to your actual test data location.

### If tests fail:
Check that:
1. The corrected UMIs CSV has a "Read Name" column
2. The CSV has UMI columns (UMI1, UMI2, etc.)
3. BAM file is properly formatted and indexed
4. samtools is available in PATH (for fast grouping tests)

---

## Next Steps

After running tests and inspecting the output:

1. **Verify UMI Sequences:** Check that all expected UMIs are present in headers
2. **Validate Grouping:** Ensure reads with same UMI are grouped together
3. **Check Sequence Quality:** Look at quality scores in FASTQ output
4. **Inspect Multiple Groups:** Review output from different UMI groups

The visual output makes it easy to verify that:
- All UMI information is correctly extracted
- Headers are properly formatted
- Sequences are complete and accurate
- Grouping is working as expected

