# BELT Extract Module - Comprehensive Test Report

**Generated:** 2024  
**Module:** `belt/scripts/extract.py`  
**Test File:** `belt/tests/test_extract.py`  
**Test Framework:** pytest

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Module Overview](#module-overview)
3. [Function Analysis](#function-analysis)
4. [Test Coverage Analysis](#test-coverage-analysis)
5. [Test Cases Breakdown](#test-cases-breakdown)
6. [Missing Test Coverage](#missing-test-coverage)
7. [Test Execution Guide](#test-execution-guide)
8. [Recommendations](#recommendations)

---

## Executive Summary

The `extract.py` module provides core functionality for extracting Unique Molecular Identifiers (UMIs) from aligned BAM files based on barcode motifs in a reference sequence. This report provides a comprehensive analysis of the current test coverage and identifies areas for improvement.

### Key Metrics

- **Total Test Functions:** 10
- **Lines of Code Tested:** ~194 lines
- **Test Coverage Areas:** Core functionality, edge cases, error handling
- **Coverage Status:** Good (with identified gaps)

### Overall Assessment

✅ **Strengths:**
- Comprehensive coverage of main functionality
- Good edge case testing
- Multiple motif scenarios tested
- Synthetic data testing included

⚠️ **Areas for Improvement:**
- Missing tests for backward compatibility (string motif input)
- No tests for empty BAM files
- Limited testing of partial UMI extraction
- Missing tests for malformed input files

---

## Module Overview

### Purpose

The `extract.py` module extracts UMIs from aligned sequencing reads by:
1. Locating barcode motifs in a reference sequence
2. Extracting corresponding sequences from aligned reads
3. Outputting results to CSV format

### Dependencies

- `pysam` - BAM file handling
- `Bio.SeqIO` (Biopython) - FASTA file parsing
- `pandas` - Data manipulation and CSV output
- `pathlib` - Path handling

### Key Features

- Supports single and multiple barcode motifs
- Sequential motif finding (each motif found after the previous)
- Minimum read length filtering
- Flexible output (directory or direct CSV file path)
- Error handling for invalid inputs

---

## Function Analysis

### 1. `extract_umis()`

**Purpose:** High-level function that orchestrates UMI extraction and CSV output generation.

**Parameters:**
- `reference` (str): Path to reference FASTA file
- `aligned_reads` (str): Path to input BAM file
- `min_read_length` (int): Minimum read length to process
- `barcode_motifs` (list): List of barcode motif sequences
- `output_path` (str): Output directory or CSV file path

**Key Functionality:**
- Determines output path type (directory vs. CSV file)
- Creates output directories if needed
- Calls `findUMI()` for actual extraction
- Handles extraction failures gracefully
- Generates CSV with UMI columns (UMI1, UMI2, etc.)

**Lines of Code:** 24-67

### 2. `findUMI()`

**Purpose:** Core function that performs UMI extraction from BAM files.

**Parameters:**
- `reference` (str): Path to reference FASTA file
- `aligned_reads` (str): Path to input BAM file
- `barcode_motifs` (list or str): Barcode motif(s) - accepts string for backward compatibility
- `min_read_length` (int): Minimum read length to process

**Returns:**
- `tuple`: `(motif_indices, raw_UMI_dict)` or `(None, None)` on failure

**Key Functionality:**
- Converts string motif to list (backward compatibility)
- Validates file extensions (.bam, .fa/.fasta)
- Reads reference sequence
- Finds motifs sequentially in reference
- Fetches reads from BAM in motif region
- Filters reads by minimum length
- Extracts UMIs using reference-to-query position mapping
- Returns dictionary with read names as keys and UMI lists as values

**Lines of Code:** 70-193

**Algorithm:**
1. Validate input files
2. Read reference sequence
3. Find each motif sequentially (starting search after previous motif)
4. Determine BAM fetch region (min to max motif positions)
5. For each read in region:
   - Skip if too short or no sequence
   - Map reference positions to query positions
   - Extract UMI sequences at motif positions
   - Store if at least one UMI found

---

## Test Coverage Analysis

### Test File Structure

The test file (`test_extract.py`) contains **10 test functions** organized into:

1. **Core Functionality Tests** (3 tests)
   - Single motif extraction
   - Multiple motifs (2, 5 motifs)
   - Full `extract_umis()` workflow

2. **Edge Case Tests** (5 tests)
   - No motifs found
   - Invalid file paths
   - Wrong file types
   - Output directory creation
   - CSV file path handling

3. **Synthetic Data Tests** (1 test)
   - Controlled test with known expected outputs

4. **Integration Tests** (1 test)
   - End-to-end extraction with CSV output validation

### Test Data

The tests use:
- **Real test data:** `test-data/test_extract/` directory
  - `extract_umi_1.fa` / `extract_umi_1.bam` (single motif)
  - `extract_umi_5.fa` / `extract_umi_5.bam` (five motifs)
- **Synthetic data:** Programmatically generated FASTA and BAM files
- **Temporary files:** Created on-the-fly for edge case testing

---

## Test Cases Breakdown

### ✅ Core Functionality Tests

#### 1. `test_findUMI` (Single Motif)
- **Purpose:** Test basic UMI extraction with one motif
- **Test Data:** `extract_umi_1.fa`, `extract_umi_1.bam`
- **Motif:** `["NNNNNNNNNNNN"]` (12 N's)
- **Validations:**
  - No duplicate read names
  - At least one UMI found
  - Correct number of UMIs per read (1)
  - Correct number of motif positions (1)
- **Status:** ✅ Covered

#### 2. `test_findUMI_five_motifs`
- **Purpose:** Test extraction with multiple motifs
- **Test Data:** `extract_umi_5.fa`, `extract_umi_5.bam`
- **Motifs:** 5 × `["NNNNNNNNNNNN"]`
- **Validations:**
  - No duplicate read names
  - At least one UMI found
  - Correct number of UMIs per read (5)
  - Correct number of motif positions (5)
- **Status:** ✅ Covered

#### 3. `test_findUMI_two_motifs`
- **Purpose:** Test extraction with two different motifs
- **Test Data:** `test-data/test-umi_map/` (optional, skipped if not available)
- **Motifs:** `["NNNNNNNNNNNNNNM", "NNNNNNNNNNNNNNM"]` (15 N's + M)
- **Validations:**
  - No duplicate read names
  - Correct UMI count per read (2)
  - Correct motif position count (2)
- **Status:** ✅ Covered (with skip if files unavailable)

#### 4. `test_extract_umis_function`
- **Purpose:** Test full `extract_umis()` workflow including CSV generation
- **Test Data:** `extract_umi_1.fa`, `extract_umi_1.bam`
- **Validations:**
  - CSV file created
  - Correct number of UMI columns
  - "Read Name" column present
  - At least one row of data
- **Status:** ✅ Covered

### ✅ Edge Case Tests

#### 5. `test_findUMI_no_motifs_found`
- **Purpose:** Test behavior when motifs don't exist in reference
- **Test Data:** Temporary FASTA without motifs
- **Expected:** Returns `(None, None)`
- **Status:** ✅ Covered

#### 6. `test_findUMI_invalid_files`
- **Purpose:** Test with non-existent file paths
- **Expected:** Returns `(None, None)`
- **Status:** ✅ Covered

#### 7. `test_findUMI_wrong_file_types`
- **Purpose:** Test with incorrect file extensions
- **Test Data:** `.txt` files instead of `.fa`/`.bam`
- **Expected:** Returns `(None, None)`
- **Status:** ✅ Covered

#### 8. `test_extract_umis_output_directory_creation`
- **Purpose:** Test automatic directory creation
- **Test Data:** Nested directory path that doesn't exist
- **Note:** Uses dummy files (will fail, but tests directory creation)
- **Status:** ⚠️ Partially covered (uses dummy files)

#### 9. `test_extract_umis_csv_output_path`
- **Purpose:** Test direct CSV file path handling
- **Test Data:** Direct CSV path instead of directory
- **Note:** Uses dummy files (will fail, but tests path handling)
- **Status:** ⚠️ Partially covered (uses dummy files)

### ✅ Synthetic Data Tests

#### 10. `test_findUMI_synthetic_data`
- **Purpose:** Test with programmatically generated data for controlled testing
- **Test Data:** Created on-the-fly with known UMI sequences
- **Reference:** Contains two 15N motifs
- **Read:** Contains specific UMI sequences (`ATCGATCGATCGAT`, `GCTAGCTAGCTAG`)
- **Validations:**
  - Motif indices found
  - UMI dictionary created
  - Exact UMI sequences match expected values
- **Status:** ✅ Covered

---

## Missing Test Coverage

### 🔴 Critical Missing Tests

#### 1. **Backward Compatibility: String Motif Input**
- **Code Location:** `extract.py:85-86`
- **Functionality:** `findUMI()` accepts a string instead of list
- **Test Needed:** Verify string motif is converted to list correctly
- **Impact:** Medium (backward compatibility feature)

#### 2. **Empty Barcode Motifs List**
- **Code Location:** `extract.py:115` (loop over empty list)
- **Functionality:** Behavior with empty `barcode_motifs` list
- **Test Needed:** Verify graceful handling or appropriate error
- **Impact:** Medium (edge case)

#### 3. **extract_umis() When findUMI() Returns None**
- **Code Location:** `extract.py:55-57`
- **Functionality:** Early return when extraction fails
- **Test Needed:** Verify no CSV file created and function returns gracefully
- **Impact:** High (error handling)

#### 4. **Reads Filtered by min_read_length**
- **Code Location:** `extract.py:148`
- **Functionality:** Reads shorter than `min_read_length` are skipped
- **Test Needed:** Verify short reads are excluded from results
- **Impact:** High (core filtering logic)

#### 5. **Reads with None query_sequence**
- **Code Location:** `extract.py:148`
- **Functionality:** Reads without sequence are skipped
- **Test Needed:** Verify reads with `query_sequence = None` are handled
- **Impact:** Medium (edge case)

### 🟡 Important Missing Tests

#### 6. **Partial UMI Extraction**
- **Code Location:** `extract.py:178`
- **Functionality:** Only adds reads if at least one UMI is found
- **Test Needed:** Test reads where some UMIs are `None` and some are found
- **Impact:** Medium (real-world scenario)

#### 7. **Empty BAM File**
- **Code Location:** `extract.py:147`
- **Functionality:** Behavior with valid but empty BAM file
- **Test Needed:** Verify returns empty dictionary, not error
- **Impact:** Medium (edge case)

#### 8. **BAM File with Reads but No Valid UMIs**
- **Code Location:** `extract.py:178`
- **Functionality:** Reads exist but don't overlap motif region or have no UMIs
- **Test Needed:** Verify empty dictionary returned
- **Impact:** Medium (real-world scenario)

#### 9. **extract_umis() with CSV File Path (Successful Case)**
- **Code Location:** `extract.py:38-43`
- **Functionality:** Direct CSV file path with successful extraction
- **Test Needed:** Complete test with real data (currently only tests with dummy files)
- **Impact:** Low (partially covered)

#### 10. **Reference File with Multiple Sequences**
- **Code Location:** `extract.py:102`
- **Functionality:** `SeqIO.read()` only reads first sequence
- **Test Needed:** Verify behavior with multi-sequence FASTA
- **Impact:** Low (documentation/clarification)

### 🟢 Nice-to-Have Tests

#### 11. **Malformed Reference File**
- **Code Location:** `extract.py:106-108`
- **Functionality:** Exception handling for invalid FASTA
- **Test Needed:** Test with malformed FASTA (missing header, invalid characters)
- **Impact:** Low (error handling)

#### 12. **BAM File Without Index**
- **Code Location:** `extract.py:147`
- **Functionality:** `bamfile.fetch()` behavior with unindexed BAM
- **Test Needed:** Verify error handling or fallback behavior
- **Impact:** Low (edge case)

#### 13. **Sequential Motif Search Behavior**
- **Code Location:** `extract.py:114-129`
- **Functionality:** Second motif found after first, not before
- **Test Needed:** Verify sequential search logic with motifs in different orders
- **Impact:** Low (algorithm validation)

#### 14. **Reads Outside Motif Region**
- **Code Location:** `extract.py:147`
- **Functionality:** Reads fetched only in motif region range
- **Test Needed:** Verify reads outside region are not processed
- **Impact:** Low (performance/validation)

#### 15. **UMI Length Validation**
- **Code Location:** `extract.py:173`
- **Functionality:** UMI sequence length must match motif length
- **Test Needed:** Test cases where extracted sequence length doesn't match
- **Impact:** Low (data validation)

---

## Test Execution Guide

### Running All Extract Tests

```bash
cd /Users/anitido/Documents/GitHub/belt/belt
pytest tests/test_extract.py -v
```

### Running Specific Test

```bash
# Single test
pytest tests/test_extract.py::test_findUMI -v

# Multiple tests
pytest tests/test_extract.py::test_findUMI tests/test_extract.py::test_findUMI_five_motifs -v
```

### Running with Output

```bash
# Show print statements
pytest tests/test_extract.py -v -s

# Show short traceback
pytest tests/test_extract.py -v --tb=short
```

### Running with Coverage

```bash
# Install coverage if needed
pip install pytest-cov

# Run with coverage report
pytest tests/test_extract.py --cov=scripts/extract --cov-report=html
```

### Test Dependencies

Ensure the following are installed:
- `pytest>=7.0.0`
- `pysam>=0.20.0`
- `biopython>=1.80`
- `pandas>=2.0.0`
- `numpy>=1.24.0`

---

## Recommendations

### Priority 1: Critical Tests (Implement First)

1. **Test `extract_umis()` failure handling**
   - Verify no CSV created when `findUMI()` returns `None`
   - Ensure graceful return without exceptions

2. **Test `min_read_length` filtering**
   - Create BAM with reads of varying lengths
   - Verify only reads ≥ `min_read_length` are processed

3. **Test backward compatibility (string motif)**
   - Pass string instead of list to `findUMI()`
   - Verify conversion and correct extraction

### Priority 2: Important Tests (Implement Next)

4. **Test partial UMI extraction**
   - Create scenario where some UMIs are found, some are `None`
   - Verify read is included if at least one UMI found

5. **Test empty BAM file**
   - Create valid but empty BAM file
   - Verify returns empty dictionary, not error

6. **Test empty `barcode_motifs` list**
   - Pass empty list
   - Verify appropriate error handling or graceful return

### Priority 3: Edge Cases (Implement When Time Permits)

7. **Test malformed reference file**
   - Invalid FASTA formats
   - Missing headers, invalid characters

8. **Test multi-sequence FASTA**
   - Verify only first sequence is used
   - Document behavior clearly

9. **Test unindexed BAM**
   - Verify error handling or fallback

### Code Quality Improvements

1. **Add type hints** for better IDE support and documentation
2. **Add docstring examples** showing usage patterns
3. **Consider logging** instead of print statements for better control
4. **Add unit tests** for helper logic (if extracted to separate functions)

### Test Infrastructure Improvements

1. **Create test fixtures** for common BAM/FASTA generation
2. **Add test data validation** to ensure test files are correct
3. **Add performance benchmarks** for large BAM files
4. **Add integration tests** with real-world data scenarios

---

## Test Statistics

### Current Test Count

- **Total Test Functions:** 10
- **Core Functionality:** 4 tests
- **Edge Cases:** 5 tests
- **Synthetic Data:** 1 test

### Code Coverage Estimate

Based on test analysis:
- **Core Functionality:** ~85% covered
- **Error Handling:** ~70% covered
- **Edge Cases:** ~60% covered
- **Overall Estimated Coverage:** ~75%

*Note: Actual coverage should be measured using `pytest-cov`*

### Test Execution Time

- **Expected Duration:** < 30 seconds for all tests
- **Slowest Tests:** Likely those using real BAM files (I/O bound)

---

## Conclusion

The `extract.py` module has **good test coverage** for core functionality and main use cases. The tests effectively validate:

✅ Multiple motif scenarios  
✅ Error handling for invalid inputs  
✅ CSV output generation  
✅ Synthetic data with known outputs  

However, there are **identifiable gaps** in coverage, particularly:

⚠️ Backward compatibility features  
⚠️ Edge cases with empty inputs  
⚠️ Failure scenarios in `extract_umis()`  
⚠️ Read filtering logic  

**Recommendation:** Implement Priority 1 and Priority 2 tests to achieve comprehensive coverage. The current test suite provides a solid foundation and can be extended incrementally.

---

## Appendix

### Test File Locations

- **Test File:** `belt/tests/test_extract.py`
- **Module Under Test:** `belt/scripts/extract.py`
- **Test Data:** `belt/test-data/test_extract/`
- **Test Reports:** `belt/tests/TestReports/`

### Related Documentation

- BELT Pipeline Overview: `docs/BELT_PIPELINE_OVERVIEW.md`
- Functions Reference: `FUNCTIONS_REFERENCE.md`
- Quick Reference: `docs/BELT_QUICK_REFERENCE.md`

---

**Report Generated:** 2024  
**Last Updated:** 2024  
**Maintainer:** BELT Development Team

