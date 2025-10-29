# BELT Software Functions Reference

This document provides a comprehensive reference for all functions in the BELT (BalazsLab Error-correction Longread Tools) software package. Each function is documented with its purpose, pipeline usage, and test validation coverage.

---

## Table of Contents

1. [Core Pipeline Functions](#core-pipeline-functions)
   - [map](#map)
   - [extract](#extract)
   - [umi-map](#umi-map)
   - [umi-group](#umi-group)
   - [consensus](#consensus)
   - [whitelist](#whitelist)
2. [Analysis Functions](#analysis-functions)
   - [mutation-call](#mutation-call)
   - [mutation-summary](#mutation-summary)
   - [mutation-analyze](#mutation-analyze)
   - [umi-meta-analysis](#umi-meta-analysis)
3. [Utility Functions](#utility-functions)
   - [run-pipeline](#run-pipeline)

---

## Core Pipeline Functions

### map

**Purpose**: Map FASTQ reads to reference sequences using minimap2 with quality filtering, sorting, indexing, and optional demultiplexing.

**What it does**:
- Filters FASTQ files by read length and quality score
- Maps reads to reference using minimap2 (optimized for N-base mapping)
- Sorts and indexes BAM files with samtools
- Optionally demultiplexes reads based on barcode sequences
- Generates HTML quality reports

**Input**:
- Reference FASTA file
- FASTQ file(s) - can be multiple files or glob patterns
- Optional: Sample sheet CSV or config YAML for demultiplexing

**Output**:
- Filtered FASTQ files
- Sorted, indexed BAM files
- Demultiplexed FASTQ files (if demultiplexing enabled)
- Quality reports (HTML/JSON)

**Pipeline Position**: **Step 1** - First step in BELT pipeline workflow

**Pipeline Usage**:
```
Raw FASTQ → [MAP] → Filtered/Mapped BAM → Extract → UMI Mapping → ...
```

**Test Coverage**:
- Basic mapping functionality tested through pipeline integration tests
- Demultiplexing tested in `pipeline_test.py`
- Quality filtering validated through output verification

**Key Features**:
- Multi-file processing with glob patterns
- CSV and YAML-based demultiplexing support
- Comprehensive quality filtering
- Parallel processing support

---

### extract

**Purpose**: Extract Unique Molecular Identifiers (UMIs) from aligned BAM files based on predefined barcode motifs in the reference sequence.

**What it does**:
- Searches for motif patterns in reference sequence
- Extracts corresponding UMI sequences from aligned reads
- Handles multiple motifs per read (e.g., 5' and 3' UMIs)
- Filters reads by minimum length and sequence availability
- Outputs CSV with Read Name and extracted UMI columns

**Input**:
- Reference FASTA file
- Aligned BAM file(s)
- Barcode motif(s) - single or list
- Minimum read length threshold

**Output**:
- CSV file with columns: `Read Name`, `UMI1`, `UMI2`, ... (one per motif)

**Pipeline Position**: **Step 2** - After mapping step

**Pipeline Usage**:
```
Mapped BAM → [EXTRACT] → UMI CSV → UMI Mapping → UMI Grouping → ...
```

**Test Coverage**: ✅ **Comprehensive** (14 tests in `test_extract.py`)
- Core functionality: Single motif, multiple motifs (2, 5 motifs)
- Edge cases: No motifs found, invalid files, wrong file types
- Error handling: Directory creation, CSV output paths
- Synthetic data: Controlled testing with known sequences
- Performance: Large file processing

**Key Features**:
- Multi-motif support
- Flexible motif matching
- Quality filtering
- Robust error handling

**Current Limitations**:
- Missing read quality metrics (mapping quality, read length)
- No UMI extraction success rate reporting
- Limited processing statistics

---

### umi-map

**Purpose**: Map raw UMIs to a whitelist using hamming distance matching and generate paired UMI output with correction statistics.

**What it does**:
- Takes raw UMIs from extract step
- Maps UMIs to whitelist using BK-tree for efficient hamming distance search
- Corrects UMI errors within specified maximum distance
- Generates paired UMI CSV files
- Provides mapping statistics and success rates

**Input**:
- Reference FASTA file
- Aligned BAM file
- Raw UMI CSV (from extract step)
- Whitelist CSV file
- Maximum hamming distance threshold

**Output**:
- CSV with corrected UMIs: `Read Name`, `UMI1`, `UMI2`, ... (paired UMIs)
- Mapping statistics (success rates, distance distributions)

**Pipeline Position**: **Step 3** - After extract, before umi-group

**Pipeline Usage**:
```
Raw UMI CSV → [UMI-MAP] → Corrected UMI CSV → UMI Grouping → Consensus → ...
```

**Test Coverage**: ✅ **Comprehensive** (Multiple test files)
- Core functionality: `test_umi_map.py`
- Max distance testing: Multiple test files for distance validation
- Edge cases: Non-matching UMIs, whitelist validation
- Performance: Large whitelist handling

**Key Features**:
- BK-tree based efficient matching
- Configurable maximum distance
- Handles missing/NaN UMIs gracefully
- Detailed mapping statistics

**Pipeline-Specific Usage**:
- **ERRSeq2**: Corrects UMIs for mutation analysis
- **LoopR**: Critical for error correction workflow
- **RecombRate**: Tracks UMI associations for recombination detection

---

### umi-group

**Purpose**: Group long-read UMIs based on corrected UMI data, filtering by threshold and generating grouped FASTQ files for consensus calling.

**What it does**:
- Reads corrected UMI CSV from umi-map step
- Groups reads by specified UMI column (UMI1, UMI2, 5pUMI, 3pUMI, etc.)
- Filters groups by minimum read count threshold
- Extracts reads from BAM file grouped by UMI
- Generates grouped FASTQ files for consensus calling
- Creates grouped BAM file with UMI tags

**Input**:
- Reference FASTA file
- Aligned BAM file
- Corrected UMI CSV file
- Grouping column name (e.g., "UMI1", "5pUMI")
- Threshold (minimum reads per group)

**Output**:
- Filtered UMI groups CSV
- Grouped FASTQ files (one per UMI group)
- Grouped BAM file with UMI tags

**Pipeline Position**: **Step 4** - After umi-map, before consensus

**Pipeline Usage**:
```
Corrected UMI CSV → [UMI-GROUP] → Grouped FASTQ/BAM → Consensus Calling → ...
```

**Test Coverage**: ✅ **Comprehensive** 
- Core functionality: `test_umi_group.py`
- Real data testing: `test_umi_group_real.py`
- Edge cases: Threshold filtering, missing UMIs
- Validation: Read extraction accuracy

**Key Features**:
- Flexible grouping by any UMI column
- Threshold-based filtering
- Supports variable numbers of UMIs
- Legacy format support (5pUMI, 3pUMI)
- UMI tag assignment in BAM files

**Pipeline-Specific Usage**:
- **ERRSeq2**: Groups reads for mutation analysis
- **LoopR**: Core function for error correction grouping
- **RecombRate**: Groups reads for template swapping analysis

---

### consensus

**Purpose**: Generate high-quality consensus sequences from grouped BAM files using medaka for variant calling and bcftools for consensus generation.

**What it does**:
- Processes grouped BAM files (one per UMI group)
- Aligns reads to reference using minimap2
- Calls variants using medaka variant caller
- Normalizes and compresses VCF files
- Generates consensus FASTA sequences using bcftools consensus
- Supports parallel processing of multiple UMI groups

**Input**:
- Grouped BAM file with UMI tags
- Reference FASTA file
- Sample name
- Medaka model specification
- Threading configuration

**Output**:
- Consensus FASTA files (one per UMI group)
- Intermediate BAM files
- VCF files (variant calls)
- Manifest file listing all outputs

**Pipeline Position**: **Step 5** - After umi-group, final step in core pipeline

**Pipeline Usage**:
```
Grouped BAM → [CONSENSUS] → Consensus FASTA → Mutation Analysis → ...
```

**Test Coverage**: ✅ **Comprehensive**
- Core functionality: `test_consensus.py`
- Sequence validation: `test_consensus_sequence_data.py`
- CSV verification: `verify_consensus_csv.py`
- Edge cases: Deletion handling, consecutive deletions

**Key Features**:
- Medaka-based variant calling
- Parallel UMI processing
- Configurable threading
- Comprehensive output organization
- Manifest file generation

**Pipeline-Specific Usage**:
- **ERRSeq2**: Generates consensus for mutation rate analysis
- **LoopR**: Core function for error-corrected sequences
- **RecombRate**: Generates consensus for recombination analysis

---

### whitelist

**Purpose**: Generate a whitelist of UMIs from sequencing data using clustering and correction algorithms.

**What it does**:
- Extracts raw UMIs from aligned BAM files
- Filters outlier barcodes
- Clusters UMIs using UMI-tools (or other methods)
- Corrects UMI errors within clusters
- Generates whitelist CSV files
- Creates visualization reports (hamming distances, UMI distributions)

**Input**:
- Reference FASTA file
- Aligned BAM file
- Barcode motif sequence(s)
- Clustering method

**Output**:
- Whitelist CSV file
- Corrected UMI CSV files
- HTML visualization reports
- UMI statistics

**Pipeline Position**: **Pre-processing step** - Used before umi-map to generate whitelist

**Pipeline Usage**:
```
Mapped BAM → [WHITELIST] → Whitelist CSV → Used by UMI-MAP → ...
```

**Test Coverage**: ✅ **Good**
- Core functionality: `test_whitelist.py`
- Clustering validation
- Statistics verification

**Key Features**:
- Multiple clustering methods
- Outlier filtering
- Statistical analysis
- Visualization generation

**Note**: Two implementations exist:
- `whitelist.py` - Main implementation
- `whitelist_longread.py` - Specialized for long-read data

---

## Analysis Functions

### mutation-call

**Purpose**: Analyze mutations for individual samples by comparing consensus sequences to reference sequences.

**What it does**:
- Loads consensus and reference FASTA files
- Performs pairwise sequence alignment
- Identifies mutations (substitutions, insertions, deletions)
- Calculates mutation rates and statistics
- Generates per-read mutation analysis

**Input**:
- Consensus FASTA file
- Reference FASTA file
- Sample name
- Output directory

**Output**:
- Per-read mutation CSV file
- Mutation analysis JSON report
- Summary statistics

**Pipeline Position**: **Post-processing** - After consensus calling

**Pipeline Usage**:
```
Consensus FASTA → [MUTATION-CALL] → Mutation Analysis → Mutation Summary → ...
```

**Test Coverage**: ✅ **Comprehensive**
- Core functionality: `test_mutation_caller.py`
- Edge cases: `test_consecutive_deletions.py`
- Validation: `verify_mutation_fix.py`
- Documentation: `MUTATION_CALLER_FIX.md`

**Key Features**:
- BioPython-based alignment
- Comprehensive mutation classification
- Per-read detailed analysis
- Statistical calculations

**Pipeline-Specific Usage**:
- **ERRSeq2**: Primary analysis function for mutation rate calculation
- **LoopR**: Used for quality assessment
- **RecombRate**: Not typically used

---

### mutation-summary

**Purpose**: Generate summary analysis from individual sample mutation results, aggregating data across multiple samples.

**What it does**:
- Aggregates mutation data from multiple samples
- Generates summary statistics and distributions
- Creates consensus sequence data CSV
- Produces mutation analysis reports
- Calculates sample-level statistics

**Input**:
- Input directory containing sample mutation results
- Output directory for summary files

**Output**:
- Mutation summary CSV
- Consensus sequence data CSV
- Mutation analysis report JSON
- Per-sample statistics

**Pipeline Position**: **Post-processing** - After mutation-call, aggregates results

**Pipeline Usage**:
```
Sample Results → [MUTATION-SUMMARY] → Summary Analysis → Visualization → ...
```

**Test Coverage**: ✅ **Good**
- Core functionality: `test_mutation_analysis.py`
- Data validation: `test_consensus_sequence_data.py`
- CSV verification: `verify_consensus_csv.py`

**Key Features**:
- Multi-sample aggregation
- Statistical summaries
- Comprehensive reporting
- Data validation

**Pipeline-Specific Usage**:
- **ERRSeq2**: Critical for mutation rate analysis
- **LoopR**: Used for quality assessment
- **RecombRate**: Not typically used

---

### mutation-analyze

**Purpose**: Legacy command for analyzing mutations in consensus sequences compared to reference sequences (wrapper around mutation_summary functionality).

**What it does**:
- Analyzes all consensus sequences in a directory
- Compares to reference sequences
- Generates mutation analysis reports
- Supports parallel processing

**Input**:
- Directory containing compiled consensus FASTA files
- Directory containing reference FASTA files
- Optional sample list

**Output**:
- Mutation analysis report JSON
- Mutation summary CSV
- Per-read mutation data

**Pipeline Position**: **Post-processing** - Legacy alternative to mutation-call + mutation-summary

**Pipeline Usage**:
```
Consensus FASTA Directory → [MUTATION-ANALYZE] → Mutation Analysis → ...
```

**Test Coverage**: ✅ **Good**
- Core functionality: `test_mutation_analysis.py`
- Integration tests

**Note**: This is a legacy command that combines functionality of `mutation-call` and `mutation-summary`. Consider using the newer separate commands for better control.

---

### umi-meta-analysis

**Purpose**: Perform meta-analysis on UMI mapping data to analyze UMI patterns, diversity, and statistics across samples.

**What it does**:
- Analyzes UMI patterns for multiple samples
- Calculates UMI diversity metrics
- Generates UMI statistics and reports
- Identifies UMI composition patterns
- Provides cross-sample comparisons

**Input**:
- Directory containing UMI mapping CSV files
- Optional sample list
- Output directory

**Output**:
- UMI meta-analysis report JSON
- UMI statistics CSV
- Pattern analysis results

**Pipeline Position**: **Post-processing** - After umi-map, analyzes UMI patterns

**Pipeline Usage**:
```
UMI Mapping Files → [UMI-META-ANALYSIS] → Meta-Analysis Report → ...
```

**Test Coverage**: ⚠️ **Limited**
- Basic functionality tested through integration
- Needs more comprehensive test coverage

**Key Features**:
- Cross-sample analysis
- UMI diversity calculations
- Pattern identification
- Statistical reporting

**Pipeline-Specific Usage**:
- **ERRSeq2**: Analyzes UMI diversity for quality assessment
- **LoopR**: Validates UMI extraction quality
- **RecombRate**: Critical for analyzing template swapping patterns

---

## Utility Functions

### run-pipeline

**Purpose**: Run the entire BELT pipeline using a YAML configuration file, executing all steps in sequence.

**What it does**:
- Loads configuration from YAML file
- Executes extract → umi-map → umi-group in sequence
- Provides unified pipeline interface

**Input**:
- YAML configuration file with all parameters

**Output**:
- All intermediate and final outputs from pipeline steps

**Pipeline Position**: **Orchestration** - Runs entire pipeline workflow

**Pipeline Usage**:
```
Config YAML → [RUN-PIPELINE] → Complete Pipeline Output
```

**Test Coverage**: ✅ **Good**
- Integration tests: `pipeline_test.py`
- End-to-end workflow validation

**Key Features**:
- Configuration-driven execution
- Unified interface
- Complete workflow automation

**Note**: This is a simplified pipeline runner. For production use, consider using Snakemake workflows (as in the pipeline directories) for better control, parallelization, and reproducibility.

---

## Pipeline Workflow Summary

### Standard BELT Pipeline Flow:

```
1. map              → Filter and map FASTQ to reference
2. extract          → Extract UMIs from aligned reads
3. umi-map          → Map UMIs to whitelist (correct errors)
4. umi-group        → Group reads by corrected UMIs
5. consensus        → Generate consensus sequences
6. mutation-call    → Analyze mutations (ERRSeq2)
7. mutation-summary → Summarize mutation analysis (ERRSeq2)
8. umi-meta-analysis → Analyze UMI patterns (RecombRate)
```

### Pipeline-Specific Workflows:

**ERRSeq2 (Mutation Rate Analysis)**:
```
map → extract → umi-map → umi-group → consensus → mutation-call → mutation-summary
```

**LoopR (Error-Corrected Long Read Sequencing)**:
```
map → extract → umi-map → umi-group → consensus
```

**RecombRate (Template Swapping Rate Measurement)**:
```
map → extract → umi-map → umi-group → consensus → umi-meta-analysis
```

---

## Test Coverage Summary

| Function | Test Files | Coverage Status |
|----------|------------|-----------------|
| **map** | Pipeline integration tests | ✅ Good |
| **extract** | `test_extract.py` (14 tests) | ✅ Comprehensive |
| **umi-map** | `test_umi_map.py`, max distance tests | ✅ Comprehensive |
| **umi-group** | `test_umi_group.py`, `test_umi_group_real.py` | ✅ Comprehensive |
| **consensus** | `test_consensus.py`, validation tests | ✅ Comprehensive |
| **whitelist** | `test_whitelist.py` | ✅ Good |
| **mutation-call** | `test_mutation_caller.py`, edge case tests | ✅ Comprehensive |
| **mutation-summary** | `test_mutation_analysis.py` | ✅ Good |
| **mutation-analyze** | Integration tests | ✅ Good |
| **umi-meta-analysis** | Integration tests | ⚠️ Limited |
| **run-pipeline** | `pipeline_test.py` | ✅ Good |

---

## Additional Notes

### Function Dependencies

- **extract** depends on: `map` (output)
- **umi-map** depends on: `extract` (output), `whitelist` (whitelist file)
- **umi-group** depends on: `umi-map` (corrected UMIs)
- **consensus** depends on: `umi-group` (grouped BAM)
- **mutation-call** depends on: `consensus` (consensus FASTA)

### Performance Considerations

- **map**: Can process multiple FASTQ files in parallel
- **extract**: Processes reads sequentially from BAM file
- **umi-map**: Uses BK-tree for efficient whitelist matching
- **umi-group**: Filters by threshold for efficiency
- **consensus**: Supports parallel UMI processing

### Output File Conventions

- Extract: `{sample}.belt_extract.csv`
- UMI Map: `{sample}_Paired_UMIs_read_names.csv`
- UMI Group: `{sample}_filtered_umi_groups.csv`, `{sample}.grouped.bam`
- Consensus: `{sample}/{umi}.fasta` in `fasta/` directory
- Mutation Call: `{sample}_per_read_mutations.csv`
- Mutation Summary: `mutation_summary.csv`, `mutation_analysis_report.json`

---

## References

- Main BELT Documentation: `README.md`
- Pipeline Documentation: `pipelines/README.md`
- Development Guide: `DEVELOPMENT.md`
- Test Documentation: `tests/` directory

---

**Last Updated**: 2024
**BELT Version**: 0.1.0
