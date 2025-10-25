# BELT Consensus Module Tests

## Overview

The consensus module (`scripts/consensus.py`) has been fully migrated from a bash script (`consensus.sh`) to a native Python implementation. This document describes the implementation and test suite.

## Consensus Pipeline

The `run_consensus_for_umi()` function implements a 5-step consensus calling pipeline:

1. **Align reads to reference** - Uses minimap2 to align reads and samtools to sort
2. **Index BAM** - Creates BAM index for efficient access
3. **Call variants** - Uses medaka_variant to identify variants
4. **Normalize VCF** - Normalizes and compresses VCF with bcftools
5. **Generate consensus** - Creates consensus FASTA from reference + variants

### Key Features

- ✅ **Pure Python implementation** - No external shell scripts needed
- ✅ **Comprehensive error handling** - Captures and logs all errors
- ✅ **Detailed logging** - Step-by-step progress tracking
- ✅ **Parallel processing** - Supports multiple UMI groups simultaneously
- ✅ **Automatic cleanup** - Removes temporary FASTQ files after processing

## Function Signature

```python
def run_consensus_for_umi(
    reads_fq: Path,        # Path to FASTQ file with reads for this UMI
    ref_fa: Path,          # Path to reference FASTA
    meta_dir: Path,        # Directory for intermediate files (BAM, VCF)
    con_dir: Path,         # Directory for final consensus output
    bam_dir: Path,         # Directory for BAM files (optional)
    umi_threads: int,      # Number of threads to use
    model: str,            # Medaka model name
    log_file: Path,        # Path to log file
    umi_safe: str          # Sanitized UMI name for output files
) -> str:                  # Returns path to consensus FASTA
```

## Test Suite

### Test Structure

The test suite (`tests/test_consensus.py`) contains 16 comprehensive tests organized into 5 test classes:

#### 1. **TestSanitizeUmiName** (7 tests)
Tests UMI name sanitization for safe file naming:
- Simple alphanumeric UMIs
- UMIs with numbers and allowed characters (`.`, `_`, `-`)
- UMIs with special characters (replaced with `_`)
- Handling of `None` and empty strings (converted to `"NA"`)

#### 2. **TestWriteGroupFastq** (2 tests)
Tests FASTQ file generation from BAM reads:
- Writing single read to FASTQ format
- Writing multiple reads to FASTQ format
- Verifying correct FASTQ formatting (header, sequence, quality)

#### 3. **TestCheckDependencies** (2 tests)
Tests external tool dependency checking:
- Verifies all required tools are available (minimap2, samtools, medaka_variant, bcftools)
- Tests error message when dependencies are missing

#### 4. **TestRunConsensusForUmi** (3 tests)
Tests the main consensus calling function:
- **Function signature test** - Verifies expected parameters
- **Basic workflow test** - Complete end-to-end consensus calling with synthetic data
- **Cleanup test** - Verifies temporary files are removed

#### 5. **TestConsensusIntegration** (1 test)
Integration test using real BELT pipeline data:
- Tests consensus calling with actual grouped BAM files
- Validates output with real sequencing data
- Skipped if test data is not available

### Running the Tests

#### Run all consensus tests:
```bash
cd belt
source ~/miniforge3/etc/profile.d/conda.sh
conda activate belt-pipeline-medaka
pytest tests/test_consensus.py -v
```

#### Run with detailed output:
```bash
pytest tests/test_consensus.py -v -s
```

#### Run specific test class:
```bash
pytest tests/test_consensus.py::TestSanitizeUmiName -v
```

#### Run single test:
```bash
pytest tests/test_consensus.py::TestRunConsensusForUmi::test_consensus_basic_workflow -v
```

## Test Results

### Current Status
✅ **15 tests passed**  
⏭️  **1 test skipped** (integration test requires real data)

### Expected Output
```
============================= test session starts ==============================
tests/test_consensus.py::TestSanitizeUmiName::test_simple_umi PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_umi_with_numbers PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_umi_with_allowed_chars PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_umi_with_special_chars PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_umi_with_spaces PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_none_umi PASSED
tests/test_consensus.py::TestSanitizeUmiName::test_empty_umi PASSED
tests/test_consensus.py::TestWriteGroupFastq::test_write_single_read PASSED
tests/test_consensus.py::TestWriteGroupFastq::test_write_multiple_reads PASSED
tests/test_consensus.py::TestCheckDependencies::test_check_dependencies_available PASSED
tests/test_consensus.py::TestCheckDependencies::test_check_dependencies_message PASSED
tests/test_consensus.py::TestRunConsensusForUmi::test_consensus_function_signature PASSED
tests/test_consensus.py::TestRunConsensusForUmi::test_consensus_basic_workflow PASSED
tests/test_consensus.py::TestRunConsensusForUmi::test_consensus_cleanup_temp_files PASSED
tests/test_consensus.py::TestConsensusIntegration::test_consensus_with_real_grouped_bam SKIPPED
tests/test_consensus.py::test_module_imports PASSED

======================== 15 passed, 1 skipped in 7.70s =========================
```

## Dependencies

The consensus module requires the following external tools:
- **minimap2** - Read alignment
- **samtools** - BAM manipulation and indexing
- **medaka** - Variant calling (specifically `medaka_variant`)
- **bcftools** - VCF manipulation and consensus generation

These are all included in the `belt-pipeline-medaka` conda environment.

## File Organization

### Module Files
- `scripts/consensus.py` - Main consensus module (Python)
- `scripts/consensus.sh` - Legacy bash script (can be removed)

### Test Files
- `tests/test_consensus.py` - Comprehensive test suite
- `tests/CONSENSUS_TESTS_README.md` - This documentation

### Output Structure
```
output_dir/
└── sample_name/
    ├── consensus/           # Final consensus FASTA files
    │   ├── UMI1_consensus.fasta
    │   ├── UMI2_consensus.fasta
    │   └── ...
    ├── bam/                 # BAM files (optional)
    └── tmp/                 # Temporary files (cleaned up)
        └── UMI1/
            ├── UMI1_reads.fastq  (deleted after processing)
            └── meta/
                ├── aligned.bam
                ├── aligned.bam.bai
                ├── medaka.vcf
                └── medaka.norm.vcf.gz
```

## Example Usage

### Basic Usage
```python
from pathlib import Path
from scripts.consensus import run_consensus_for_umi

# Generate consensus for one UMI group
consensus_fasta = run_consensus_for_umi(
    reads_fq=Path("reads/UMI1_reads.fastq"),
    ref_fa=Path("reference.fa"),
    meta_dir=Path("tmp/UMI1/meta"),
    con_dir=Path("consensus"),
    bam_dir=Path("bam"),
    umi_threads=2,
    model="r1041_e82_400bps_sup_variant_v5.0.0",
    log_file=Path("consensus.log"),
    umi_safe="UMI1"
)

print(f"Consensus created: {consensus_fasta}")
```

### Process Grouped BAM
```python
from scripts.consensus import process_grouped_bam

# Process all UMI groups in a BAM file
consensus_files = process_grouped_bam(
    grouped_bam="grouped.bam",
    reference="reference.fa",
    output_dir="output",
    sample="sample1",
    tag="ZU",
    procs=4,              # Process 4 UMIs in parallel
    umi_threads=2,        # Use 2 threads per UMI
    model="r1041_e82_400bps_sup_variant_v5.0.0"
)

print(f"Generated {len(consensus_files)} consensus sequences")
```

## Migration Notes

### Differences from consensus.sh

The Python implementation maintains **100% functional compatibility** with the original bash script while adding:

1. **Better error handling** - All errors are caught and logged
2. **Progress tracking** - Step-by-step logging with timestamps
3. **Parallel processing** - Multiple UMIs processed simultaneously
4. **Cleaner code** - Easier to maintain and extend
5. **Testable** - Comprehensive test suite ensures reliability

### Breaking Changes
**None** - The function signatures and behavior are identical to the bash script interface.

### Performance
The Python implementation has **equivalent performance** to the bash script, as it calls the same external tools (minimap2, samtools, medaka, bcftools).

## Troubleshooting

### Common Issues

#### 1. Missing Dependencies
```
RuntimeError: Missing required tools: medaka_variant
```
**Solution:** Activate the conda environment:
```bash
conda activate belt-pipeline-medaka
```

#### 2. VCF File Not Found
```
FileNotFoundError: VCF file not found at meta/medaka.vcf
```
**Solution:** This usually means medaka_variant failed. Check the log file for medaka error messages.

#### 3. Empty Consensus FASTA
If the consensus FASTA exists but is very short, check:
- Reference FASTA is correct
- Reads align to the reference
- Sufficient read coverage (at least 3-5 reads recommended)

## Future Enhancements

Potential improvements for the consensus module:

- [ ] Add support for different consensus callers (racon, canu, etc.)
- [ ] Implement quality filtering for input reads
- [ ] Add coverage depth statistics to log
- [ ] Support for paired-end reads
- [ ] Parallel processing within each UMI (currently sequential steps)
- [ ] Option to retain intermediate files for debugging

## Contributing

When modifying the consensus module:

1. **Run the test suite** before committing:
   ```bash
   pytest tests/test_consensus.py -v
   ```

2. **Add tests** for new features

3. **Update documentation** if function signatures change

4. **Check backwards compatibility** with existing pipelines

## License

This module is part of the BELT pipeline and follows the same license.

