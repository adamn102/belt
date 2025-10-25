# BELT Documentation

This directory contains all documentation for the BELT (Barcode Extraction and Lineage Tracing) pipeline.

## 📚 Documentation Structure

```
docs/
├── README.md (you are here)
│
├── Main Documentation/
│   ├── BELT_PIPELINE_OVERVIEW.md        # Complete architecture guide (1000+ lines)
│   └── BELT_QUICK_REFERENCE.md          # One-page quick reference
│
├── reports/                              # Validation and evaluation reports
│   ├── PIPELINE_VALIDATION_REPORT.md    # Full pipeline validation (threshold=5)
│   ├── EVALUATION_SUMMARY.md            # Filtered run summary (threshold=1000)
│   ├── FILTERED_PIPELINE_EVALUATION.md  # Detailed threshold analysis
│   └── SNAKEMAKE_UPDATE_SUMMARY.md      # Snakemake update notes
│
└── test-documentation/                   # Testing documentation
    ├── CONSENSUS_TESTS_README.md        # Consensus testing guide
    ├── PARTIAL_UMI_HANDLING.md          # Partial UMI implementation
    └── UMI_GROUP_TESTS_README.md        # UMI grouping tests
```

---

## 🚀 Quick Start

**New to BELT?** Start here:
1. Read [`BELT_QUICK_REFERENCE.md`](BELT_QUICK_REFERENCE.md) for a one-page overview
2. Review [`BELT_PIPELINE_OVERVIEW.md`](BELT_PIPELINE_OVERVIEW.md) for complete details

**Looking for validation results?**
- See [`reports/EVALUATION_SUMMARY.md`](reports/EVALUATION_SUMMARY.md)

**Setting up tests?**
- See [`test-documentation/`](test-documentation/)

---

## 📖 Document Descriptions

### Main Documentation

#### `BELT_PIPELINE_OVERVIEW.md`
**Complete architecture guide** covering:
- Full pipeline architecture with diagrams
- All 8 pipeline steps explained in detail
- BELT software components (extract, umi-map, umi-group, consensus)
- Configuration guide and parameters
- Data flow mapping
- Performance optimization
- Troubleshooting guide
- External dependencies
- Future development roadmap

**Length:** ~1000 lines  
**Audience:** Developers, advanced users, bioinformaticians  
**When to read:** For comprehensive understanding of the pipeline

#### `BELT_QUICK_REFERENCE.md`
**One-page quick reference** with:
- Quick start commands
- Key parameters (threshold, threads)
- Common workflows
- Troubleshooting tips
- Performance recommendations

**Length:** 1 page  
**Audience:** All users  
**When to read:** Need quick answers, running the pipeline

---

### Validation Reports

#### `reports/PIPELINE_VALIDATION_REPORT.md`
**Full pipeline validation** from initial development:
- Complete end-to-end test with threshold=5
- 237 UMI groups processed successfully
- 939 KB output generated
- Performance metrics and timing
- Quality validation checks

**Test Date:** October 24, 2025  
**Threshold:** 5 reads (comprehensive coverage)  
**Status:** ✅ Passed - 237/237 consensus sequences

#### `reports/EVALUATION_SUMMARY.md`
**Filtered pipeline validation summary:**
- Test with threshold=1000 (high-confidence mode)
- 5 UMI groups processed
- 20 KB output generated
- 32x speedup compared to threshold=5
- Quality assessment and recommendations

**Test Date:** October 25, 2025  
**Threshold:** 1000 reads (clinical/diagnostic mode)  
**Status:** ✅ Passed - 5/5 consensus sequences

#### `reports/FILTERED_PIPELINE_EVALUATION.md`
**Detailed threshold analysis:**
- Filtering impact analysis
- Sequence validation
- Performance comparison (threshold=5 vs 1000)
- Threshold selection guide for different use cases
- Sample consensus sequences
- Biological validity checks

**Length:** 235 lines  
**Focus:** Threshold optimization

#### `reports/SNAKEMAKE_UPDATE_SUMMARY.md`
**Snakemake pipeline updates:**
- Fixes to consensus output paths
- Compile consensus rule updates
- Demo summary syntax fixes
- Testing results with updated paths

**Date:** October 24, 2025  
**Status:** All updates validated

---

### Test Documentation

#### `test-documentation/CONSENSUS_TESTS_README.md`
**Consensus module testing guide:**
- Test suite structure (16 tests)
- Test coverage breakdown
- Running tests (commands and examples)
- Expected output
- Test validation results (15/15 passed)

**Test Suite:** `tests/test_consensus.py`  
**Status:** ✅ All tests passing

#### `test-documentation/PARTIAL_UMI_HANDLING.md`
**Partial UMI implementation documentation:**
- Problem description (missing UMIs in real data)
- Solution architecture (BK_scan refactoring)
- Test coverage for partial UMIs
- Real data validation
- Error handling strategies

**Feature:** Graceful degradation for missing UMIs  
**Status:** ✅ Implemented and tested

#### `test-documentation/UMI_GROUP_TESTS_README.md`
**UMI grouping tests documentation:**
- Dynamic UMI column support
- Fast BAM-based grouping tests
- FASTQ output validation
- Multi-UMI testing

**Test Suite:** `tests/test_umi_group.py`  
**Status:** ✅ All tests passing

---

## 🔑 Key Concepts

### Pipeline Steps
1. **Quality Filtering** - Filter reads by length and quality
2. **Read Mapping** - Align reads to reference
3. **BAM Processing** - Sort and index alignments
4. **UMI Extraction** - Extract and correct UMIs (`belt umi-map`)
5. **UMI Grouping** - Group reads by UMI (`belt umi-group`)
6. **Consensus Calling** - Generate consensus sequences (`belt consensus`)
7. **Compile** - Concatenate all consensus FASTA files

### BELT Software Components
- **`belt extract`** - Basic UMI extraction
- **`belt umi-map`** - UMI extraction + whitelist correction 🔧
- **`belt umi-group`** - Read grouping by UMI 🔧
- **`belt consensus`** - Consensus sequence generation 🔧 ⭐

🔧 = Used in main pipeline  
⭐ = Recently updated to Python-native implementation

### Threshold Parameter
The most important parameter for quality control:
```python
# Snakefile line 212
threshold = 1000  # Minimum reads per UMI group
```

**Guidelines:**
- `≥5000` - Ultra-high confidence (1-2 groups)
- `≥1000` - Clinical/diagnostic (5-10 groups) ⭐ **Recommended**
- `≥100` - Research balanced (50-100 groups)
- `≥5` - Maximum coverage (200-300 groups)

---

## 📊 Validation Summary

### Python-Native Consensus
- ✅ **15/15 unit tests passed**
- ✅ **100% success rate** with real data
- ✅ **Identical performance** to bash version
- ✅ **Better error handling** and logging
- ✅ **Production ready**

### Pipeline Performance
| Metric | Threshold=5 | Threshold=1000 | Improvement |
|--------|-------------|----------------|-------------|
| UMI groups | 237 | 5 | 47x fewer |
| Runtime | 32 min | 1 min | 32x faster |
| Output size | 939 KB | 20 KB | 47x smaller |

---

## 🛠 Common Tasks

### Read Documentation
```bash
# Quick reference
cat docs/BELT_QUICK_REFERENCE.md

# Full overview
less docs/BELT_PIPELINE_OVERVIEW.md

# Validation results
cat docs/reports/EVALUATION_SUMMARY.md
```

### Run Pipeline
```bash
cd "Belt Pipeline"
conda activate belt-pipeline-medaka
snakemake --cores 8
```

### Run Tests
```bash
cd belt
pytest tests/test_consensus.py -v      # Consensus tests
pytest tests/test_umi_map.py -v        # UMI mapping tests
pytest tests/test_umi_group.py -v      # UMI grouping tests
```

---

## 📝 Contributing

When adding new documentation:
1. Place general docs in `docs/`
2. Place validation reports in `docs/reports/`
3. Place test documentation in `docs/test-documentation/`
4. Update this README with links and descriptions

---

## 🔗 Quick Links

### Essential Reading
- [Quick Reference](BELT_QUICK_REFERENCE.md) - Start here
- [Pipeline Overview](BELT_PIPELINE_OVERVIEW.md) - Complete guide
- [Evaluation Summary](reports/EVALUATION_SUMMARY.md) - Validation results

### Advanced Topics
- [Threshold Analysis](reports/FILTERED_PIPELINE_EVALUATION.md)
- [Partial UMI Handling](test-documentation/PARTIAL_UMI_HANDLING.md)
- [Consensus Testing](test-documentation/CONSENSUS_TESTS_README.md)

### Source Code
- [`../scripts/belt.py`](../scripts/belt.py) - CLI entry point
- [`../scripts/consensus.py`](../scripts/consensus.py) - Consensus module
- [`../Belt Pipeline/Snakefile`](../Belt%20Pipeline/Snakefile) - Pipeline definition

---

**Documentation Version:** 2.0  
**Last Updated:** October 25, 2025  
**Status:** ✅ Current and Comprehensive

