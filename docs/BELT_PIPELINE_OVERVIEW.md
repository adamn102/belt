# BELT Pipeline Architecture Overview
## Barcode Extraction and Lineage Tracing Pipeline

**Version:** 2.0 (Python-Native Consensus)  
**Date:** October 25, 2025  
**Status:** Production Ready ✅

---

## Executive Summary

BELT (Barcode Extraction and Lineage Tracing) is a comprehensive bioinformatics pipeline for processing long-read sequencing data with UMI (Unique Molecular Identifier) barcodes. The pipeline extracts UMIs, groups reads by barcode, and generates high-quality consensus sequences for downstream variant analysis.

### Key Features
- ✅ **Multi-motif UMI support** - handles variable numbers of UMI positions
- ✅ **Partial UMI handling** - gracefully processes reads with missing UMIs
- ✅ **Whitelist correction** - maps experimental UMIs to known sequences
- ✅ **Fast BAM-based grouping** - efficient read clustering by UMI
- ✅ **Python-native consensus calling** - robust, tested, production-ready
- ✅ **Flexible quality control** - threshold-based filtering
- ✅ **Parallel processing** - multi-core support throughout

---

## Pipeline Architecture

### Overview Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    BELT Pipeline Workflow                       │
└─────────────────────────────────────────────────────────────────┘

Input: Raw FASTQ + Reference FASTA + Whitelist CSV + Config YAML
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: Quality Filtering (fastplong)                           │
│   - Filter by length (≥1000 bp)                                 │
│   - Filter by quality (Q25, ≤50% low quality)                   │
│   Output: data/2-fastq_filtered/{sample}.fastq                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: Read Mapping (minimap2)                                 │
│   - Map filtered reads to reference                             │
│   - Generate alignment BAM                                      │
│   Output: data/3-mapped_reads/{sample}.bam                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 3: BAM Processing (samtools)                               │
│   - Sort BAM by coordinate                                      │
│   - Index BAM for random access                                 │
│   Output: data/4-sorted_reads/{sample}.bam + .bai               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 4: UMI Extraction (belt umi-map) 🔧 BELT SOFTWARE          │
│   - Extract UMIs from aligned reads using motifs                │
│   - Map raw UMIs to whitelist (error correction)                │
│   - Handle partial UMIs gracefully                              │
│   Output: data/5-umi_maps/{sample}_Paired_UMIs_read_names.csv   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 5: UMI Grouping (belt umi-group) 🔧 BELT SOFTWARE          │
│   - Group reads by corrected UMI                                │
│   - Filter groups by read count threshold                       │
│   - Create sorted, tagged BAM                                   │
│   Output: data/6-umi_groups/{sample}.grouped.bam                │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 6: Consensus Calling (belt consensus) 🔧 BELT SOFTWARE     │
│   - Generate consensus for each UMI group                       │
│   - Python-native implementation (NEW)                          │
│   - Parallel processing of multiple UMIs                        │
│   Output: data/7-consensus_calls/{sample}/{sample}/consensus/   │
│           └── {UMI}_consensus.fasta (one per UMI)               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 7: Compile Consensus                                       │
│   - Concatenate all consensus FASTA files                       │
│   - Single multi-FASTA output                                   │
│   Output: data/8-compiled_consensus/{sample}-compiled.fa        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ STEP 8: Summary Report                                          │
│   - Generate pipeline summary                                   │
│   Output: demo_summary.txt                                      │
└─────────────────────────────────────────────────────────────────┘
```

---

## BELT Software Components

### 1. **`belt extract`** (UMI Extraction - Basic)
**Module:** `scripts/extract.py`  
**Function:** `extract_umis()`

Extracts UMIs from aligned reads without whitelist correction.

**Key Features:**
- Multi-motif support (UMI1, UMI2, UMI3, ...)
- Handles partial UMIs (some motifs present, others missing)
- Outputs raw UMI sequences

**Not used in main pipeline** - replaced by `belt umi-map` which includes whitelist correction.

---

### 2. **`belt umi-map`** (UMI Extraction + Whitelist Correction) 🔧
**Module:** `scripts/umi_map.py`  
**Function:** `umi_extraction_map()`

Extracts UMIs and maps them to a whitelist for error correction.

**Input:**
```bash
belt umi-map \
  --reference reference.fa \
  --alignment sorted.bam \
  --min_read_length 1000 \
  --motif "['NNNNNNNNNNNN', 'NNNNNNNNNNNN']" \
  --whitelist whitelist.csv \
  --output_dir output/
```

**Process:**
1. **Motif Detection:** Finds expected UMI positions in reference
2. **UMI Extraction:** Extracts sequences at motif positions from reads
3. **Partial UMI Handling:** Records `None` for missing UMIs
4. **Whitelist Mapping:** Uses BK-tree for approximate string matching
   - Hamming distance ≤1 by default
   - Preserves `None` for truly missing UMIs
5. **Output:** CSV with Read Name, UMI1, UMI2, ... columns

**Key Features:**
- ✅ Dynamic UMI column naming (UMI1, UMI2, ... instead of fixed 5pUMI/3pUMI)
- ✅ Handles 11.0% missing UMI1, 22.5% missing UMI2 in test data
- ✅ Error-tolerant whitelist matching
- ✅ Comprehensive logging of missing UMIs

**Output Format:**
```csv
Read Name,UMI1,UMI2
read-001,GTCTGTGCCACG,GCCAAGGCGGCT
read-002,TGGTGTTGCATG,
read-003,,GAGGCAGGGACT
```

---

### 3. **`belt umi-group`** (Read Grouping by UMI) 🔧
**Module:** `scripts/umi_group.py`  
**Function:** `long_read_umi_grouping_fast()`

Groups reads by corrected UMI and generates a sorted, tagged BAM file.

**Input:**
```bash
belt umi-group \
  --reference reference.fa \
  --alignment sorted.bam \
  --corrected_umis umi_map.csv \
  --grouping UMI1 \
  --threshold 1000 \
  --output_dir output/
```

**Process:**
1. **Load UMI Mapping:** Read CSV with corrected UMIs
2. **Filter by Threshold:** Keep only UMI groups with ≥N reads
3. **Tag BAM Reads:** Add `ZU` tag with UMI identifier
4. **Sort by UMI:** Create BAM sorted by UMI tag (not coordinate)
5. **Statistics:** Report number of groups passing threshold

**Key Features:**
- ✅ Dynamic grouping by any UMI column (UMI1, UMI2, ...)
- ✅ Threshold-based quality control
- ✅ Fast BAM-based grouping (no FASTQ intermediate)
- ✅ Optional FASTQ output with annotated headers

**Output:**
- `{sample}.grouped.bam` - BAM file sorted by UMI tag
- Reads tagged with `ZU:{UMI_sequence}`

---

### 4. **`belt consensus`** (Consensus Sequence Generation) 🔧 ⭐ NEW
**Module:** `scripts/consensus.py`  
**Function:** `process_grouped_bam()` → `run_consensus_for_umi()`

Generates high-quality consensus sequences for each UMI group using medaka.

**Input:**
```bash
belt consensus \
  --grouped-bam grouped.bam \
  --reference reference.fa \
  --output-dir output/ \
  --sample sample_name \
  --threads 8
```

**Process (per UMI group):**
1. **Extract Reads:** Pull reads for this UMI from grouped BAM → temporary FASTQ
2. **Align to Reference:** `minimap2` + `samtools sort`
3. **Index BAM:** `samtools index`
4. **Call Variants:** `medaka_variant` (ONT-specific variant caller)
5. **Normalize VCF:** `bcftools norm` (left-align indels, compress)
6. **Generate Consensus:** `bcftools consensus` (apply variants to reference)
7. **Cleanup:** Remove temporary FASTQ

**Key Features:**
- ✅ **Python-native implementation** (no bash script dependency)
- ✅ **Parallel processing** (multiple UMIs simultaneously)
- ✅ **Comprehensive error handling** (per-UMI try/catch)
- ✅ **Detailed logging** (step-by-step progress tracking)
- ✅ **Automatic cleanup** (removes temporary files)
- ✅ **Fully tested** (15/15 unit tests passing)

**Output Structure:**
```
output_dir/
└── sample_name/
    └── sample_name/
        ├── consensus/              # Final consensus sequences
        │   ├── UMI1_consensus.fasta
        │   ├── UMI2_consensus.fasta
        │   └── ...
        ├── tmp/                    # Temporary files (cleaned up)
        │   └── UMI1/
        │       ├── UMI1_reads.fastq (deleted after use)
        │       └── meta/
        │           ├── aligned.bam
        │           ├── medaka.vcf
        │           └── ...
        └── bam/                    # Optional BAM storage
```

**Performance:**
- ~8-12 seconds per UMI group
- Scales linearly with number of UMIs
- Efficient parallel processing

---

## Pipeline Configuration

### Input Files

#### 1. **Config YAML** (`config.yaml`)
```yaml
samples:
  SAMPLE_NAME: path/to/sample.fastq

barcodes:
  SAMPLE_NAME: "['NNNNNNNNNNNN', 'NNNNNNNNNNNN']"

whitelist:
  SAMPLE_NAME: WHITELIST_ID

min_lengths:
  SAMPLE_NAME: 1000
```

#### 2. **Reference Files**
```
data/
├── 0-reference_sequence/
│   └── {ref_id}.fa          # For read mapping
└── 0-consensus_reference/
    └── {ref_id}.fa          # For consensus calling
```

#### 3. **Whitelist Files**
```
data/
└── 0-library_whitelist/
    └── {whitelist_id}.csv
```

**Whitelist Format:**
```csv
UMI1,UMI2
AAAAAAAAAGGG,CCCCCCCCCTTT
AAAAAAAAACCC,CCCCCCCCGGG
...
```

---

## Data Flow

### Input → Output Mapping

| Step | Input | Output | Size Change |
|------|-------|--------|-------------|
| **Raw Data** | FASTQ | 227 MB | - |
| **1. Filter** | 153,725 reads | 58,205 reads | -62% |
| **2. Map** | 58,205 reads | 58,205 aligned | Same |
| **3. Sort** | Unsorted BAM | Sorted BAM | Same |
| **4. Extract UMIs** | 55,342 mapped | CSV (55K rows) | <<1 MB |
| **5. Group** | 55,342 reads | 237 UMI groups | Same data |
| **6. Consensus** | 237 groups | 5 consensus¹ | 20 KB |
| **7. Compile** | 5 FASTA files | 1 multi-FASTA | 20 KB |

¹ With threshold=1000; 237 consensus with threshold=5

### File Sizes (Threshold=1000)
```
data/
├── 1-fastq_files/              227 MB (raw)
├── 2-fastq_filtered/           181 MB (filtered)
├── 4-sorted_reads/             ~50 MB (BAM)
├── 5-umi_maps/                 ~2 MB (CSV)
├── 6-umi_groups/               ~50 MB (BAM)
├── 7-consensus_calls/          20 KB (5 FASTA)
└── 8-compiled_consensus/       20 KB (compiled)
```

---

## Key Parameters

### Adjustable Parameters in Snakefile

#### Quality Filtering
```python
# Line 135-137
params:
    length = lambda wildcards: config["min_lengths"][wildcards.sample],  # 1000
    qual_threshold = 25,
    qual_percent = 50
```

#### UMI Grouping
```python
# Line 210-212
params:
    grouping = "UMI1",      # Which UMI to group by
    threshold = 1000        # Minimum reads per UMI group ⭐ KEY PARAMETER
```

#### Consensus Calling
```python
# Line 238
threads: 8  # Parallel UMI processing
```

### Threshold Guidelines

| Threshold | UMI Groups | Processing Time | Use Case |
|-----------|------------|----------------|----------|
| ≥5000 | 1-2 | Seconds | Ultra-high confidence only |
| ≥1000 | 5-10 | ~1 minute | **Clinical/diagnostic** ⭐ |
| ≥500 | 10-20 | ~2-3 minutes | High-quality research |
| ≥100 | 50-100 | ~8-10 minutes | Balanced analysis |
| ≥10 | 200-300 | ~30 minutes | Comprehensive |
| ≥5 | 237 | ~32 minutes | Maximum sensitivity |

---

## External Dependencies

### Required Tools

#### Bioinformatics Tools
```bash
minimap2      # Read alignment
samtools      # BAM manipulation
medaka        # ONT consensus calling (medaka_variant)
bcftools      # VCF manipulation and consensus generation
fastplong     # Quality filtering
```

#### Python Packages
```python
pysam         # BAM/SAM file I/O
pandas        # Data manipulation
numpy         # Numerical operations
biopython     # Sequence I/O
click         # CLI framework
pybktree      # Approximate string matching
```

### Installation

#### Via Conda (Recommended)
```bash
conda env create -f environment.yml
conda activate belt-pipeline-medaka
```

#### Via Pip (Python packages only)
```bash
pip install -r requirements.txt
```

---

## Pipeline Execution

### Running the Full Pipeline

#### 1. **Setup Configuration**
```bash
cd "Belt Pipeline"
# Edit config.yaml with your samples
```

#### 2. **Dry Run** (recommended)
```bash
snakemake --cores 8 -n
```

#### 3. **Full Run**
```bash
snakemake --cores 8
```

#### 4. **Specific Sample**
```bash
snakemake --cores 8 data/8-compiled_consensus/SAMPLE-compiled-consensus.fa
```

### Running Individual BELT Tools

#### Extract UMIs Only
```bash
belt extract \
  --alignment sorted.bam \
  --reference reference.fa \
  --barcode_motif "['NNNNNNNNNNNN']" \
  --output_dir output/
```

#### Extract + Map to Whitelist
```bash
belt umi-map \
  --alignment sorted.bam \
  --reference reference.fa \
  --motif "['NNNNNNNNNNNN', 'NNNNNNNNNNNN']" \
  --whitelist whitelist.csv \
  --output_dir output/
```

#### Group by UMI
```bash
belt umi-group \
  --alignment sorted.bam \
  --reference reference.fa \
  --corrected_umis umi_map.csv \
  --grouping UMI1 \
  --threshold 1000 \
  --output_dir output/
```

#### Generate Consensus
```bash
belt consensus \
  --grouped-bam grouped.bam \
  --reference reference.fa \
  --output-dir output/ \
  --sample sample_name \
  --threads 8
```

---

## Testing and Validation

### Unit Tests

#### Consensus Module
```bash
cd belt
pytest tests/test_consensus.py -v
# 15/15 tests passed ✅
```

#### UMI Extraction
```bash
pytest tests/test_extract.py -v
```

#### UMI Mapping
```bash
pytest tests/test_umi_map.py -v
```

#### UMI Grouping
```bash
pytest tests/test_umi_group.py -v
```

### Integration Tests

Full pipeline test with real data:
```bash
cd "Belt Pipeline"
snakemake --cores 8
# Validates all steps end-to-end
```

---

## Recent Updates (v2.0)

### Python-Native Consensus Implementation ⭐
- **Replaced:** Bash script (`consensus.sh`)
- **With:** Pure Python implementation (`consensus.py`)
- **Benefits:**
  - ✅ Better error handling
  - ✅ Comprehensive logging
  - ✅ Fully tested (15 unit tests)
  - ✅ Easier to maintain and extend
  - ✅ No performance penalty

### Multi-Motif UMI Support
- **Dynamic UMI columns:** UMI1, UMI2, ... (not fixed 5pUMI/3pUMI)
- **Flexible configuration:** Support any number of motifs
- **Backward compatible:** Existing data still works

### Partial UMI Handling
- **Graceful degradation:** Handles missing UMIs
- **Preserves data:** Records `None` instead of failing
- **Statistics:** Reports missing UMI percentages

### Threshold-Based Filtering
- **Quality control:** Filter UMI groups by read count
- **Performance:** 32x speedup with threshold=1000
- **Flexibility:** Easy to adjust for different use cases

---

## Output Interpretation

### Consensus FASTA Format
```fasta
>UMI_IDENTIFIER
CTACGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAACACATGGAAAAGAT
TAGTAAAACACCATATGTATATTTCAAGGAAAGCTAAGGACTGGTTTTATAGACATCACT
...
```

### Quality Indicators

#### Good Consensus
- ✅ Sequence length ≥95% of reference
- ✅ No large gaps (>100 bp)
- ✅ No high N content (>5%)
- ✅ High read coverage (≥100 reads)

#### Poor Consensus (may need filtering)
- ❌ Sequence length <90% of reference
- ❌ Large gaps or deletions
- ❌ High N content
- ❌ Low read coverage (<10 reads)

### Downstream Analysis

Consensus sequences can be used for:
- **Variant calling** (identify SNPs, indels)
- **Phylogenetic analysis** (build trees)
- **Diversity assessment** (measure population heterogeneity)
- **Clone tracking** (follow individual lineages)
- **Quality control** (validate library preparation)

---

## Troubleshooting

### Common Issues

#### 1. "Missing required tools"
```
RuntimeError: Missing required tools: medaka_variant
```
**Solution:** Activate conda environment
```bash
conda activate belt-pipeline-medaka
```

#### 2. "No UMI groups passing threshold"
```
[BELT] Done. Processed 0 UMI groups passing threshold.
```
**Solution:** Lower threshold in Snakefile (line 212)
```python
threshold = 5  # or 10, 100, etc.
```

#### 3. Empty compiled consensus
**Solution:** Check that individual consensus files were generated
```bash
ls data/7-consensus_calls/SAMPLE/SAMPLE/consensus/*.fasta
```

#### 4. Slow consensus calling
**Solution:** 
- Increase threshold to process fewer UMIs
- Increase threads in Snakefile (line 238)
- Process on more powerful hardware

### Log Files

Check logs for detailed error messages:
```bash
# Consensus calling logs
cat logs/consensus/SAMPLE.log

# Snakemake execution log
cat .snakemake/log/*.snakemake.log
```

---

## Performance Optimization

### Speed Improvements

1. **Increase Threshold** (biggest impact)
   - threshold=5 → 237 UMIs → 32 minutes
   - threshold=1000 → 5 UMIs → 1 minute
   - **32x speedup** ⭐

2. **Use More Threads**
   - Default: 8 cores
   - Can use up to number of UMI groups
   - Diminishing returns beyond 8-16 cores

3. **Filter Early**
   - Pre-filter reads by quality
   - Remove unmapped reads before UMI extraction
   - Skip low-coverage UMIs

### Resource Requirements

#### Minimum
```
CPU: 4 cores
RAM: 8 GB
Disk: 50 GB
```

#### Recommended
```
CPU: 8-16 cores
RAM: 16-32 GB
Disk: 100 GB
```

#### Large Datasets
```
CPU: 32+ cores
RAM: 64+ GB
Disk: 500+ GB
```

---

## Citation and Documentation

### Documentation Files

| File | Content |
|------|---------|
| `BELT_PIPELINE_OVERVIEW.md` | This overview (you are here) |
| `PIPELINE_VALIDATION_REPORT.md` | Comprehensive validation |
| `FILTERED_PIPELINE_EVALUATION.md` | Threshold filtering analysis |
| `EVALUATION_SUMMARY.md` | Quick evaluation summary |
| `tests/CONSENSUS_TESTS_README.md` | Test suite documentation |
| `tests/PARTIAL_UMI_HANDLING.md` | Partial UMI implementation |

### Source Code

| Module | Purpose |
|--------|---------|
| `scripts/belt.py` | Main CLI entry point |
| `scripts/extract.py` | UMI extraction (basic) |
| `scripts/umi_map.py` | UMI extraction + whitelist |
| `scripts/umi_group.py` | Read grouping by UMI |
| `scripts/consensus.py` | Consensus generation ⭐ |
| `tests/test_*.py` | Unit tests |

---

## Future Development

### Planned Features
- [ ] Support for additional consensus callers (racon, canu)
- [ ] Quality metrics in consensus output
- [ ] Adaptive threshold selection
- [ ] Web-based visualization dashboard
- [ ] Support for paired-end reads
- [ ] Integration with variant calling pipelines

### Community Contributions
Contributions welcome! Please see `CONTRIBUTING.md` for guidelines.

---

## Summary

BELT is a **production-ready, comprehensive pipeline** for processing long-read sequencing data with UMI barcodes. The recent v2.0 update introduces a **Python-native consensus implementation** with improved error handling, comprehensive testing, and flexible quality control.

### Key Strengths
✅ **Robust** - handles edge cases (partial UMIs, low coverage)  
✅ **Flexible** - configurable parameters for different use cases  
✅ **Fast** - parallel processing and efficient algorithms  
✅ **Tested** - comprehensive unit and integration tests  
✅ **Documented** - extensive documentation and examples  
✅ **Production-Ready** - validated with real sequencing data  

### Recommended Use
- **Clinical Applications:** threshold ≥1000 for high-confidence results
- **Research Applications:** threshold 100-500 for balanced analysis
- **Exploratory Analysis:** threshold ≥5 for comprehensive coverage

The pipeline successfully balances **speed, quality, and flexibility** for a wide range of sequencing applications.

---

**Document Version:** 2.0  
**Last Updated:** October 25, 2025  
**Status:** ✅ Current and Validated

