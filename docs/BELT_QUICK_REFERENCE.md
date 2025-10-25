# BELT Pipeline - Quick Reference Guide

## One-Page Overview

### Pipeline Steps
```
1. Filter   → Quality & length filtering (fastplong)
2. Map      → Align to reference (minimap2)
3. Sort     → Sort & index BAM (samtools)
4. Extract  → UMI extraction + whitelist correction (belt umi-map) 🔧
5. Group    → Group reads by UMI (belt umi-group) 🔧
6. Consensus→ Generate consensus sequences (belt consensus) 🔧
7. Compile  → Concatenate all consensus FASTA
```

### Quick Start
```bash
# 1. Setup
cd "Belt Pipeline"
conda activate belt-pipeline-medaka

# 2. Edit config
vim config.yaml

# 3. Run
snakemake --cores 8
```

### Key Parameters (Snakefile)

**Line 212 - UMI Threshold** ⭐ MOST IMPORTANT
```python
threshold = 1000  # Minimum reads per UMI group
```

| Threshold | Groups | Time | Use Case |
|-----------|--------|------|----------|
| ≥1000 | 5-10 | ~1 min | Clinical ⭐ |
| ≥100 | 50-100 | ~10 min | Research |
| ≥5 | 200-300 | ~30 min | Comprehensive |

**Line 238 - Parallel Threads**
```python
threads: 8  # More cores = faster consensus
```

### Output Structure
```
data/
├── 5-umi_maps/                    # UMI extraction results
│   └── {sample}_Paired_UMIs_read_names.csv
├── 6-umi_groups/                  # Grouped BAM files
│   └── {sample}.grouped.bam
├── 7-consensus_calls/             # Individual consensus
│   └── {sample}/{sample}/consensus/
│       └── {UMI}_consensus.fasta
└── 8-compiled_consensus/          # Final output ⭐
    └── {sample}-compiled-consensus.fa
```

### BELT Commands

#### Extract UMIs with Whitelist
```bash
belt umi-map \
  --alignment sorted.bam \
  --reference reference.fa \
  --motif "['NNNNNNNNNNNN', 'NNNNNNNNNNNN']" \
  --whitelist whitelist.csv \
  --output_dir output/
```

#### Group Reads by UMI
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

### Troubleshooting

#### Pipeline Fails
```bash
# Check logs
cat logs/consensus/*.log
cat .snakemake/log/*.snakemake.log

# Unlock if interrupted
snakemake --unlock

# Dry run to validate
snakemake --cores 8 -n
```

#### No Consensus Generated
```bash
# Check threshold - may be too high
# Line 212 in Snakefile
threshold = 5  # Lower threshold

# Check if UMI groups exist
samtools view -H data/6-umi_groups/*.bam | grep @PG
```

#### Missing Tools
```bash
# Activate conda environment
conda activate belt-pipeline-medaka

# Verify tools
which minimap2 samtools medaka_variant bcftools
```

### Testing

#### Run Unit Tests
```bash
cd belt
pytest tests/test_consensus.py -v      # 15 tests
pytest tests/test_extract.py -v
pytest tests/test_umi_map.py -v
pytest tests/test_umi_group.py -v
```

#### Validate Pipeline
```bash
cd "Belt Pipeline"
snakemake --cores 8
# Check demo_summary.txt for results
```

### Performance Tips

1. **Adjust threshold** - biggest impact on speed
2. **Use more threads** - line 238 in Snakefile
3. **Pre-filter reads** - remove low-quality early
4. **Use SSD** - faster I/O for temp files

### Key Features

✅ **Multi-motif UMI support** - any number of UMI positions  
✅ **Partial UMI handling** - gracefully handles missing UMIs  
✅ **Python-native consensus** - robust, tested, production-ready  
✅ **Flexible threshold** - tune speed vs coverage  
✅ **Parallel processing** - efficient multi-core utilization  

### Common Workflows

#### High-Confidence Clinical
```python
# Snakefile line 212
threshold = 1000
# Fast processing, highest confidence
```

#### Balanced Research
```python
# Snakefile line 212
threshold = 100
# Good balance of coverage and quality
```

#### Comprehensive Discovery
```python
# Snakefile line 212
threshold = 5
# Maximum UMI recovery
```

### Documentation

| File | Purpose |
|------|---------|
| `BELT_PIPELINE_OVERVIEW.md` | Complete architecture guide |
| `BELT_QUICK_REFERENCE.md` | This quick reference |
| `EVALUATION_SUMMARY.md` | Validation results |
| `tests/CONSENSUS_TESTS_README.md` | Test documentation |

### Getting Help

1. Check documentation in `docs/` and `tests/`
2. Review log files in `logs/`
3. Run tests to validate installation
4. Check Snakemake logs for detailed errors

---

**Version:** 2.0 (Python-Native Consensus)  
**Status:** ✅ Production Ready  
**Last Updated:** October 25, 2025

