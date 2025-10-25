# Output Structure Update - Option 4

**Date:** October 25, 2025  
**Change:** Simplified consensus output directory structure  
**Status:** ✅ Implemented and Tested

---

## Summary

Updated the consensus calling output structure from a messy nested format to a clean, type-based organization (Option 4).

---

## Old Structure (Before) ❌

```
data/7-consensus_calls/
└── SAMPLE/
    └── SAMPLE/                          ❌ Redundant nesting
        ├── consensus/
        │   └── UMI_consensus.fasta
        ├── tmp/                         ❌ Files hard to find
        │   └── UMI/
        │       └── meta/
        │           └── aligned.bam
        └── bam/                         ❌ Empty/unused
```

**Problems:**
- Redundant sample name nesting
- BAM files buried in tmp/meta directories
- Inconsistent naming (UMI_consensus.fasta vs aligned.bam)
- tmp directory mixed with final outputs

---

## New Structure (Option 4) ✅

```
data/7-consensus_calls/
└── SAMPLE/
    ├── bam/                  📂 Aligned BAM files
    │   ├── UMI1.bam
    │   ├── UMI1.bam.bai
    │   ├── UMI2.bam
    │   ├── UMI2.bam.bai
    │   └── ...
    │
    ├── fasta/                📂 Consensus sequences
    │   ├── UMI1.fasta
    │   ├── UMI2.fasta
    │   └── ...
    │
    └── tmp/                  📂 Temporary files (optional cleanup)
        └── UMI/
            ├── UMI_reads.fastq
            └── meta/
                ├── medaka.vcf
                └── medaka.norm.vcf.gz
```

**Benefits:**
- ✅ Clean, flat organization
- ✅ Files grouped by type (bam/ vs fasta/)
- ✅ Easy to find UMI pairs
- ✅ Consistent naming scheme
- ✅ Simple to add more output types later

---

## Real Example

### Sample: LP249-NL43-A3G-1-A (5 UMI groups, threshold=1000)

```
data/7-consensus_calls/LP249-NL43-A3G-1-A/
├── bam/                              Total: 5 files (17 MB)
│   ├── GCGTAAGTCTAA.bam              3.7 MB (2,200 reads)
│   ├── GCGTAAGTCTAA.bam.bai          24 KB
│   ├── TAGTGGTTACTC.bam              4.3 MB
│   ├── TAGTGGTTACTC.bam.bai          24 KB
│   ├── TATTTGGACCTA.bam              1.4 MB
│   ├── TATTTGGACCTA.bam.bai          24 KB
│   ├── TGTATTTGTCAT.bam              3.5 MB
│   ├── TGTATTTGTCAT.bam.bai          24 KB
│   ├── TTGCGCTGGGGC.bam              3.8 MB
│   └── TTGCGCTGGGGC.bam.bai          24 KB
│
├── fasta/                            Total: 5 files (20 KB)
│   ├── GCGTAAGTCTAA.fasta            4.0 KB
│   ├── TAGTGGTTACTC.fasta            3.9 KB
│   ├── TATTTGGACCTA.fasta            3.7 KB
│   ├── TGTATTTGTCAT.fasta            4.0 KB
│   ├── TTGCGCTGGGGC.fasta            4.0 KB
│   └── runner.log                    Processing log
│
└── tmp/                              Intermediate files
    └── {UMI}/
        ├── {UMI}_reads.fastq         Temporary (deleted after processing)
        └── meta/
            ├── aligned.bam           Intermediate alignment
            ├── medaka.vcf            Raw variant calls
            └── medaka.norm.vcf.gz    Normalized variants
```

---

## File Naming Convention

### BAM Files
```
{UMI}.bam       - Aligned reads for this UMI
{UMI}.bam.bai   - BAM index for random access
```

**Example:** `GCGTAAGTCTAA.bam`

### FASTA Files
```
{UMI}.fasta     - Consensus sequence for this UMI
```

**Example:** `GCGTAAGTCTAA.fasta`

**Header Format:**
```fasta
>GCGTAAGTCTAA
CTACGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAACACATGGAAAAGAT
...
```

---

## Implementation Details

### Files Modified

1. **`scripts/consensus.py`**
   - Changed output directories: `consensus` → `fasta`
   - Added Step 6: Copy aligned BAM to final location
   - Updated file naming: `{UMI}_consensus.fasta` → `{UMI}.fasta`
   - BAM saved to `bam/{UMI}.bam` instead of `tmp/{UMI}/meta/aligned.bam`

2. **`Belt Pipeline/Snakefile`**
   - Updated `consensus_call` rule:
     - Old: `output: directory("data/7-consensus_calls/{sample}/{sample}/consensus")`
     - New: `output: fasta_dir, bam_dir`
   - Updated `compile_consensus` rule:
     - Old glob: `{sample}/{sample}/consensus/*_consensus.fasta`
     - New glob: `{sample}/fasta/*.fasta`

### Code Changes

#### consensus.py (Line 260-264)
```python
# OLD
con_dir = out_root / "consensus"
bam_dir = out_root / "bam"

# NEW
con_dir = out_root / "fasta"          # Consensus FASTA files
bam_dir = out_root / "bam"            # Aligned BAM files
```

#### consensus.py (Line 68-73)
```python
# OLD
align_bam = meta_dir / "aligned.bam"
consensus_fasta = con_dir / f"{umi_safe}_consensus.fasta"

# NEW
align_bam = meta_dir / "aligned.bam"                    # Temporary
final_bam = bam_dir / f"{umi_safe}.bam"                # Final: bam/UMI.bam
consensus_fasta = con_dir / f"{umi_safe}.fasta"        # Final: fasta/UMI.fasta
```

#### consensus.py (Line 199-210) - New Step 6
```python
# Step 6: Copy aligned BAM to final location
if logfh is not subprocess.DEVNULL:
    logfh.write(f"✅ Step 6: Saving aligned BAM to {final_bam}\n")
    logfh.flush()

import shutil
shutil.copy2(str(align_bam), str(final_bam))
shutil.copy2(str(align_bam) + ".bai", str(final_bam) + ".bai")
```

---

## Testing

### Test Run
- **Sample:** LP249-NL43-A3G-1-A
- **Threshold:** 1000 reads
- **UMI Groups:** 5
- **Status:** ✅ All passed

### Validation
```bash
# BAM files created
✅ 5 BAM files in bam/
✅ 5 BAM indices (.bai)
✅ All BAM files readable with samtools

# FASTA files created  
✅ 5 FASTA files in fasta/
✅ All FASTA files valid format
✅ Headers match UMI identifiers

# Compiled output
✅ Combined FASTA created successfully
✅ All 5 sequences included
```

---

## Usage

### Access BAM Files
```bash
# List all BAM files for a sample
ls data/7-consensus_calls/SAMPLE/bam/*.bam

# View specific UMI alignment
samtools view data/7-consensus_calls/SAMPLE/bam/UMI.bam | less

# Get read count
samtools view -c data/7-consensus_calls/SAMPLE/bam/UMI.bam
```

### Access FASTA Files
```bash
# List all consensus sequences
ls data/7-consensus_calls/SAMPLE/fasta/*.fasta

# View specific UMI consensus
cat data/7-consensus_calls/SAMPLE/fasta/UMI.fasta

# Get sequence length
grep -v "^>" data/7-consensus_calls/SAMPLE/fasta/UMI.fasta | tr -d '\n' | wc -c
```

### Downstream Analysis
```bash
# Compare UMI to reference
minimap2 -ax map-ont reference.fa data/7-consensus_calls/SAMPLE/fasta/UMI.fasta

# Call variants
bcftools mpileup -f reference.fa data/7-consensus_calls/SAMPLE/bam/UMI.bam | \
  bcftools call -mv -Oz -o variants.vcf.gz

# Extract coverage
samtools depth data/7-consensus_calls/SAMPLE/bam/UMI.bam > coverage.txt
```

---

## Cleanup

### Remove Temporary Files
```bash
# Safe to delete tmp/ after consensus generation
rm -rf data/7-consensus_calls/SAMPLE/tmp/

# Keep bam/ and fasta/ directories - these are your final outputs!
```

### Archive Old Runs
```bash
# Archive a complete run
tar -czf SAMPLE_consensus_$(date +%Y%m%d).tar.gz \
  data/7-consensus_calls/SAMPLE/bam/ \
  data/7-consensus_calls/SAMPLE/fasta/
```

---

## Migration from Old Structure

If you have data in the old structure:

```bash
# For each sample, reorganize files
SAMPLE="your_sample_name"
OLD="data/7-consensus_calls/$SAMPLE/$SAMPLE"
NEW="data/7-consensus_calls/$SAMPLE"

# Move FASTA files
mkdir -p $NEW/fasta
mv $OLD/consensus/*_consensus.fasta $NEW/fasta/
rename 's/_consensus\.fasta/.fasta/' $NEW/fasta/*.fasta

# Move BAM files if needed
mkdir -p $NEW/bam
find $OLD/tmp -name "aligned.bam*" -exec mv {} $NEW/bam/ \;
rename 's/aligned\.bam/APPROPRIATE_UMI.bam/' $NEW/bam/*.bam

# Remove old nested directory
rm -rf $OLD
```

---

## Future Enhancements

The new structure makes it easy to add more outputs per UMI:

### Potential Additions
```
data/7-consensus_calls/SAMPLE/
├── bam/
├── fasta/
├── vcf/              🆕 Variant calls per UMI
│   ├── UMI1.vcf.gz
│   └── ...
├── coverage/         🆕 Coverage statistics
│   ├── UMI1.depth
│   └── ...
└── qc/              🆕 Quality metrics
    ├── UMI1_stats.txt
    └── ...
```

---

## Summary

✅ **Implemented Option 4** - Clean type-based structure  
✅ **Tested with real data** - 5 UMI groups processed successfully  
✅ **Documentation updated** - All changes documented  
✅ **Backward compatible** - Old data can be migrated  
✅ **Future-proof** - Easy to extend with new output types  

The new structure provides a professional, maintainable organization that keeps related files together while maintaining clarity through type-based separation.

---

**Implementation Date:** October 25, 2025  
**Tested By:** BELT Pipeline v2.0  
**Status:** ✅ Production Ready

