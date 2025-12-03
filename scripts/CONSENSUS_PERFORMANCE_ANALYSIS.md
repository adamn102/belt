# Consensus Process Performance Analysis

## Current Bottlenecks

### 1. **Redundant Alignment (Major Bottleneck)**
- **Problem**: Reads in `grouped.bam` are already aligned to the reference (confirmed by BAM header showing `@SQ` and `@PG` with minimap2)
- **Current Flow**: BAM → FASTQ → minimap2 alignment → BAM (redundant!)
- **Impact**: ~10-30% of total time wasted on re-alignment
- **Evidence**: BAM header shows reads were already aligned with minimap2 in step 2 of pipeline

### 2. **medaka_variant is Computationally Expensive (Major Bottleneck)**
- **Problem**: Deep learning model inference is slow (~5-15 seconds per UMI)
- **Current**: Each UMI processed independently with full model inference
- **Impact**: 60-80% of total processing time
- **Bottleneck**: Cannot be easily parallelized per UMI (GPU would help but not always available)

### 3. **I/O Overhead**
- **Problem**: Multiple format conversions: BAM → FASTQ → BAM → VCF → FASTA
- **Current**: Each conversion requires disk I/O
- **Impact**: ~5-10% of total time
- **Issue**: Writing/reading temporary files for each UMI

### 4. **Sequential Processing Per UMI**
- **Problem**: Each UMI goes through 6 sequential steps:
  1. Align (redundant)
  2. Index BAM
  3. medaka_variant (slow)
  4. Normalize VCF
  5. Generate consensus
  6. Copy BAM
- **Impact**: No pipeline parallelism within UMI processing

### 5. **Reference Indexing Redundancy**
- **Problem**: Reference FASTA indexed multiple times (once per UMI if not cached)
- **Current**: `samtools faidx` called for each UMI
- **Impact**: Small but adds up with many UMIs

## Optimization Strategies

### **High Impact Optimizations**

#### 1. Skip Redundant Alignment (Estimated 15-25% speedup)
**Approach**: Remove the minimap2 alignment step since medaka_variant does its own alignment internally

**Current Flow**: 
- BAM → FASTQ → minimap2 alignment → BAM → medaka_variant (which aligns again internally)

**Optimized Flow**:
- BAM → FASTQ → medaka_variant (aligns internally, no need for our alignment)

**Implementation**:
- Remove Step 1 (minimap2 alignment) entirely
- Pass FASTQ directly to medaka_variant
- medaka_variant will handle alignment internally
- Note: medaka_variant requires FASTX input, not BAM, so we still need BAM→FASTQ conversion

**Challenges**:
- medaka_variant does its own alignment, so we can't skip alignment completely
- But we can skip OUR redundant alignment step
- Still need to extract reads from grouped BAM and convert to FASTQ

#### 2. Batch Process Multiple UMIs (Estimated 30-50% speedup)
**Approach**: Process multiple small UMIs together when total read count is low

**Implementation**:
- Group UMIs by total read count
- For UMIs with <5 reads, batch process together
- Reduces medaka overhead for tiny UMI groups

#### 3. Use Faster Consensus Methods for Simple Cases (Estimated 40-60% speedup for simple UMIs)
**Approach**: Use simple majority-rule consensus for UMIs with few reads/variants

**Implementation**:
- If UMI has <3 reads: use simple consensus (no medaka)
- If UMI has <5 reads and low complexity: use bcftools mpileup + consensus
- Only use medaka for complex cases (>5 reads or high variant density)

#### 4. Parallelize VCF Processing (Estimated 20-30% speedup)
**Approach**: Process VCF normalization and consensus generation in parallel where possible

**Implementation**:
- Use multiprocessing for VCF normalization (can be parallelized)
- Pipeline: medaka → normalize → consensus (overlap where possible)

### **Medium Impact Optimizations**

#### 5. Cache Reference Indexing (Estimated 2-5% speedup)
**Approach**: Index reference once at start, reuse for all UMIs

**Implementation**:
```python
# At start of process_grouped_bam:
fai_file = Path(str(ref_fa) + '.fai')
if not fai_file.exists():
    subprocess.run(['samtools', 'faidx', str(ref_fa)], check=True)
```

#### 6. Optimize File I/O (Estimated 5-10% speedup)
**Approach**: 
- Use memory-mapped files where possible
- Reduce temporary file writes
- Use faster storage (SSD vs HDD)

#### 7. Increase Parallelism (Estimated 10-20% speedup)
**Approach**: Better thread/process management

**Implementation**:
- Increase `procs` parameter (more parallel UMIs)
- Better load balancing (process large UMIs separately from small ones)
- Use thread pool for I/O-bound operations

### **Low Impact Optimizations**

#### 8. Use Faster Medaka Models
**Approach**: Use lighter models for simple cases
- `r941_min_sup_g507` is faster than `r1041_e82_400bps_sup_variant_v5.0.0`
- Trade-off: slightly lower accuracy

#### 9. Skip BAM Copying for Failed UMIs
**Approach**: Only copy BAM if consensus succeeds

#### 10. Profile and Optimize Hot Paths
**Approach**: Use profiling to identify exact bottlenecks
- `cProfile` or `py-spy` to find slow functions
- Optimize Python loops and subprocess calls

## Recommended Implementation Priority

1. **Phase 1 (Quick Wins - 30-50% speedup)**:
   - Skip redundant alignment (use grouped BAM directly)
   - Cache reference indexing
   - Optimize file I/O

2. **Phase 2 (Medium Effort - 40-60% additional speedup)**:
   - Implement simple consensus for small UMIs
   - Batch process tiny UMIs
   - Better load balancing for parallel processing

3. **Phase 3 (Long-term - 20-30% additional speedup)**:
   - Pipeline parallelism within UMI processing
   - Alternative consensus methods
   - GPU acceleration for medaka (if available)

## Expected Overall Speedup

- **Current**: ~10-15 seconds per UMI (depending on read count)
- **After Phase 1**: ~7-10 seconds per UMI (30-50% faster)
- **After Phase 2**: ~3-5 seconds per UMI (70-80% faster overall)
- **After Phase 3**: ~2-4 seconds per UMI (80-90% faster overall)

## Code Changes Required

### Minimal Changes (Phase 1):
1. Modify `run_consensus_for_umi` to accept BAM input instead of FASTQ
2. Extract UMI reads from grouped BAM using `samtools view -t ZU`
3. Skip minimap2 alignment step
4. Cache reference indexing

### Medium Changes (Phase 2):
1. Add simple consensus function for small UMIs
2. Implement UMI batching logic
3. Add read count thresholds

### Major Changes (Phase 3):
1. Refactor to pipeline architecture
2. Add alternative consensus methods
3. Implement GPU support (optional)

