# BELT Demo Pipeline

This is a **minimal demo pipeline** with synthetic data designed to be small enough for GitHub hosting. It demonstrates the complete BELT workflow for UMI-based consensus sequence generation from long-read data.

## 📊 Demo Data Overview

All data files are **synthetic** and optimized for size:

- **Reference sequence**: 2.5 kb synthetic DNA sequence
- **FASTQ reads**: 150 synthetic reads (~734 KB)
- **Whitelist**: 12 UMI pairs
- **Total data size**: < 1 MB (GitHub-friendly!)

### Data Distribution

The synthetic reads contain:
- **6 highly abundant UMI pairs** (used in ~80% of reads)
- **6 less abundant UMI pairs** (used in ~20% of reads)
- Reads are ~2.5 kb with UMI barcodes at positions ~500 bp and ~2000 bp
- Quality scores simulate realistic long-read data (mean Q30)

## 🚀 Quick Start

### 1. Install Dependencies

```bash
# Create conda environment
mamba env create -f environment.yml

# Activate environment
conda activate belt-pipeline-medaka

# Install BELT (if not already installed)
pip install -e ../..
```

### 2. Run the Demo Pipeline

```bash
# Run with 4 cores
snakemake --cores 4

# Or run with more cores for faster execution
snakemake --cores 8
```

The complete demo should finish in **< 2 minutes** on a modern laptop.

### 3. Check Results

```bash
# View summary
cat demo_summary.txt

# Check compiled consensus
cat data/8-compiled_consensus/demo_sample-compiled-consensus.fa

# View individual UMI consensus sequences
ls data/7-consensus_calls/demo_sample/fasta/

# View aligned reads per UMI
ls data/7-consensus_calls/demo_sample/bam/
```

## 📁 Output Structure

```
Belt Pipeline (demo)/
├── data/
│   ├── 0-reference_sequence/      # Reference for mapping
│   ├── 0-consensus_reference/     # Reference for consensus calling
│   ├── 0-library_whitelist/       # UMI whitelists
│   ├── 1-fastq_files/             # Raw synthetic reads
│   ├── 2-fastq_filtered/          # Quality filtered reads
│   ├── 3-mapped_reads/            # Mapped reads (temp)
│   ├── 4-sorted_reads/            # Sorted & indexed BAM files
│   ├── 5-umi_maps/                # Extracted UMIs mapped to whitelist
│   ├── 6-umi_groups/              # Reads grouped by UMI
│   ├── 7-consensus_calls/
│   │   └── demo_sample/
│   │       ├── bam/               # Aligned reads per UMI
│   │       ├── fasta/             # Consensus sequences per UMI
│   │       └── tmp/               # Intermediate files
│   └── 8-compiled_consensus/      # All consensus sequences compiled
├── configs/
│   └── demo.yaml                  # Pipeline configuration
├── logs/                          # Pipeline logs
├── Snakefile                      # Pipeline workflow
├── environment.yml                # Conda environment
└── demo_summary.txt               # Results summary
```

## 🔧 Pipeline Steps

1. **Quality Filtering** (`fastplong`)
   - Filter reads ≥ 1500 bp
   - Quality threshold: Q25
   - Maximum 50% low-quality bases

2. **Read Mapping** (`minimap2`)
   - Map filtered reads to reference
   - Sort and index BAM files

3. **UMI Extraction** (`belt umi-map`)
   - Extract UMIs from aligned reads
   - Map extracted UMIs to whitelist
   - Handle missing/partial UMIs gracefully

4. **UMI Grouping** (`belt umi-group`)
   - Group reads by UMI1
   - Filter groups with ≥ 10 reads
   - Create tagged & sorted BAM

5. **Consensus Calling** (`belt consensus`)
   - Align reads within each UMI group
   - Call variants with Medaka
   - Generate consensus FASTA sequences
   - Save aligned BAM and consensus FASTA per UMI

6. **Compile Results**
   - Combine all consensus sequences
   - Generate summary report

## 📈 Expected Results

With the demo data (threshold = 10 reads per UMI):
- **~3-5 UMI groups** will pass the threshold
- **~3-5 consensus sequences** will be generated
- Total runtime: **< 2 minutes**

### Sample Output

```
BELT Demo Pipeline Results
=========================

Samples processed: 1
Samples: demo_sample

Pipeline Steps Completed:
1. Quality filtering: 1 files
2. Read mapping: 1 files
3. UMI extraction: 1 files
4. UMI grouping: 1 files
5. Consensus calling: 1 files

Output Files:
  - data/8-compiled_consensus/demo_sample-compiled-consensus.fa
```

## 🎯 Customization

### Adjust UMI Threshold

Edit `Snakefile` line ~210:

```python
params:
    threshold = 10  # Change this value
```

Lower values = more UMI groups (but lower quality)
Higher values = fewer UMI groups (but higher quality)

### Change Configuration

Edit `configs/demo.yaml`:

```yaml
samples:
  demo_sample: data/1-fastq_files/demo_sample.fastq

barcodes:
  demo_sample: "['NNNNNNNNNNNN', 'NNNNNNNNNNNN']"

whitelist:
  demo_sample: demo_whitelist

min_lengths:
  demo_sample: 1500  # Minimum read length
```

## 🧪 Testing & Validation

This demo pipeline is ideal for:
- Testing BELT installation
- Learning the BELT workflow
- Validating pipeline changes
- CI/CD integration testing
- Teaching/tutorial purposes

## 📝 Notes

- **Data is synthetic**: Generated for demonstration purposes only
- **Not for production**: Use real experimental data for actual analysis
- **GitHub-friendly**: All files are small enough for version control
- **Fast execution**: Complete pipeline runs in < 2 minutes

## 🔗 Full Pipeline

For the complete pipeline with real data, see:
```
../Belt Pipeline/
```

## 📚 Documentation

For complete BELT documentation, see:
```
../../docs/
```

## 🐛 Troubleshooting

### Pipeline fails to start
```bash
# Check conda environment
conda activate belt-pipeline-medaka

# Verify BELT is installed
belt --help

# Check dependencies
which samtools minimap2 medaka
```

### "No files found" error
```bash
# Regenerate synthetic data
cd "Belt Pipeline (demo)"
python3 << 'EOF'
# [Run the data generation script from earlier]
EOF
```

### Permission errors
```bash
# Unlock Snakemake
snakemake --unlock

# Clean and restart
snakemake --cores 4 --forceall
```

## 📄 License

This demo pipeline is part of the BELT software package.

## 🙋 Questions?

For issues or questions about BELT:
- Check the main documentation in `../../docs/`
- Open an issue on GitHub
- Contact the BELT development team

