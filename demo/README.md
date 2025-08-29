# BELT Demo Pipeline

This directory contains a complete, runnable demo of the BELT pipeline using Snakemake. This demo demonstrates the full power of BELT by processing real test data through the entire workflow.

## 🎯 What This Demo Does

The demo pipeline processes long-read sequencing data through the complete BELT workflow:

1. **Quality Filtering** - Filter raw FASTQ files for quality and length
2. **Read Mapping** - Map reads to reference sequences using minimap2
3. **UMI Extraction** - Extract UMIs from mapped reads
4. **UMI Mapping** - Map UMIs to whitelists and correct errors
5. **UMI Grouping** - Group reads by corrected UMIs
6. **Consensus Calling** - Generate consensus sequences from UMI groups
7. **Analysis** - Compile results and generate statistics

## 📁 Directory Structure

```
demo/
├── README.md                 # This file
├── Snakefile                 # Main Snakemake pipeline
├── config.yaml              # Pipeline configuration
├── requirements.txt          # Python dependencies
├── run_demo.sh              # Demo execution script
├── data/                    # Input and output data
│   ├── 0-reference_sequence/    # Reference genomes
│   ├── 0-library_whitelist/     # UMI whitelists
│   ├── 1-fastq_files/           # Raw FASTQ files
│   ├── 2-fastq_filtered/        # Quality-filtered FASTQ
│   ├── 3-mapped_reads/          # BAM files after mapping
│   ├── 4-sorted_reads/          # Sorted and indexed BAM
│   ├── 5-umi_maps/              # UMI extraction results
│   ├── 6-umi_groups/            # UMI grouping results
│   ├── 7-consensus_calls/       # Consensus sequences
│   ├── 8-compiled_consensus/    # Compiled consensus files
│   └── logs/                    # Pipeline logs
└── scripts/                 # Helper scripts
    └── demo_setup.py        # Demo data setup script
```

## 🚀 Quick Start

### Prerequisites

1. **BELT installed** - Follow the main README installation instructions
2. **Snakemake** - Install with `pip install snakemake`
3. **External tools**:
   - `minimap2` - For read mapping
   - `samtools` - For BAM file manipulation
   - `fastplong` - For quality filtering

### Running the Demo

1. **Setup the demo environment:**
   ```bash
   cd demo
   python scripts/demo_setup.py
   ```

2. **Run the pipeline:**
   ```bash
   ./run_demo.sh
   ```

3. **Or run manually:**
   ```bash
   snakemake --cores 4 --use-conda
   ```

## 📊 Demo Data

The demo uses real test data from the BELT test suite:

- **Sample**: B145-XTMMF9-12
- **Reference**: B145.fa (HIV reference sequence)
- **UMI Pattern**: NNNNNNNNNNNNN (13 N's)
- **Reads**: ~3,600 long reads with UMIs

## ⚙️ Configuration

The pipeline is configured via `config.yaml`:

```yaml
samples:
  B145-XTMMF9-12: data/1-fastq_files/B145-XTMMF9-12.fastq

barcodes:
  B145-XTMMF9-12: NNNNNNNNNNNNN

min_lengths:
  B145-XTMMF9-12: 2000

whitelist:
  B145-XTMMF9-12: false  # Use clustering instead of whitelist
```

## 🔧 Customization

### Adding Your Own Data

1. **Place your FASTQ files** in `data/1-fastq_files/`
2. **Add your reference** to `data/0-reference_sequence/`
3. **Update config.yaml** with your sample information
4. **Run the pipeline**

### Modifying Parameters

Edit `config.yaml` to adjust:
- **Barcode motifs** for UMI extraction
- **Minimum read lengths** for filtering
- **UMI whitelists** for error correction
- **Grouping thresholds** for consensus calling

## 📈 Expected Results

After running the demo, you should see:

- **Quality filtered reads**: ~3,000 reads passing filters
- **UMI extraction**: ~276 reads with extracted UMIs
- **UMI groups**: ~5 groups meeting threshold
- **Consensus sequences**: ~5 consensus FASTA files
- **Analysis reports**: Quality metrics and statistics

## 🐛 Troubleshooting

### Common Issues

1. **Missing dependencies**: Install required tools
2. **Memory issues**: Reduce `--cores` parameter
3. **File permissions**: Ensure write access to demo directory

### Getting Help

- Check the logs in `data/logs/`
- Review the Snakemake dry-run: `snakemake -n`
- Consult the main BELT documentation

## 🎓 Learning Resources

- **Snakemake Tutorial**: https://snakemake.readthedocs.io/
- **BELT Documentation**: See main README
- **Pipeline Components**: Each rule is documented in the Snakefile

## 📝 License

This demo is part of the BELT project and follows the same license terms.
