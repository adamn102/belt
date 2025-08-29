# BELT Demo Pipeline Summary

## 🎯 Overview

The BELT Demo Pipeline is a complete, runnable example that demonstrates the full power of the BELT toolkit. It processes real test data through the entire workflow from raw FASTQ files to consensus sequences.

## 📊 What the Demo Demonstrates

### 1. **Complete Workflow Integration**
- Shows how all BELT components work together
- Demonstrates proper data flow between steps
- Illustrates best practices for pipeline design

### 2. **Real Data Processing**
- Uses actual test data from the BELT test suite
- Processes ~3,600 long reads with UMIs
- Generates real consensus sequences

### 3. **Professional Pipeline Design**
- Snakemake-based workflow for reproducibility
- Proper error handling and logging
- Configurable parameters via YAML
- Modular rule design for easy customization

## 🔧 Pipeline Components

### **Quality Filtering**
- Filters raw FASTQ files for quality and length
- Uses `fastplong` for comprehensive quality control
- Generates HTML reports for quality assessment

### **Read Mapping**
- Maps reads to reference sequences using `minimap2`
- Sorts and indexes BAM files with `samtools`
- Prepares data for UMI extraction

### **UMI Processing**
- **Extraction**: Extracts UMIs from mapped reads
- **Mapping**: Corrects UMI errors using clustering
- **Grouping**: Groups reads by corrected UMIs

### **Consensus Calling**
- Generates consensus sequences from UMI groups
- Uses external tools for high-quality consensus
- Compiles results into final FASTA files

## 📈 Expected Results

When run successfully, the demo produces:

```
data/
├── 2-fastq_filtered/          # Quality-filtered reads
├── 4-sorted_reads/            # Mapped and sorted BAM files
├── 5-umi_maps/                # UMI extraction results
├── 6-umi_groups/              # UMI grouping results
├── 7-consensus_calls/         # Consensus sequences
└── 8-compiled_consensus/      # Final compiled results
```

**Key Metrics:**
- **Input**: ~3,600 raw reads
- **Filtered**: ~3,000 quality reads
- **UMI Groups**: ~5 groups meeting threshold
- **Consensus**: ~5 high-quality consensus sequences

## 🎓 Using as a Template

### **For Your Own Data**

1. **Replace Input Data**:
   ```bash
   # Place your FASTQ files
   cp your_data.fastq data/1-fastq_files/
   
   # Add your reference
   cp your_reference.fa data/0-reference_sequence/
   ```

2. **Update Configuration**:
   ```yaml
   # Edit config.yaml
   samples:
     your_sample: data/1-fastq_files/your_data.fastq
   
   barcodes:
     your_sample: YOUR_UMI_PATTERN
   ```

3. **Run the Pipeline**:
   ```bash
   snakemake --cores 8
   ```

### **Customization Options**

- **Add new samples**: Update `config.yaml`
- **Modify parameters**: Edit quality thresholds, UMI patterns
- **Extend workflow**: Add new rules for additional analysis
- **Change tools**: Replace minimap2, fastplong, etc.

## 🔍 Key Learning Points

### **Snakemake Best Practices**
- Use wildcards for sample-specific processing
- Implement proper input/output dependencies
- Add logging and error handling
- Use configuration files for parameters

### **BELT Integration**
- Proper command-line usage for each BELT component
- Data format requirements and expectations
- Error handling and validation
- Output file organization

### **Bioinformatics Workflow Design**
- Modular rule design
- Reproducible execution
- Scalable architecture
- Clear documentation

## 🚀 Next Steps

After running the demo:

1. **Examine Results**: Check output files and logs
2. **Modify Parameters**: Experiment with different settings
3. **Add Your Data**: Replace test data with your own
4. **Extend Workflow**: Add additional analysis steps
5. **Deploy**: Use as template for production pipelines

## 📚 Resources

- **BELT Documentation**: Main README
- **Snakemake Tutorial**: https://snakemake.readthedocs.io/
- **Test Data**: Available in `scripts/test-data/`
- **Real Examples**: Check the main Snakefile for production examples

---

*This demo demonstrates the complete BELT workflow and serves as both a learning tool and a template for your own pipelines.*
