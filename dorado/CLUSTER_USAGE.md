# Running Dorado Demultiplexing on a Cluster

This directory contains scripts for running Dorado demultiplexing jobs on SLURM-based clusters.

## Files Overview

- `run_dorado.sh` - Main Dorado demultiplexing script
- `run_dorado.sbatch` - Basic SLURM batch script
- `run_dorado_template.sbatch` - Flexible template with more options
- `barcode_sequences.fasta` - Example barcode sequences file
- `ED15_barcodes.fasta` - ED15-specific barcode sequences

## Quick Start

### Direct sbatch submission

```bash
# Using the basic sbatch script
sbatch run_dorado.sbatch input.fastq output_dir

# Using the template with custom resources
sbatch --time=48:00:00 --mem=64G --cpus-per-task=32 run_dorado_template.sbatch
```

## Customizing for Your Cluster

### 1. Module Loading

Edit the sbatch files to load the correct modules for your cluster:

```bash
# Common patterns:
module load dorado
module load python/3.9
module load gcc/9.3.0
module load conda
source activate dorado_env
```

### 2. Resource Requirements

Adjust the SLURM directives based on your data size and cluster:

```bash
#SBATCH --time=24:00:00        # Time limit
#SBATCH --mem=32G              # Memory limit
#SBATCH --cpus-per-task=16     # Number of threads
#SBATCH --partition=normal     # Partition name
```

### 3. Email Notifications

Update the email address in the sbatch files:

```bash
#SBATCH --mail-user=your_email@domain.com
```

## Resource Guidelines

### Memory Requirements
- **Small datasets** (< 1GB FASTQ): 16-32GB RAM
- **Medium datasets** (1-10GB FASTQ): 32-64GB RAM  
- **Large datasets** (> 10GB FASTQ): 64GB+ RAM

### Time Requirements
- **Small datasets**: 2-6 hours
- **Medium datasets**: 6-24 hours
- **Large datasets**: 24-48+ hours

### CPU Requirements
- Use 8-32 threads depending on cluster availability
- More threads = faster processing but higher resource usage

## Creating Barcode Files

If you need to create a new barcode sequences file:

```bash
# Create a template barcode file
./run_dorado.sh create my_barcodes.fasta

# Edit the file to add your sequences
nano my_barcodes.fasta
```

Barcode file format:
```
>barcode_name_1
SEQUENCE1
>barcode_name_2
SEQUENCE2
```

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Check job details
scontrol show job <job_id>

# Cancel a job
scancel <job_id>

# View job output
tail -f dorado_demux_<job_id>.out
```

## Troubleshooting

### Common Issues

1. **Dorado not found**: Ensure the dorado module is loaded
2. **Out of memory**: Increase `--mem` parameter
3. **Job timeout**: Increase `--time` parameter
4. **Permission denied**: Check file permissions and paths

### Debug Mode

Check the sbatch script parameters before submitting:

```bash
# Review the sbatch script
cat run_dorado.sbatch

# Submit with custom parameters
sbatch --time=48:00:00 --mem=64G run_dorado.sbatch input.fastq output_dir
```

## Example Workflows

### Single Sample Processing
```bash
sbatch run_dorado.sbatch sample1.fastq sample1_demux
```

### Multiple Samples (Array Job)
```bash
# Create a sample list
echo "sample1.fastq" > samples.txt
echo "sample2.fastq" >> samples.txt

# Submit array job
sbatch --array=1-2 run_dorado_template.sbatch
```

### High-Resource Job
```bash
sbatch --time=48:00:00 --mem=128G --cpus-per-task=64 run_dorado_template.sbatch
```

## Output Files

After successful completion, you'll find:
- Demultiplexed FASTQ files in the output directory
- Job logs: `dorado_demux_<job_id>.out` and `dorado_demux_<job_id>.err`
- Summary statistics in the output directory
