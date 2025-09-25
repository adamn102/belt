# Dorado Basecalling SLURM Scripts

Simple SLURM sbatch scripts for running Dorado basecalling on HPC clusters.

## Scripts

### 1. `dorado_basecall.sbatch` - Basic Basecalling
Simple basecalling without demultiplexing.

**Usage:**
```bash
sbatch dorado_basecall.sbatch <input_dir> <output_dir> [model]
```

**Example:**
```bash
sbatch dorado_basecall.sbatch /path/to/raw_data /path/to/output sup
```

### 2. `dorado_demux.sbatch` - Basecalling + Demultiplexing
Basecalling with barcode demultiplexing.

**Usage:**
```bash
sbatch dorado_demux.sbatch <input_dir> <output_dir> [model] [barcode_kit]
```

**Example:**
```bash
sbatch dorado_demux.sbatch /path/to/raw_data /path/to/output sup SQK-RBK114
```

## Resource Requirements

Both scripts request:
- **Time:** 24 hours
- **CPUs:** 16 cores
- **Memory:** 32 GB
- **Partition:** general

Adjust these in the script headers as needed for your cluster.

## Common Dorado Models

- `sup` - Super accuracy model (recommended)
- `hac` - High accuracy model
- `fast` - Fast model

## Common Barcode Kits

- `SQK-RBK114` - Rapid barcoding kit
- `SQK-RBK004` - Rapid barcoding kit v4
- `SQK-RBK110-96` - 96-sample barcoding kit

## Before Running

1. **Load modules** (uncomment and adjust in scripts):
   ```bash
   module load dorado
   module load minimap2
   module load samtools
   ```

2. **Check Dorado installation:**
   ```bash
   dorado --version
   ```

3. **Verify input directory** contains raw nanopore data files

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Monitor job output
tail -f dorado_basecall_<job_id>.out

# Cancel job if needed
scancel <job_id>
```

## Output

- Basecalled FASTQ files in the specified output directory
- Job logs: `dorado_basecall_<job_id>.out` and `dorado_basecall_<job_id>.err`
- For demultiplexing: separate FASTQ files for each barcode
