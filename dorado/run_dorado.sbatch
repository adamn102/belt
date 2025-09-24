#!/bin/bash
#SBATCH --job-name=dorado_demux
#SBATCH --output=dorado_demux_%j.out
#SBATCH --error=dorado_demux_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nomur@mgh.harvard.edu

# Dorado demultiplexing job
# Usage: sbatch run_dorado.sbatch <input.fastq> <output_dir> [kit_name] [barcode_sequences.fasta] [threads]

# Load required modules (adjust based on your cluster)
module load dorado
# module load python/3.9  # if needed
# module load gcc/9.3.0    # if needed

# Set working directory
cd $SLURM_SUBMIT_DIR

# Get job parameters from command line or set defaults
INPUT_FASTQ=${1:-"input.fastq"}
OUTPUT_DIR=${2:-"demultiplexed_output"}
KIT_NAME=${3:-"custom_kit"}
BARCODE_FILE=${4:-"barcode_sequences.fasta"}
THREADS=${5:-16}

# Validate inputs
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ not found!"
    echo "Usage: sbatch run_dorado.sbatch <input.fastq> <output_dir> [kit_name] [barcode_sequences.fasta] [threads]"
    exit 1
fi

if [ ! -f "$BARCODE_FILE" ]; then
    echo "ERROR: Barcode file $BARCODE_FILE not found!"
    echo "Please ensure the barcode sequences file exists or create one using:"
    echo "./run_dorado.sh create my_barcodes.fasta"
    exit 1
fi

# Print job information
echo "=========================================="
echo "Dorado Demultiplexing Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Input: $INPUT_FASTQ"
echo "Output: $OUTPUT_DIR"
echo "Kit name: $KIT_NAME"
echo "Barcode file: $BARCODE_FILE"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo "=========================================="

# Create output directory
mkdir -p $OUTPUT_DIR

# Check if Dorado is available
if ! command -v dorado &> /dev/null; then
    echo "ERROR: Dorado is not available!"
    echo "Please load the dorado module or ensure it's in your PATH"
    exit 1
fi

# Run the dorado demultiplexing script
echo "Starting Dorado demultiplexing..."
$SLURM_SUBMIT_DIR/run_dorado.sh $INPUT_FASTQ $OUTPUT_DIR $KIT_NAME $BARCODE_FILE $THREADS
DORADO_EXIT_CODE=$?

# Check if the job completed successfully
if [ $DORADO_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Job completed successfully!"
    echo "End time: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo "Output files:"
    ls -la $OUTPUT_DIR/
    echo "=========================================="
else
    echo "=========================================="
    echo "Job failed with exit code $DORADO_EXIT_CODE"
    echo "End time: $(date)"
    echo "Check the error log for details"
    echo "=========================================="
    exit 1
fi
