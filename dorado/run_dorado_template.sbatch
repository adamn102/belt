#!/bin/bash
#SBATCH --job-name=dorado_demux
#SBATCH --output=dorado_demux_%j.out
#SBATCH --error=dorado_demux_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@domain.com

# =============================================================================
# Dorado Demultiplexing Job Template
# =============================================================================
# 
# This template provides a flexible sbatch script for running Dorado demultiplexing
# on a cluster. Customize the parameters below based on your cluster and requirements.
#
# Usage Examples:
#   sbatch run_dorado_template.sbatch
#   sbatch --time=48:00:00 --mem=64G run_dorado_template.sbatch
#   sbatch --array=1-10 run_dorado_template.sbatch  # for array jobs
#
# =============================================================================

# Load required modules (customize for your cluster)
# Common module loading patterns:
module load dorado
# module load python/3.9
# module load gcc/9.3.0
# module load conda
# source activate dorado_env

# Set working directory
cd $SLURM_SUBMIT_DIR

# =============================================================================
# JOB PARAMETERS - CUSTOMIZE THESE
# =============================================================================

# Input parameters (modify these for your specific run)
INPUT_FASTQ="input.fastq"                    # Path to your input FASTQ file
OUTPUT_DIR="demultiplexed_output"            # Output directory
KIT_NAME="custom_kit"                        # Custom kit name
BARCODE_FILE="barcode_sequences.fasta"       # Barcode sequences file
THREADS=16                                   # Number of threads (should match --cpus-per-task)

# Optional: Override with command line arguments
if [ $# -ge 1 ]; then INPUT_FASTQ=$1; fi
if [ $# -ge 2 ]; then OUTPUT_DIR=$2; fi
if [ $# -ge 3 ]; then KIT_NAME=$3; fi
if [ $# -ge 4 ]; then BARCODE_FILE=$4; fi
if [ $# -ge 5 ]; then THREADS=$5; fi

# =============================================================================
# VALIDATION AND SETUP
# =============================================================================

echo "=========================================="
echo "Dorado Demultiplexing Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID (if applicable)"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID (if applicable)"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "=========================================="

# Validate input file
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ not found!"
    echo "Please check the file path and try again."
    exit 1
fi

# Validate barcode file
if [ ! -f "$BARCODE_FILE" ]; then
    echo "ERROR: Barcode file $BARCODE_FILE not found!"
    echo "Please ensure the barcode sequences file exists."
    echo "You can create one using: ./run_dorado.sh create my_barcodes.fasta"
    exit 1
fi

# Check if Dorado is available
if ! command -v dorado &> /dev/null; then
    echo "ERROR: Dorado is not available!"
    echo "Please load the dorado module or ensure it's in your PATH"
    echo "Try: module load dorado"
    exit 1
fi

# Create output directory
mkdir -p $OUTPUT_DIR

# =============================================================================
# RUN DORADO DEMULTIPLEXING
# =============================================================================

echo "Starting Dorado demultiplexing..."
echo "Input: $INPUT_FASTQ"
echo "Output: $OUTPUT_DIR"
echo "Kit name: $KIT_NAME"
echo "Barcode file: $BARCODE_FILE"
echo "Threads: $THREADS"
echo ""

# Run the dorado demultiplexing script
./run_dorado.sh $INPUT_FASTQ $OUTPUT_DIR $KIT_NAME $BARCODE_FILE $THREADS

# =============================================================================
# POST-PROCESSING AND REPORTING
# =============================================================================

# Check if the job completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Job completed successfully!"
    echo "End time: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo ""
    echo "Output files:"
    ls -la $OUTPUT_DIR/
    echo ""
    echo "File sizes:"
    du -sh $OUTPUT_DIR/*
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "Job failed with exit code $?"
    echo "End time: $(date)"
    echo "Check the error log for details"
    echo "=========================================="
    exit 1
fi
