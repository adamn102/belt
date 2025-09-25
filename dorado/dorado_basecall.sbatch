#!/bin/bash
#SBATCH --job-name=dorado_basecall
#SBATCH --output=dorado_basecall_%j.out
#SBATCH --error=dorado_basecall_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=compute

# Simple Dorado Basecalling SLURM Script
# Usage: sbatch dorado_basecall.sbatch <input_dir> <output_dir> [model]

set -e

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: sbatch dorado_basecall.sbatch <input_dir> <output_dir> [model]"
    echo "Example: sbatch dorado_basecall.sbatch /path/to/raw_data /path/to/output sup"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
MODEL=${3:-"sup"}  # Default to sup model if not specified

# Create output directory
mkdir -p $OUTPUT_DIR

# Load modules (adjust based on your cluster)
# module load dorado
# module load minimap2
# module load samtools

# Check if Dorado is available
if ! command -v dorado &> /dev/null; then
    echo "Error: Dorado is not installed or not in PATH"
    echo "Please load the dorado module or install Dorado"
    exit 1
fi

echo "=========================================="
echo "Dorado Basecalling Job Started"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Model: $MODEL"
echo "=========================================="

# Start time
start_time=$(date)
echo "Start time: $start_time"

# Run Dorado basecalling
echo "Starting Dorado basecalling..."
dorado basecaller \
    $MODEL \
    $INPUT_DIR \
    --threads $SLURM_CPUS_PER_TASK \
    --recursive \
    --output-dir $OUTPUT_DIR

# Check if basecalling was successful
if [ $? -eq 0 ]; then
    echo "✅ Dorado basecalling completed successfully!"
else
    echo "❌ Dorado basecalling failed!"
    exit 1
fi

# End time
end_time=$(date)
echo "End time: $end_time"

# Show output files
echo "=========================================="
echo "Output files:"
ls -lh $OUTPUT_DIR/
echo "=========================================="

# Show job statistics
echo "Job completed successfully!"
echo "Total runtime: $(date -d@$SECONDS -u +%H:%M:%S)"
