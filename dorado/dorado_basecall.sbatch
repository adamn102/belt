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
#SBATCH --gres=gpu:h200_3g.71gb:1 

# Simple Dorado Basecalling SLURM Script with Barcoding
# Usage: sbatch dorado_basecall.sbatch <input_dir> <output_dir> [model] [kit_name]

set -e

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: sbatch dorado_basecall.sbatch <input_dir> <output_dir> [model] [kit_name]"
    echo "Example: sbatch dorado_basecall.sbatch /path/to/raw_data /path/to/output sup SQK-NBD114-24"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
MODEL=${3:-"/project/balazslab/models/dorado/sup"}  # Default to local model path if not specified
KIT_NAME=${4:-"SQK-NBD114-24"}  # Default to SQK-NBD114-24 kit

# Create output directory
mkdir -p $OUTPUT_DIR

# Load modules
module load dorado
# module load minimap2
# module load samtools

# Use default model (no download)
MODEL_DIR="/project/balazslab/models/dorado"

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
echo "Kit: $KIT_NAME"
echo "=========================================="

# Check GPU usage if available
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Information:"
    nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv
    echo "=========================================="
fi

# Start time
start_time=$(date)
echo "Start time: $start_time"

# Run Dorado basecalling with barcoding
echo "Starting Dorado basecalling with barcoding..."
dorado basecaller \
    $MODEL \
    $INPUT_DIR \
    --recursive \
    --output-dir $OUTPUT_DIR \
    --device cuda:all \
    --kit-name $KIT_NAME \
    --emit-fastq

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
