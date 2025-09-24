#!/bin/bash
# Dorado demultiplexing script for custom barcodes
# Usage: ./run_dorado.sh input.fastq output_dir [kit_name] [barcode_sequences.fasta] [threads]

set -e

# Default values
DEFAULT_KIT_NAME="custom_kit"
DEFAULT_BARCODE_FILE="barcode_sequences.fasta"
DEFAULT_THREADS=8

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.fastq> <output_dir> [kit_name] [barcode_sequences.fasta] [threads]"
    echo "Example: $0 input.fastq demultiplexed_output"
    echo "Example: $0 input.fastq demultiplexed_output my_kit_name my_barcodes.fasta 16"
    echo ""
    echo "Parameters:"
    echo "  input.fastq              - Input FASTQ file"
    echo "  output_dir               - Output directory"
    echo "  kit_name                 - Custom kit name (default: $DEFAULT_KIT_NAME)"
    echo "  barcode_sequences.fasta   - Barcode sequences file (default: $DEFAULT_BARCODE_FILE)"
    echo "  threads                  - Number of threads (default: $DEFAULT_THREADS)"
    exit 1
fi

INPUT_FASTQ=$1
OUTPUT_DIR=$2
KIT_NAME=${3:-$DEFAULT_KIT_NAME}
BARCODE_FILE=${4:-$DEFAULT_BARCODE_FILE}
THREADS=${5:-$DEFAULT_THREADS}

# Helper function to create new barcode files
create_barcode_file() {
    local output_file=$1
    cat > "$output_file" << 'EOF'
# Add your barcode sequences below
# Format: >barcode_name
#         SEQUENCE

>barcode_1
CCTACGTTCTCGTTGCCGTT
>barcode_2
GAGGTTGCAATCGTACTCCC
>barcode_3
GCCGAACATAGTGATATCAG
EOF
}

# If the script is called with 'create' as first argument, create a new barcode file
if [ "$1" = "create" ]; then
    if [ $# -lt 2 ]; then
        echo "Usage: $0 create <output_barcode_file.fasta>"
        echo "Example: $0 create my_new_barcodes.fasta"
        exit 1
    fi
    create_barcode_file $2
    echo "✅ Created new barcode file: $2"
    echo "Edit the file to add your custom barcodes"
    exit 0
fi

# Create output directory
mkdir -p $OUTPUT_DIR

# Check if Dorado is available
if ! command -v dorado &> /dev/null; then
    echo "Error: Dorado is not installed or not in PATH"
    echo "Please install Dorado from: https://github.com/nanoporetech/dorado"
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input file $INPUT_FASTQ does not exist"
    exit 1
fi

# Check if barcode file exists
if [ ! -f "$BARCODE_FILE" ]; then
    echo "Error: Barcode sequences file $BARCODE_FILE does not exist"
    echo "Please create a FASTA file with your barcode sequences"
    echo "Example format:"
    echo ">barcode_1"
    echo "CCTACGTTCTCGTTGCCGTT"
    echo ">barcode_2"
    echo "GAGGTTGCAATCGTACTCCC"
    exit 1
fi

echo "Starting Dorado demultiplexing..."
echo "Input: $INPUT_FASTQ"
echo "Output: $OUTPUT_DIR"
echo "Kit name: $KIT_NAME"
echo "Barcode file: $BARCODE_FILE"
echo "Threads: $THREADS"

# Run Dorado with custom barcodes
dorado basecaller \
    --barcode-kits custom \
    --kit-name $KIT_NAME \
    --barcode-sequences $BARCODE_FILE \
    --output-dir $OUTPUT_DIR \
    --threads $THREADS \
    $INPUT_FASTQ

echo "Dorado demultiplexing completed!"
echo "Check the output directory: $OUTPUT_DIR"

# List the output files
echo "Output files:"
ls -la $OUTPUT_DIR/
