#!/bin/bash
# Dorado demultiplexing script for ED15 custom barcodes
# Usage: ./run_dorado_demux.sh input.fastq output_dir

set -e

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.fastq> <output_dir>"
    echo "Example: $0 input.fastq demultiplexed_output"
    exit 1
fi

INPUT_FASTQ=$1
OUTPUT_DIR=$2

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

echo "Starting Dorado demultiplexing..."
echo "Input: $INPUT_FASTQ"
echo "Output: $OUTPUT_DIR"

# Run Dorado with custom barcodes
# Note: You may need to adjust the model and barcode kit parameters
# based on your specific Dorado installation and requirements

dorado basecaller \
    --barcode-kits custom \
    --kit-name ED15_custom_kit \
    --barcode-sequences dorado_custom_barcodes.fasta \
    --output-dir $OUTPUT_DIR \
    --threads 8 \
    $INPUT_FASTQ

echo "Dorado demultiplexing completed!"
echo "Check the output directory: $OUTPUT_DIR"

# List the output files
echo "Output files:"
ls -la $OUTPUT_DIR/
