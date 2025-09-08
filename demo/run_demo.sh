#!/bin/bash
# BELT Demo Pipeline Runner
# This script runs the complete BELT demo pipeline using conda

set -e  # Exit on any error

echo "🧪 BELT Demo Pipeline"
echo "===================="

# Check if we're in the demo directory
if [ ! -f "Snakefile" ]; then
    echo "❌ Error: Please run this script from the demo directory"
    echo "   cd demo && ./run_demo.sh"
    exit 1
else
    echo "✅ Running demo from demo directory"
fi

# Check if conda/mamba is available
if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "❌ Error: Neither conda nor mamba found"
    echo "   Please install conda or mamba first:"
    echo "   https://docs.conda.io/en/latest/miniconda.html"
    echo "   or https://mamba.readthedocs.io/en/latest/installation.html"
    exit 1
fi

# Use mamba if available, otherwise conda
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "✅ Using mamba for environment management"
else
    CONDA_CMD="conda"
    echo "✅ Using conda for environment management"
fi

# Check if conda environment exists
ENV_NAME="belt_demo"
if ! $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
    echo "📋 Creating conda environment..."
    $CONDA_CMD env create -f environment.yml
    if [ $? -ne 0 ]; then
        echo "❌ Failed to create conda environment"
        exit 1
    fi
    echo "✅ Conda environment created successfully"
else
    echo "✅ Conda environment already exists"
fi

# Activate conda environment
echo "🔧 Activating conda environment..."
source $($CONDA_CMD info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Check if setup has been run
if [ ! -d "data/1-fastq_files" ] || [ ! -f "data/1-fastq_files/B145-XTMMF9-12.fastq" ]; then
    echo "📋 Setting up demo environment..."
    python scripts/demo_setup.py
    if [ $? -ne 0 ]; then
        echo "❌ Setup failed"
        exit 1
    fi
fi

# Check for required tools in conda environment
echo "🔍 Checking required tools in conda environment..."
tools_missing=0

for tool in snakemake minimap2 samtools fastplong; do
    if ! command -v $tool &> /dev/null; then
        echo "   ❌ $tool not found in conda environment"
        tools_missing=1
    else
        echo "   ✅ $tool found"
    fi
done

if [ $tools_missing -eq 1 ]; then
    echo ""
    echo "⚠️  Some required tools are missing from the conda environment."
    echo "   Try updating the environment:"
    echo "   $CONDA_CMD env update -f environment.yml"
    exit 1
fi

# Check for BELT
if ! command -v belt &> /dev/null && [ ! -f "../belt.py" ]; then
    echo "   ❌ BELT not found"
    echo "   Installing BELT in conda environment..."
    pip install -e ..
    if [ $? -ne 0 ]; then
        echo "   ❌ Failed to install BELT"
        exit 1
    fi
    echo "   ✅ BELT installed successfully"
else
    echo "   ✅ BELT found"
fi

echo ""
echo "🚀 Starting BELT Demo Pipeline..."
echo "   This will process the test data through the complete workflow"
echo ""

# Run Snakemake
echo "📊 Running pipeline with 4 cores..."
snakemake --cores 4 --printshellcmds

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Demo pipeline completed successfully!"
    echo ""
    echo "📈 Results:"
    echo "   - Quality filtered reads: data/2-fastq_filtered/"
    echo "   - Mapped reads: data/4-sorted_reads/"
    echo "   - UMI extraction: data/5-umi_maps/"
    echo "   - UMI groups: data/6-umi_groups/"
    echo "   - Consensus sequences: data/7-consensus_calls/"
    echo "   - Compiled consensus: data/8-compiled_consensus/"
    echo "   - Pipeline summary: demo_summary.txt"
    echo ""
    echo "📝 Check demo_summary.txt for detailed results"
    echo ""
    echo "🔧 To deactivate the conda environment:"
    echo "   conda deactivate"
else
    echo ""
    echo "❌ Pipeline failed. Check the logs in data/logs/"
    exit 1
fi
