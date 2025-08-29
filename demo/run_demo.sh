#!/bin/bash
# BELT Demo Pipeline Runner
# This script runs the complete BELT demo pipeline

set -e  # Exit on any error

echo "🧪 BELT Demo Pipeline"
echo "===================="

# Check if we're in the demo directory
if [ ! -f "Snakefile" ]; then
    echo "❌ Error: Please run this script from the demo directory"
    echo "   cd demo && ./run_demo.sh"
    exit 1
fi

# Check if setup has been run
if [ ! -d "data/1-fastq_files" ] || [ ! -f "data/1-fastq_files/B145-XTMMF9-12.fastq" ]; then
    echo "📋 Setting up demo environment..."
    python scripts/demo_setup.py
    if [ $? -ne 0 ]; then
        echo "❌ Setup failed"
        exit 1
    fi
fi

# Check for required tools
echo "🔍 Checking required tools..."
tools_missing=0

for tool in snakemake minimap2 samtools; do
    if ! command -v $tool &> /dev/null; then
        echo "   ❌ $tool not found"
        tools_missing=1
    else
        echo "   ✅ $tool found"
    fi
done

if [ $tools_missing -eq 1 ]; then
    echo ""
    echo "⚠️  Some required tools are missing. Please install:"
    echo "   pip install snakemake"
    echo "   conda install -c bioconda minimap2 samtools"
    echo "   # or install via your system package manager"
    exit 1
fi

# Check for BELT
if ! command -v belt &> /dev/null && [ ! -f "../belt.py" ]; then
    echo "   ❌ BELT not found"
    echo "   Please install BELT first: ../install.sh"
    exit 1
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
else
    echo ""
    echo "❌ Pipeline failed. Check the logs in data/logs/"
    exit 1
fi
