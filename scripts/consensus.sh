#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./make_consensus.sh <reads.fastq> <reference.fasta> <output_dir> [threads] [medaka_model]

READS=$1           # e.g., AAAAATTTCCCAAAC_reads.fastq
REF=$2             # e.g., reference.fasta
OUTDIR=$3          # e.g., output_folder
THREADS=${4:-2}    # default 2 threads
MODEL=${5:-r1041_e82_400bps_sup_variant_v5.0.0}  # default model for R10.4.1 sup

# Extract UMI from reads filename
# e.g., AAAAATTTCCCAAAC_reads.fastq --> AAAAATTTCCCAAAC
BASENAME=$(basename "$READS")
UMI="${BASENAME%%_*}"

echo "✅ UMI detected: $UMI"

# Create or clean output dir
if [ -d "$OUTDIR" ]; then
    echo "⚠️  Output directory exists, clearing old files..."
    rm -rf "$OUTDIR"/*
else
    mkdir -p "$OUTDIR"
fi

# Filenames
ALIGN_BAM="$OUTDIR/aligned.bam"
VCF_RAW="$OUTDIR/medaka.vcf"
VCF_NORM="$OUTDIR/medaka.norm.vcf.gz"
CONSENSUS_FASTA="$OUTDIR/${UMI}_consensus.fasta"

echo "✅ Step 1: Aligning reads to reference..."
minimap2 -t "$THREADS" -ax map-ont "$REF" "$READS" | samtools sort -o "$ALIGN_BAM"

echo "✅ Step 2: Indexing BAM..."
samtools index "$ALIGN_BAM"

echo "✅ Step 3: Calling variants with medaka_variant..."
medaka_variant \
    -i "$READS" \
    -r "$REF" \
    -o "$OUTDIR" \
    -m "$MODEL"

# Detect output VCF
if [ -f "$VCF_RAW" ]; then
    echo "✅ Step 4: Normalizing and compressing VCF..."
    samtools faidx "$REF"
    bcftools norm \
        -f "$REF" \
        -Oz \
        -o "$VCF_NORM" \
        "$VCF_RAW"

    bcftools index "$VCF_NORM"

    VCF_TO_USE="$VCF_NORM"
else
    echo "❌ Error: VCF file not found at $VCF_RAW"
    exit 1
fi

echo "✅ Step 5: Generating consensus FASTA..."
# We pipe the output of bcftools consensus into sed to replace the FASTA header with the UMI
bcftools consensus \
    -f "$REF" \
    "$VCF_TO_USE" | \
    sed "1s|^>.*|>${UMI}|" > "$CONSENSUS_FASTA"

echo "🎉 Consensus FASTA created: $CONSENSUS_FASTA"
