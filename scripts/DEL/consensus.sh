#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./consensus.sh <reads.fastq> <reference.fasta> <meta_dir> <consensus_dir> <threads> <medaka_model> [bam_dir]

READS=$1           # e.g., AAAAATTTCCCAAAC_reads.fastq
REF=$2             # e.g., reference.fasta
METADIR=$3         # e.g., meta directory for intermediate files
OUTDIR=$4          # e.g., consensus output directory
THREADS=${5:-2}    # number of threads
MODEL=${6:-r1041_e82_400bps_sup_variant_v5.0.0}  # medaka model
# $7 is optional bam_dir (not used currently)

# Extract UMI from reads filename
# e.g., AAAAATTTCCCAAAC_reads.fastq --> AAAAATTTCCCAAAC
BASENAME=$(basename "$READS")
UMI="${BASENAME%%_*}"

echo "✅ UMI detected: $UMI"

# Create meta and output directories
mkdir -p "$METADIR"
mkdir -p "$OUTDIR"

# Filenames (intermediate files in METADIR)
ALIGN_BAM="$METADIR/aligned.bam"
VCF_RAW="$METADIR/medaka.vcf"
VCF_NORM="$METADIR/medaka.norm.vcf.gz"
# Final consensus in OUTDIR
CONSENSUS_FASTA="$OUTDIR/${UMI}_consensus.fasta"

echo "✅ Step 1: Aligning reads to reference..."
minimap2 -t "$THREADS" -ax map-ont "$REF" "$READS" | samtools sort -o "$ALIGN_BAM"

echo "✅ Step 2: Indexing BAM..."
samtools index "$ALIGN_BAM"

echo "✅ Step 3: Calling variants with medaka_variant..."
medaka_variant \
    -i "$READS" \
    -r "$REF" \
    -o "$METADIR" \
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
