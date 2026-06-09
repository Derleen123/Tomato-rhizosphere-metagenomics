#!/bin/bash
# =============================================================================
# Script 01: Quality Control of Raw Sequencing Reads
# =============================================================================
# Description:
#   Evaluates the accuracy and reliability of raw paired-end sequencing reads
#   using FastQC. Outputs HTML quality reports for visual inspection.
#
# Tool: FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# =============================================================================

### USER-DEFINED PATHS — edit before running
READ1="A1_S3_L001_R1_001.fastq.gz"   # Forward reads
READ2="A1_S3_L001_R2_001.fastq.gz"   # Reverse reads
OUTDIR="fastqc_results"               # Output directory
###

set -euo pipefail

echo "=== Step 1: Quality Control ==="

# Install FastQC if not available
if ! command -v fastqc &>/dev/null; then
  echo "FastQC not found. Installing..."
  sudo apt-get install -y fastqc
fi

# Create output directory
mkdir -p "$OUTDIR"

# Run FastQC on both reads
echo "Running FastQC on forward reads..."
fastqc "$READ1" -o "$OUTDIR"

echo "Running FastQC on reverse reads..."
fastqc "$READ2" -o "$OUTDIR"

echo "Quality control complete. Reports saved to: $OUTDIR"
echo "Open the .html files in a browser to inspect read quality."
