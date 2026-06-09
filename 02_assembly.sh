#!/bin/bash
# =============================================================================
# Script 02: Metagenomic Assembly
# =============================================================================
# Description:
#   Reconstructs genomic sequences from mixed microbial communities using
#   either MEGAHIT (recommended for large datasets) or MetaSPAdes.
#   Both options are provided; comment/uncomment as needed.
#
# Tools:
#   MEGAHIT   - https://github.com/voutcn/megahit
#   MetaSPAdes - https://github.com/ablab/spades
# =============================================================================

### USER-DEFINED PATHS — edit before running
READ1="A1_S3_L001_R1_001.fastq.gz"   # Forward reads
READ2="A1_S3_L001_R2_001.fastq.gz"   # Reverse reads
MEGAHIT_OUT="MEGAHIT-OUT"            # MEGAHIT output directory
SPADES_OUT="SPADES-OUT"              # MetaSPAdes output directory
ASSEMBLER="megahit"                  # Choose: "megahit" or "metaspades"
###

set -euo pipefail

echo "=== Step 2: Metagenomic Assembly ==="
echo "Selected assembler: $ASSEMBLER"

# ---- Option A: MEGAHIT ----
if [ "$ASSEMBLER" == "megahit" ]; then
  echo "Downloading and extracting MEGAHIT v1.2.9..."
  wget -q https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
  tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

  echo "Running MEGAHIT assembly..."
  MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
    -1 "$READ1" \
    -2 "$READ2" \
    -o "$MEGAHIT_OUT"

  echo "Assembly complete. Output files in: $MEGAHIT_OUT"
  echo "  contigs.fa    — longest assembled contigs"
  echo "  contigs.fasta — all assembled contigs"
  echo "  log           — assembly summary log"

# ---- Option B: MetaSPAdes ----
elif [ "$ASSEMBLER" == "metaspades" ]; then
  echo "Downloading and extracting MetaSPAdes v3.15.5..."
  wget -q http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz
  tar -xzf SPAdes-3.15.5-Linux.tar.gz

  echo "Running MetaSPAdes assembly..."
  SPAdes-3.15.5-Linux/bin/metaspades.py \
    -1 "$READ1" \
    -2 "$READ2" \
    -o "$SPADES_OUT"

  echo "Assembly complete. Output files in: $SPADES_OUT"
  echo "  contigs.fasta    — assembled contigs"
  echo "  scaffolds.fasta  — assembled scaffolds"
  echo "  metaSPAdes.log   — assembly log"

else
  echo "ERROR: Unknown assembler '$ASSEMBLER'. Choose 'megahit' or 'metaspades'."
  exit 1
fi

echo "=== Assembly step complete ==="
