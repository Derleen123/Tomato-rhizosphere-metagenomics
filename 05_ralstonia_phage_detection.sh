#!/bin/bash
# =============================================================================
# Script 05: Ralstonia-Specific Phage Detection
# =============================================================================
# Description:
#   Applies stringent filtering to identify high-confidence Ralstonia phage
#   sequences from Kraken2 output, then confirms them using BLAST against
#   the NCBI viral RefSeq nucleotide database.
#   Filters: contig length >500 bp AND Kraken2 confidence score >20.
#
# Tools: BLAST+, seqtk
# Input:  ralstonia_confidence.txt — Kraken2 output filtered for Ralstonia hits
#         ralstonia_phages.fna     — FASTA sequences of Ralstonia phage candidates
# =============================================================================

### USER-DEFINED PATHS — edit before running
PROJECT_DIR=~/derleen/kraken_direct
VIRAL_DB_FASTA=~/databases/viral_refseq/viral.1.1.genomic.fna   # Raw FASTA for DB build
VIRAL_DB=~/databases/viral_refseq/viral_refseq_nucl              # BLAST DB path
THREADS=4
MIN_LENGTH=500       # Minimum contig length (bp)
MIN_CONFIDENCE=20    # Minimum Kraken2 confidence score
###

set -euo pipefail
cd "$PROJECT_DIR"

echo "=== Script 05: Ralstonia-Specific Phage Detection ==="

# ---- Check conda environment ----
echo "--- Checking environment ---"
echo "Ensure your conda environment is active: conda activate virsorter2env"

# ---- Check BLAST installation ----
echo "--- BLAST version ---"
blastn -version || { echo "ERROR: blastn not found in PATH. Install with: conda install -c bioconda blast"; exit 1; }

# ---- Build BLAST database if not present ----
echo "--- Checking BLAST database ---"
if ! blastdbcmd -info -db "$VIRAL_DB" >/dev/null 2>&1; then
  echo "BLAST database not found or invalid. Building from: $VIRAL_DB_FASTA"
  makeblastdb -in "$VIRAL_DB_FASTA" -dbtype nucl -out "$VIRAL_DB"
  echo "Database built: $VIRAL_DB"
else
  echo "BLAST database OK: $VIRAL_DB"
fi

# ---- Check required input files ----
echo "--- Checking input files ---"
[ -f ralstonia_confidence.txt ] || { echo "ERROR: ralstonia_confidence.txt not found in $PROJECT_DIR"; exit 1; }
[ -f ralstonia_phages.fna ] || { echo "ERROR: ralstonia_phages.fna not found in $PROJECT_DIR"; exit 1; }
[ -s ralstonia_confidence.txt ] || { echo "seqtk not installed. Install: conda install -c bioconda seqtk"; exit 1; }

# ---- Filter high-confidence sequences ----
echo "--- Filtering: length >$MIN_LENGTH bp AND confidence >$MIN_CONFIDENCE ---"
awk -v minlen="$MIN_LENGTH" -v mincf="$MIN_CONFIDENCE" \
  '$2 > minlen && ($4+0) > mincf {print $1}' \
  ralstonia_confidence.txt > high_confidence_ids.txt

n_hc=$(wc -l < high_confidence_ids.txt)
if [ "$n_hc" -eq 0 ]; then
  echo "No sequences passed the high-confidence filter. Try relaxing thresholds."
  exit 0
fi
echo "High-confidence Ralstonia phage candidates: $n_hc"

# ---- Extract high-confidence sequences ----
seqtk subseq ralstonia_phages.fna high_confidence_ids.txt \
  > ralstonia_phages_highconf.fna
echo "Sequences extracted: $(grep -c '^>' ralstonia_phages_highconf.fna)"

# ---- Run BLAST confirmation ----
echo "--- Running BLAST confirmation (top 3 hits per sequence) ---"
blastn \
  -query ralstonia_phages_highconf.fna \
  -db "$VIRAL_DB" \
  -out ralstonia_highconf_blast.tsv \
  -outfmt "6 qseqid stitle pident length evalue bitscore" \
  -max_target_seqs 3 \
  -evalue 1e-5 \
  -num_threads "$THREADS"

echo "BLAST complete. Total hits: $(wc -l < ralstonia_highconf_blast.tsv)"

# ---- Extract Ralstonia-confirmed hits ----
echo "--- Ralstonia-confirmed BLAST hits ---"
RALSTONIA_HITS=$(grep -i "ralstonia" ralstonia_highconf_blast.tsv || true)
if [ -z "$RALSTONIA_HITS" ]; then
  echo "No Ralstonia-specific hits found in BLAST output."
  echo "Review ralstonia_highconf_blast.tsv for all hits."
else
  echo "Confirmed Ralstonia phage hits:"
  echo "$RALSTONIA_HITS"
  echo ""
  echo "Total Ralstonia-confirmed hits: $(echo "$RALSTONIA_HITS" | wc -l)"
fi

echo "=== Script 05 complete ==="
echo "Full BLAST results: $PROJECT_DIR/ralstonia_highconf_blast.tsv"
