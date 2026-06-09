#!/bin/bash
# =============================================================================
# Script 03: Taxonomic Classification using DIAMOND + MEGAN
# =============================================================================
# Description:
#   Identifies the taxonomic origin of assembled contigs using DIAMOND
#   for high-speed protein alignment against the NCBI-nr database,
#   followed by MEGAN for LCA-based taxonomic and functional classification.
#
# Tools:
#   DIAMOND - https://github.com/bbuchfink/diamond
#   MEGAN CE - https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/
# =============================================================================

### USER-DEFINED PATHS — edit before running
CONTIGS="contigs.fasta"                        # Assembled contigs (from Script 02)
NR_DB="nr.gz"                                  # NCBI-nr database (downloaded below)
DIAMOND_DB="nr.dmnd"                           # DIAMOND database index
DAA_OUT="alignment_result.daa"                 # DIAMOND output
MEGAN_MAP="megan-map-Feb2022.db"               # MEGAN mapping database
THREADS=4                                      # Number of CPU threads
TARGET_TAXON="Siphoviridae"                    # Taxon for read extraction (edit as needed)
TAXON_OUT="Siphoviridae.fasta.gz"              # Output file for extracted reads
###

set -euo pipefail

echo "=== Step 3: Taxonomic Classification (DIAMOND + MEGAN) ==="

# ---- 3a. Install DIAMOND ----
echo "--- 3a. Installing DIAMOND ---"
if ! command -v diamond &>/dev/null; then
  wget -q https://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz
  tar xzf diamond-linux64.tar.gz
  echo "DIAMOND installed."
else
  echo "DIAMOND already available: $(diamond --version)"
fi

# ---- 3b. Download and index NCBI-nr database ----
echo "--- 3b. Downloading NCBI-nr database ---"
echo "WARNING: This download is ~100 GB and may take several hours."
if [ ! -f "$NR_DB" ]; then
  wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
else
  echo "nr.gz already present, skipping download."
fi

echo "Building DIAMOND database index (this may take 30-60 minutes)..."
diamond makedb --in "$NR_DB" --db nr

# ---- 3c. Run DIAMOND alignment ----
echo "--- 3c. Running DIAMOND blastx alignment ---"
# blastx = translated BLAST (nucleotide query vs protein database)
# -f 100 = DAA output format (required for MEGAN)
diamond blastx \
  -d "$DIAMOND_DB" \
  -q "$CONTIGS" \
  -o "$DAA_OUT" \
  -f 100 \
  -p "$THREADS"

echo "Alignment complete: $DAA_OUT"

# ---- 3d. Install MEGAN CE ----
echo "--- 3d. MEGAN CE installation ---"
echo "Download MEGAN CE from: https://uni-tuebingen.de/...megan6/"
echo "Then run: chmod +x MEGAN_Community_unix_6_24_20.sh"
echo "          ./MEGAN_Community_unix_6_24_20.sh"
echo "Unzip mapping databases:"
echo "  unzip megan-map-Feb2022.db.zip"
echo "  unzip megan-nucl-Feb2022.db.zip"

# ---- 3e. MEGANize the DAA file ----
echo "--- 3e. MEGANizing DAA file ---"
if command -v daa-meganizer &>/dev/null; then
  daa-meganizer \
    -i "$DAA_OUT" \
    -mdb "$MEGAN_MAP"
  echo "MEGANization complete."
else
  echo "daa-meganizer not found in PATH. Add MEGAN bin/ to your PATH and rerun."
fi

# ---- 3f. Extract reads for a specific taxon ----
echo "--- 3f. Extracting reads for taxon: $TARGET_TAXON ---"
if command -v read-extractor &>/dev/null; then
  read-extractor \
    -b \
    -i "$DAA_OUT" \
    -o "$TAXON_OUT" \
    -c Taxonomy \
    -n "$TARGET_TAXON"
  echo "Reads extracted to: $TAXON_OUT"
  echo "  -b: includes all reads below the target node in the taxonomy tree"
else
  echo "read-extractor not found. Ensure MEGAN bin/ is in your PATH."
fi

echo "=== Taxonomic classification step complete ==="
echo "Open the .daa file in MEGAN GUI for interactive exploration."
echo "Navigate to KEGG pathway viewer for functional analysis."
