#!/bin/bash
# =============================================================================
# Script 06: Taxonomic Abundance Extraction at All Levels
# =============================================================================
# Description:
#   Extracts taxonomic abundance counts from the Kraken2 report at all
#   hierarchical levels (Phylum, Class, Order, Family, Genus, Species).
#   Also creates plotting-ready files (Top 15 taxa + Others) and a
#   BLAST-confirmed species abundance table.
#
# Input:  all_contigs_report.txt             — Kraken2 summary report
#         final_results/blast_species_counts.txt — from Script 04
#         final_results/blast_confirmed_ids.txt  — from Script 04
# Output: Abundance TSV files for all taxonomic levels + plotting-ready files
# =============================================================================

### USER-DEFINED PATHS — edit before running
PROJECT_DIR=~/derleen/kraken_direct
RESULTS_DIR=~/derleen/kraken_direct/final_results
###

set -euo pipefail
cd "$RESULTS_DIR"

echo "=== Script 06: Taxonomic Abundance Extraction ==="

# ---- Extract abundance at each taxonomic level from Kraken2 report ----
# Kraken2 report columns: %, seqs_rooted, seqs_direct, tax_level, taxid, name

echo "--- Phylum level ---"
awk '$4=="P" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > phylum_abundances.tsv
echo "Phyla detected: $(wc -l < phylum_abundances.tsv)"

echo "--- Class level ---"
awk '$4=="C" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > class_abundances.tsv
echo "Classes detected: $(wc -l < class_abundances.tsv)"

echo "--- Order level ---"
awk '$4=="O" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > order_abundances.tsv
echo "Orders detected: $(wc -l < order_abundances.tsv)"

echo "--- Family level ---"
awk '$4=="F" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > family_abundances.tsv
echo "Families detected: $(wc -l < family_abundances.tsv)"

echo "--- Genus level ---"
awk '$4=="G" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > genus_abundances.tsv
echo "Genera detected: $(wc -l < genus_abundances.tsv)"

echo "--- Species level (Kraken2) ---"
awk '$4=="S" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' \
  "$PROJECT_DIR/all_contigs_report.txt" > species_abundances_kraken.tsv
echo "Species detected (Kraken2): $(wc -l < species_abundances_kraken.tsv)"

echo "--- Species level (BLAST-confirmed) ---"
if [ -f blast_species_counts.txt ] && [ -f blast_confirmed_ids.txt ]; then
  total=$(wc -l < blast_confirmed_ids.txt)
  awk -v total="$total" '{
    pct = ($1 / total) * 100;
    printf "%.2f\t%d\t%s\n", pct, $1, substr($0, index($0,$2))
  }' blast_species_counts.txt > species_abundances_blast.tsv
  echo "Species detected (BLAST-confirmed): $(wc -l < species_abundances_blast.tsv)"
else
  echo "BLAST results not found. Run Script 04 first."
fi

# ---- Create Top N + Others plotting files ----
echo ""
echo "--- Creating Top 15 + Others plotting files ---"

create_top_n_table() {
  local input_file=$1
  local output_file=$2
  local n=${3:-15}
  head -n "$n" "$input_file" > "$output_file"
  tail -n +"$((n+1))" "$input_file" | \
    awk '{sum+=$2} END { if(sum>0) printf "N/A\t%d\tOthers\n", sum }' \
    >> "$output_file"
  echo "Created: $output_file (top $n + Others)"
}

create_top_n_table class_abundances.tsv    class_top15_plot.tsv    15
create_top_n_table family_abundances.tsv   family_top15_plot.tsv   15
create_top_n_table order_abundances.tsv    order_top15_plot.tsv    15
create_top_n_table genus_abundances.tsv    genus_top20_plot.tsv    20

# ---- Taxonomic richness summary ----
echo ""
echo "=== Taxonomic Richness Summary ==="
printf "%-20s %s\n" "Taxonomic Level" "Groups Detected"
printf "%-20s %s\n" "---------------" "---------------"
printf "%-20s %d\n" "Phylum"   "$(wc -l < phylum_abundances.tsv)"
printf "%-20s %d\n" "Class"    "$(wc -l < class_abundances.tsv)"
printf "%-20s %d\n" "Order"    "$(wc -l < order_abundances.tsv)"
printf "%-20s %d\n" "Family"   "$(wc -l < family_abundances.tsv)"
printf "%-20s %d\n" "Genus"    "$(wc -l < genus_abundances.tsv)"
printf "%-20s %d\n" "Species (Kraken2)"   "$(wc -l < species_abundances_kraken.tsv)"
[ -f species_abundances_blast.tsv ] && \
  printf "%-20s %d\n" "Species (BLAST)"  "$(wc -l < species_abundances_blast.tsv)"

echo ""
echo "=== Script 06 complete ==="
echo "Plotting-ready files: *_top15_plot.tsv, genus_top20_plot.tsv"
ls -lh *_abundances.tsv *_plot.tsv 2>/dev/null
