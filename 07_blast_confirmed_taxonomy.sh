#!/bin/bash
# =============================================================================
# Script 07: BLAST-Confirmed Taxonomy — Final Summary for Publication
# =============================================================================
# Description:
#   Organises all Kraken2 and BLAST-confirmed results into separate folders,
#   parses BLAST hits to infer taxonomy at multiple levels, and produces a
#   publication-ready comparison of Kraken2-only vs BLAST-confirmed diversity.
#
#   This script produces the final abundance tables used for manuscript
#   figures and the data reported in Results Section 3.3.
#
# Input:  final_results/ folder from Scripts 04 and 06
# Output: kraken2_results/, blast_confirmed_results/, comparison/
# =============================================================================

### USER-DEFINED PATHS — edit before running
RESULTS_DIR=~/derleen/kraken_direct/final_results
###

set -euo pipefail
cd "$RESULTS_DIR"

echo "=== Script 07: BLAST-Confirmed Taxonomy — Final Publication Summary ==="

# ---- Organise folder structure ----
echo "--- Organising results into subfolders ---"
mkdir -p kraken2_results blast_confirmed_results comparison

# Move Kraken2-derived files
for f in phylum_abundances.tsv class_abundances.tsv order_abundances.tsv \
          family_abundances.tsv genus_abundances.tsv \
          species_abundances_kraken.tsv *_plot.tsv \
          all_viral_classified.txt all_viral_sequences.fna; do
  [ -f "$f" ] && mv "$f" kraken2_results/ && echo "Moved: $f"
done

# Copy BLAST files for processing
for f in all_viral_blast.tsv blast_confirmed_ids.txt \
          blast_confirmed_viruses.fna blast_species_counts.txt; do
  [ -f "$f" ] && cp "$f" blast_confirmed_results/
done
[ -f combined_taxonomy_table.tsv ] && cp combined_taxonomy_table.tsv comparison/

echo "Folder structure organised."

# ---- Parse BLAST results for taxonomy ----
echo "--- Parsing BLAST hits for taxonomic levels ---"
cd blast_confirmed_results

awk -F'\t' '{
  species = $2;
  gsub(/ \(taxid.*/, "", species);
  gsub(/, complete.*/, "", species);
  gsub(/, partial.*/, "", species);
  split(species, parts, " ");
  genus = parts[1];
  # Infer host-based family grouping
  if (genus ~ /Ralstonia/)    family = "Ralstonia_phages";
  else if (genus ~ /Streptomyces/) family = "Streptomyces_phages";
  else if (genus ~ /Mycobacterium/) family = "Mycobacterium_phages";
  else if (genus ~ /Salmonella/)  family = "Salmonella_phages";
  else if (genus ~ /Escherichia/) family = "Escherichia_phages";
  else if (genus ~ /Pseudomonas/) family = "Pseudomonas_phages";
  else family = genus "_phages";
  # Output: seqid, species, genus, family, % identity, e-value
  print $1"\t"species"\t"genus"\t"family"\t"$3"\t"$5;
}' all_viral_blast.tsv > blast_taxonomy_parsed.tsv

echo "Parsed taxonomy saved: blast_taxonomy_parsed.tsv"

# ---- Count unique taxa ----
total=$(wc -l < blast_taxonomy_parsed.tsv)

cut -f2 blast_taxonomy_parsed.tsv | sort | uniq -c | sort -rn > species_counts.txt
cut -f3 blast_taxonomy_parsed.tsv | sort | uniq -c | sort -rn > genus_counts.txt
cut -f4 blast_taxonomy_parsed.tsv | sort | uniq -c | sort -rn > family_counts.txt

echo "BLAST-confirmed taxa: Species=$(wc -l < species_counts.txt) | Genera=$(wc -l < genus_counts.txt) | Families=$(wc -l < family_counts.txt)"

# ---- Create formatted abundance tables (% + count + name) ----
echo "--- Creating abundance tables ---"
awk -v total="$total" '{
  pct = ($1/total)*100;
  name = substr($0, index($0,$2));
  printf "%.2f\t%d\t%s\n", pct, $1, name;
}' species_counts.txt > species_abundances.tsv

awk -v total="$total" '{printf "%.2f\t%d\t%s\n", ($1/total)*100, $1, $2}' \
  genus_counts.txt > genus_abundances.tsv

awk -v total="$total" '{printf "%.2f\t%d\t%s\n", ($1/total)*100, $1, $2}' \
  family_counts.txt > family_abundances.tsv

# ---- Create Top 15 + Others plotting files ----
for level in species genus family; do
  head -15 "${level}_abundances.tsv" > "${level}_top15_plot.tsv"
  tail -n +16 "${level}_abundances.tsv" | \
    awk '{sum+=$2; pct+=$1} END { if(sum>0) printf "%.2f\t%d\tOthers\n", pct, sum }' \
    >> "${level}_top15_plot.tsv"
  echo "Plotting file: ${level}_top15_plot.tsv"
done

# ---- Publication summary ----
echo "--- Creating publication summary ---"
k2_total=$(wc -l < ../kraken2_results/all_viral_classified.txt 2>/dev/null || echo "N/A")
blast_total=$(wc -l < blast_confirmed_ids.txt)

cat > BLAST_CONFIRMED_SUMMARY.txt << SUMMARY
=== BLAST-CONFIRMED VIRAL TAXONOMY — PUBLICATION SUMMARY ===

Total contigs analyzed:          3,353,930
Kraken2 viral classifications:   $k2_total
BLAST-confirmed viruses:         $blast_total

BLAST-Confirmed Taxonomic Richness:
  Families: $(wc -l < family_abundances.tsv)
  Genera:   $(wc -l < genus_abundances.tsv)
  Species:  $(wc -l < species_abundances.tsv)

Top 10 BLAST-Confirmed Species:
$(head -10 species_abundances.tsv | awk '{printf "  %d sequences (%.1f%%): %s\n", $2, $1, $3}')

Ralstonia Phages:
$(grep -i "ralstonia" species_abundances.tsv | \
  awk '{printf "  %d sequences (%.1f%%): %s\n", $2, $1, $3}')
Total Ralstonia sequences: $(grep -c -i "ralstonia" blast_taxonomy_parsed.tsv || echo 0)

Files for publication figures:
  species_top15_plot.tsv  — species-level bar chart
  genus_top15_plot.tsv    — genus-level bar chart
  family_top15_plot.tsv   — family-level bar chart
  blast_taxonomy_parsed.tsv — full per-sequence taxonomy
SUMMARY

cat BLAST_CONFIRMED_SUMMARY.txt

# ---- Kraken2 vs BLAST comparison ----
echo "--- Creating Kraken2 vs BLAST comparison ---"
cd ../comparison

cat > kraken_vs_blast_comparison.txt << COMPARE
=== KRAKEN2 vs BLAST-CONFIRMED COMPARISON ===

                        Kraken2-only   BLAST-confirmed
Total viral sequences:  $k2_total          $blast_total
Families detected:      $(wc -l < ../kraken2_results/family_abundances.tsv 2>/dev/null || echo "N/A")              $(wc -l < ../blast_confirmed_results/family_abundances.tsv)
Genera detected:        $(wc -l < ../kraken2_results/genus_abundances.tsv 2>/dev/null || echo "N/A")              $(wc -l < ../blast_confirmed_results/genus_abundances.tsv)
Species detected:       $(wc -l < ../kraken2_results/species_abundances_kraken.tsv 2>/dev/null || echo "N/A")              $(wc -l < ../blast_confirmed_results/species_abundances.tsv)

Reporting decision: BLAST-confirmed results used for publication
Rationale: More conservative and scientifically rigorous; every
           reported taxon has direct sequence-level evidence.

Manuscript text: "Kraken2 identified X viral sequences, of which
Y (Z%) were confirmed by BLAST against NCBI viral RefSeq."
COMPARE

cat kraken_vs_blast_comparison.txt

echo ""
echo "=== Script 07 complete ==="
echo ""
echo "Final folder structure:"
echo "  kraken2_results/           — Kraken2-derived abundances"
echo "  blast_confirmed_results/   — BLAST-confirmed taxa (USE FOR PUBLICATION)"
echo "  comparison/                — Kraken2 vs BLAST comparison"
