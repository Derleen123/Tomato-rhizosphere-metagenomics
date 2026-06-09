#!/bin/bash
# =============================================================================
# Script 04: Kraken2 Viral Classification + BLAST Confirmation
# =============================================================================
# Description:
#   Extracts all viral sequences classified by Kraken2 from assembled contigs,
#   then confirms classifications using BLAST against the NCBI viral RefSeq
#   nucleotide database. Creates abundance tables at multiple taxonomic levels
#   for plotting and reporting.
#
# Tools: Kraken2, BLAST+, seqtk
# Input:  all_contigs_kraken.txt  — Kraken2 per-sequence output
#         all_contigs_report.txt  — Kraken2 summary report
#         combined_contigs.fa     — assembled contigs FASTA
# Output: final_results/          — all output files
# =============================================================================

### USER-DEFINED PATHS — edit before running
PROJECT_DIR=~/derleen/kraken_direct
CONTIGS_FA=~/derleen/combined_contigs.fa
VIRAL_DB=~/databases/viral_refseq/viral_refseq_nucl   # BLAST viral RefSeq DB
THREADS=2
###

set -euo pipefail
cd "$PROJECT_DIR"

echo "=== Script 04: Kraken2 Viral Classification + BLAST Confirmation ==="

# ---- Clean up previous runs ----
echo "--- Cleaning up previous results ---"
rm -rf combined_results blast_results final_results
rm -f all_viral_classified.txt all_viral_confidence.txt \
      high_quality_viral_ids.txt high_quality_viruses.fna \
      genus_abundances.tsv family_abundances.tsv order_abundances.tsv \
      class_abundances.tsv blast_confirmed_ids.txt confirmed_viral_kraken.txt \
      combined_taxonomy.txt genus_summary.txt all_viruses_blast.tsv \
      blast_species_counts.txt
mkdir -p final_results

# ---- Step 1: Extract all classified viral sequences ----
echo "--- Step 1: Extracting all Kraken2-classified viral sequences ---"
awk '$1=="C"' all_contigs_kraken.txt > final_results/all_viral_classified.txt
echo "Total viral sequences classified by Kraken2: $(wc -l < final_results/all_viral_classified.txt)"

cut -f2 final_results/all_viral_classified.txt > final_results/all_viral_ids.txt
seqtk subseq "$CONTIGS_FA" final_results/all_viral_ids.txt > final_results/all_viral_sequences.fna
echo "Sequences extracted: $(grep -c '>' final_results/all_viral_sequences.fna)"

# ---- Step 2: BLAST all viral sequences against viral RefSeq ----
echo "--- Step 2: BLAST confirmation against NCBI viral RefSeq ---"
# Output format 6 columns: query ID, subject title, % identity,
# alignment length, e-value, bitscore, query length
blastn \
  -query final_results/all_viral_sequences.fna \
  -db "$VIRAL_DB" \
  -out final_results/all_viral_blast.tsv \
  -outfmt "6 qseqid stitle pident length evalue bitscore qlen" \
  -max_target_seqs 1 \
  -evalue 1e-5 \
  -num_threads "$THREADS"

echo "Total BLAST hits: $(wc -l < final_results/all_viral_blast.tsv)"

# ---- Step 3: Count BLAST-confirmed sequences by species ----
echo "--- Step 3: Counting BLAST-confirmed sequences by species ---"
cut -f2 final_results/all_viral_blast.tsv | \
  sed 's/ (taxid.*//' | \
  sed 's/, complete.*//' | \
  sed 's/, partial.*//' | \
  sort | uniq -c | sort -rn > final_results/blast_species_counts.txt

echo "Top 30 BLAST-confirmed species:"
head -30 final_results/blast_species_counts.txt
echo "Total unique species detected: $(wc -l < final_results/blast_species_counts.txt)"

# ---- Step 4: Extract BLAST-confirmed sequences ----
echo "--- Step 4: Extracting BLAST-confirmed sequences ---"
cut -f1 final_results/all_viral_blast.tsv | sort -u > final_results/blast_confirmed_ids.txt
seqtk subseq final_results/all_viral_sequences.fna \
  final_results/blast_confirmed_ids.txt > final_results/blast_confirmed_viruses.fna
echo "BLAST-confirmed viral sequences: $(grep -c '>' final_results/blast_confirmed_viruses.fna)"

# ---- Step 5: Create combined Kraken2 + BLAST taxonomy table ----
echo "--- Step 5: Creating combined taxonomy table ---"
grep -f final_results/blast_confirmed_ids.txt \
  final_results/all_viral_classified.txt > final_results/confirmed_kraken_taxonomy.txt

echo -e "SeqID\tLength\tKraken_Classification\tBLAST_Match\tIdentity\tEvalue" \
  > final_results/combined_taxonomy_table.tsv

join -1 1 -2 1 \
  <(awk '{print $2, $4, $3}' final_results/confirmed_kraken_taxonomy.txt | sort -k1,1) \
  <(awk '{print $1, $2, $3, $5}' final_results/all_viral_blast.tsv | sort -k1,1) | \
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  >> final_results/combined_taxonomy_table.tsv

echo "Combined taxonomy table: final_results/combined_taxonomy_table.tsv"

# ---- Step 6: Abundance tables at multiple taxonomic levels ----
echo "--- Step 6: Creating abundance tables at multiple taxonomic levels ---"
# Extract from Kraken2 report (columns: %, seqs, tax_level, taxon_name)
awk '$4=="C" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' all_contigs_report.txt \
  > final_results/class_abundances.tsv
awk '$4=="F" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' all_contigs_report.txt \
  > final_results/family_abundances.tsv
awk '$4=="G" && $2>0 {printf "%s\t%d\t%s\n",$1,$2,$6}' all_contigs_report.txt \
  > final_results/genus_abundances.tsv

echo "Abundance tables created for Class, Family, Genus levels."

# ---- Step 7: Summary statistics ----
echo "--- Step 7: Generating summary ---"
cat > final_results/SUMMARY.txt << SUMMARY
=== VIRAL DIVERSITY ANALYSIS SUMMARY ===
Total contigs analyzed: 3,353,930
Kraken2 viral classifications: $(wc -l < final_results/all_viral_classified.txt)
BLAST-confirmed viruses:        $(grep -c '>' final_results/blast_confirmed_viruses.fna)
Unique viral species (BLAST):   $(wc -l < final_results/blast_species_counts.txt)

Taxonomic levels detected (Kraken2 report):
  Classes:  $(wc -l < final_results/class_abundances.tsv)
  Families: $(wc -l < final_results/family_abundances.tsv)
  Genera:   $(wc -l < final_results/genus_abundances.tsv)

Key output files:
  all_viral_sequences.fna        — All Kraken2 viral hits (FASTA)
  blast_confirmed_viruses.fna    — BLAST-confirmed subset (FASTA)
  combined_taxonomy_table.tsv    — Full taxonomy for confirmed viruses
  blast_species_counts.txt       — Species abundance counts
  class/family/genus_abundances.tsv — For plotting
SUMMARY

cat final_results/SUMMARY.txt
echo "=== Script 04 complete. Results in: $PROJECT_DIR/final_results/ ==="
