# =============================================================================
# Script 08: Bacterial Classification (Kraken2) + Ralstonia Species
#            Identification and BLASTn Confirmation
# =============================================================================
# Description:
#   DIAMOND/MEGAN (Script 03) assigned most Ralstonia sequences only to genus
#   level. This script gets species-level calls by classifying the assembled
#   contigs with Kraken2 and then confirms every Ralstonia contig by BLASTn
#   against a local nucleotide database (NCBI bacterial RefSeq subset plus all
#   RefSeq Ralstonia assemblies; see 08a_build_ralstonia_blast_db.sh).
#
#   1. Kraken2 classification of all contigs (full species names kept)
#   2. Abundance tables for each taxonomic rank (bacteria only, full names)
#   3. Ralstonia taxa from the Kraken2 report (species + strain level)
#   4. Extraction of Ralstonia-assigned contigs (length >= MIN_LEN)
#   5. BLASTn of those contigs against bacterial RefSeq
#   6. Per-contig and per-species confirmation tables
#
# NOTE ON NAMES: Kraken2 report/output files are TAB-separated but taxon names
#   contain spaces ("Ralstonia pickettii"). Always split on tabs (awk -F'\t').
#   Splitting on whitespace keeps only the first word (the genus).
#
# Tools:  Kraken2, BLAST+, seqtk
# Input:  combined_contigs.fa  — MEGAHIT and metaSPAdes contigs (Script 02)
#         concatenated into one file, with no length/coverage filtering and no
#         dereplication. This combined set is used for Kraken2 and BLASTn only
#         (not for abundance estimates); DIAMOND/MEGAN used the MEGAHIT assembly.
# Output: kraken_bacteria/     — all output files (see end of script)
# =============================================================================

### USER-DEFINED PATHS — edit before running
WORKING_DIR=~/derleen
CONTIGS_FA=$WORKING_DIR/combined_contigs.fa
OUT_DIR=$WORKING_DIR/kraken_bacteria
KRAKEN_DB=~/databases/kraken2_standard                        # Kraken2 Standard-8 (k2_standard_08_GB_20250714)
BLAST_DB=~/ralstonia_all/blastdb/bpr   # local BLAST DB: NCBI bacterial RefSeq subset + all RefSeq Ralstonia assemblies (built with 08a_build_ralstonia_blast_db.sh)
THREADS=4
TARGET_GENUS="Ralstonia"
MIN_PIDENT=97    # minimum % identity for a "confirmed_species" call (0 = no cut-off)
MIN_QCOV=50      # minimum % query coverage for a "confirmed_species" call (0 = no cut-off)
MIN_LEN=200      # minimum contig length (bp) of a Ralstonia contig to be BLASTed/confirmed
###

set -euo pipefail

echo "=== Script 08: Bacterial Kraken2 + ${TARGET_GENUS} BLAST confirmation ==="

# ---- Record software and database versions (for Methods) ----
echo "--- Versions ---"
mkdir -p "$OUT_DIR"
{
  echo "Date: $(date)"
  echo "Kraken2 database: $KRAKEN_DB"
  echo "BLAST database:   $BLAST_DB"
  if command -v kraken2 >/dev/null; then kraken2 --version | head -1; fi
  if command -v blastn  >/dev/null; then blastn -version | head -2; fi
  if command -v seqtk   >/dev/null; then echo "seqtk: $(seqtk 2>&1 | grep -i '^Version' || true)"; fi
} > "$OUT_DIR/software_versions.txt"
cat "$OUT_DIR/software_versions.txt"

cd "$OUT_DIR"

# ---- Step 1: Kraken2 classification ----
echo "--- Step 1: Kraken2 classification ---"
if [ -s all_contigs_kraken.txt ] && [ -s all_contigs_report.txt ]; then
  echo "Kraken2 output already present, skipping (delete files to rerun)."
else
  [ -f "$KRAKEN_DB/hash.k2d" ] || { echo "ERROR: Kraken2 database not found: $KRAKEN_DB"; exit 1; }
  [ -f "$CONTIGS_FA" ]         || { echo "ERROR: contigs not found: $CONTIGS_FA"; exit 1; }
  # --use-names is required: column 3 of the output is then "Name (taxid N)"
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --use-names \
    --output all_contigs_kraken.txt \
    --report all_contigs_report.txt \
    "$CONTIGS_FA"
fi

if ! awk -F'\t' '$1=="C"{print $3; exit}' all_contigs_kraken.txt | grep -q '(taxid '; then
  echo "ERROR: all_contigs_kraken.txt has no taxon names. Rerun Kraken2 with --use-names."
  exit 1
fi
echo "Contigs classified: $(awk -F'\t' '$1=="C"' all_contigs_kraken.txt | wc -l) of $(wc -l < all_contigs_kraken.txt)"

# ---- Step 2: Abundance tables (bacteria only, full names, tab-aware) ----
# Report columns: % | clade contigs | direct contigs | rank | taxid | name (indented)
# Output columns: % | contigs | name
echo "--- Step 2: Abundance tables per rank (Bacteria only) ---"
for pair in "P:phylum" "C:class" "O:order" "F:family" "G:genus" "S:species"; do
  rank=${pair%%:*}; level=${pair##*:}
  awk -F'\t' -v OFS='\t' -v r="$rank" '
    { name=$6; gsub(/^ +/, "", name) }
    $4=="U" || $4=="R"  { inb=0 }
    $4=="D"             { inb=(name=="Bacteria") }
    inb && $4==r && $2>0 { print $1, $2, name }
  ' all_contigs_report.txt | sort -t$'\t' -k2,2nr > "bacterial_${level}_abundances.tsv"
  echo "  ${level}: $(wc -l < "bacterial_${level}_abundances.tsv")"
done

# ---- Step 3: Target genus taxa from the Kraken2 report ----
echo "--- Step 3: ${TARGET_GENUS} taxa in the Kraken2 report ---"
# Rank codes: G genus, S species, S1/S2.. strain-level nodes below a species.
# Strain counts are already included in their parent species (clade counts).
printf 'percent\tclade_contigs\tdirect_contigs\trank\ttaxid\tname\n' > ralstonia_kraken2_taxa.tsv
awk -F'\t' -v OFS='\t' -v g="$TARGET_GENUS" '
  { name=$6; gsub(/^ +/, "", name) }
  name ~ ("^" g "( |$)") && ($4=="G" || $4 ~ /^S/) { print $1, $2, $3, $4, $5, name }
' all_contigs_report.txt >> ralstonia_kraken2_taxa.tsv
cat ralstonia_kraken2_taxa.tsv

# ---- Step 4: Extract contigs assigned to the target genus ----
echo "--- Step 4: Extracting ${TARGET_GENUS} contigs ---"
# Output columns: contig_id | taxid | assigned name | contig length
# Only contigs of at least MIN_LEN bp are kept.
awk -F'\t' -v OFS='\t' -v g="$TARGET_GENUS" -v minlen="$MIN_LEN" '
  $1=="C" && $3 ~ ("^" g "( |$)") && ($4+0 >= minlen) {
    name=$3; taxid=$3
    sub(/ \(taxid [0-9]+\)$/, "", name)
    sub(/^.*\(taxid /, "", taxid); sub(/\)$/, "", taxid)
    print $2, taxid, name, $4
  }
' all_contigs_kraken.txt > ralstonia_kraken_contigs.tsv
N_TARGET=$(wc -l < ralstonia_kraken_contigs.tsv)
echo "${TARGET_GENUS} contigs assigned by Kraken2 (>= ${MIN_LEN} bp): $N_TARGET"
[ "$N_TARGET" -gt 0 ] || { echo "No ${TARGET_GENUS} contigs found. Nothing to confirm."; exit 0; }

cut -f1 ralstonia_kraken_contigs.tsv > ralstonia_contig_ids.txt
if [ ! -s ralstonia_contigs.fna ]; then
  command -v seqtk >/dev/null || { echo "ERROR: seqtk not found (conda install -c bioconda seqtk)"; exit 1; }
  seqtk subseq "$CONTIGS_FA" ralstonia_contig_ids.txt > ralstonia_contigs.fna
fi
echo "Sequences extracted: $(grep -c '>' ralstonia_contigs.fna)"

# ---- Step 5: BLASTn confirmation ----
echo "--- Step 5: BLASTn against the local bacterial + Ralstonia RefSeq database ---"
if [ -s ralstonia_blast_all_hits.tsv ]; then
  echo "BLAST output already present, skipping (delete file to rerun)."
else
  command -v blastn >/dev/null || { echo "ERROR: blastn not found (conda install -c bioconda blast)"; exit 1; }
  [ -f "${BLAST_DB}.nhr" ] || [ -f "${BLAST_DB}.00.nhr" ] || { echo "ERROR: BLAST database not found: $BLAST_DB"; exit 1; }
  # Several hits are kept per contig and the best one (highest bitscore) is
  # chosen below. -max_target_seqs 1 does not guarantee the best hit.
  blastn \
    -query ralstonia_contigs.fna \
    -db "$BLAST_DB" \
    -out ralstonia_blast_all_hits.tsv \
    -outfmt "6 qseqid sseqid stitle pident length qlen evalue bitscore qcovs" \
    -max_target_seqs 10 \
    -evalue 1e-5 \
    -num_threads "$THREADS"
fi
# Columns: 1 qseqid 2 sseqid 3 stitle 4 pident 5 length 6 qlen 7 evalue 8 bitscore 9 qcovs
sort -t$'\t' -k1,1 -k8,8gr ralstonia_blast_all_hits.tsv | awk -F'\t' '!seen[$1]++' > ralstonia_blast_best_hit.tsv
echo "Contigs with a BLAST hit: $(wc -l < ralstonia_blast_best_hit.tsv) of $N_TARGET"

# ---- Step 6: Per-contig confirmation table ----
echo "--- Step 6: Confirmation tables ---"
# Status codes:
#   confirmed_species     best BLAST hit is the target genus AND the same species as Kraken2
#                         (and passes MIN_PIDENT / MIN_QCOV)
#   confirmed_genus_only  best BLAST hit is the target genus but a different species,
#                         Kraken2 stopped at genus, or the hit is below the cut-offs
#   not_confirmed         best BLAST hit is a different genus
#   no_blast_hit          no hit at e-value 1e-5
# Species names are the first two words after the accession ("Ralstonia pickettii");
# "Ralstonia sp. R-29" keeps three words.
awk -F'\t' -v OFS='\t' -v g="$TARGET_GENUS" -v minid="$MIN_PIDENT" -v mincov="$MIN_QCOV" '
  function species(name,   w, n) {
    n = split(name, w, " ")
    if (n < 2) return name
    if (w[2] == "sp." && n >= 3) return w[1] " " w[2] " " w[3]
    return w[1] " " w[2]
  }
  FNR==NR { n++; id[n]=$1; ktax[$1]=$2; kname[$1]=$3; klen[$1]=$4; next }
  { hit[$1]=$0 }
  END {
    print "contig_id","length","kraken_taxid","kraken_assignment","kraken_species", \
          "blast_best_hit","blast_species","pident","query_coverage","evalue","status"
    for (i=1; i<=n; i++) {
      c = id[i]; ks = species(kname[c])
      if (c in hit) {
        split(hit[c], f, "\t")
        title = f[3]; sub(/^[^ ]+ /, "", title)       # drop accession
        bs = species(title); split(bs, bg, " ")
        if (bg[1] != g)                                        st = "not_confirmed"
        else if (bs == ks && f[4]+0 >= minid && f[9]+0 >= mincov) st = "confirmed_species"
        else                                                   st = "confirmed_genus_only"
        print c, klen[c], ktax[c], kname[c], ks, title, bs, f[4], f[9], f[7], st
      } else {
        print c, klen[c], ktax[c], kname[c], ks, "NA", "NA", "NA", "NA", "NA", "no_blast_hit"
      }
    }
  }
' ralstonia_kraken_contigs.tsv ralstonia_blast_best_hit.tsv > ralstonia_contig_confirmation.tsv

# Per Kraken2 species: how many contigs were confirmed
awk -F'\t' -v OFS='\t' '
  NR==1 { next }
  { s=$5; tot[s]++; cnt[s SUBSEP $11]++ }
  END {
    print "kraken_species","contigs","confirmed_species","confirmed_genus_only","not_confirmed","no_blast_hit"
    for (s in tot) print s, tot[s], cnt[s,"confirmed_species"]+0, cnt[s,"confirmed_genus_only"]+0, \
                          cnt[s,"not_confirmed"]+0, cnt[s,"no_blast_hit"]+0
  }
' ralstonia_contig_confirmation.tsv | { read -r header; echo "$header"; sort -t$'\t' -k2,2nr; } \
  > ralstonia_species_confirmation.tsv

# Species according to BLAST (all contigs whose best hit is the target genus)
awk -F'\t' -v OFS='\t' -v g="$TARGET_GENUS" '
  NR==1 { next }
  $11=="confirmed_species" || $11=="confirmed_genus_only" {
    b=$7; tot[b]++; if ($11=="confirmed_species") agree[b]++
  }
  END {
    print "blast_species","contigs","agree_with_kraken2"
    for (b in tot) print b, tot[b], agree[b]+0
  }
' ralstonia_contig_confirmation.tsv | { read -r header; echo "$header"; sort -t$'\t' -k2,2nr; } \
  > ralstonia_blast_species_counts.tsv

# ---- Summary ----
echo "--- Summary ---"
N_SP=$(awk -F'\t' 'NR>1 && $11=="confirmed_species"' ralstonia_contig_confirmation.tsv | wc -l)
N_GN=$(awk -F'\t' 'NR>1 && $11=="confirmed_genus_only"' ralstonia_contig_confirmation.tsv | wc -l)
N_NC=$(awk -F'\t' 'NR>1 && $11=="not_confirmed"' ralstonia_contig_confirmation.tsv | wc -l)
N_NH=$(awk -F'\t' 'NR>1 && $11=="no_blast_hit"' ralstonia_contig_confirmation.tsv | wc -l)
cat > RALSTONIA_CONFIRMATION_SUMMARY.txt << SUMMARY
=== ${TARGET_GENUS} SPECIES IDENTIFICATION AND BLAST CONFIRMATION ===
Contigs assigned to ${TARGET_GENUS} by Kraken2 (>= ${MIN_LEN} bp):   $N_TARGET
  confirmed at species level:        $N_SP
  confirmed at genus level only:     $N_GN
  not confirmed (other genus):       $N_NC
  no BLAST hit:                      $N_NH
Cut-offs for species confirmation:   identity >= ${MIN_PIDENT}%, query coverage >= ${MIN_QCOV}%, contig length >= ${MIN_LEN} bp

Per Kraken2 species:
$(cat ralstonia_species_confirmation.tsv)

Species according to BLAST:
$(cat ralstonia_blast_species_counts.tsv)
SUMMARY
cat RALSTONIA_CONFIRMATION_SUMMARY.txt

# ---- Optional: redundancy check between the two assemblies (off by default) ----
# Set RUN_CDHIT_CHECK=1 to cluster the Ralstonia contigs (95% identity, 80% shorter-
# sequence coverage). This only reports how many clusters there are; it does not
# change any result above. (In this study: 269 contigs -> 264 clusters.)
if [ "${RUN_CDHIT_CHECK:-0}" = "1" ]; then
  command -v cd-hit-est >/dev/null || { echo "ERROR: cd-hit-est not found (conda install -c bioconda cd-hit)"; exit 1; }
  cd-hit-est -i ralstonia_contigs.fna -o ralstonia_contigs_cdhit95.fna -c 0.95 -aS 0.8 -T "$THREADS" > cdhit.log
  echo "Ralstonia contigs: $N_TARGET; CD-HIT-EST clusters (95% identity): $(grep -c '>' ralstonia_contigs_cdhit95.fna)"
fi

echo ""
echo "=== Script 08 complete. Results in: $OUT_DIR ==="
echo "  bacterial_<rank>_abundances.tsv      — Kraken2 abundances, full names (percent, contigs, name)"
echo "  ralstonia_kraken2_taxa.tsv           — ${TARGET_GENUS} taxa from the Kraken2 report"
echo "  ralstonia_contig_confirmation.tsv    — per-contig Kraken2 call, best BLAST hit, status"
echo "  ralstonia_species_confirmation.tsv   — confirmation status per Kraken2 species"
echo "  ralstonia_blast_species_counts.tsv   — species according to BLAST"
echo "  software_versions.txt                — Kraken2 / BLAST+ / seqtk versions and databases"