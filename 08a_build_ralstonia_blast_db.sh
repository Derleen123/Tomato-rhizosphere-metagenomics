#!/bin/bash
# =============================================================================
# Script 08a: Build the local BLAST database used by Script 08
# =============================================================================
# The database combines (1) a local NCBI bacterial RefSeq nucleotide subset and
# (2) all RefSeq Ralstonia assemblies, so that every Ralstonia species (including
# R. insidiosa, R. pseudosolanacearum, R. syzygii) has a reference to hit.
#
# Tools: ncbi-genome-download, BLAST+ (makeblastdb)   (conda install -c bioconda ncbi-genome-download blast)
#
# EDIT these paths before running. BASE_FASTA must be the FASTA file(s) of the
# bacterial RefSeq subset you already have (uncompressed or .gz; see note below).
# =============================================================================
set -euo pipefail
WORK=~/ralstonia_all
BASE_FASTA=(~/databases/bacterial_refseq/bacterial_refseq_subset.fna)   # <-- edit: your existing bacterial RefSeq FASTA file(s)
THREADS=4
mkdir -p "$WORK" && cd "$WORK"

# 1. Download all RefSeq Ralstonia assemblies (FASTA)
ncbi-genome-download bacteria --genera Ralstonia --section refseq \
  --assembly-levels all --formats fasta --output-folder genomes --parallel "$THREADS"

# 2. Concatenate the Ralstonia genomes into one FASTA (files are .fna.gz)
find genomes -name '*_genomic.fna.gz' -print0 | xargs -0 zcat > ralstonia_all.fna
echo "Ralstonia sequences: $(grep -c '>' ralstonia_all.fna)"

# 3. Build one database from the bacterial subset + Ralstonia sequences.
#    makeblastdb needs plain (uncompressed) FASTA: do not mix .gz files in -in.
mkdir -p blastdb
cat "${BASE_FASTA[@]}" ralstonia_all.fna > blastdb/bpr_input.fna
nohup makeblastdb -in blastdb/bpr_input.fna -dbtype nucl -parse_seqids \
  -out blastdb/bpr -title "bacterial RefSeq subset + Ralstonia RefSeq" > makeblastdb.log 2>&1 &
echo "makeblastdb started in the background; check makeblastdb.log, then run: blastdbcmd -db blastdb/bpr -info"
