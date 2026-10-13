# Metagenomic Analysis Pipeline
## Metagenomic Analysis of Bacterial-Wilt Tomato Rhizosphere Microbiomes

**Authors:** Derleen Mogire¹, Johnstone Neondo¹, Christabel Muhonja², Consolate Awuor², Daniel Mureithi², Stephen Ogada¹, Tofick Barasa Wekesa¹
 
**Affiliations:** 
¹ Institute for Biotechnology Research, Jomo Kenyatta University of Agriculture and Technology (JKUAT), Nairobi, Kenya
² Machakos University, Machakos, Kenya

**Corresponding author:** Derleen Mogire (contact via journal correspondence)

**Manuscript:** Mogire et al. (2026). Metagenomic analysis of bacterial-wilt tomato rhizosphere microbiomes (under review)

---

## Overview

This repository contains all bioinformatics scripts used for the shotgun metagenomic analysis of tomato rhizosphere soil samples from Kirinyaga, Kiambu, and Laikipia counties, Kenya (August–October 2022).

The pipeline covers:
1. Quality control of raw sequencing reads
2. Metagenomic assembly
3. Taxonomic classification (Kraken2 + DIAMOND/MEGAN)
4. BLAST confirmation of viral sequences
5. Ralstonia-specific phage detection
6. Taxonomic abundance extraction and summary
7. BLAST-confirmed viral taxonomy summary
8. Bacterial classification (Kraken2) and species-level identification of *Ralstonia*, confirmed by BLASTn

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | ≥0.11.9 | Read quality control |
| MEGAHIT | 1.2.9 | Metagenomic assembly |
| MetaSPAdes | 3.15.5 | Alternative assembler |
| Kraken2 | ≥2.1.2 | Taxonomic classification |
| DIAMOND | ≥0.9.36 | Protein alignment |
| MEGAN | 6.24.20 | Taxonomic/functional classification |
| BLAST+ | ≥2.12 | Sequence confirmation |
| seqtk | ≥1.3 | Sequence extraction |

### Databases

| Database | Used in | Description |
|----------|---------|-------------|
| NCBI nr + MEGAN mapping `megan-map-Feb2022.db` | Script 03 | DIAMOND/MEGAN classification |
| Kraken2 viral (`k2_viral_20240605`) | Scripts 04, 05, 06 | Prebuilt Kraken2 viral index (release 5 June 2024; downloaded October 2025), run with Kraken2 v2.1.3 on all assembled contigs |
| Kraken2 Standard-8 (`k2_standard_08_GB_20250714`) | Script 08 | Prebuilt Kraken2 index capped at 8 GB, built 14 July 2025 |
| `bacterial_refseq_nucl` | Script 08 | Local BLAST nucleotide database built with `makeblastdb` from NCBI bacterial RefSeq sequences |
| `viral_refseq_nucl` | Scripts 04, 05 | Local BLAST nucleotide database built with `makeblastdb` from NCBI viral RefSeq (`viral.1.1.genomic.fna`, downloaded 26 October 2025) |


---

## Pipeline Structure

```
├── README.md                          # This file
├── 01_quality_control.sh              # FastQC quality assessment
├── 02_assembly.sh                     # Metagenomic assembly (MEGAHIT/MetaSPAdes)
├── 03_taxonomic_classification.sh     # DIAMOND + MEGAN classification
├── 04_kraken2_viral_blast.sh          # Kraken2 viral classification + BLAST confirmation
├── 05_ralstonia_phage_detection.sh    # Ralstonia-specific phage BLAST
├── 06_taxonomy_extraction.sh          # Taxonomic abundance tables
├── 07_blast_confirmed_taxonomy.sh     # BLAST-confirmed viral taxonomy summary
└── 08_bacterial_kraken2_ralstonia_blast.sh  # Bacterial Kraken2 + Ralstonia species ID + BLASTn confirmation
```

---

## Usage

### Step 1 — Set up your environment

```bash
conda create -n metagenomics_env
conda activate metagenomics_env
conda install -c bioconda fastqc megahit spades kraken2 blast seqtk diamond
```

### Step 2 — Edit paths in each script

Each script contains a `### USER-DEFINED PATHS` section at the top.  
Edit these to match your directory structure before running.

### Step 3 — Run scripts in order

```bash
bash 01_quality_control.sh
bash 02_assembly.sh
bash 03_taxonomic_classification.sh
bash 04_kraken2_viral_blast.sh
bash 05_ralstonia_phage_detection.sh
bash 06_taxonomy_extraction.sh
bash 07_blast_confirmed_taxonomy.sh
bash 08_bacterial_kraken2_ralstonia_blast.sh
```

---

### Ralstonia species identification (Script 08)

DIAMOND/MEGAN assigned most *Ralstonia* sequences only to genus level. Script 08 classifies the assembled contigs with Kraken2 to obtain species-level calls, extracts every contig assigned to *Ralstonia*, and confirms each one by BLASTn against bacterial RefSeq. Each contig receives one status:

| Status | Meaning |
|--------|---------|
| `confirmed_species` | Best BLAST hit is *Ralstonia* and the same species as the Kraken2 call |
| `confirmed_genus_only` | Best BLAST hit is *Ralstonia* but a different species, or Kraken2 stopped at genus |
| `not_confirmed` | Best BLAST hit is a different genus |
| `no_blast_hit` | No hit at e-value ≤ 1e-5 |

Main outputs: `ralstonia_contig_confirmation.tsv` (per contig), `ralstonia_species_confirmation.tsv` (per Kraken2 species) and `ralstonia_blast_species_counts.tsv` (species according to BLAST). Counts are numbers of **contigs**, not reads.

Note: species within the *R. solanacearum* species complex (*R. solanacearum*, *R. pseudosolanacearum*, *R. syzygii*) are very closely related, so Kraken2 and BLAST can disagree between them. Such contigs are reported as `confirmed_genus_only`.

---

## Input Data

Raw sequencing reads (paired-end, Illumina NovaSeq 2×150 bp) are deposited at NCBI SRA under BioProject accession PRJNA1273494: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1273494

---

## License

This work is released under the Creative Commons Attribution 4.0 International (CC BY 4.0) licence. See the LICENSE file for details.

---
