# Metagenomic Analysis Pipeline
## Metagenomic Analysis of Bacterial-Wilt Tomato Rhizosphere Microbiomes

**Authors:** Derleen Mogire, Johnstone Neondo, Christabel Muhonja, Consolate Awuor, Daniel Mureithi, Stephen Ogada, Tofick Barasa Wekesa
**Institution:** Machakos University / JKUAT, Kenya 

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

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | ≥0.11.9 | Read quality control |
| MEGAHIT | 1.2.9 | Metagenomic assembly |
| MetaSPAdes | 3.15.5 | Alternative assembler |
| Kraken2 | ≥2.1.2 | Taxonomic classification |
| DIAMOND | ≥0.9.36 | Protein alignment |
| MEGAN CE | 6.24.20 | Taxonomic/functional classification |
| BLAST+ | ≥2.12 | Sequence confirmation |
| seqtk | ≥1.3 | Sequence extraction |

---

## Pipeline Structure

```
pipeline/
├── README.md                          # This file
├── 01_quality_control.sh              # FastQC quality assessment
├── 02_assembly.sh                     # Metagenomic assembly (MEGAHIT/MetaSPAdes)
├── 03_taxonomic_classification.sh     # DIAMOND + MEGAN classification
├── 04_kraken2_viral_blast.sh          # Kraken2 viral classification + BLAST confirmation
├── 05_ralstonia_phage_detection.sh    # Ralstonia-specific phage BLAST
├── 06_taxonomy_extraction.sh          # Taxonomic abundance tables
└── 07_blast_confirmed_taxonomy.sh     # BLAST-confirmed taxonomy summary
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
```

---

## Input Data

Raw sequencing reads (paired-end, Illumina NovaSeq 2×150 bp) are deposited at NCBI SRA under BioProject accession PRJNA1273494: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1273494

---

## License

This code is released under the MIT License. See LICENSE file for details.

---

