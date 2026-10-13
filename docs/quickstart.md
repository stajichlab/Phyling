# Quick Start

This guide walks through a complete example pipeline execution using sample data included in the Phyling repository repository (`example/` folder).

---

## Example Data Overview

Navigate to the `example` directory:

```sh
cd example
```

The example folder contains:
- `pep/`: 5 gzipped protein FASTA files (`.aa.fasta.gz`).
- `cds/`: 5 gzipped DNA coding sequence FASTA files (`.cds.fasta.gz`).

> **Note:** A minimum of **4 sample genomes** is required to construct a phylogenetic tree.

---

## Pipeline Workflow (Protein Input)

### Step 1: Download Marker Set

First, list available online and local marker sets from the BUSCO database:

```sh
phyling download list
```

Download the fungal lineage marker set (`fungi_odb10`):

```sh
phyling download fungi_odb10
```

---

### Step 2: Align Orthologs

Run `phyling align` on all peptide FASTA files in `example/pep`:

```sh
phyling align -I pep -o align -m fungi_odb10 -t 8
```

Alternatively, pass explicit input files:

```sh
phyling align -i pep/*.fasta.gz -o align -m fungi_odb10 -t 8
```

**Output:** Individual Multiple Sequence Alignment (`.mfa`) files created in the `align/` folder, alongside a checkpoint file for quick re-runs.

---

### Step 3: Filter Alignments by Informativeness

Calculate treeness/RCV (`toverr`) scores and select the top 20 most informative ortholog markers:

```sh
phyling filter -I align -o filtered_align -n 20 -t 8
```

**Output:**
- `treeness.tsv`: Summary of calculated treeness, RCV, and `toverr` scores for each marker.
- `filtered_align/`: Directory containing symlinks to the top 20 `.mfa` alignment files.

---

### Step 4: Construct Species Tree

Build a consensus species tree using the default tree algorithm (FastTree):

```sh
phyling tree -I filtered_align -o tree_output -f -t 8
```

Adding `-f` generates a graphical tree image (`.png`/`.pdf`) using Matplotlib.

**Output:**
- `final_tree.nwk`: Newick formatted species tree file.
- `tree_figure.png`: Tree visualization plot (if `-f` is specified).

---

## Alternative Workflow (DNA CDS Input)

For closely related species, DNA coding sequences (CDS) provide enhanced evolutionary signal compared to conserved protein sequences.

```sh
# 1. Align using CDS inputs (auto-translates to peptide for HMM search/alignment, then back-translates to DNA)
phyling align -I cds -o align_cds -m fungi_odb10 -t 8 --seqtype dna

# 2. Filter top 20 DNA alignments
phyling filter -I align_cds -o filtered_cds -n 20 -t 8 --seqtype dna

# 3. Build tree with concatenated DNA alignments using IQ-TREE
phyling tree -I filtered_cds -o tree_cds_output -c -M iqtree -t 8 --seqtype dna
```

---

## Next Steps

- Explore detailed command options for each module:
  - **[Download](usage/download.md)**
  - **[Align](usage/align.md)**
  - **[Filter](usage/filter.md)**
  - **[Tree](usage/tree.md)**
