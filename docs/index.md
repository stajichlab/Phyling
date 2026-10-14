# Phyling

[![CI/pip](https://github.com/stajichlab/Phyling/actions/workflows/ci-pip.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/ci-pip.yml)
[![Conda](https://anaconda.org/bioconda/phyling/badges/version.svg)](https://anaconda.org/bioconda/phyling)
[![Container](https://ghcr-badge.egpl.dev/stajichlab/phyling/latest_tag?trim=major&label=Image)](https://github.com/stajichlab/Phyling/pkgs/container/phyling)
[![DOI](https://img.shields.io/badge/DOI-10.1093/g3journal/jkag062-blue)](https://doi.org/10.1093/g3journal/jkag062)

**Phyling** is a fast, scalable, and user-friendly software pipeline for phylogenomic reconstruction of species phylogenies directly from protein-encoded or nucleotide genomic data.

Published in *G3: Genes&#124;Genomes&#124;Genetics*: [jkag062 (2026)](https://doi.org/10.1093/g3journal/jkag062).

---

## Key Features

- **Profile HMM Marker Identification**: Directly searches input protein/CDS sequences against pre-curated Hidden Markov Model (HMM) marker sets from the [BUSCO](https://busco-data.ezlab.org/v5/data/lineages/) database, bypassing expensive all-vs-all reciprocal best-hit searches.
- **Support for Protein & DNA Sequences**: Supports protein peptide sequences as well as DNA coding sequences (CDS). For CDS inputs, Phyling automatically translates DNA to protein for alignment and back-translates back to DNA to preserve nucleotide signal among closely related taxa.
- **Multithreaded & High Performance**: Powered by [PyHMMER](https://pyhmmer.readthedocs.io/), achieving massive speedups and significantly reduced memory usage when scaling to thousands of genomes.
- **Incremental Checkpointing**: Saves search and alignment checkpoints to avoid repeating HMM searches when adding or removing samples in iterative analyses.
- **Informativeness Filtering**: Integrates [PhyKIT](https://jlsteenwyk.com/PhyKIT/) to rank and select orthologs based on treeness/RCV (`toverr`) scores, reducing compositional bias and noise.
- **Flexible Tree Inferences**: Supports both **Consensus** (coalescence-based via [ASTER](https://github.com/chaoszhang/ASTER) / [ASTRAL](https://github.com/smirarab/ASTRAL)) and **Concatenation** (supermatrix-based) strategies, with partition model support ([IQ-TREE](http://www.iqtree.org/), [RAxML-NG](https://github.com/amkozlov/raxml-ng), [FastTree](https://github.com/morgannprice/fasttree)).

---

## Pipeline Flowchart

```{image} _static/images/phyling_flowchart-light.svg
:alt: Phyling flowchart
:class: only-light
:align: center
```

```{image} _static/images/phyling_flowchart-dark.svg
:alt: Phyling flowchart
:class: only-dark
:align: center
```

---

## Contents

- **[Installation](installation.md)** — Conda, Pixi, containers, pip, source, developer setup, and `$PHYLING_DB`.
- **[Quick Start](quickstart.md)** — Walkthrough using the bundled example peptide and CDS datasets.
- **Command modules**: [download](usage/download.md) · [align](usage/align.md) · [filter](usage/filter.md) · [tree](usage/tree.md)
- **[Pipeline Architecture & Benchmarks](architecture.md)** — Methodology and scaling performance.
- **[Nextflow Workflow](nextflow.md)** — HPC/cloud execution with [nf_phyling](https://github.com/stajichlab/nf_phyling).
- **[Citation and References](citation.md)** — How to cite Phyling and its dependencies.

```{toctree}
:caption: Getting started
:hidden:

installation
quickstart
```

```{toctree}
:caption: Command modules
:hidden:

usage/download
usage/align
usage/filter
usage/tree
```

```{toctree}
:caption: Background
:hidden:

architecture
nextflow
```

```{toctree}
:caption: Resources
:hidden:

citation
changelog
GitHub repository <https://github.com/stajichlab/Phyling>
Bioconda package <https://anaconda.org/bioconda/phyling>
```
