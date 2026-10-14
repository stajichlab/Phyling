[![CI/pip](https://github.com/stajichlab/Phyling/actions/workflows/ci-pip.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/ci-pip.yml)
[![CI/conda](https://github.com/stajichlab/Phyling/actions/workflows/ci-conda.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/ci-conda.yml)
[![CI/pixi](https://github.com/stajichlab/Phyling/actions/workflows/ci-pixi.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/ci-pixi.yml)
[![Python](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12_%7C_3.13_%7C_3.14-blue?logo=python)](https://github.com/stajichlab/Phyling/actions/workflows/ci-pip.yml)
[![codecov](https://codecov.io/gh/stajichlab/Phyling/graph/badge.svg?token=ZH5GBQYKZ6)](https://codecov.io/gh/stajichlab/Phyling)
[![License](https://img.shields.io/github/license/stajichlab/Phyling?label=license)](https://github.com/stajichlab/Phyling/blob/main/LICENSE)
[![Conda](https://anaconda.org/bioconda/phyling/badges/version.svg)](https://anaconda.org/bioconda/phyling)
[![Container](https://ghcr-badge.egpl.dev/stajichlab/phyling/latest_tag?trim=major&label=Image)](https://github.com/stajichlab/Phyling/pkgs/container/phyling)
[![Documentation Status](https://app.readthedocs.org/projects/phyling/badge/?version=latest)](https://phyling.readthedocs.io)
[![DOI](https://img.shields.io/badge/DOI-10.1093/g3journal/jkag062-blue)](https://doi.org/10.1093/g3journal/jkag062)

# Phyling

**Phyling** is a fast, scalable, and user-friendly software pipeline for phylogenomic reconstruction of species phylogenies directly from protein-encoded or nucleotide genomic data.

Published in *G3: Genes&#124;Genomes&#124;Genetics*: [jkag062 (2026)](https://doi.org/10.1093/g3journal/jkag062).

---

## Pipeline Flowchart

<p align="center">
  <img src="https://raw.githubusercontent.com/stajichlab/Phyling/main/docs/_static/images/phyling_flowchart-light.svg#gh-light-mode-only" alt="Phyling flowchart" width="800">
  <img src="https://raw.githubusercontent.com/stajichlab/Phyling/main/docs/_static/images/phyling_flowchart-dark.svg#gh-dark-mode-only" alt="Phyling flowchart" width="800">
</p>

---

## Quick Installation

### Via Conda (Recommended)

```sh
conda install bioconda::phyling
```

### Via Pixi

```sh
pixi global install -c bioconda phyling
```

*For alternative installation options (pip, Git clone, developer setup), see the [Installation Guide](https://phyling.readthedocs.io/page/installation.html).*

---

## Quick Start Example

```sh
cd example

# 1. Download HMM marker set
phyling download fungi_odb10

# 2. Identify orthologs and align
phyling align -I pep -o align -m fungi_odb10

# 3. Filter top 20 informative markers (via PhyKIT treeness/RCV)
phyling filter -I align -o filtered_align -n 20

# 4. Construct species tree (Consensus via ASTER/FastTree)
phyling tree -I filtered_align -o tree_out -f
```

---

## Documentation

Detailed documentation, full CLI argument references, and technical deep-dives are available on **[Read the Docs](https://phyling.readthedocs.io)**:

- **[Installation Guide](https://phyling.readthedocs.io/page/installation.html)** — Conda, Pixi, PyPI, and Developer Setup (`pyproject.toml`).
- **[Quick Start Guide](https://phyling.readthedocs.io/page/quickstart.html)** — Complete tutorial on protein (`pep`) and DNA (`cds`) inputs.
- **Command Modules**:
  - **[Download](https://phyling.readthedocs.io/page/usage/download.html)** — Download & manage BUSCO HMM sets.
  - **[Align](https://phyling.readthedocs.io/page/usage/align.html)** — HMM search, CDS translation/back-translation, PyHMMER multithreading, and checkpoints.
  - **[Filter](https://phyling.readthedocs.io/page/usage/filter.html)** — PhyKIT treeness/RCV (`toverr`) ranking and marker selection.
  - **[Tree](https://phyling.readthedocs.io/page/usage/tree.html)** — Consensus (ASTER) vs Concatenation, partition analysis, and tree engines (FastTree, IQ-TREE, RAxML-NG).
- **[Pipeline Architecture & Benchmarks](https://phyling.readthedocs.io/page/architecture.html)** — Algorithmic design & scaling benchmarks.
- **[Nextflow Workflow](https://phyling.readthedocs.io/page/nextflow.html)** — Scalable HPC (SLURM) & Cloud pipeline ([nf_phyling](https://github.com/stajichlab/nf_phyling)).
- **[Citation & References](https://phyling.readthedocs.io/page/citation.html)** — Publication metadata & BibTeX.

---

## Citation

If you use Phyling in your research, please cite:

> **Phyling: Fast, scalable, and user-friendly phylogenomic reconstruction of species phylogenies.**
> *G3: Genes&#124;Genomes&#124;Genetics*, Volume 16, Issue 5, 2026, jkag062.
> DOI: [10.1093/g3journal/jkag062](https://doi.org/10.1093/g3journal/jkag062)
