[![CI/build and test](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml)
[![CI/Conda build and test](https://github.com/stajichlab/Phyling/actions/workflows/conda_build_and_test.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/conda_build_and_test.yml)
[![Python](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12_%7C_3.13_%7C_3.14-blue?logo=python)](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml)
[![codecov](https://codecov.io/gh/stajichlab/Phyling/graph/badge.svg?token=ZH5GBQYKZ6)](https://codecov.io/gh/stajichlab/Phyling)
[![License](https://img.shields.io/github/license/stajichlab/Phyling?label=license)](https://github.com/stajichlab/Phyling/blob/main/LICENSE)
[![Conda](https://anaconda.org/bioconda/phyling/badges/version.svg)](https://anaconda.org/bioconda/phyling)
[![Container](https://ghcr-badge.egpl.dev/stajichlab/phyling/latest_tag?trim=major&label=image)](https://github.com/stajichlab/Phyling/pkgs/container/phyling/versions?filters%5Bversion_type%5D=tagged)
[![DOI](https://img.shields.io/badge/DOI-10.1093/g3journal/jkag062-blue)](https://doi.org/10.1093/g3journal/jkag062)

# Phyling

**Phyling** is a fast, scalable, and user-friendly software pipeline for phylogenomic reconstruction of species phylogenies directly from protein-encoded or nucleotide genomic data.

Published in *G3: Genes&#124;Genomes&#124;Genetics*: [jkag062 (2026)](https://doi.org/10.1093/g3journal/jkag062).

---

## Pipeline Flowchart

<p align="center">
  <img src="misc/phyling_flowchart-light.svg#gh-light-mode-only" alt="Phyling flowchart" width="800">
  <img src="misc/phyling_flowchart-dark.svg#gh-dark-mode-only" alt="Phyling flowchart" width="800">
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

*For alternative installation options (pip, Git clone, developer setup), see the [Installation Guide](https://github.com/stajichlab/Phyling/wiki/Installation).*

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

## Documentation & Wiki

Detailed documentation, full CLI argument references, and technical deep-dives are available on the **[Phyling GitHub Wiki](https://github.com/stajichlab/Phyling/wiki)**:

- **[Installation Guide](https://github.com/stajichlab/Phyling/wiki/Installation)** — Conda, Pixi, PyPI, and Developer Setup (`pyproject.toml`).
- **[Quick Start Guide](https://github.com/stajichlab/Phyling/wiki/Quick-Start)** — Complete tutorial on protein (`pep`) and DNA (`cds`) inputs.
- **Command Modules**:
  - **[Download](https://github.com/stajichlab/Phyling/wiki/Download)** — Download & manage BUSCO HMM sets.
  - **[Align](https://github.com/stajichlab/Phyling/wiki/Align)** — HMM search, CDS translation/back-translation, PyHMMER multithreading, and checkpoints.
  - **[Filter](https://github.com/stajichlab/Phyling/wiki/Filter)** — PhyKIT treeness/RCV (`toverr`) ranking and marker selection.
  - **[Tree](https://github.com/stajichlab/Phyling/wiki/Tree)** — Consensus (ASTER) vs Concatenation, partition analysis, and tree engines (FastTree, IQ-TREE, RAxML-NG).
- **[Pipeline Architecture & Benchmarks](https://github.com/stajichlab/Phyling/wiki/Pipeline-Architecture)** — Algorithmic design & scaling benchmarks.
- **[Nextflow Workflow](https://github.com/stajichlab/Phyling/wiki/Nextflow-Workflow)** — Scalable HPC (SLURM) & Cloud pipeline ([nf_phyling](https://github.com/stajichlab/nf_phyling)).
- **[Citation & References](https://github.com/stajichlab/Phyling/wiki/Citation-and-References)** — Publication metadata & BibTeX.

---

## Citation

If you use Phyling in your research, please cite:

> **Phyling: Fast, scalable, and user-friendly phylogenomic reconstruction of species phylogenies.**
> *G3: Genes&#124;Genomes&#124;Genetics*, Volume 16, Issue 5, 2026, jkag062.
> DOI: [10.1093/g3journal/jkag062](https://doi.org/10.1093/g3journal/jkag062)
