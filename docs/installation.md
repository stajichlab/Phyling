# Installation

Phyling can be installed via **Conda** (recommended), **Pixi**, a **container image** (Docker/Apptainer), from the **latest release**, or cloned directly from **GitHub**.

---

## 1. Conda Installation (Recommended)

The simplest way to install Phyling and all its system dependencies (including HMMER, IQ-TREE, FastTree, RAxML-NG, and ASTER) is through [Bioconda](https://anaconda.org/bioconda/phyling):

```sh
conda install bioconda::phyling
```

---

## 2. Pixi Installation (Fast Alternative)

[Pixi](https://pixi.sh) is a high-performance package manager built on top of the Conda ecosystem.

### Global CLI Installation

To install `phyling` globally using Pixi:

```sh
pixi global install -c bioconda phyling
```

### Adding to a Pixi Project

To use `phyling` within a local Pixi project environment:

```sh
pixi init my_phyling_project
cd my_phyling_project
pixi add -c bioconda phyling

# Run phyling commands inside the pixi environment
pixi run phyling --help
```

---

## 3. Container Installation (Docker / Apptainer)

A pre-built container image with Phyling and all its dependencies is published to the GitHub Container Registry at [ghcr.io/stajichlab/phyling](https://github.com/stajichlab/Phyling/pkgs/container/phyling). This is a convenient option for HPC clusters or environments where Conda/Pixi is not available.

### Docker

```sh
docker pull ghcr.io/stajichlab/phyling:latest

docker run --rm -v "$PWD:/data" -w /data ghcr.io/stajichlab/phyling:latest phyling --help
```

### Apptainer (Singularity)

```sh
apptainer pull phyling.sif docker://ghcr.io/stajichlab/phyling:latest

apptainer exec phyling.sif phyling --help
```

> **Tip:** Replace `latest` with a specific version tag (e.g. `2.X.X`) for reproducible pipeline runs. Available tags are listed on the [package page](https://github.com/stajichlab/Phyling/pkgs/container/phyling/versions?filters%5Bversion_type%5D=tagged).

---

## 4. Install from Latest Stable Release

You can download the source code from the [latest stable release](https://github.com/stajichlab/Phyling/releases/latest) on GitHub, decompress it, and install via pip:

```sh
tar -zxvf Phyling-2.X.X.tar.gz
cd Phyling-2.X.X
pip install .
```

> **Note:** Some command-line binary dependencies (e.g. IQ-TREE, RAxML-NG, ASTER, FastTree) are not available on PyPI and must be installed separately via Conda or system package managers if installing via pip.

---

## 5. Install from GitHub (Git Clone)

To get the latest development version:

```sh
# Via SSH
git clone git@github.com:stajichlab/Phyling.git

# Via HTTPS
git clone https://github.com/stajichlab/Phyling.git

cd Phyling
pip install -e .
```

---

## 6. Setup Marker Set Storage (`$PHYLING_DB`)

By default, Phyling downloads and saves BUSCO HMM marker sets to `~/.phyling` in the user's home directory.

To specify custom database paths or share marker sets across a team/cluster, set the `$PHYLING_DB` environment variable (colon-separated list of directories):

```sh
export PHYLING_DB=/path/user/phyling_db:/path/group/phyling_db:/path/sys/phyling_db
```

### Search & Storage Behavior
1. Phyling searches for marker sets sequentially across all paths listed in `$PHYLING_DB`.
2. If the marker set is not found locally, Phyling downloads it into the **first path** specified in `$PHYLING_DB` (provided it is writable).
3. If no `$PHYLING_DB` is set or the paths are not writable, Phyling falls back to `~/.phyling`.

> **Conda Tip:** Refer to the [Conda environment variable documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables) to set `$PHYLING_DB` automatically upon environment activation.

---

## 7. Developer Setup

Developers can set up a local development environment using either **Pixi** (recommended, configured via `pyproject.toml`) or **Conda** (`dev_environment.yml`).

### Option A: Using Pixi (Recommended)

Phyling includes native Pixi workspace configuration inside `pyproject.toml` (`[tool.pixi.*]`), defining channels (`conda-forge`, `bioconda`), dependencies, dev environments, editable package installation, and test tasks.

```sh
git clone https://github.com/stajichlab/Phyling.git
cd Phyling

# Install development environment defined in pyproject.toml
pixi install --environment dev

# Activate the development environment shell
pixi shell -e dev

# Set up pre-commit hooks
pre-commit install

# Run automated tests defined in pyproject.toml
pixi run test
```

#### Key Pixi Dev Commands (via `pyproject.toml`):
- **Run Tests**: `pixi run test` (executes `pytest`)
- **Run Linting**: `pixi run -e dev ruff check .` (uses `ruff` installed in dev feature)
- **Interactive Shell**: `pixi shell -e dev` (spawns shell with `phyling` installed in editable mode)

---

### Option B: Using Conda

Alternatively, use `dev_environment.yml` to set up a Conda environment:

```sh
git clone https://github.com/stajichlab/Phyling.git
cd Phyling

# Create development conda environment
conda env create --file dev_environment.yml
conda activate phyling

# Install pre-commit hooks before committing changes
pre-commit install
```

To update an existing Conda development environment:

```sh
conda env update -f dev_environment.yml --prune
```

---

## System Requirements & Software Stack

- **[Python](https://www.python.org/)**: `>= 3.10,< 3.15`
- **[Biopython](https://biopython.org/)**: Sequence parsing and file format support.
- **[PyHMMER](https://pyhmmer.readthedocs.io/)**: Python bindings for HMMER3 for fast profile HMM searches.
- **[pyfaidx](https://github.com/mdshw5/pyfaidx)**: Efficient random-access indexing of FASTA files.
- **[Muscle](https://drive5.com/muscle5/)**: Alternative multiple sequence alignment engine.
- **[ClipKIT](https://jlsteenwyk.com/ClipKIT/)**: Alignment trimming for uninformative/noisy sites.
- **[PhyKIT](https://jlsteenwyk.com/PhyKIT/)**: Computes treeness/RCV (`toverr`) scores for marker filtering.
- **[FastTree](https://github.com/morgannprice/fasttree)**: Fast approximately maximum-likelihood phylogenetic tree inference.
- **[RAxML-NG](https://github.com/amkozlov/raxml-ng)**: Modern maximum-likelihood tree inference tool.
- **[IQ-TREE](http://www.iqtree.org/)**: Advanced maximum-likelihood tree builder with ModelFinder and concordance analysis.
- **[ASTER](https://github.com/chaoszhang/ASTER)** / **[ASTRAL](https://github.com/smirarab/ASTRAL)**: Coalescence-based consensus species tree inference.
