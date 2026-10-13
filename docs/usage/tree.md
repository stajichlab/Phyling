# phyling tree

The `tree` module constructs the final species phylogenetic tree from selected multiple sequence alignment (`.mfa`) files, supporting both **consensus (coalescence-based)** and **concatenation (supermatrix-based)** strategies.

---

## Command Syntax

```sh
phyling tree -I FILTERED_ALIGN_DIR [OPTIONS]
# or
phyling tree -i FILE1.mfa FILE2.mfa ... [OPTIONS]
```

---

## Arguments & Options Reference

### Required Arguments

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| `-i` | `--inputs` | Selected Multiple Sequence Alignment (`.mfa`) files. |
| `-I` | `--input_dir` | Directory containing selected alignment files (e.g. output of `phyling filter`). |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-o` | `--output` | `phyling-tree-[TIMESTAMP]` | Output directory path for tree files and figures. |
| | `--seqtype` | `AUTO` | Input data sequence type (`pep`, `dna`, or `AUTO`). |
| `-M` | `--method` | `ft` | Tree inference engine: `ft` (FastTree), `raxml` (RAxML-NG), or `iqtree` (IQ-TREE). |
| `-c` | `--concat` | `False` | Use concatenated supermatrix strategy instead of consensus trees. |
| `-p` | `--partition` | `False` | Enable partitioned analysis by sequence marker (only valid with `-c` in IQ-TREE or RAxML-NG). |
| `-f` | `--figure` | `False` | Render a simple tree visualization (`.png`/`.pdf`) via Matplotlib. |
| | `--seed` | `-1` | Random seed integer for reproducible inferences (`-1` for random seed). |
| `-t` | `--threads` | `All CPUs` | CPU threads for parallel tree construction (default: all available CPUs). |
| `-v` | `--verbose` | `False` | Verbose debug output logging. |
| `-h` | `--help` | | Display help message and exit. |

---

## Tree Construction Strategies

### 1. Consensus Tree Strategy (Default)

In consensus mode, single-gene phylogenetic trees are inferred individually for each ortholog marker, and a coalescence-based species tree is computed across all gene trees using **ASTER / ASTRAL**.

```sh
phyling tree -I filtered_align -o consensus_tree -t 16
```

> **Note on Substitution Models:**
> All tree-building tools offer more than a single substitution model for phylogenetic inference. To simplify model selection, Phyling applies the `GTR` model for DNA sequences and the `LG` model for peptide data in consensus mode. In concatenation mode, Phyling follows a series of steps to select the best-suited model(s) via ModelFinder (with IQ-TREE).

---

### 2. Concatenation Strategy (`-c / --concat`)

In concatenation mode, all selected marker alignments are concatenated end-to-end into a single supermatrix alignment (`concat_alignments.mfa`), and a single species tree is inferred from the combined matrix.

```sh
phyling tree -I filtered_align -o concat_tree -c -t 16
```

#### Automatic Model Selection & Support Calculation:
When `--concat` is enabled with IQ-TREE:
1. Phyling invokes **ModelFinder** to determine the optimal substitution model for the dataset.
2. Infers the species tree using the selected optimal model.
3. Computes branch support using **Ultrafast Bootstrap (UFBoot)** and **Site Concordance Factor (sCF)**.

---

### 3. Partitioned Analysis Mode (`-p / --partition`)

When using `--concat` with **IQ-TREE** or **RAxML-NG**, enabling `-p` performs a partitioned analysis where substitution models and rates are estimated independently for each marker partition.

```sh
phyling tree -I filtered_align -M iqtree -o partition_tree -c -p -t 16
```

> **Note:** Partition mode is supported in **IQ-TREE** and **RAxML-NG**, but is **not** supported by FastTree (`ft`).

---

## Supported Tree Inference Engines

| Code | Engine | Description | Partition Support |
| :--- | :--- | :--- | :--- |
| `ft` *(default)* | **FastTree** | Fast approximately maximum-likelihood inference engine. Ideal for rapid exploratory analyses. | No |
| `raxml` | **RAxML-NG** | Highly accurate maximum-likelihood tree engine with model evaluation. | Yes |
| `iqtree` | **IQ-TREE** | Advanced maximum-likelihood tree engine featuring ModelFinder, UFBoot, and sCF concordance metrics. | Yes |

---

## Output Visualization & Seed Reproducibility

### Graphical Tree Figures (`-f / --figure`)

Specify `-f` to render a simple tree visualization using Matplotlib:

```sh
phyling tree -I filtered_align -o tree_output -f -t 8
```

### Reproducible Inferences (`--seed`)

To ensure identical tree topology across multiple runs, supply a fixed integer seed:

```sh
phyling tree -I filtered_align -o tree_output --seed 42 -t 8
```
