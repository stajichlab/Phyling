# phyling filter

The `filter` module evaluates the phylogenetic informativeness of individual ortholog alignment files using **PhyKIT**, ranks them by their **treeness over RCV (`toverr`)** score, and selects the top N markers for species tree inference.

---

## Command Syntax

```sh
phyling filter -I ALIGN_DIR -n TOP_N [OPTIONS]
# or
phyling filter -i FILE1.mfa FILE2.mfa ... -n TOP_N [OPTIONS]
```

---

## Arguments & Options Reference

### Required Arguments

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| `-i` | `--inputs` | Input Multiple Sequence Alignment (`.mfa`) FASTA files. |
| `-I` | `--input_dir` | Input directory containing alignment FASTA files. |
| `-n` | `--top_n_toverr` | Integer specifying the top N ranked markers to select for tree building. |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-o` | `--output` | `phyling-tree-[TIMESTAMP]` | Output directory for `treeness.tsv` and selected alignment symlinks. |
| | `--seqtype` | `AUTO` | Sequence type (`pep`, `dna`, or `AUTO`). |
| | `--ml` | `False` | Enable maximum-likelihood estimation during initial single-gene tree building. |
| `-t` | `--threads` | `All CPUs` | CPU threads allocated for parallel marker evaluation (default: all available CPUs). |
| `-v` | `--verbose` | `False` | Verbose logging mode. |
| `-h` | `--help` | | Display help message and exit. |

---

## How Treeness/RCV (`toverr`) Filtering Works

In large phylogenomic datasets, not all orthologous markers are equally informative. Some markers suffer from high sequence saturation, missing data, or strong compositional bias across taxa, which can introduce noise or misleading topologies into species tree construction.

Phyling integrates **PhyKIT** to measure two key metrics:
1. **Treeness**: The proportion of overall tree length corresponding to internal branches (higher treeness indicates stronger phylogenetic signal relative to terminal branch noise).
2. **Relative Composition Variability (RCV)**: Measures compositional heterogeneity across taxa (lower RCV indicates lower susceptibility to compositional bias).

```
toverr = Treeness / RCV
```

Markers with higher `toverr` values are ranked at the top and selected for final tree inference.

---

## Output Files

Executing `phyling filter` generates:

1. **`treeness.tsv`**: A tab-separated file listing calculated treeness, RCV, and `toverr` metrics for all evaluated input markers, sorted in descending order of informativeness.
2. **Selected Alignments (`filtered_align/`)**: Symbolic links pointing to the top N selected `.mfa` alignment files.

---

## Usage Examples

### Select Top 20 Most Informative Markers

```sh
phyling filter -I align -o filtered_align -n 20 -t 16
```

### Enable Maximum-Likelihood Estimation (`--ml`)

By default, single-gene trees built during filtering disable maximum-likelihood optimization to maximize speed. To enable ML tree estimation for finer score calculation:

```sh
phyling filter -I align -o filtered_align -n 20 --ml -t 16
```

---

## Handling Missing Samples Warnings

When filtering for the top N markers, certain sample genomes with high rates of missing genes may lack ortholog representations in the top N set.

If any sample is completely absent from the top N alignments, `phyling filter` will emit a **warning message**:

```
[WARNING] Sample 'Genome_X' has no ortholog markers in the selected top 20 alignments!
```

> **Tip:** If maintaining representation for all sample taxa is critical, consider increasing `-n / --top_n_toverr` until all taxa are represented in at least one marker.
