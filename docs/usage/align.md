# phyling align

The `align` module identifies single-copy orthologous markers across sample genomes using profile HMM searches (`hmmsearch`), extracts matching sequences, and aligns them into Multiple Sequence Alignment (MSA) files.

---

## Command Syntax

```sh
phyling align -i FILE [FILE ...] -m MARKERSET [OPTIONS]
# or
phyling align -I DIRECTORY -m MARKERSET [OPTIONS]
```

---

## Arguments & Options Reference

### Required Arguments

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| `-i` | `--inputs` | Input peptide/CDS FASTA files (plain or `.gz`). Space-separated list. |
| `-I` | `--input_dir` | Directory containing query peptide/CDS FASTA files (plain or `.gz`). |
| `-m` | `--markerset` | Name or directory path of the HMM marker set (e.g., `fungi_odb10`). |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-o` | `--output` | `phyling-align-[TIMESTAMP]` | Output directory path for alignment files and checkpoint data. |
| | `--seqtype` | `AUTO` | Input sequence type (`pep`, `dna`, or `AUTO`). Manually set to `dna` if IUPAC ambiguity codes break autodetection. |
| `-E` | `--evalue` | `1e-10` | `hmmsearch` E-value threshold (used when bitscore cutoff file is not available in HMM set). |
| `-M` | `--method` | `hmmalign` | Alignment engine (`hmmalign` or `muscle`). |
| | `--non_trim` | `False` | Retain non-trimmed alignment columns (do not trim unaligned regions). |
| `-t` | `--threads` | `All CPUs` | Total CPU threads to allocate for parallel HMM search and alignment (default: all available CPUs). |
| `-v` | `--verbose` | `False` | Verbose debug output logging. |
| `-h` | `--help` | | Display help message and exit. |

---

## Technical Details & Features

### 1. Ortholog Identification Thresholds

- **Sample Frequency Cutoff**: HMM profiles matching at least **4 sample genomes** are identified as orthologs and saved for downstream alignment.
- **Bitscore vs E-value Cutoffs**: Phyling checks for a bitscore cutoff file in the HMM marker directory. If present, bitscore cutoffs are applied. If absent, the threshold falls back to `-E/--evalue` (default: `1e-10`).

---

### 2. Protein vs DNA (CDS) Processing

Phyling natively supports both protein peptide sequences and DNA coding sequences (CDS):

```
[DNA CDS Input] ──► Auto-Translate to Protein ──► Profile HMM Search ──► MSA Alignment ──► Back-Translate to DNA ──► [DNA MSA Output]
```

- **Autodetection (`--seqtype AUTO`)**: Phyling inspects initial sequence characters to determine whether inputs are DNA or protein.
- **Explicit Override (`--seqtype dna`)**: DNA sequences containing degenerate [IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html) can cause autodetection failures. Override manually by specifying `--seqtype dna`.

```sh
phyling align -I cds -o align_cds -m fungi_odb10 -t 16 --seqtype dna
```

---

### 3. Multithreading Strategy

Phyling utilizes **PyHMMER** to parallelize C-level HMMER routines:

- **Threads < 8**: `hmmsearch` runs sequentially per sample, allocating all available CPU threads to each search job.
- **Threads >= 8**: `hmmsearch` runs in parallel blocks. Every parallel job receives **4 threads**. For example, `-t 16` launches 4 parallel HMM search jobs running simultaneously across 4 threads each.
- **MSA Alignment Phase**: Launches parallel single-threaded alignment tasks matching total thread allocation (`-t`).

---

### 4. Incremental Checkpointing System

When `phyling align` completes, it writes a checkpoint file inside the designated output folder (`-o/--output`).

#### Key Checkpoint Features:
- **Resuming Interrupted Runs**: If an alignment run is interrupted, re-running with the same output directory reads the checkpoint and skips already completed samples.
- **Adding/Removing Genomes**:
  - If new genome FASTA files are added to `-i/-I`, HMM search runs *only* on the newly added samples.
  - If genomes are omitted, Phyling removes them from the ortholog matrix without re-searching remaining samples.
- **Automatic Invalidation**: Changing critical parameters (HMM marker set, `--seqtype`, or `--evalue`) invalidates the checkpoint and triggers a full rerun.

#### Example: Adding a Sample via Checkpoint

```sh
# Initial run with 4 samples
phyling align -i sample1.aa.fa sample2.aa.fa sample3.aa.fa sample4.aa.fa -o align_out -m fungi_odb10 -t 16

# Add sample5 to the same output directory (sample1-4 are skipped automatically)
phyling align -i sample1.aa.fa sample2.aa.fa sample3.aa.fa sample4.aa.fa sample5.aa.fa -o align_out -m fungi_odb10 -t 16
```
