[![CI/build and test](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml)
[![CI/Conda build and test](https://github.com/stajichlab/Phyling/actions/workflows/conda_build_and_test.yml/badge.svg?branch=main)](https://github.com/stajichlab/Phyling/actions/workflows/conda_build_and_test.yml)
[![Python](https://img.shields.io/badge/python-3.9_%7C_3.10_%7C_3.11_%7C_3.12_%7C_3.13-blue?logo=python)](https://github.com/stajichlab/Phyling/actions/workflows/build_and_test.yml)
[![codecov](https://codecov.io/gh/stajichlab/Phyling/graph/badge.svg?token=ZH5GBQYKZ6)](https://codecov.io/gh/stajichlab/Phyling)
[![License](https://img.shields.io/github/license/stajichlab/Phyling?label=license)](https://github.com/stajichlab/Phyling/blob/main/LICENSE)
[![Conda](https://anaconda.org/bioconda/phyling/badges/version.svg)][Bioconda]
[![Last updated](https://anaconda.org/bioconda/phyling/badges/latest_release_date.svg)][Bioconda]
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.07.30.666921-blue)](https://www.biorxiv.org/content/10.1101/2025.07.30.666921)


# Phyling tool

Phyling is a fast, scalable, and user-friendly tool supporting phylogenomic reconstruction of species phylogenies directly from
protein-encoded genomic data. It identifies orthologous genes by searching a sample's protein sequences against a Hidden Markov
Models marker set, containing single-copy orthologs, retrieved from the [BUSCO database][Busco]. In the final step, users can
choose between consensus and concatenation strategies to construct the species tree from the aligned orthologs.

### Flow chart

<img src= "misc/phyling_flowchart-light.svg#gh-light-mode-only" alt="Phyling flowchart" width="800">
<img src= "misc/phyling_flowchart-dark.svg#gh-dark-mode-only" alt="Phyling flowchart" width="800">

## Usage

First of all, install the package following the [instruction](#install) below.

Phyling is a package to extract phylogenomic markers and build a phylogenetic tree upon them. It comprises 4 modules - download,
align, filter and tree. Use `phyling --help` to see more details.

```
Modules:
  {download,align,filter,tree}
    download            Download HMM markers
    align               Run multiple sequence alignments against orthologs found among samples
    filter              Filter the multiple sequence alignment results
    tree                Build a phylogenetic tree based on selected multiple sequence alignment results

Options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
```

To test run on the example files, please `cd` into the folder `example`.

```sh
cd example
```

In general, Phyling takes fasta as input. The bgzipped fasta is also valid.

The folder `example/pep` includes 5 example peptide fasta which can be used for test run.

In addition to the peptide sequences, Phyling can also takes DNA coding sequences as inputs to more accurately estimate the
phylogeny of closely related species. When taking DNA coding sequences as inputs, DNA sequences will be translated into peptide
sequences and all the hmmsearch/align are done on the peptide version. The final MSA results will be back-translated into DNA at
the final stage. The example DNA sequences are placed under the folder `example/cds`.

### Download HMM marker sets

The download module is used to download HMM marker sets from BUSCO website. (Currently is updated to [v5][Busco]) See all options
with `phyling download --help`.

```
positional arguments:
  HMM markerset or "list"
                        Name of the HMM markerset

Options:
  -v, --verbose         Verbose mode for debug
  -h, --help            show this help message and exit
```

Firstly, use `download list` to show the available BUSCO marker sets.

```sh
phyling download list
```

By default the downloaded marker sets will be saved to the `~/.phyling` or the first path in `$PHYLING_DB` (if have been set) if
it is writable. The **Datasets available online** section lists all marker sets that are available on the [BUSCO][Busco] website.

And the **Datasets available on local** section lists all marker sets that have already been downloaded.

To download the marker set, copy the name from the list and paste it to the download module directly. Here we use `fungi_odb10` as
example.

```sh
phyling download fungi_odb10
```

The download module will automatically check for updates to marker sets online each time it runs. Local marker sets which have
available updates online will be marked as _\[Outdated\]_. You can rerun `phyling download [marker set]` to update the local files.

### Find the orthologs and align them

The align module identify the orthologs among all the samples using _hmmsearch_. HMM profiles that have matches on more than 4
samples are considered **orthologs**.

Before conducting _hmmsearch_, the module will first search for the bitscore cutoff file within the root HMM folder. If the cutoff
file is not found, the reporting threshold for _hmmsearch_ will be determined based on the `-E/-evalue` (default is 1e-10).

Once the orthologs are identified, the sequences extracted from each sample undergo multiple sequence alignment. By default, the
alignment is performed using the _hmmalign_ method. However, users have the option to switch to _muscle_ by specifying the
`-M/--method muscle` flag.

Finally, each alignment result is output separately. You can decide whether you want to filter it with [treeness/RCV] score or use
them all for tree building. Please check out the filter command through `phyling filter --help`

```
Required arguments:
  -i file [files ...], --inputs file [files ...]
                        Query pepetide/cds fasta or gzipped fasta
  -I directory, --input_dir directory
                        Directory containing query pepetide/cds fasta or gzipped fasta
  -m directory, --markerset directory
                        Directory of the HMM markerset

Options:
  -o directory, --output directory
                        Output directory of the alignment results (default: phyling-align-[YYYYMMDD-HHMMSS] (UTC timestamp))
  --seqtype {dna,pep,AUTO}
                        Input data sequence type (default: AUTO)
  -E float, --evalue float
                        Hmmsearch reporting threshold (default: 1e-10, only being used when bitscore cutoff file is not available)
  -M {hmmalign,muscle}, --method {hmmalign,muscle}
                        Program used for multiple sequence alignment (default: hmmalign)
  --non_trim            Report non-trimmed alignment results
  -t THREADS, --threads THREADS
                        Threads for hmmsearch and the number of parallelized jobs in MSA step. Better be multiple of 4 if using
                        more than 8 threads (default: 8)
  -v, --verbose         Verbose mode for debug
  -h, --help            show this help message and exit
```

Run the align module with all the fasta files under folder `pep`.

```sh
phyling align -I pep -o align -m fungi_odb10
```

An equivalent way to send inputs.

```sh
phyling align -i pep/*.fasta.gz -o align -m fungi_odb10
```

Or if you're just interested in part of the fasta, you can specify the inputs one-by-one.

```sh
phyling align -i pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  pep/Zygorhynchus_heterogamous_NRRL_1489.aa.fasta.gz \
  -o align \
  -m fungi_odb10
```

**Note: Required at least 4 samples to build a tree!**

Accelerate by using 16 cpus.

```sh
phyling align -I pep -o align -m fungi_odb10 -t 16
```

#### Multithreading strategy

According to [pyhmmer benchmark](https://pyhmmer.readthedocs.io/en/stable/benchmarks.html), the acceleration benefits from
multithreading drop significantly as more CPUs are utilized. When less then 8 cpus are given, the hmmsearch step will run on
single-thread manner and all cpus will be used for each round of hmmsearch. When 8 or more cpus are given, the hmmsearch step will
use 4 cpus for each parallel job. In the example above, 4 hmmsearch jobs will run parallelly and each job utilize 4 cpus. For the
alignment step, 16 parallel jobs will be launched and each parallel job is running on single-thread manner.

In addition to the hmmsearch step, the assigned cpus will also be used to accelerate the alignment steps. Highly recommended to
enable the multithreading option.

#### Use coding sequences instead of peptide sequences

In some circumstances, the highly shared peptide sequences make it difficult to resolve the relationship among closely related
species. To address the issue, one can use DNA coding sequences (CDS), which can contain more recent evolutionary changes, instead of peptide
sequences for phylogeny analysis.

Run the align module with cds fasta files under folder `cds`.

```sh
phyling align -I cds -o align_cds -m fungi_odb10 -t 16
```

The CDS inputs will be translated into peptide sequences in the first steps. The translated peptide sequences will be used for
hmmsearch and the alignment steps. The peptide alignment results will then being back-translated according to the original CDS
inputs. And the back-translated DNA version alignments will be output.

By default, Phyling automatically detects the sequence type based on the input data. However, this autodetection can occasionally
fail, particularly with DNA sequences that contain [IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html). If the
autodetection is incorrect or fails, users can override it by specifying the sequence type manually using the `--seqtype` option.

```sh
phyling align -I cds -o align_cds -m fungi_odb10 -t 16 --seqtype dna
```

#### Checkpoint for quick rerun

Once the align module complete, a checkpoint file will be generated to the output folder. This checkpoint file stores the
parameters, samples and identified orthologs, which will be loaded to the pipeline when rerun with the same output folder. Then
the align module will determine whether to skip the hmmsearch on some of the samples that were already completed in the previous
run.

For example, we first run the align module by:

```sh
phyling align -i pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  pep/Zygorhynchus_heterogamous_NRRL_1489.aa.fasta.gz \
  -o align \
  -m fungi_odb10 \
  -t 16
```

And later we want to add `Actinomucor elegans` to the analysis but kick out `Zygorhynchus heterogamous`. We can start another run
and specifying the same output folder:

```sh
phyling align -i pep/Actinomucor_elegans_CBS_100.09.aa.fasta.gz \
  pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  -o align \
  -m fungi_odb10 \
  -t 16
```

In this case, `Pilobolus umbonatus`, `Rhizopus homothallicus` and `Rhizopus rhizopodiformis` will be skipped from the hmmsearch
process since they have already been searched in the previous run. The `Actinomucor elegans` is the only sample need to be
hmmsearched. On the other hand, the `Zygorhynchus heterogamous` will be removed from the current run.

Note that if the content of input files have changes, they will also being detected by align module and trigger the rerun. If the
hmmsearch evalue/bitscore cutoff is changed, all the samples will need to rerun the hmmsearch step. Also, the changes on HMM
marker sets or input samples with different seqtype will terminate the align module. (since this case should be considered an
entirely different analysis)

### Filter the alignment results

The align module reports a lot of orthologs depending on the size of BUSCO dataset and the gene set homogeneity among samples.
However, not every marker is equally informative in resolving phylogeny; some substitutions may not contribute significantly to
the branches in the phylogenetic tree. To address this issue, Phyling incorporates [PhyKIT] to compute the [treeness/RCV] scores
(toverr) and ranks the markers accordingly. Higher ranks indicate greater informativeness and lower susceptibility to composition
bias and will be selected for final tree building (thru _consensus_ or _concatenate_ strategy).

```
Required arguments:
  -i file [files ...], --inputs file [files ...]
                        Multiple sequence alignment fasta of the markers
  -I directory, --input_dir directory
                        Directory containing multiple sequence alignment fasta of the markers
  -n TOP_N_TOVERR, --top_n_toverr TOP_N_TOVERR
                        Select the top n markers based on their treeness/RCV for final tree building

Options:
  -o directory, --output directory
                        Output directory of the treeness.tsv and selected MSAs (default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))
  --seqtype {pep,dna,AUTO}
                        Input data sequence type (default: AUTO)
  --ml                  Use maximum-likelihood estimation during tree building
  -t THREADS, --threads THREADS
                        Threads for filtering (default: 6)
  -v, --verbose         Verbose mode for debug
  -h, --help            show this help message and exit
```

For example, we pick only the top 20 markers by specifying `-n/--top_n_toverr`:

```sh
phyling filter -I align -o filtered_align -n 20 -t 16
```

[FastTree] will be used to build a tree for all markers and the toverr will be computed on these trees. To accelerate the process,
the maximum-likelihood estimation is turned off by default and users can enable it by specifying `--ml`. A file `treeness.tsv`
recording the toverr of all markers and the symlinks to the mfa of these selected markers will be output.

Note that some of the samples may be entirely excluded after filtering by markers (e.g., if they lack orthologs in the top n
selected markers). When missing samples are detected, the filter module will issue a warning. If the presence of all samples is
important, consider adjusting the `-n/--top_n_toverr` argument accordingly.

Again, the filter module automatically detects the sequence type of the input data. Try manually assign it by `--seqtype` if it
fails. (The case below are peptide fasta)

```sh
phyling filter -I align -o filtered_align -n 20 -t 16 --seqtype pep
```

### Build tree from multiple sequence alignment results

Finally, we can run the tree module, use the multiple sequence alignment results to build a phylogenetic tree. By default, it uses
the _consensus tree_ strategy (conclude the majority of trees which was built upon each single gene) but you can choose to use
_concatenated alignment_ strategy by specifying `-c/--concat`. Currently, 3 methods ([FastTree], [RAxML-NG] and [IQ-TREE]) are
available for tree building. You can choose your own preferred method by specifying `-M/--method`. (default is FastTree) See all
the options with `phyling tree --help`.

```
Required arguments:
  -i file [files ...], --inputs file [files ...]
                        Multiple sequence alignment fasta of the markers
  -I directory, --input_dir directory
                        Directory containing multiple sequence alignment fasta of the markers

Options:
  -o directory, --output directory
                        Output directory of the newick treefile (default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))
  --seqtype {pep,dna,AUTO}
                        Input data sequence type (default: AUTO)
  -M {ft,raxml,iqtree}, --method {ft,raxml,iqtree}
                        Algorithm used for tree building. (default: ft)
                        Available options:
                        ft: FastTree
                        raxml: RAxML-NG
                        iqtree: IQTree
  -c, --concat          Concatenated alignment results
  -p, --partition       Partitioned analysis by sequence. Only works when --concat enabled.
  -f, --figure          Generate a matplotlib tree figure
  --seed SEED           Seed number for stochastic elements during inferences. (default: -1 to generate randomly)
  -t THREADS, --threads THREADS
                        Threads for tree construction (default: 6)
  -v, --verbose         Verbose mode for debug
  -h, --help            show this help message and exit
```

Run the tree module with all the alignment results under folder `filtered_align` got from the filter module.

```sh
phyling tree -I filtered_align
```

You can also use only a part of the alignment results to build the tree.

```sh
phyling tree -i filtered_align/100957at4751.aa.mfa filtered_align/174653at4751.aa.mfa filtered_align/255412at4751.aa.mfa
```

By default, the tree module creates an output folder named using the current UTC timestamp. You can specify a output path using
the `--output/-o` option.

```sh
phyling tree -I filtered_align -o tree
```

Use [IQ-TREE] instead of the default [FastTree] method for tree building and run with 16 threads.

```sh
phyling tree -I filtered_align -o tree -m iqtree -t 16
```

Manually assign sequence type.

```sh
phyling tree -I filtered_align -o tree -m iqtree -t 16 --seqtype pep
```

Use matplotlib to generate a tree figure.

```sh
phyling tree -I filtered_align -o tree -f -t 16
```

Many tree-building algorithms include stochastic steps that rely on random number generation to determine initial states. To
ensure reproducibility and minimize randomness, users can specify a fixed seed using the `--seed` option.

```sh
phyling tree -I filtered_align -o tree --seed 12345 -t 16
```

#### Use concatenate strategy

The consensus strategy is used for tree construction by default but users can choose to concatenate the markers and generate a
single tree on it.

```sh
phyling tree -I filtered_align -o concat_tree -c -t 16
```

When concatenation mode is enabled, Phyling concatenates all the inputs to a `concat_alignments.mfa` file and use that to generate
a single final tree.

Meanwhile, users can construct tree with a more sophisticated **partition mode** when using [RAxML-NG] and [IQ-TREE]. In general,
the partition mode expects different genes exhibit different evolutionary rates, which should be estimated with different models.

The example below concatenates all the markers and run tree building with partitioning through IQ-TREE.

```sh
phyling tree -I align -M iqtree -o partition_tree -c -p -t 16
```

**Note: Partition mode is not supported in FastTree.**

#### Model selection

All tree-building tools offer more than a single substitution model for phylogenetic inference. To simplify model selection, we
use the GTR model for DNA and the LG model for peptide data in consensus mode. In concatenation mode, Phyling follows a series of
steps to achieve the best possible results:

1. Find the best suited model(s) through the ModelFinder module from [IQ-TREE].
2. Use the best model(s) to build tree with either [FastTree], [RAxML-NG] or [IQ-TREE].
3. Calculate the branch support by both Ultrafast bootstrap and site concordance factor by [IQ-TREE], using the parameters
   retrieve from the previous step.

## Requirements

- Python >= 3.9
- [Biopython](https://biopython.org/)
- [PyHMMER], a HMMER3 implementation on python3.
- [pyfaidx](https://github.com/mdshw5/pyfaidx) for quickly retrieving fasta sequences.
- [Muscle] for alternative method for multiple sequence alignment. (Optional)
- [ClipKIT] for removing sites that are poor of phylogenetic signal.
- [PhyKIT] for calculating treeness/RCV to filter uninformative orthologs.
- [FastTree], a faster version of FastTree, use approximately maximum-likelihood to build trees.
- [RAxML-NG], a more sophisticated maximum-likelihood-based tree building tool. (Optional)
- [IQ-TREE], a modern maximum-likelihood-based tool for tree building.
- [ASTER], a C++ re-implementation of [ASTRAL] to resolve consensus among trees.

## Install

### through Conda (recommended)

The most easiest way to install Phyling is through Conda.

```sh
conda install bioconda::phyling
```

### from latest stable release

Alternatively, you can download the source code from the [latest stable
release](https://github.com/stajichlab/Phyling/releases/latest) and decompress it.

Go into the Phyling folder and install the package through pip.

```sh
cd Phyling-2.X.X
pip install .
```

However, some of the Phyling dependencies do not available on pip, so you need to install them manually or through Conda as well.

### git clone the project

You can also directly `git clone` the main or other branches.

ssh
```sh
git clone git@github.com:stajichlab/Phyling.git
```

or https
```sh
git clone https://github.com/stajichlab/Phyling.git
```

Note that the source code on main and other branches may still under development. Please use it cautiously.

### Setup folder for marker sets

By default, Phyling creates a `.phyling` folder under users' home directory to save marker sets downloaded from [BUSCO] on the
first launch. To avoid cluttering the home directory, you can use symlink or assign a different location by setting the
`$PHYLING_DB` environment variable.

```sh
export PHYLING_DB=/path/user/phyling_db:/path/group/phyling_db:/path/sys/phyling_db
```

If you installed Phyling through Conda, please refer to [this section in Conda
documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables)
to see how to set environment variables during activating a Conda environment.

Similar to `$PATH`, `$PHYLING_DB` allows Phyling to search for marker sets across multiple specified paths when running align
module. If the specified marker set is not found in the first path, phyling will keep looking for the rest of paths listed in
`$PHYLING_DB`. However, the first path should be writable by user; otherwise, Phyling will still create a `.phyling` folder in the
home directory to store marker sets.

### Install Phyling for developing (developer only)

Developer should clone the GitHub project directly instead of downloading from the releases. Some of the files for developing
purpose only are not included in the releases.

For convenience we provide a Conda env file. Please use the dev_environment.yml to setup a conda environment and install the
required packages and Phyling.

```sh
conda env create --file dev_environment.yml
```

or use the command to update the existing environment.

```sh
conda env update -f dev_environment.yml --prune
```

Before making commit. Please go into the Phyling folder and set up the [pre-commit] first.

```sh
cd Phyling
pre-commit install
```

[Busco]: https://busco-data.ezlab.org/v5/data/lineages/
[PyHMMER]: https://pyhmmer.readthedocs.io/en/stable/index.html
[Muscle]: https://drive5.com/muscle5/
[ClipKIT]: https://jlsteenwyk.com/ClipKIT/
[PhyKIT]: https://jlsteenwyk.com/PhyKIT/
[Treeness/RCV]: https://jlsteenwyk.com/PhyKIT/usage/index.html#treeness-over-rcv
[FastTree]: https://github.com/morgannprice/fasttree
[RAxML-NG]: https://github.com/amkozlov/raxml-ng
[IQ-TREE]: http://www.iqtree.org/
[ASTER]: https://github.com/chaoszhang/ASTER
[ASTRAL]: https://github.com/smirarab/ASTRAL
[pre-commit]: https://pre-commit.com/
[Bioconda]: https://anaconda.org/bioconda/phyling
