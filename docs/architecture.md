# Pipeline Architecture & Benchmarks

This page details the underlying methodology, algorithmic design, and benchmark performance of **Phyling**, as published in *G3: Genes&#124;Genomes&#124;Genetics* (DOI: [10.1093/g3journal/jkag062](https://doi.org/10.1093/g3journal/jkag062)).

---

## Background & Motivation

Phylogenomic species tree reconstruction traditionally relies on identifying orthologous gene families across target genomes. Established tools like **OrthoFinder** rely on all-vs-all reciprocal best-hit (RBH) sequence comparisons.

While accurate, all-vs-all blast/diamond searches scale **quadratically (O(N²))** with the number of input species (N). As genomic datasets grow to hundreds or thousands of species, RBH approaches become computationally prohibitive in terms of CPU hours and memory usage.

**Phyling** resolves this scalability bottleneck by adopting a **Profile HMM Marker Search** paradigm:

```
Traditional RBH (OrthoFinder):   Species A ◄── All-vs-All ──► Species B   [O(N²) Complexity]

Phyling Profile HMM Search:       Species A/B/C/D ──► BUSCO HMM Set        [O(N) Complexity]
```

By searching input proteomes directly against pre-curated single-copy ortholog HMM profiles from the [BUSCO database](https://busco-data.ezlab.org/v5/data/lineages/), search time scales **linearly (O(N))** with sample count.

---

## Benchmark Highlights & Performance

In extensive comparative benchmarks against OrthoFinder and GToTree across fungal and eukaryotic genomic datasets:

1. **Runtime Efficiency**: Phyling completes ortholog extraction and tree construction substantially faster than OrthoFinder, achieving an order-of-magnitude reduction in total CPU time.
2. **Memory Footprint**: By utilizing PyHMMER in-memory search structures, memory usage remains lightweight even when scaling to large proteomic datasets.
3. **Topology Accuracy & Congruence**: Phylogenetic trees constructed by Phyling exhibit high topological congruence with trees derived from full OrthoFinder orthogroups, confirming that curated BUSCO profile HMMs effectively capture core phylogenomic signal.

---

## Modular Pipeline Breakdown

```{graphviz}
:align: center
:caption: Phyling module workflow

digraph phyling {
    graph [rankdir=TB, fontname="Helvetica", fontsize=12, nodesep=0.4, ranksep=0.45, compound=true];
    node  [shape=box, style="rounded,filled", fillcolor="#f5f7fa", color="#5b6b7f", fontname="Helvetica", fontsize=11];
    edge  [color="#5b6b7f", fontname="Helvetica", fontsize=10];

    subgraph cluster_align {
        label="Align module"; style="rounded,dashed"; color="#8a99ad";
        InputData     [label="Input data type", shape=diamond, fillcolor="#e8eef7"];
        Translate     [label="Translate CDS to protein"];
        HMMSearch     [label="Profile HMM search via PyHMMER"];
        MSA           [label="Multiple sequence alignment\n(hmmalign / muscle)"];
        IsCDS         [label="Was input CDS?", shape=diamond, fillcolor="#e8eef7"];
        BackTranslate [label="Back-translate protein MSA to DNA"];
        OutputMSA     [label="Generate MSA (.mfa)\n& checkpoint file"];

        InputData -> Translate [label="DNA CDS"];
        InputData -> HMMSearch [label="Protein"];
        Translate -> HMMSearch;
        HMMSearch -> MSA;
        MSA -> IsCDS;
        IsCDS -> BackTranslate [label="Yes"];
        IsCDS -> OutputMSA [label="No"];
        BackTranslate -> OutputMSA;
    }

    subgraph cluster_filter {
        label="Filter module"; style="rounded,dashed"; color="#8a99ad";
        F1 [label="Compute treeness & RCV via PhyKIT"];
        F2 [label="Rank markers by\ntoverr = treeness / RCV"];
        F3 [label="Select top N informative markers"];
        F1 -> F2 -> F3;
    }

    subgraph cluster_tree {
        label="Tree module"; style="rounded,dashed"; color="#8a99ad";
        Strategy   [label="Tree strategy", shape=diamond, fillcolor="#e8eef7"];
        GeneTrees  [label="Build individual gene trees\n(FastTree)"];
        Consensus  [label="Consensus analysis\n(ASTER / ASTRAL)"];
        Concat     [label="Concatenate alignments\nto supermatrix"];
        ConcatTree [label="Infer supermatrix tree\n(FastTree / IQ-TREE / RAxML-NG)"];
        Partition  [label="Partitioned analysis\n& ModelFinder"];

        Strategy -> GeneTrees [label="default"];
        GeneTrees -> Consensus;
        Strategy -> Concat [label="-c / --concat"];
        Concat -> ConcatTree;
        ConcatTree -> Partition [label="-p / --partition"];
    }

    FinalTree [label="Final species tree\n(.nwk & figure)", fillcolor="#dff0e3", color="#4f8a5b"];

    OutputMSA -> F1;
    F3 -> Strategy;
    Consensus -> FinalTree;
    ConcatTree -> FinalTree;
    Partition -> FinalTree;
}
```
