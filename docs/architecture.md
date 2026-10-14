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

<!-- Both variants come from _diagrams/phyling_workflow.dot.in; Furo shows the one matching the active theme. -->
:::{container} only-light
```{graphviz} _generated/phyling_workflow-light.dot
:align: center
:caption: Phyling module workflow
```
:::

:::{container} only-dark
```{graphviz} _generated/phyling_workflow-dark.dot
:align: center
:caption: Phyling module workflow
```
:::
