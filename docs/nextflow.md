# Nextflow Workflow

For high-throughput execution across large genomic cohorts, high-performance computing (HPC) clusters, or cloud infrastructure, the Phyling pipeline is wrapped in a dedicated Nextflow workflow: **[nf_phyling](https://github.com/stajichlab/nf_phyling)**.

---

## Overview

**`nf_phyling`** automates the end-to-end execution of Phyling modules (`download`, `align`, `filter`, `tree`) and more dedicated
final tree inference using [Nextflow](https://www.nextflow.io/). It handles data parallelization, resource management, dependency
containerization, and HPC scheduler submission automatically.

### Key Features

- **HPC Ready**: Out-of-the-box support for SLURM and UCR HPCC environment modules.
- **Containerized Execution**: Native support for **Singularity/Apptainer** and **pixi** profiles for 100% reproducible environments.
- **Parallel Processing**: Automatically shards HMM search and alignment tasks across cluster compute nodes to maximize throughput.
- **Resumable Pipelines**: Nextflow's `-resume` capability ensures failed or interrupted cluster runs pick up right where they left off.

---

## Quick Start

### Prerequisites

Ensure [Nextflow](https://www.nextflow.io/) (v21.04+) and either Docker, Singularity, or Conda are installed on your system.

### Running the Workflow

To run `nf_phyling` directly from GitHub:

```sh
nextflow run stajichlab/nf_phyling \
  -profile singularity_slurm,ucr_hpcc \
  --seq_type protein \
  --input /path/to/pep \
  --prefix my_project \
  --markerset fungi_odb12,mucoromycota_odb12 \
  --outdir results/pep
```

---

## Repository Links

- **GitHub Repository**: [stajichlab/nf_phyling](https://github.com/stajichlab/nf_phyling)
- **Issue Tracker**: [stajichlab/nf_phyling/issues](https://github.com/stajichlab/nf_phyling/issues)
