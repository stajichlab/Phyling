# phyling download

The `download` module fetches and manages Hidden Markov Model (HMM) marker sets from the [BUSCO database](https://busco-data.ezlab.org/v5/data/lineages/) (v5 database supported).

---

## Command Syntax

```sh
phyling download [MARKER_SET | list] [OPTIONS]
```

### Positional Arguments

| Argument | Description |
| :--- | :--- |
| `list` | Displays all available online BUSCO datasets and currently cached local datasets. |
| `<markerset>` | Downloads the specified BUSCO HMM dataset (e.g. `fungi_odb10`, `eukaryota_odb10`, `metazoa_odb10`). |

### Options

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| `-v` | `--verbose` | Enable verbose output for debugging. |
| `-h` | `--help` | Show command help message and exit. |

---

## Listing Datasets

To view all available online datasets and local downloads, run:

```sh
phyling download list
```

### Output Breakdown

- **Datasets available online**: Displays all lineage marker datasets available on the official BUSCO server.
- **Datasets available on local**: Displays marker sets stored in your `$PHYLING_DB` paths or `~/.phyling`.

---

## Downloading a Dataset

To download a dataset, supply its exact name from `phyling download list`:

```sh
phyling download fungi_odb10
```

```sh
phyling download metazoa_odb10
```

---

## Automatic Update Checking

Each time `phyling download` runs, Phyling checks online metadata for dataset updates.

- If a newer version is available online, the local dataset will be marked as **[Outdated]** in `phyling download list`.
- To update an outdated dataset, re-run `phyling download [markerset]`.

---

## Storage & Database Fallback Logic

1. Phyling checks `$PHYLING_DB` for writable storage paths.
2. If `$PHYLING_DB` is set (e.g., `export PHYLING_DB=/data/busco_db`), the dataset is saved into the first path listed.
3. If `$PHYLING_DB` is unset or not writable, downloads default to `~/.phyling/`.
