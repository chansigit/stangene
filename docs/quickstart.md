# Quick Start

## Installation

```bash
pip install -e .

# With development dependencies (pytest, scanpy):
pip install -e ".[dev]"
```

**Dependencies:** pandas, anndata, pyarrow. Downloads use stdlib `urllib` (no `requests`).

## Step 1: Build reference data

Before harmonizing, build the reference annotation databases. This is a one-time step per species.

```bash
stangene build-refs --species human   # downloads HGNC (~15 MB)
stangene build-refs --species mouse   # downloads MGI + BioMart (~10 MB)
```

Or from Python:

```python
from stangene import build_reference

build_reference("human")
build_reference("mouse")
```

References are stored in a local `references/` directory (gitignored by default). Re-run with `--force` to update from the latest upstream sources.

## Step 2: Harmonize a dataset

### Python API

```python
import stangene

result = stangene.run(
    path="my_data.h5ad",       # or .tsv / .csv
    species="human",            # or "mouse"
    output_dir="results/",      # where to write reports
    dataset_name="pbmc_10k",    # optional label
)

# Inspect results programmatically
print(result.stats)
print(result.mapping_table.head())
```

### CLI

```bash
stangene harmonize --input my_data.h5ad --species human --output-dir results/
```

## Step 3: Review outputs

After running, the output directory contains:

| File | Contents |
|---|---|
| `harmonization_table.tsv` | Full mapping table, one row per original feature |
| `summary.json` | Dataset-level statistics as JSON |
| `report.md` | Human-readable markdown report |
| `conflicts.tsv` | Many-to-one collisions, ambiguities, outdated names |
| `unmapped.tsv` | Unmapped features for manual review |
| `*_harmonized.h5ad` | Enriched h5ad with harmonization columns in `adata.var` |

The markdown report (`report.md`) is the best starting point for understanding your results. It includes summary tables, tier breakdowns, conflict details, and warnings about potential issues like Excel-corrupted gene names.

## Example output

Running on the 10x pbmc3k dataset (32,738 human genes):

```
Stats: {
    'exact_id': 24260,        # 74.1% - matched by Ensembl ID
    'exact_symbol': 411,      #  1.3% - matched by approved symbol
    'previous_symbol': 172,   #  0.5% - matched by old gene name
    'alias_symbol': 33,       #  0.1% - matched by alias
    'unmapped': 7859,          # 24.0% - GENCODE novel loci not in HGNC
    'ambiguous': 3
}
```
