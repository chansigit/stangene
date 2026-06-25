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
stangene build-refs --species human       # downloads HGNC (~15 MB)
stangene build-refs --species mouse       # downloads MGI + BioMart (~10 MB)
stangene build-refs --species rat         # downloads RGD (~5 MB)
stangene build-refs --species zebrafish   # downloads ZFIN
stangene build-refs --species fruit_fly   # downloads FlyBase
stangene build-refs --species c_elegans   # downloads WormBase
# Ensembl-only (BioMart) species:
stangene build-refs --species cynomolgus
stangene build-refs --species rhesus
stangene build-refs --species marmoset
stangene build-refs --species mouse_lemur
```

Or from Python:

```python
from stangene import build_reference

# Dedicated nomenclature authorities
for sp in ["human", "mouse", "rat", "zebrafish", "fruit_fly", "c_elegans"]:
    build_reference(sp)

# Ensembl BioMart (no dedicated authority)
for sp in ["cynomolgus", "rhesus", "marmoset", "mouse_lemur"]:
    build_reference(sp)
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

## Step 4: Single-cell QC preparation

After harmonization, stangene can generate the boolean masks that feed directly into scanpy/AnnData QC metrics — no manual gene-list curation needed.

```python
import stangene
import scanpy as sc

adata = sc.read_h5ad("my_data.h5ad")
result = stangene.run(path="my_data.h5ad", species="human")

# Canonical symbols from the harmonization result
symbols = result.mapping_table["gene_symbol_harmonized"].fillna(
    result.mapping_table["original_feature_name"]
)

# One line per QC metric
adata.var["is_mito"]  = stangene.mito_mask(symbols, "human")   # MT- prefix genes
adata.var["is_hb"]    = stangene.hb_mask(symbols, "human")     # HBA1, HBB, ...
adata.var["biotype"]  = result.mapping_table["canonical_biotype"].values

# Feed into scanpy QC
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["is_mito", "is_hb"],
    inplace=True,
)
```

`canonical_biotype` uses a unified 13-category vocabulary (`protein_coding`, `lncRNA`, `pseudogene`, `miRNA`, `rRNA`, ...) normalised from HGNC / MGI / RGD / Ensembl / WormBase / ZFIN / FlyBase — the same vocabulary across all 10 supported species.

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
