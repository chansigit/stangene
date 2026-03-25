# Reference Data

stangene harmonizes gene names against official annotation databases. References must be built once per species before harmonization.

## Building references

```bash
# CLI
stangene build-refs --species human
stangene build-refs --species mouse
stangene build-refs --species human --force  # re-download and rebuild

# Python
from stangene import build_reference
build_reference("human")
build_reference("mouse", force=True)
```

## Human (HGNC)

**Source:** [HGNC complete gene set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) (~15 MB)

Provides:
- Approved gene symbols
- Alias (alternative) symbols
- Previous (old) symbols
- Ensembl gene IDs
- HGNC IDs
- Gene types (protein-coding, lncRNA, pseudogene, etc.)
- Approval status (Approved, Entry Withdrawn)

## Mouse (MGI + Ensembl BioMart)

**Sources:**
- [MGI marker list](https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt) (~7 MB) — approved symbols, synonyms, feature types
- [MGI-to-Ensembl mapping](https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt) — MGI ID to Ensembl ID links
- Ensembl BioMart (supplementary, non-fatal if unavailable) — fills Ensembl ID gaps for mouse genes not covered by MGI mapping

## Internal format

Built references are stored as parquet files:

```
references/<species>/
├── gene_table.parquet       # one row per gene
├── symbol_lookup.parquet    # flattened symbol → gene index
└── build_metadata.json      # source URLs, timestamps, checksums
```

### gene_table columns

| Column | Description |
|---|---|
| `ensembl_id` | Ensembl gene ID (nullable for some mouse genes) |
| `symbol` | Approved gene symbol |
| `alias_symbols` | Pipe-delimited alias symbols |
| `prev_symbols` | Pipe-delimited previous symbols |
| `gene_type` | Gene biotype (protein-coding, lncRNA, etc.) |
| `status` | Approval status (Approved / Entry Withdrawn) |
| `source` | Reference authority (HGNC / MGI) |
| `source_id` | Authority-specific ID (HGNC:12345 / MGI:12345) |

### symbol_lookup columns

| Column | Description |
|---|---|
| `lookup_string` | The symbol/alias/prev string (original case) |
| `lookup_string_upper` | Uppercased for case-insensitive matching |
| `ensembl_id` | Target Ensembl gene ID (nullable) |
| `source_id` | Target authority ID (always present) |
| `lookup_type` | `approved_symbol`, `alias_symbol`, or `prev_symbol` |
| `source` | Reference authority |

### build_metadata.json

Records exactly what was downloaded and when, for reproducibility:
- Source URLs
- SHA-256 checksums of downloaded files
- Download timestamps
- Row counts

## Custom reference directory

By default, references are stored in a `references/` directory relative to the package. To use a custom location:

```python
build_reference("human", reference_dir="/path/to/my/refs")
result = stangene.run("data.h5ad", species="human", reference_dir="/path/to/my/refs")
```

## Versioning references

The `references/` directory is gitignored by default. If you want to version-control your references (recommended for reproducibility), you can:

1. Commit the parquet files to a separate git repo or GitHub release
2. Or remove `references/` from `.gitignore` and commit them directly
