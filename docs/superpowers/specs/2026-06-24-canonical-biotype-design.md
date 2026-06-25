# Design: canonical_biotype — Unified Gene Biotype Normalization

**Date:** 2026-06-24
**Status:** Approved

## Problem

The `gene_type` field in `gene_table.parquet` is populated directly from upstream
databases, each with its own vocabulary:

| Source | Example values |
|--------|---------------|
| HGNC | `"protein-coding gene"`, `"non-coding RNA"`, `"pseudogene"` |
| MGI | `"protein coding gene"`, `"lncRNA gene"`, `"miRNA gene"` |
| RGD | `"protein-coding"`, `"ncrna"`, `"lincrna"`, `"snorna"` |
| Ensembl BioMart | `"protein_coding"`, `"lincRNA"`, `"miRNA"` |
| WormBase | `"protein_coding_gene"`, `"pseudogene"` |
| ZFIN | *(empty — SO IDs discarded during ETL)* |
| FlyBase | *(empty — feature type not extracted during ETL)* |

This prevents any cross-species biotype-based filtering without per-species special-casing.
It also blocks implementing QC masks (e.g., `ribo_mask`) cleanly.

## Goal

Add a `canonical_biotype: str` column to `gene_table.parquet` for all species, with
values drawn from a fixed 13-category vocabulary. The column is written at
`build_reference()` time and stored in the parquet file.

## Canonical Vocabulary

```python
CANONICAL_BIOTYPES = frozenset({
    "protein_coding",   # protein-coding genes
    "lncRNA",           # long non-coding RNA
    "pseudogene",       # all pseudogene subtypes
    "miRNA",            # microRNA
    "snoRNA",           # small nucleolar RNA
    "snRNA",            # small nuclear RNA
    "rRNA",             # ribosomal RNA
    "tRNA",             # transfer RNA
    "IG_gene",          # immunoglobulin gene segments
    "TR_gene",          # T-cell receptor gene segments
    "other_ncrna",      # ncRNA not covered above (piRNA, circRNA, etc.)
    "other",            # non-RNA non-coding (enhancers, etc.)
    "unknown",          # no biotype information available
})
```

Unmapped raw values fall through to `"unknown"`. Matching is case-insensitive after
stripping whitespace.

## Architecture

### New module: `src/stangene/biotype.py`

Single responsibility: raw `gene_type` string → `canonical_biotype`.

Public API:
- `CANONICAL_BIOTYPES: frozenset[str]`
- `normalize_biotype(raw: str, source: str) -> str`

Internals: one mapping dict per source (`_HGNC_MAP`, `_MGI_MAP`, `_RGD_MAP`,
`_ENSEMBL_MAP`, `_WORMBASE_MAP`, `_ZFIN_SO_MAP`, `_FLYBASE_MAP`). All dicts map
lowercase-stripped raw values to canonical strings.

Exposed in `__init__.py` alongside `hb_mask` and `mito_mask` for advanced users.

### Changes to `references.py`

**Standard species** (human, mouse, rat, cynomolgus, rhesus, marmoset, mouse_lemur):

After constructing `gene_table`, add:

```python
from stangene.biotype import normalize_biotype
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source=r["source"]), axis=1
)
```

**zebrafish** — ETL change required:

`genes.txt` has columns `zfin_id, so_id, symbol, ensembl_id`. The `so_id` column
(Sequence Ontology ID, e.g. `SO:0001217`) is currently discarded. Change: retain
`so_id` during parsing, pass it to `normalize_biotype(..., source="ZFIN")` which
uses `_ZFIN_SO_MAP`. Do not store `so_id` in the final `gene_table`.

Key SO ID mappings:
```
SO:0001217 → protein_coding
SO:0001263 → miRNA
SO:0001900 → lncRNA
SO:0000336 → pseudogene
SO:0001984 → snoRNA
SO:0001268 → snRNA
SO:0001637 → rRNA
SO:0001272 → tRNA
```

**fruit_fly** — ETL change required:

`fbgn_annotation_ID.tsv.gz` does not have a biotype column. Use the FlyBase
`precomputed_files/genes/` annotation which includes feature type (`gene`, `mRNA`,
`ncRNA`, `pseudogene`, `tRNA`, `rRNA`, etc.). The `_build_fruitfly_reference()`
function needs to extract this column and map via `_FLYBASE_MAP`.

### Changes to `harmonize.py`

Replace the raw `gene_type` string check (lines 153–155) with a `canonical_biotype`
check:

```python
# before
if gene_type and "protein" not in str(gene_type).lower():
    notes.append(f"Non-protein-coding gene type: {gene_type}")

# after
if gene_info.get("canonical_biotype", "unknown") not in ("protein_coding", "unknown"):
    notes.append(f"Non-protein-coding: {gene_info['canonical_biotype']}")
```

### parquet schema change

`gene_table.parquet` gains one column: `canonical_biotype: str`, positioned after
`gene_type`. All 10 pre-built species files must be regenerated as part of this PR
(requires `sh_dev` or a small Slurm job — cannot run on the login node).

## Testing

New file: `tests/test_biotype.py`

| Test | What it checks |
|------|---------------|
| per-source typical mappings | `normalize_biotype("protein-coding gene", "HGNC") == "protein_coding"` for each source |
| unknown passthrough | unmapped raw values return `"unknown"` |
| case-insensitive | `"Protein-Coding Gene"` maps same as `"protein-coding gene"` |
| integration | load a built species ref, assert `canonical_biotype` column exists, no NaN, all values in `CANONICAL_BIOTYPES` |

No changes to existing tests. The `harmonize.py` note text changes but the
output schema (columns, types) does not.

## Out of Scope

- SO ID column stored in `gene_table` (only the mapped `canonical_biotype` is stored)
- `hb_mask()` / `mito_mask()` — these use curated symbol sets, not biotype
- New public mask functions (e.g., `ribo_mask`) — separate feature, enabled by this work
- Storing user-overridable mapping as a JSON/YAML file
