# stangene

Gene identifier harmonization for single-cell transcriptomics datasets.

## What it does

Maps gene features from individual datasets into a shared canonical gene identity
system using a tiered matching cascade, while preserving all original information
and tracking mapping provenance.

## Install

```bash
pip install -e .
```

## Quick start

### 1. Build reference data (one-time)

```bash
stangene build-refs --species human
stangene build-refs --species mouse
```

### 2. Harmonize a dataset

```python
import stangene

result = stangene.run(
    path="my_data.h5ad",
    species="human",
    output_dir="results/",
)
```

Or via CLI:

```bash
stangene harmonize --input my_data.h5ad --species human --output-dir results/
```

### 3. Review outputs

- `results/harmonization_table.tsv` — full mapping table
- `results/summary.json` — statistics
- `results/conflicts.tsv` — conflicts and ambiguities
- `results/unmapped.tsv` — unmapped features

## Matching cascade

Features are matched in priority order:

1. **Tier 1 — Exact stable ID:** Ensembl gene ID exact match (high confidence)
2. **Tier 2 — Version-stripped ID:** Match after removing `.N` version suffix (high confidence)
3. **Tier 3 — Exact approved symbol:** Official gene symbol match (high confidence)
4. **Tier 4 — Alias/previous symbol:** Match via synonyms or old names (medium confidence)
5. **Tier 5 — Unmapped:** No confident match found

Non-gene features (antibody capture, CRISPR guides, spike-ins, peaks) are
classified and excluded from gene matching.

## Design principles

- **Never destroy original information.** Original identifiers are always preserved.
- **Stable IDs over symbols.** Ensembl gene IDs are the canonical key.
- **Conservative by default.** Ambiguous > incorrect.
- **Full traceability.** Every mapping records its tier, confidence, and source.

## Supported species

- Human (via HGNC)
- Mouse (via MGI + Ensembl)

## Reference data

References are built from official sources:
- **Human:** [HGNC complete gene set](https://www.genenames.org/download/statistics-and-files/)
- **Mouse:** [MGI marker files](https://www.informatics.jax.org/downloads/reports/) + Ensembl BioMart

Build with `stangene build-refs --species <name>`. References are stored locally
in `references/` and can be version-controlled or hosted on GitHub.
