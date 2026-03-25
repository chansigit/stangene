# Conflicts and Merge Policy

## Understanding conflicts

After harmonization, some features may have conflicts that require attention.

### Many-to-one collisions

Multiple original features in a dataset map to the same canonical gene. This happens when:

- A gene and its alias both appear in the input (e.g., `TP53` and `p53`)
- A current name and a previous name coexist (e.g., `BRCA1` and `RNF53`)
- A version-specific and version-stripped ID both appear

stangene flags these in the conflict report but **never merges them automatically**.

### Ambiguous features

A feature name matches multiple candidate genes in the reference. For example, an alias that is shared by two different genes. These are marked `ambiguous` with all candidates recorded in `mapping_notes`.

### Unmapped features

Features that couldn't be matched to any reference gene. Common reasons:

- Novel loci from GENCODE not registered in HGNC (e.g., `RP11-*`, `AL*`)
- Non-coding RNAs absent from the reference
- Genuinely unknown or custom features
- Excel-corrupted gene names (date formats like `1-Mar`)

## Conservative merge

Merging is **never automatic**. You must explicitly opt in.

```python
from stangene import merge_features

# Strict: only merge rows sharing the same harmonized ID via Tier 1 or 2
merge_result = merge_features(result, policy="strict")

# Symbol: also merge Tier 3 (exact approved symbol) matches
merge_result = merge_features(result, policy="symbol")
```

### What gets merged

| Policy | Eligible statuses |
|---|---|
| `strict` | `exact_id`, `id_no_version` only |
| `symbol` | `exact_id`, `id_no_version`, `exact_symbol` |

### What never gets merged

- Tier 4 matches (`alias_symbol`, `previous_symbol`)
- `ambiguous` features
- `unmapped` features
- `non_gene_feature` entries

### MergeResult

```python
merge_result.merged_table    # collapsed DataFrame
merge_result.provenance      # which originals contributed to each merged row
merge_result.merge_log       # human-readable merge decisions
```

Every merge is recorded in provenance, making it fully reversible.

## Pitfalls to avoid

1. **Don't force symbols to uppercase.** This erases species-specific naming conventions (mouse uses capitalized form like `Trp53`, not `TRP53`).

2. **Don't merge on symbol collision alone.** Two features mapping to the same display symbol doesn't mean they should be summed.

3. **Don't discard unmapped genes.** They may be rescuable with updated annotation resources.

4. **Don't confuse gene and transcript IDs.** `ENST*` identifiers are transcript-level, not gene-level.

5. **Don't assume all datasets use the same annotation release.** Different GENCODE/Ensembl versions can change symbols.
