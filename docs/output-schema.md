# Output Schema

## Harmonization table

The harmonization table (`harmonization_table.tsv`) contains one row per original feature with all mapping metadata.

| Column | Description | Always present |
|---|---|---|
| `original_feature_name` | Raw feature name from input | Yes |
| `original_feature_id` | Raw feature ID from input (Ensembl ID if available) | No |
| `original_feature_type` | Classified feature type (`gene`, `spike_in`, etc.) | Yes |
| `feature_id_no_version` | Ensembl ID with version suffix stripped | No |
| `species` | Species identifier | Yes |
| `dataset` | Dataset name | Yes |
| `reference_source` | Inferred annotation source of input (e.g., Ensembl/GENCODE) | No |
| `reference_release` | Timestamp of reference data used for harmonization | Yes |
| `gene_id_harmonized` | Canonical Ensembl gene ID (or source_id fallback) | Mapped only |
| `gene_symbol_harmonized` | Official approved symbol | Mapped only |
| `mapping_status` | Which tier resolved the mapping | Yes |
| `mapping_confidence` | `high`, `medium`, `low`, or null | Yes |
| `mapping_source` | Lookup that resolved this feature | Mapped only |
| `mapping_notes` | Warnings, candidates, version mismatches | When applicable |
| `stangene_version` | Version of stangene that produced this mapping | Yes |

## Mapping status values

| Status | Tier | Meaning |
|---|---|---|
| `exact_id` | 1 | Ensembl gene ID matched exactly |
| `id_no_version` | 2 | Ensembl ID matched after stripping version suffix |
| `exact_symbol` | 3 | Official approved gene symbol matched |
| `alias_symbol` | 4 | Matched via an alias (alternative) name |
| `previous_symbol` | 4 | Matched via a previous (old) name |
| `ambiguous` | - | Multiple candidate genes; not resolved |
| `unmapped` | 5 | No confident match found |
| `non_gene_feature` | - | Not a gene; excluded from matching |

## Confidence levels

| Level | When assigned |
|---|---|
| `high` | Tier 1, 2, or 3 match (exact ID or approved symbol) |
| `medium` | Tier 4 match (alias/previous symbol), or Tier 3 match to a withdrawn gene |
| `low` | Ambiguous — multiple candidates found |
| null | Unmapped or non-gene feature |

## Summary JSON

`summary.json` contains:

```json
{
    "total_features": 32738,
    "gene_features": 32738,
    "non_gene_features": 0,
    "status_counts": {
        "exact_id": 24260,
        "unmapped": 7859,
        "exact_symbol": 411,
        "previous_symbol": 172,
        "alias_symbol": 33,
        "ambiguous": 3
    },
    "duplicate_harmonized_ids": 165,
    "duplicate_harmonized_symbols": 165
}
```

## Conflict report

`conflicts.tsv` lists issues that may require manual review:

| Conflict type | Description |
|---|---|
| `many_to_one` | Multiple original features map to the same canonical gene |
| `unmapped` | Feature could not be matched |
| `outdated_name` | Feature resolved via a previous symbol (gene was renamed) |
| `ambiguous` | Feature matched multiple candidate genes |

## Markdown report

`report.md` is a human-readable summary containing all of the above in formatted tables, plus warnings about Excel-corrupted gene names and detailed conflict listings.
