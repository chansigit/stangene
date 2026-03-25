# Input Formats

stangene supports two input formats. In both cases, only feature metadata is extracted — the expression matrix is never loaded into memory.

## h5ad (AnnData)

The primary format. stangene reads `adata.var` and `adata.var_names` to extract feature metadata.

Recognized `adata.var` columns:

| Column | Maps to |
|---|---|
| `gene_ids` | `original_feature_id` (Ensembl ID) |
| `feature_types` | `original_feature_type` (e.g., Gene Expression, Antibody Capture) |

When writing results, harmonization columns are added to `adata.var` in a new `*_harmonized.h5ad` file. The original `var_names` are **never overwritten**.

## TSV / CSV

stangene auto-detects common column names:

| Detected column name | Maps to |
|---|---|
| `gene`, `gene_name`, `feature_name`, `gene_symbol`, `symbol` | `original_feature_name` |
| `gene_id`, `gene_ids`, `ensembl_id`, `ensembl_gene_id`, `feature_id` | `original_feature_id` |
| `feature_types`, `feature_type` | `original_feature_type` |

If your columns have different names, pass an explicit `column_map`:

```python
ft = stangene.load_features(
    "features.tsv",
    species="human",
    column_map={
        "my_gene_col": "original_feature_name",
        "my_id_col": "original_feature_id",
    },
)
```

File extension determines the delimiter:
- `.tsv`, `.txt` → tab-separated
- `.csv` → comma-separated
