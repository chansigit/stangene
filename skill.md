# Skill: Harmonize Gene Identifiers

Use this skill when the user asks to "harmonize genes", "standardize gene names",
"map gene identifiers", "gene name mapping", or when working with single-cell
transcriptomics datasets that need cross-dataset gene alignment.

## Prerequisites

1. Check if stangene is installed:
   ```python
   python -c "import stangene; print(stangene.__version__)"
   ```
   If not installed, install from the project directory:
   ```bash
   pip install -e /path/to/stangene
   ```

2. Check if references are built for the target species:
   ```python
   from stangene.references import load_reference
   try:
       load_reference("human")  # or "mouse"
   except Exception:
       from stangene.references import build_reference
       build_reference("human")  # downloads ~15MB for human, ~7MB for mouse
   ```

## Usage

Run the harmonization pipeline on a single dataset:

```python
import stangene

result = stangene.run(
    path="path/to/data.h5ad",   # or .tsv/.csv
    species="human",             # or "mouse"
    output_dir="results/",       # where to write reports
    dataset_name="my_dataset",   # optional label
)
```

## Interpreting Results

After running, read the summary and report to the user:

```python
import json

# Read summary
with open("results/summary.json") as f:
    summary = json.load(f)

# Key stats to report:
# - summary["total_features"]: total features in dataset
# - summary["gene_features"]: features classified as genes
# - summary["status_counts"]: breakdown by mapping tier
# - summary["duplicate_harmonized_ids"]: potential collisions

# Read conflicts if any
import pandas as pd
conflicts = pd.read_csv("results/conflicts.tsv", sep="\t")
```

Report to the user:
- How many features were mapped at each tier (exact_id, id_no_version, exact_symbol, alias_symbol, previous_symbol)
- How many are ambiguous or unmapped
- Any notable conflicts (many-to-one mappings, Excel-corrupted names)
- If there are unmapped/ambiguous features, offer to show the conflict table

## Important

- Do NOT auto-resolve ambiguities. Present them to the user for decisions.
- The pipeline never overwrites original identifiers.
- Species must be specified explicitly (human or mouse).
- For cross-species work, harmonize each species separately first.

## Optional: Conservative Merge

Only if the user explicitly requests merging duplicate features:

```python
from stangene import merge_features
merge_result = merge_features(result, policy="strict")  # or "symbol"
```

## Output Files

- `harmonization_table.tsv` — full mapping, one row per original feature
- `summary.json` — dataset-level statistics
- `conflicts.tsv` — conflict report
- `unmapped.tsv` — unmapped features for manual review
- `*_harmonized.h5ad` — enriched h5ad (if input was h5ad)
