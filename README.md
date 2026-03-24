# stangene

Gene identifier harmonization for single-cell transcriptomics datasets.

`stangene` maps gene features from individual datasets into a shared canonical gene identity system using a tiered matching cascade. It preserves all original information, tracks mapping provenance, and never silently forces ambiguous mappings.

Designed to be invoked as a Claude Code / Codex skill or used as a standalone Python library and CLI tool.

---

## Workflow

```mermaid
flowchart TD
    Input["Input Dataset\n(.h5ad / .tsv / .csv)"]
    Load["load_features()\nExtract feature metadata"]
    Classify["classify_features()\nTriage: gene vs non-gene"]
    Gene["Gene features"]
    NonGene["Non-gene features\n(antibody, CRISPR, spike-in,\npeaks, transcripts)"]
    NonGeneLabel["Labeled non_gene_feature\nand passed through"]
    LoadRef["load_reference()\nHGNC / MGI + BioMart"]
    Harmonize["harmonize()\n5-tier matching cascade"]

    T1["Tier 1\nExact ID\n(high)"]
    T2["Tier 2\nID no version\n(high)"]
    T3["Tier 3\nExact symbol\n(high)"]
    T4["Tier 4\nAlias / prev sym\n(medium)"]
    T5["Tier 5\nUnmapped"]
    Excel["Excel date\ndetection"]

    Result["HarmonizationResult"]
    WriteResults["write_results()\nharmonization_table.tsv\n*_harmonized.h5ad"]
    WriteReports["write_reports()\nsummary.json\nconflicts.tsv\nunmapped.tsv"]
    Merge["merge_features()\nstrict / symbol policy"]
    MergeOut["MergeResult\n(merged_table + provenance)"]

    Input --> Load --> Classify
    Classify --> Gene
    Classify --> NonGene --> NonGeneLabel

    Gene --> LoadRef --> Harmonize
    Harmonize --> T1
    Harmonize --> T2
    Harmonize --> T3
    Harmonize --> T4
    Harmonize --> T5
    Harmonize --> Excel

    T1 --> Result
    T2 --> Result
    T3 --> Result
    T4 --> Result
    T5 --> Result
    Excel --> Result

    Result --> WriteResults
    Result --> WriteReports
    Result -.->|"optional\nexplicit opt-in"| Merge --> MergeOut

    style Input fill:#e8f4f8,stroke:#2196F3
    style Result fill:#e8f5e9,stroke:#4CAF50
    style MergeOut fill:#e8f5e9,stroke:#4CAF50
    style NonGeneLabel fill:#fff3e0,stroke:#FF9800
    style T1 fill:#c8e6c9,stroke:#388E3C
    style T2 fill:#c8e6c9,stroke:#388E3C
    style T3 fill:#c8e6c9,stroke:#388E3C
    style T4 fill:#fff9c4,stroke:#FBC02D
    style T5 fill:#ffcdd2,stroke:#D32F2F
    style Excel fill:#ffcdd2,stroke:#D32F2F
```

---

## Install

```bash
pip install -e .

# With development dependencies (pytest, scanpy):
pip install -e ".[dev]"
```

**Dependencies:** pandas, anndata, pyarrow. Downloads use stdlib `urllib` (no `requests`).

---

## Quick start

### 1. Build reference data (one-time per species)

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

References are stored in a local `references/` directory (gitignored by default). Re-run with `--force` to update.

### 2. Harmonize a dataset

```python
import stangene

result = stangene.run(
    path="my_data.h5ad",       # or .tsv / .csv
    species="human",            # or "mouse"
    output_dir="results/",      # where to write reports
    dataset_name="pbmc_10k",    # optional label
)
```

Or via CLI:

```bash
stangene harmonize --input my_data.h5ad --species human --output-dir results/
```

### 3. Review outputs

| File | Contents |
|---|---|
| `harmonization_table.tsv` | Full mapping table, one row per original feature |
| `summary.json` | Dataset-level statistics |
| `conflicts.tsv` | Many-to-one collisions, ambiguities, outdated names |
| `unmapped.tsv` | Unmapped features for manual review |
| `*_harmonized.h5ad` | Enriched h5ad with harmonization columns in `adata.var` (if input was h5ad) |

---

## How it works

### Feature classification

Before harmonization, every feature is classified:

| Pattern | Type |
|---|---|
| `ENSG*` / `ENSMUSG*` | `gene` |
| `ENST*` / `ENSMUST*` | `transcript` |
| `*_ADT`, `*_HTO`, `*TotalSeq*` | `antibody_capture` |
| `sg-*`, `gRNA-*` | `crispr_guide` |
| `ERCC-*` | `spike_in` |
| `chrN:start-end` | `peak` |
| 10x Cell Ranger labels | Mapped directly (`Gene Expression` -> `gene`, etc.) |
| Anything else | `gene` (default, flagged) |

Non-gene features are labeled `non_gene_feature` and excluded from the matching cascade. They are preserved in the output.

### Matching cascade

Gene features are matched in strict priority order. Once a feature is resolved at a tier, lower tiers are skipped.

| Tier | Strategy | Confidence | Example |
|---|---|---|---|
| 1 | Exact Ensembl gene ID | high | `ENSG00000141510` matches directly |
| 2 | Version-stripped Ensembl ID | high | `ENSG00000141510.18` -> `ENSG00000141510` |
| 3 | Exact approved symbol | high | `TP53` matches HGNC approved symbol |
| 4 | Alias or previous symbol | medium | `p53` matches via HGNC alias for TP53 |
| 5 | Unmapped | -- | No confident match; preserved as-is |

**Additional checks at each tier:**
- Withdrawn genes matched at Tier 3 are downgraded to `medium` confidence
- Non-protein-coding gene types are noted in `mapping_notes`
- Excel-corrupted date formats (`1-Mar`, `2-Sep`) are detected and marked unmapped
- Known HGNC-renamed genes (MARCH1->MARCHF1, SEPT2->SEPTIN2, etc.) are flagged

### Output schema

Every feature receives these harmonization columns:

| Column | Description |
|---|---|
| `gene_id_harmonized` | Canonical Ensembl gene ID (or MGI source_id fallback) |
| `gene_symbol_harmonized` | Official approved symbol |
| `mapping_status` | `exact_id`, `id_no_version`, `exact_symbol`, `alias_symbol`, `previous_symbol`, `ambiguous`, `unmapped`, or `non_gene_feature` |
| `mapping_confidence` | `high`, `medium`, `low`, or null |
| `mapping_source` | Which lookup resolved this feature (e.g., `HGNC:approved_symbol`) |
| `mapping_notes` | Warnings, candidate lists, version mismatches |
| `stangene_version` | Version of stangene that produced the mapping |
| `reference_release` | Timestamp of the reference data used |

Original columns (`original_feature_name`, `original_feature_id`, `original_feature_type`, `feature_id_no_version`) are always preserved.

---

## Detailed API

### Top-level

```python
import stangene

# Full pipeline
result = stangene.run(path, species, output_dir=None, dataset_name=None, reference_dir=None)
# Returns: HarmonizationResult

# Individual stages
stangene.build_reference(species, reference_dir=None, force=False)
ref = stangene.load_reference(species, reference_dir=None)
ft = stangene.load_features(path, species, dataset_name=None, column_map=None)
ft = stangene.classify_features(ft)
result = stangene.harmonize(ft, ref)

# Reporting
stats = stangene.summary(result)
conflicts = stangene.conflict_report(result)
stangene.write_reports(result, output_dir)

# Optional merge
merge_result = stangene.merge_features(result, policy="strict")
```

### `stangene.run()`

Convenience function that chains the full pipeline:

1. `load_features()` -- extract feature metadata from h5ad or TSV
2. `classify_features()` -- triage gene vs. non-gene features
3. `load_reference()` -- load pre-built annotation tables (raises `ReferenceNotFoundError` if not built)
4. `harmonize()` -- run the 5-tier matching cascade
5. `write_results()` + `write_reports()` -- write outputs (if `output_dir` provided)

### `HarmonizationResult`

```python
@dataclass
class HarmonizationResult:
    mapping_table: pd.DataFrame    # full FeatureTable with harmonization columns
    conflicts: pd.DataFrame        # many-to-one collision rows
    stats: dict                    # counts per mapping_status
```

### `MergeResult`

```python
@dataclass
class MergeResult:
    merged_table: pd.DataFrame     # collapsed rows
    provenance: pd.DataFrame       # which originals contributed to each merge
    merge_log: list[str]           # human-readable merge decisions
```

---

## Conservative merge

Merge is **never automatic**. It must be explicitly requested.

```python
from stangene import merge_features

# Strict: only merge rows sharing gene_id_harmonized via Tier 1 or 2 (ID-based)
merge_result = merge_features(result, policy="strict")

# Symbol: also merge Tier 3 (exact approved symbol) matches
merge_result = merge_features(result, policy="symbol")
```

**Never merges:** Tier 4 (alias/previous symbol), ambiguous, unmapped, or non-gene features.

Every merge is recorded in `merge_result.provenance` with the original rows and merge reason, making it fully reversible.

---

## Reference data

### Human (HGNC)

Source: [HGNC complete gene set](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt) (~15 MB)

Provides: approved symbols, alias symbols, previous symbols, Ensembl gene IDs, HGNC IDs, gene types, approval status.

### Mouse (MGI + Ensembl BioMart)

Sources:
- [MGI marker list](https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt) (~7 MB) -- symbols, synonyms, feature types
- [MGI-to-Ensembl mapping](https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt) -- MGI ID to Ensembl ID links
- Ensembl BioMart (supplementary) -- fills Ensembl ID gaps for mouse genes not covered by MGI mapping

### Reference internals

Built references are stored as:

```
references/<species>/
  gene_table.parquet       # one row per gene
  symbol_lookup.parquet    # flattened symbol->gene index
  build_metadata.json      # source URLs, timestamps, checksums
```

`gene_table` columns: `ensembl_id`, `symbol`, `alias_symbols`, `prev_symbols`, `gene_type`, `status`, `source`, `source_id`

`symbol_lookup` columns: `lookup_string`, `lookup_string_upper`, `ensembl_id`, `source_id`, `lookup_type`, `source`

---

## Conflict detection

The conflict report (`conflicts.tsv`) flags:

| Conflict type | Description |
|---|---|
| `many_to_one` | Multiple original features map to the same canonical gene ID |
| `unmapped` | Feature could not be matched to any reference gene |
| `ambiguous` | Feature matched multiple candidate genes |
| `outdated_name` | Feature resolved via a previous symbol (likely renamed) |

Additionally, the report detects:
- **Excel date artifacts:** Gene names corrupted by spreadsheet auto-formatting (`1-Mar`, `2-Sep`)
- **Known Excel-renamed genes:** Symbols renamed by HGNC in 2020 (e.g., MARCH1->MARCHF1, SEPT2->SEPTIN2)

---

## Input formats

### h5ad (primary)

Reads `adata.var` metadata. Extracts `gene_ids`, `feature_types`, and `genome` columns if present. Does **not** load the expression matrix.

When writing results, harmonization columns are added to `adata.var` in a new `*_harmonized.h5ad` file. Original `var_names` are never overwritten.

### TSV / CSV

Reads the table and auto-detects common column names:

| Detected column name | Maps to |
|---|---|
| `gene`, `gene_name`, `feature_name`, `gene_symbol`, `symbol` | `original_feature_name` |
| `gene_id`, `gene_ids`, `ensembl_id`, `ensembl_gene_id`, `feature_id` | `original_feature_id` |
| `feature_types`, `feature_type` | `original_feature_type` |

Or pass an explicit `column_map`:

```python
ft = stangene.load_features(
    "features.tsv",
    species="human",
    column_map={"my_gene_col": "original_feature_name", "my_id_col": "original_feature_id"},
)
```

---

## Design principles

1. **Never destroy original information.** All original identifiers are preserved. Harmonized identifiers live in separate columns.
2. **Stable IDs over symbols.** Ensembl gene IDs are the canonical key. Gene symbols are a display layer.
3. **Separate identity from display.** `gene_id_harmonized` (identity) and `gene_symbol_harmonized` (display) are distinct.
4. **Within-species only.** Cross-species linkage via orthology is a separate concern, not conflated with name harmonization.
5. **Conservative by default.** An explicit `ambiguous` or `unmapped` is always better than an incorrect forced mapping.
6. **Full traceability.** Every mapping records its tier, confidence, source, and notes.

---

## Adding a new species

1. Add a `SpeciesConfig` entry in `src/stangene/species.py` with the species' Ensembl prefix, naming convention, and reference source URLs.
2. Add classification patterns if needed (e.g., new Ensembl prefix patterns).
3. Implement a `_build_<species>_reference()` function in `references.py` that downloads source files, builds `gene_table` and `symbol_lookup` DataFrames, and calls `_save_reference()`.
4. Add the dispatch to `build_reference()`.

---

## Testing

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific module tests
python -m pytest tests/test_harmonize.py -v
python -m pytest tests/test_references.py -v
```

70 tests covering: species config, feature classification, I/O adapters, human/mouse reference building, all 5 harmonization tiers, Excel corruption detection, withdrawn gene handling, conservative merge, reporting, and end-to-end integration.

---

## Project structure

```
stangene/
├── pyproject.toml
├── README.md
├── skill.md                     # Claude Code skill definition
├── src/stangene/
│   ├── __init__.py              # public API: run(), re-exports
│   ├── __main__.py              # CLI: harmonize, build-refs
│   ├── species.py               # SpeciesConfig, classification patterns, Excel lists
│   ├── classify.py              # classify_features()
│   ├── io.py                    # load_features(), write_results()
│   ├── references.py            # build_reference(), load_reference()
│   ├── harmonize.py             # harmonize(), HarmonizationResult
│   ├── merge.py                 # merge_features(), MergeResult
│   ├── report.py                # summary(), conflict_report(), write_reports()
│   └── _logging.py              # structured logging
├── references/                  # built by build_reference(), gitignored
├── tests/
│   ├── conftest.py              # shared mock fixtures
│   ├── test_species.py
│   ├── test_classify.py
│   ├── test_io.py
│   ├── test_references.py
│   ├── test_harmonize.py
│   ├── test_merge.py
│   ├── test_report.py
│   └── test_run.py
└── docs/superpowers/
    ├── specs/                   # design spec
    └── plans/                   # implementation plan
```
