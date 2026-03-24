# Gene Identifier Harmonization Pipeline — Design Spec

**Date:** 2026-03-24
**Project:** stangene
**Status:** Reviewed

---

## 1. Goal

Build a Python package (`stangene`) that harmonizes gene identifiers in single-cell transcriptomics datasets. The pipeline maps features from individual datasets into a shared canonical gene identity system while preserving all original information and tracking mapping provenance.

This is designed to be invoked as a Claude Code / Codex skill, processing one dataset (one matrix) at a time.

### Non-goals (v1)

- Cross-dataset matrix concatenation.
- Cross-species ortholog mapping.
- Expression normalization.
- Interactive resolution of ambiguities (the tool presents them; the user decides).

---

## 2. Core Principles

1. **Never destroy original information.** All original identifiers are preserved in dedicated columns. Harmonized identifiers live in separate columns.
2. **Stable gene IDs over symbols.** Ensembl gene IDs are the preferred canonical key. Symbols are a display layer, not the identity key.
3. **Separate identity from display.** `gene_id_harmonized` (canonical key) and `gene_symbol_harmonized` (display name) are distinct outputs.
4. **Within-species only.** Each species is harmonized independently. Cross-species linkage is a separate concern for a future version.
5. **Conservative by default.** An explicit `ambiguous` is better than an incorrect forced mapping.

---

## 3. Scope

### Supported species (v1)

- Human
- Mouse

Extensible to other model organisms by adding a `SpeciesConfig` entry and a reference builder.

### Supported input formats

- AnnData `.h5ad` (primary)
- TSV/CSV count matrices or feature tables

### Typical scale

~10,000 features per dataset. Performance is not a bottleneck at this scale.

---

## 4. Package Structure

```
stangene/
├── pyproject.toml
├── README.md
├── src/
│   └── stangene/
│       ├── __init__.py          # public API: run(), version
│       ├── io.py                # load_features(), write_results()
│       ├── classify.py          # classify_features()
│       ├── references.py        # build_reference(), load_reference()
│       ├── harmonize.py         # harmonize() — matching cascade
│       ├── merge.py             # conservative merge logic
│       ├── report.py            # summary(), conflict_report(), write_reports()
│       ├── species.py           # species-specific config
│       └── _logging.py          # structured logging setup
├── references/                  # gitignored; built by build_reference()
│   ├── human/
│   │   ├── gene_table.parquet
│   │   ├── symbol_lookup.parquet
│   │   └── build_metadata.json
│   └── mouse/
│       ├── gene_table.parquet
│       ├── symbol_lookup.parquet
│       └── build_metadata.json
├── tests/
│   ├── test_classify.py
│   ├── test_harmonize.py
│   └── fixtures/
└── docs/
    └── superpowers/
        └── specs/
```

### Dependencies

- `anndata` — h5ad I/O
- `scanpy` — optional, for broader scRNA-seq workflows
- `pandas` — tabular operations
- `pyarrow` — parquet I/O (explicit dependency, not assumed via pandas)
- Standard library `urllib` for downloads (no `requests` dependency)

---

## 5. Module Design

### 5.1 Input Adapters (`io.py`)

```python
def load_features(
    path: str,
    species: str,
    dataset_name: str = None,
    column_map: dict = None,  # optional TSV column name overrides
) -> pd.DataFrame
```

Extracts feature metadata only (not the expression matrix) into a standardized FeatureTable DataFrame.

**For h5ad:** Reads `adata.var` and `adata.var_names`. Extracts `gene_ids`, `feature_types`, `genome` columns if present. Infers `reference_source` and `reference_release` from `adata.uns` metadata or ID patterns.

**For TSV/CSV:** Reads the table. Auto-detects common column names (`gene`, `gene_id`, `gene_name`, `feature_name`) or accepts explicit `column_map`. The `column_map` maps source column names to FeatureTable column names, e.g.: `column_map={"gene_name": "original_feature_name", "ensembl_id": "original_feature_id"}`.

**FeatureTable schema:**

| Column | Description | Always present |
|---|---|---|
| `original_feature_name` | From var_names or name column | Yes |
| `original_feature_id` | From gene_ids column if present | No |
| `original_feature_type` | From feature_types column if present | No |
| `feature_id_no_version` | Derived by stripping `.N` suffix from ID | No |
| `species` | User-provided | Yes |
| `dataset` | User-provided or inferred from filename | Yes |
| `reference_source` | Inferred from ID patterns if possible | No |
| `reference_release` | Inferred if available in metadata | No |

```python
def write_results(
    result: HarmonizationResult,
    output_dir: str,
    input_path: str = None,      # if h5ad, enrich adata.var and save
    overwrite_h5ad: bool = False, # if True, modify in place; else save copy
) -> None
```

- Always writes `harmonization_table.tsv` to `output_dir`.
- If input was h5ad, adds harmonization columns to `adata.var` and saves. Never overwrites `var_names`.

---

### 5.2 Feature Classification (`classify.py`)

```python
def classify_features(ft: pd.DataFrame) -> pd.DataFrame
```

Adds/updates `original_feature_type` column.

**Classification cascade:**

1. **Trust explicit labels.** If the input has a `feature_types` column with values like `Gene Expression`, `Antibody Capture`, etc., use them directly.
2. **Pattern-based heuristics** for unlabeled rows:

| Pattern | Classification |
|---|---|
| Starts with `ENSG`, `ENSMUSG`, etc. | `gene` |
| Starts with `ENST`, `ENSMUST` | `transcript` |
| Ends with `_ADT`, `_HTO`, or matches `TotalSeq` patterns | `antibody_capture` |
| Starts with `sg-`, `gRNA-` | `crispr_guide` |
| Starts with `ERCC-` | `spike_in` |
| Matches `chrN:start-end` | `peak` |
| None of the above | `gene` (default, flagged in `mapping_notes`) |

3. **Transcript features** are classified as `transcript` and treated as non-gene features in v1. They skip the harmonization cascade. A future version could add a Tier 0 that maps transcript IDs to parent gene IDs via an Ensembl transcript-to-gene table.

4. All non-gene features (including transcripts) receive `mapping_status = non_gene_feature` immediately and skip the harmonization cascade.

Patterns are defined as `(compiled_regex, feature_type)` tuples in `species.py` for extensibility.

---

### 5.3 Reference Data (`references.py` + `species.py`)

#### Species Config (`species.py`)

```python
@dataclass
class SpeciesConfig:
    name: str                    # "human", "mouse"
    ensembl_prefix: str          # "ENSG", "ENSMUSG"
    transcript_prefix: str       # "ENST", "ENSMUST"
    naming_convention: str       # "uppercase" (human), "capitalized" (mouse)
    reference_sources: dict      # source_name → {url, parser, description}
```

**Human references:**
- HGNC complete set (`hgnc_complete_set.txt`, ~15MB) — the authority for approved symbols, aliases, previous symbols, Ensembl gene IDs, HGNC IDs, Entrez IDs, status.

**Mouse references:**
- MGI marker list (`MRK_List2.rpt`, ~7MB) — approved symbols, synonyms, feature types.
- MGI-to-Ensembl mapping (`MRK_ENSEMBL.rpt`) — links MGI IDs to Ensembl IDs.
- Ensembl BioMart mouse gene table (supplementary) — fills Ensembl ID coverage gaps for mouse genes that lack an MGI-to-Ensembl mapping.

**Fallback canonical key:** `ensembl_id` can be null in `gene_table.parquet`. For genes without an Ensembl ID (primarily mouse), `source_id` (e.g., `MGI:1234567`) serves as the fallback canonical key. In `gene_id_harmonized`, Ensembl ID is preferred; if unavailable, `source_id` is used instead, and `mapping_notes` records that the canonical key is a source-specific ID rather than an Ensembl ID.

#### Reference build (`references.py`)

```python
def build_reference(
    species: str,
    reference_dir: str = None,  # defaults to package-relative references/
    force: bool = False,
) -> None
```

Downloads source files and produces normalized parquet tables:

**`gene_table.parquet`** — one row per gene:
- `ensembl_id` (str)
- `symbol` (str, approved symbol)
- `alias_symbols` (str, pipe-delimited)
- `prev_symbols` (str, pipe-delimited)
- `gene_type` (str, e.g. protein-coding, lncRNA)
- `status` (str, approved/withdrawn)
- `source` (str, HGNC/MGI)
- `source_id` (str, HGNC:1234 or MGI:1234)

**`symbol_lookup.parquet`** — flattened index, one row per (string → gene) mapping:
- `lookup_string` (str, stored in original case from source)
- `lookup_string_upper` (str, uppercased for case-insensitive joins)
- `ensembl_id` (str, nullable — null for genes without Ensembl mapping)
- `source_id` (str, e.g. HGNC:1234 or MGI:1234 — always present)
- `lookup_type` (str: `approved_symbol`, `alias_symbol`, `prev_symbol`)
- `source` (str)

The `lookup_string` preserves original case; `lookup_string_upper` enables efficient case-insensitive matching without full-table scans. Lookups join on `lookup_string` by default (case-sensitive) or `lookup_string_upper` when species config enables case-insensitive matching.

This design makes the matching cascade a series of DataFrame joins.

**`build_metadata.json`** — records source URLs, download timestamps, file checksums.

```python
def load_reference(
    species: str,
    reference_dir: str = None,
) -> dict  # {"gene_table": pd.DataFrame, "symbol_lookup": pd.DataFrame, "metadata": dict}
```

Reads parquet files into memory. Raises clear error if references not built yet.

**Hybrid update:** `build_reference(force=True)` re-downloads. Users can also place updated raw source files in `references/{species}/raw/` and rebuild from local copies.

---

### 5.4 Harmonization Cascade (`harmonize.py`)

```python
def harmonize(ft: pd.DataFrame, ref: dict) -> HarmonizationResult
```

Runs only on rows where `original_feature_type == "gene"` (or equivalent). Non-gene rows pass through with `mapping_status = non_gene_feature`.

#### Matching tiers

**Tier 1 — Exact stable ID match:**
If `original_feature_id` exactly matches an `ensembl_id` in `gene_table` → `mapping_status = exact_id`, `mapping_confidence = high`.

**Tier 2 — Version-stripped ID match:**
If `feature_id_no_version` matches an `ensembl_id` in `gene_table` → `mapping_status = id_no_version`, `mapping_confidence = high`. Version mismatch noted in `mapping_notes`.

**Tier 3 — Exact approved symbol match:**
Look up `original_feature_name` in `symbol_lookup` where `lookup_type == "approved_symbol"`.
- Exactly one match → `mapping_status = exact_symbol`, `mapping_confidence = high`.
- Multiple matches → `mapping_status = ambiguous`.
- If the match is against a withdrawn gene (status == "withdrawn" in `gene_table`), set `mapping_status = exact_symbol` but add `mapping_confidence = medium` and note "matched withdrawn gene" in `mapping_notes`.

**Tier 4 — Alias / previous symbol match:**
Look up `original_feature_name` in `symbol_lookup` where `lookup_type in ("alias_symbol", "prev_symbol")`.
- Exactly one match → `mapping_status = alias_symbol` or `previous_symbol`, `mapping_confidence = medium`.
- Multiple matches → `mapping_status = ambiguous`, `mapping_confidence = low`. All candidates recorded in `mapping_notes`.

**Tier 5 — Unmapped:**
`mapping_status = unmapped`, `mapping_confidence` = null.

**Gene type filtering:** Tiers 3 and 4 do NOT filter by `gene_type` — all gene types (protein-coding, lncRNA, pseudogene, etc.) are eligible candidates. However, when a symbol match hits a non-protein-coding gene type or a withdrawn gene, this is recorded in `mapping_notes` as additional context. This avoids false negatives while preserving information for the user to review.

#### Rules enforced

- **Early exit:** Once resolved at a tier, skip lower tiers.
- **No case coercion:** Matching is case-sensitive by default. Species-specific case-insensitive matching can be enabled for Tiers 3-4 only via `SpeciesConfig`, but original casing is always preserved.
- **One-to-many:** If one original feature maps to multiple canonical genes → `ambiguous`.
- **Many-to-one:** After the cascade, multiple original features mapping to the same `gene_id_harmonized` are flagged in the conflict report but NOT merged.
- **Duplicate input features:** If the same `original_feature_name` appears multiple times in the input (common in CITE-seq / multi-modal data), each row is harmonized independently. If two duplicates resolve to the same `gene_id_harmonized`, they appear in the conflict report as many-to-one collisions.

#### Output columns added to FeatureTable

| Column | Description |
|---|---|
| `gene_id_harmonized` | Canonical Ensembl gene ID (or null) |
| `gene_symbol_harmonized` | Official approved symbol (or null) |
| `mapping_status` | One of: `exact_id`, `id_no_version`, `exact_symbol`, `alias_symbol`, `previous_symbol`, `ambiguous`, `unmapped`, `non_gene_feature` |
| `mapping_confidence` | `high`, `medium`, `low`, or null |
| `mapping_source` | Which reference resource provided the mapping |
| `mapping_notes` | Free text: version mismatches, multiple candidates, warnings |

#### HarmonizationResult

```python
@dataclass
class HarmonizationResult:
    mapping_table: pd.DataFrame    # full FeatureTable with harmonization columns
    conflicts: pd.DataFrame        # many-to-one and one-to-many cases
    stats: dict                    # counts per mapping_status
```

---

### 5.5 Conservative Merge (`merge.py`)

```python
def merge_features(
    result: HarmonizationResult,
    policy: str = "strict",
) -> MergeResult
```

**Never called automatically.** Explicit opt-in only.

**`policy = "strict"` (default):** Merge only rows sharing the same `gene_id_harmonized` where both were resolved via Tier 1 or Tier 2 (ID-based). Covers version-only differences.

**`policy = "symbol"`:** Also merges rows sharing `gene_id_harmonized` via Tier 3 (exact symbol). Refuses to merge alias/previous symbol matches.

**Never merges:** Tier 4 matches, ambiguous rows, unmapped rows, non-gene features.

```python
@dataclass
class MergeResult:
    merged_table: pd.DataFrame     # collapsed DataFrame
    provenance: pd.DataFrame       # original rows per merged row
    merge_log: list[str]           # human-readable decisions
```

`write_reports()` also accepts `MergeResult` (via a `merge_result` kwarg). When provided, it additionally writes `merged_table.tsv` and `merge_provenance.tsv` to the output directory.

---

### 5.6 Reporting (`report.py`)

```python
def summary(result: HarmonizationResult) -> dict
```

Returns per-dataset statistics:
- Total features, gene features, non-gene features
- Counts per `mapping_status`
- Count of duplicate `gene_id_harmonized` values
- Count of duplicate `gene_symbol_harmonized` values

```python
def conflict_report(result: HarmonizationResult) -> pd.DataFrame
```

Flat table listing:
- Many-to-one collisions (multiple originals → same harmonized ID)
- One-to-many ambiguities (one original → multiple candidates)
- Unmapped features
- Suspicious symbols, detected by two checks:
  - **Old renamed symbols:** Gene names that HGNC renamed in 2020 due to Excel corruption (e.g., MARCH1→MARCHF1, SEPT2→SEPTIN2). These are caught by Tier 4 as `previous_symbol` matches, but additionally flagged as Excel-related.
  - **Date-like strings:** Regex detection of Excel-converted date formats (e.g., `1-Mar`, `2-Sep`, `1-Sep-2023`, `Sep-01`) that indicate the input was corrupted by spreadsheet software. These are marked `unmapped` with a note.
- Likely outdated names (features resolved only via `previous_symbol`)

```python
def write_reports(result: HarmonizationResult, output_dir: str) -> None
```

Writes:
- `harmonization_table.tsv` — full mapping table, one row per original feature
- `summary.json` — stats dict
- `conflicts.tsv` — conflict report
- `unmapped.tsv` — unmapped rows for easy review

All filenames are deterministic.

---

## 6. Top-level API

```python
# stangene/__init__.py

def run(
    path: str,
    species: str,
    output_dir: str = None,
    dataset_name: str = None,
    reference_dir: str = None,
) -> HarmonizationResult
```

Wraps the full pipeline:
1. `load_features(path, species, dataset_name)`
2. `classify_features(ft)`
3. `load_reference(species, reference_dir)` — raises `ReferenceNotFoundError` if not built
4. `harmonize(ft, ref)`
5. `write_reports(result, output_dir)` — if `output_dir` provided
6. Returns `HarmonizationResult`

`run()` does NOT auto-build references. It raises a clear error directing the user to call `build_reference(species)` first. The Claude Code skill (Section 7) handles the build-if-missing logic as a separate explicit step.

### CLI Entry Point

A `python -m stangene` entry point using `argparse` (no extra dependency):

```
python -m stangene harmonize --input data.h5ad --species human --output-dir results/
python -m stangene build-refs --species human [--force]
```

Defined in `src/stangene/__main__.py`. Also exposed via `pyproject.toml` `[project.scripts]` as:
```
stangene = "stangene.__main__:main"
```

---

## 7. Claude Code Skill Integration

The skill is a markdown file that instructs the agent when and how to invoke `stangene`.

**Trigger phrases:** "harmonize genes", "gene name mapping", "standardize gene identifiers", "integrate datasets with different gene names", working with h5ad files needing cross-dataset gene alignment.

**Agent workflow:**
1. Check if `stangene` is installed; if not, `pip install -e /path/to/stangene`.
2. Check if references are built for the target species; if not, run `build_reference(species)`.
3. Run `stangene.run(path, species, output_dir)`.
4. Read `summary.json` and `conflicts.tsv`.
5. Report key stats: features per tier, ambiguous/unmapped counts, notable conflicts.
6. If ambiguous/unmapped features exist, offer to show the conflict table for manual review.

The skill does NOT auto-resolve ambiguities. It presents them to the user for decisions.

---

## 8. Output Schema (Complete)

One row per original feature:

| Column | Description |
|---|---|
| `dataset` | Source dataset identifier |
| `species` | Species identifier |
| `original_feature_name` | Raw feature name from input |
| `original_feature_id` | Raw feature ID from input (if available) |
| `original_feature_type` | Classified feature type |
| `feature_id_no_version` | Ensembl ID with version suffix stripped |
| `gene_id_harmonized` | Canonical Ensembl gene ID |
| `gene_symbol_harmonized` | Official approved symbol |
| `mapping_status` | Tier that resolved the mapping |
| `mapping_confidence` | high / medium / low / null |
| `mapping_source` | Which lookup resolved this feature (e.g., `HGNC:approved_symbol`, `MGI:alias_symbol`) |
| `mapping_notes` | Free text: warnings, candidates, evidence |
| `reference_source` | Original annotation of the input dataset (e.g., `GENCODE v32`, `Cell Ranger 2020-A`, inferred from input metadata) |
| `reference_release` | Version/date of reference used for harmonization (from `build_metadata.json`) |
| `stangene_version` | Version of stangene that produced this mapping |

---

## 9. Pitfall Mitigations

| Pitfall | Mitigation |
|---|---|
| Forcing symbols to uppercase | Case-sensitive matching by default; species naming convention preserved |
| Merging on symbol collision | Merge is opt-in, never automatic, requires shared `gene_id_harmonized` |
| Discarding unmapped genes | Always preserved with `mapping_status = unmapped` |
| Confusing transcript and gene IDs | Classification step detects transcript ID patterns |
| Assuming same annotation release | `reference_release` tracked per dataset; `build_metadata.json` records reference version |
| Excel-corrupted gene names | Suspicious symbol list flags known problem names |

---

## 10. Testing Strategy

- **Unit tests** for each module: classify, harmonize, merge, report.
- **Fixture datasets** in `tests/fixtures/`: small h5ad and TSV files with known features covering all tiers (exact ID, version-stripped, symbol, alias, ambiguous, unmapped, non-gene).
- **Round-trip tests:** load → classify → harmonize → write → reload and verify all columns.
- **Edge cases:** empty datasets, all-unmapped, mixed species error, duplicate features in input.
