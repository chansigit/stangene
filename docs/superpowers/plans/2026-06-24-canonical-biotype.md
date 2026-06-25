# canonical_biotype Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `canonical_biotype` column to every species' `gene_table.parquet`, normalizing raw upstream `gene_type` strings into a 13-category vocabulary via a new `biotype.py` module.

**Architecture:** A new `src/stangene/biotype.py` holds all per-source mapping dicts and `normalize_biotype(raw, source) -> str`. Each `_build_*_reference()` function in `references.py` calls this function and writes the result into the parquet at build time. Two species (zebrafish, fruit_fly) require minor ETL changes to surface biotype information that is already present in the raw source files but currently discarded.

**Tech Stack:** Python 3.10+, pandas, pyarrow (parquet), pytest

## Global Constraints

- All values of `canonical_biotype` must be members of `CANONICAL_BIOTYPES` (defined in Task 1). No freeform strings allowed.
- Mapping is case-insensitive (lowercase-strip before lookup).
- Unmapped raw values → `"unknown"`. Never raise.
- Applied at `build_reference()` time, persisted to parquet. Not recomputed at `load_reference()`.
- No new network dependencies for any species.
- Must pass `pytest` after each task.

---

## File Map

| Action | Path | Responsibility |
|--------|------|---------------|
| CREATE | `src/stangene/biotype.py` | Vocabulary constant + per-source mapping dicts + `normalize_biotype()` |
| CREATE | `tests/test_biotype.py` | Unit + integration tests for biotype module |
| MODIFY | `src/stangene/references.py` | All `_build_*_reference()` add `canonical_biotype` col; human uses `locus_type`; zebrafish preserves `so_id`; fruit_fly extracts annotation type |
| MODIFY | `src/stangene/harmonize.py` | Lines 153-155: use `canonical_biotype` instead of raw `gene_type` string |
| MODIFY | `src/stangene/__init__.py` | Expose `CANONICAL_BIOTYPES` and `normalize_biotype` |

---

## Task 1: `biotype.py` — vocabulary + mapping + tests

**Files:**
- Create: `src/stangene/biotype.py`
- Create: `tests/test_biotype.py`

**Interfaces:**
- Produces:
  - `CANONICAL_BIOTYPES: frozenset[str]`
  - `normalize_biotype(raw: str, source: str) -> str`

---

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_biotype.py
import pytest
from stangene.biotype import CANONICAL_BIOTYPES, normalize_biotype


# --- vocabulary ---

def test_canonical_biotypes_is_frozenset():
    assert isinstance(CANONICAL_BIOTYPES, frozenset)


def test_canonical_biotypes_contains_expected():
    expected = {
        "protein_coding", "lncRNA", "pseudogene",
        "miRNA", "snoRNA", "snRNA", "rRNA", "tRNA",
        "IG_gene", "TR_gene", "other_ncrna", "other", "unknown",
    }
    assert expected == CANONICAL_BIOTYPES


# --- normalize_biotype: HGNC (locus_type values) ---

@pytest.mark.parametrize("raw,expected", [
    ("gene with protein product", "protein_coding"),
    ("RNA, micro",                "miRNA"),
    ("RNA, small nucleolar",      "snoRNA"),
    ("RNA, small nuclear",        "snRNA"),
    ("RNA, ribosomal",            "rRNA"),
    ("RNA, transfer",             "tRNA"),
    ("RNA, long non-coding",      "lncRNA"),
    ("immunoglobulin gene",       "IG_gene"),
    ("T cell receptor gene",      "TR_gene"),
    ("pseudogene",                "pseudogene"),
    ("endogenous retrovirus",     "other"),
    ("unknown",                   "unknown"),
])
def test_hgnc_mappings(raw, expected):
    assert normalize_biotype(raw, "HGNC") == expected


def test_hgnc_case_insensitive():
    assert normalize_biotype("Gene With Protein Product", "HGNC") == "protein_coding"
    assert normalize_biotype("  RNA, Micro  ", "HGNC") == "miRNA"


# --- normalize_biotype: MGI ---

@pytest.mark.parametrize("raw,expected", [
    ("protein coding gene",              "protein_coding"),
    ("lncRNA gene",                      "lncRNA"),
    ("pseudogene",                       "pseudogene"),
    ("miRNA gene",                       "miRNA"),
    ("snoRNA gene",                      "snoRNA"),
    ("snRNA gene",                       "snRNA"),
    ("rRNA gene",                        "rRNA"),
    ("tRNA gene",                        "tRNA"),
    ("unclassified non-coding RNA gene", "other_ncrna"),
    ("transgene",                        "other"),
    ("unclassified gene",                "other"),
    ("gene segment",                     "other"),
])
def test_mgi_mappings(raw, expected):
    assert normalize_biotype(raw, "MGI") == expected


# --- normalize_biotype: RGD ---

@pytest.mark.parametrize("raw,expected", [
    ("protein-coding", "protein_coding"),
    ("ncrna",          "other_ncrna"),
    ("lincrna",        "lncRNA"),
    ("lncrna",         "lncRNA"),
    ("mirna",          "miRNA"),
    ("snorna",         "snoRNA"),
    ("snrna",          "snRNA"),
    ("rrna",           "rRNA"),
    ("trna",           "tRNA"),
    ("pseudo",         "pseudogene"),
    ("pseudogene",     "pseudogene"),
])
def test_rgd_mappings(raw, expected):
    assert normalize_biotype(raw, "RGD") == expected


# --- normalize_biotype: Ensembl (BioMart gene_biotype) ---

@pytest.mark.parametrize("raw,expected", [
    ("protein_coding",                          "protein_coding"),
    ("lincRNA",                                 "lncRNA"),
    ("lncRNA",                                  "lncRNA"),
    ("miRNA",                                   "miRNA"),
    ("snoRNA",                                  "snoRNA"),
    ("snRNA",                                   "snRNA"),
    ("rRNA",                                    "rRNA"),
    ("Mt_rRNA",                                 "rRNA"),
    ("tRNA",                                    "tRNA"),
    ("Mt_tRNA",                                 "tRNA"),
    ("pseudogene",                              "pseudogene"),
    ("processed_pseudogene",                    "pseudogene"),
    ("transcribed_unprocessed_pseudogene",      "pseudogene"),
    ("IG_C_gene",                               "IG_gene"),
    ("IG_V_gene",                               "IG_gene"),
    ("TR_C_gene",                               "TR_gene"),
    ("TR_J_gene",                               "TR_gene"),
    ("misc_RNA",                                "other_ncrna"),
    ("ribozyme",                                "other_ncrna"),
    ("vault_RNA",                               "other_ncrna"),
    ("TEC",                                     "other"),
])
def test_ensembl_mappings(raw, expected):
    assert normalize_biotype(raw, "Ensembl") == expected


# --- normalize_biotype: WormBase ---

@pytest.mark.parametrize("raw,expected", [
    ("protein_coding_gene", "protein_coding"),
    ("pseudogene",          "pseudogene"),
    ("lncRNA_gene",         "lncRNA"),
    ("miRNA_gene",          "miRNA"),
    ("snoRNA_gene",         "snoRNA"),
    ("snRNA_gene",          "snRNA"),
    ("rRNA_gene",           "rRNA"),
    ("tRNA_gene",           "tRNA"),
    ("ncRNA_gene",          "other_ncrna"),
    ("transposable_element_gene", "other"),
])
def test_wormbase_mappings(raw, expected):
    assert normalize_biotype(raw, "WormBase") == expected


# --- normalize_biotype: ZFIN (SO IDs) ---

@pytest.mark.parametrize("raw,expected", [
    ("SO:0001217", "protein_coding"),
    ("SO:0000336", "pseudogene"),
    ("SO:0001263", "miRNA"),
    ("SO:0001267", "snoRNA"),
    ("SO:0001268", "snRNA"),
    ("SO:0001637", "rRNA"),
    ("SO:0001272", "tRNA"),
    ("SO:0001877", "lncRNA"),
])
def test_zfin_mappings(raw, expected):
    assert normalize_biotype(raw, "ZFIN") == expected


# --- normalize_biotype: FlyBase (annotation_id prefix) ---

@pytest.mark.parametrize("raw,expected", [
    ("CG12345",  "protein_coding"),
    ("CR12345",  "other_ncrna"),
    ("FBtr0000", "other"),
    ("",         "unknown"),
])
def test_flybase_mappings(raw, expected):
    assert normalize_biotype(raw, "FlyBase") == expected


# --- unknown source / unknown value ---

def test_unknown_source_returns_unknown():
    assert normalize_biotype("protein-coding gene", "UnknownDB") == "unknown"


def test_empty_raw_returns_unknown():
    assert normalize_biotype("", "HGNC") == "unknown"


def test_whitespace_raw_returns_unknown():
    assert normalize_biotype("   ", "MGI") == "unknown"


# --- all return values are in vocabulary ---

@pytest.mark.parametrize("source,raw", [
    ("HGNC",      "gene with protein product"),
    ("MGI",       "miRNA gene"),
    ("RGD",       "lincrna"),
    ("Ensembl",   "processed_pseudogene"),
    ("WormBase",  "rRNA_gene"),
    ("ZFIN",      "SO:0001217"),
    ("FlyBase",   "CG9999"),
    ("HGNC",      "something_totally_new"),
])
def test_return_value_always_in_vocabulary(source, raw):
    result = normalize_biotype(raw, source)
    assert result in CANONICAL_BIOTYPES
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /scratch/users/chensj16/projects/stangene
python -m pytest tests/test_biotype.py -v 2>&1 | head -30
```

Expected: `ImportError` or `ModuleNotFoundError` for `stangene.biotype`.

- [ ] **Step 3: Implement `src/stangene/biotype.py`**

```python
"""Canonical gene biotype vocabulary and per-source normalization."""

from __future__ import annotations

CANONICAL_BIOTYPES: frozenset[str] = frozenset({
    "protein_coding",
    "lncRNA",
    "pseudogene",
    "miRNA",
    "snoRNA",
    "snRNA",
    "rRNA",
    "tRNA",
    "IG_gene",
    "TR_gene",
    "other_ncrna",
    "other",
    "unknown",
})

# HGNC: locus_type field (more granular than locus_group)
_HGNC_MAP: dict[str, str] = {
    "gene with protein product": "protein_coding",
    "rna, micro":                "miRNA",
    "rna, small nucleolar":      "snoRNA",
    "rna, small nuclear":        "snRNA",
    "rna, ribosomal":            "rRNA",
    "rna, transfer":             "tRNA",
    "rna, long non-coding":      "lncRNA",
    "rna, y":                    "other_ncrna",
    "rna, vault":                "other_ncrna",
    "rna, small cytoplasmic":    "other_ncrna",
    "immunoglobulin gene":       "IG_gene",
    "t cell receptor gene":      "TR_gene",
    "pseudogene":                "pseudogene",
    "endogenous retrovirus":     "other",
    "complex locus constituent": "other",
    "fragile site":              "other",
    "readthrough":               "other",
    "region":                    "other",
    "unknown":                   "unknown",
}

# MGI: feature_type field
_MGI_MAP: dict[str, str] = {
    "protein coding gene":               "protein_coding",
    "lncrna gene":                       "lncRNA",
    "pseudogene":                        "pseudogene",
    "mirna gene":                        "miRNA",
    "snorna gene":                       "snoRNA",
    "snrna gene":                        "snRNA",
    "rrna gene":                         "rRNA",
    "trna gene":                         "tRNA",
    "ig gene":                           "IG_gene",
    "tr gene":                           "TR_gene",
    "unclassified non-coding rna gene":  "other_ncrna",
    "pirna gene":                        "other_ncrna",
    "scrna gene":                        "other_ncrna",
    "rnase mrp rna gene":                "other_ncrna",
    "rnase p rna gene":                  "other_ncrna",
    "srp rna gene":                      "other_ncrna",
    "telomerase rna gene":               "other_ncrna",
    "antisense lncrna gene":             "lncRNA",
    "bidirectional promoter lncrna":     "lncRNA",
    "transgene":                         "other",
    "unclassified gene":                 "other",
    "gene segment":                      "other",
    "other feature types":               "other",
}

# RGD: GENE_TYPE field (lowercase in data)
_RGD_MAP: dict[str, str] = {
    "protein-coding":  "protein_coding",
    "ncrna":           "other_ncrna",
    "lincrna":         "lncRNA",
    "lncrna":          "lncRNA",
    "mirna":           "miRNA",
    "snorna":          "snoRNA",
    "snrna":           "snRNA",
    "rrna":            "rRNA",
    "trna":            "tRNA",
    "pseudo":          "pseudogene",
    "pseudogene":      "pseudogene",
    "gene":            "other",
}

# Ensembl BioMart: gene_biotype field
# Pseudogene subtypes all end in "pseudogene" — handled by suffix fallback below.
# IG/TR subtypes all match IG_*_gene / TR_*_gene — handled by prefix fallback.
_ENSEMBL_MAP: dict[str, str] = {
    "protein_coding":   "protein_coding",
    "lincrna":          "lncRNA",
    "lncrna":           "lncRNA",
    "mirna":            "miRNA",
    "snorna":           "snoRNA",
    "snrna":            "snRNA",
    "rrna":             "rRNA",
    "mt_rrna":          "rRNA",
    "trna":             "tRNA",
    "mt_trna":          "tRNA",
    "pseudogene":       "pseudogene",
    "misc_rna":         "other_ncrna",
    "ribozyme":         "other_ncrna",
    "vault_rna":        "other_ncrna",
    "scrna":            "other_ncrna",
    "pirna":            "other_ncrna",
    "scarna":           "other_ncrna",
    "known_ncrna":      "other_ncrna",
    "processed_transcript": "other_ncrna",
    "tec":              "other",
    "artifact":         "other",
    "disrupted_domain": "other",
}

# WormBase: biotype field
_WORMBASE_MAP: dict[str, str] = {
    "protein_coding_gene":       "protein_coding",
    "pseudogene":                "pseudogene",
    "lncrna_gene":               "lncRNA",
    "mirna_gene":                "miRNA",
    "snorna_gene":               "snoRNA",
    "snrna_gene":                "snRNA",
    "rrna_gene":                 "rRNA",
    "trna_gene":                 "tRNA",
    "ncrna_gene":                "other_ncrna",
    "pirna_gene":                "other_ncrna",
    "scrna_gene":                "other_ncrna",
    "transposable_element_gene": "other",
    "operon":                    "other",
}

# ZFIN: SO (Sequence Ontology) ID, e.g. "SO:0001217"
_ZFIN_SO_MAP: dict[str, str] = {
    "so:0001217": "protein_coding",   # protein_coding_gene
    "so:0000336": "pseudogene",       # pseudogene
    "so:0001263": "miRNA",            # miRNA_gene
    "so:0001267": "snoRNA",           # snoRNA_gene
    "so:0001268": "snRNA",            # snRNA_gene
    "so:0001637": "rRNA",             # rRNA_gene
    "so:0001272": "tRNA",             # tRNA_gene
    "so:0001877": "lncRNA",           # lncRNA_gene
    "so:0001900": "other_ncrna",      # ncRNA_gene
    "so:0001265": "other_ncrna",      # miRNA_primary_transcript (pre-miRNA)
    "so:0001032": "other_ncrna",      # other ncRNA
    "so:0000704": "other",            # gene (unspecified)
    "so:0000233": "other",            # processed_transcript
}

_SOURCE_MAPS: dict[str, dict[str, str]] = {
    "HGNC":      _HGNC_MAP,
    "MGI":       _MGI_MAP,
    "RGD":       _RGD_MAP,
    "Ensembl":   _ENSEMBL_MAP,
    "WormBase":  _WORMBASE_MAP,
    "ZFIN":      _ZFIN_SO_MAP,
}


def normalize_biotype(raw: str, source: str) -> str:
    """Map a raw gene_type string to a canonical biotype.

    Args:
        raw: The raw gene_type string from the upstream database.
        source: One of "HGNC", "MGI", "RGD", "Ensembl", "WormBase",
                "ZFIN", or "FlyBase".

    Returns:
        A member of CANONICAL_BIOTYPES. Unmapped values return "unknown".
    """
    key = raw.strip().lower()
    if not key:
        return "unknown"

    # FlyBase uses annotation_id prefix convention (no dict lookup needed)
    if source == "FlyBase":
        if key.startswith("cg"):
            return "protein_coding"
        if key.startswith("cr"):
            return "other_ncrna"
        if key:
            return "other"
        return "unknown"

    source_map = _SOURCE_MAPS.get(source)
    if source_map is None:
        return "unknown"

    result = source_map.get(key)
    if result is not None:
        return result

    # Ensembl fallbacks: many pseudogene subtypes end in "pseudogene";
    # IG/TR gene subtypes start with ig_ or tr_ and end with _gene.
    if source == "Ensembl":
        if key.endswith("pseudogene"):
            return "pseudogene"
        if key.startswith("ig_") and key.endswith("_gene"):
            return "IG_gene"
        if key.startswith("tr_") and key.endswith("_gene"):
            return "TR_gene"

    return "unknown"
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_biotype.py -v
```

Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/stangene/biotype.py tests/test_biotype.py
git commit -m "feat(biotype): add canonical biotype vocabulary and normalize_biotype()"
```

---

## Task 2: Update mammal + primate + worm builders

Add `canonical_biotype` to human, mouse, rat, c_elegans, cynomolgus, rhesus, marmoset, mouse_lemur.

**Files:**
- Modify: `src/stangene/references.py`
  - `_build_human_reference()`: switch `gene_type` to use `locus_type` (not `locus_group`); add `canonical_biotype`
  - `_build_mouse_reference()`: add `canonical_biotype` after gene_table build
  - `_build_rat_reference()`: add `canonical_biotype` after gene_table build
  - `_build_celegans_reference()`: add `canonical_biotype` after gene_table build
  - `_build_ensembl_biomart_reference()`: add `canonical_biotype` after gene_table build

**Interfaces:**
- Consumes: `normalize_biotype(raw: str, source: str) -> str` from Task 1

---

- [ ] **Step 1: Write integration test for canonical_biotype presence**

Add to `tests/test_biotype.py`:

```python
import os
import pandas as pd
from stangene.biotype import CANONICAL_BIOTYPES

_REFS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "src", "stangene", "data", "refs"
)


def _gene_table(species: str) -> pd.DataFrame:
    path = os.path.join(_REFS_DIR, species, "gene_table.parquet")
    if not os.path.exists(path):
        pytest.skip(f"Reference not built for {species}")
    return pd.read_parquet(path)


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_column_exists(species):
    gt = _gene_table(species)
    assert "canonical_biotype" in gt.columns, f"{species}: missing canonical_biotype"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_no_nulls(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    assert gt["canonical_biotype"].isna().sum() == 0, f"{species}: NaN in canonical_biotype"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_in_vocabulary(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    bad = set(gt["canonical_biotype"].unique()) - CANONICAL_BIOTYPES
    assert not bad, f"{species}: out-of-vocabulary values: {bad}"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat",
])
def test_canonical_biotype_has_protein_coding(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    n = (gt["canonical_biotype"] == "protein_coding").sum()
    assert n > 10_000, f"{species}: unexpectedly few protein_coding genes ({n})"
```

- [ ] **Step 2: Run new integration tests to verify they fail**

```bash
python -m pytest tests/test_biotype.py -v -k "canonical_biotype_column"
```

Expected: FAIL — `canonical_biotype` column absent from parquets.

- [ ] **Step 3: Update `_build_human_reference()` in `references.py`**

Find the block that builds `gene_table` in `_build_human_reference()` (around line 185). Replace the `gene_type` line and add canonical_biotype:

```python
# At top of references.py, add import alongside existing imports:
from stangene.biotype import normalize_biotype

# Inside _build_human_reference(), replace:
#   "gene_type": hgnc.get("locus_group", pd.Series(dtype=str)).fillna(""),
# with:
#   "gene_type": hgnc.get("locus_type", pd.Series(dtype=str)).fillna(""),
# Then after the gene_table = pd.DataFrame({...}) block, add:

gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="HGNC"), axis=1
)
```

The full updated `gene_table` construction in `_build_human_reference()`:

```python
gene_table = pd.DataFrame({
    "ensembl_id": hgnc["ensembl_gene_id"].where(hgnc["ensembl_gene_id"].notna(), None),
    "symbol": hgnc["symbol"],
    "alias_symbols": hgnc["alias_symbol"].fillna(""),
    "prev_symbols": hgnc["prev_symbol"].fillna(""),
    "gene_type": hgnc.get("locus_type", pd.Series(dtype=str)).fillna(""),
    "status": hgnc["status"].fillna(""),
    "source": "HGNC",
    "source_id": hgnc["hgnc_id"],
})
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="HGNC"), axis=1
)
```

- [ ] **Step 4: Update `_build_mouse_reference()` in `references.py`**

After the `gene_table = pd.DataFrame(rows)` line and before the MGI filter, add:

```python
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="MGI"), axis=1
)
```

- [ ] **Step 5: Update `_build_rat_reference()` in `references.py`**

After `gene_table = pd.DataFrame(rows)`, add:

```python
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="RGD"), axis=1
)
```

- [ ] **Step 6: Update `_build_celegans_reference()` in `references.py`**

After `gene_table = pd.DataFrame(rows)`, add:

```python
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="WormBase"), axis=1
)
```

- [ ] **Step 7: Update `_build_ensembl_biomart_reference()` in `references.py`**

After `gene_table = pd.DataFrame(rows)`, add:

```python
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="Ensembl"), axis=1
)
```

- [ ] **Step 8: Rebuild the 8 affected species' parquets**

Run from an interactive session (not the login node):

```bash
sh_dev -c 2 --mem 8GB -t 0:30:00
cd /scratch/users/chensj16/projects/stangene
pip install -e . -q
python - <<'EOF'
import stangene
for sp in ["human", "mouse", "rat", "c_elegans",
           "cynomolgus", "rhesus", "marmoset", "mouse_lemur"]:
    print(f"Building {sp}...")
    stangene.build_reference(sp, force=True)
    print(f"  done")
EOF
```

Expected: no exceptions, logs show gene counts similar to before.

- [ ] **Step 9: Run integration tests**

```bash
python -m pytest tests/test_biotype.py -v -k "canonical_biotype"
```

Expected: all PASS for the 8 species.

- [ ] **Step 10: Run full test suite**

```bash
python -m pytest --tb=short -q
```

Expected: all existing tests still PASS.

- [ ] **Step 11: Commit**

```bash
git add src/stangene/references.py
git commit -m "feat(biotype): add canonical_biotype to mammal + primate + worm builders"
```

---

## Task 3: Update zebrafish builder (preserve SO ID)

**Files:**
- Modify: `src/stangene/references.py` — `_build_zebrafish_reference()`

**Interfaces:**
- Consumes: `normalize_biotype(raw: str, source: str) -> str` from Task 1
- The `so_id` column is already read into `genes_df` (column named `"so_id"`) — it is just not used. No change to the CSV parsing is needed.

---

- [ ] **Step 1: Add integration tests for zebrafish**

Add to `tests/test_biotype.py`:

```python
def test_zebrafish_canonical_biotype_column_exists():
    gt = _gene_table("zebrafish")
    assert "canonical_biotype" in gt.columns


def test_zebrafish_canonical_biotype_no_nulls():
    gt = _gene_table("zebrafish")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    assert gt["canonical_biotype"].isna().sum() == 0


def test_zebrafish_canonical_biotype_in_vocabulary():
    gt = _gene_table("zebrafish")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    bad = set(gt["canonical_biotype"].unique()) - CANONICAL_BIOTYPES
    assert not bad, f"out-of-vocabulary: {bad}"


def test_zebrafish_has_protein_coding():
    gt = _gene_table("zebrafish")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    n = (gt["canonical_biotype"] == "protein_coding").sum()
    assert n > 5_000, f"unexpectedly few protein_coding genes: {n}"
```

- [ ] **Step 2: Verify tests fail**

```bash
python -m pytest tests/test_biotype.py -v -k "zebrafish"
```

Expected: FAIL — column absent.

- [ ] **Step 3: Update `_build_zebrafish_reference()` in `references.py`**

In the loop `for _, g in genes_df.iterrows()`, `so_id` is already available as `g["so_id"]`. Pass it to `normalize_biotype`:

```python
# Replace the existing rows.append({...}) block in _build_zebrafish_reference()
# with this version that captures so_id:

so_id_raw = str(g["so_id"]).strip() if pd.notna(g["so_id"]) else ""

rows.append({
    "ensembl_id": ensembl_id,
    "symbol": symbol,
    "alias_symbols": alias_symbols,
    "prev_symbols": prev_symbols,
    "gene_type": so_id_raw,       # store the SO ID as gene_type for traceability
    "status": "approved",
    "source": "ZFIN",
    "source_id": f"ZFIN:{zid}",
})

# After gene_table = pd.DataFrame(rows), add:
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="ZFIN"), axis=1
)
```

- [ ] **Step 4: Rebuild zebrafish parquet**

```bash
# From sh_dev interactive session:
python -c "
import stangene
stangene.build_reference('zebrafish', force=True)
print('done')
"
```

- [ ] **Step 5: Run zebrafish tests**

```bash
python -m pytest tests/test_biotype.py -v -k "zebrafish"
```

Expected: all PASS.

- [ ] **Step 6: Run full suite**

```bash
python -m pytest --tb=short -q
```

Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add src/stangene/references.py
git commit -m "feat(biotype): zebrafish builder preserves SO ID for canonical_biotype mapping"
```

---

## Task 4: Update fruit_fly builder (annotation_id prefix)

FlyBase protein-coding genes have annotation IDs starting with `CG`; non-coding RNA genes start with `CR`. This convention is stable and documented by FlyBase. No additional download is needed.

**Files:**
- Modify: `src/stangene/references.py` — `_build_fruitfly_reference()`

---

- [ ] **Step 1: Add integration tests for fruit_fly**

Add to `tests/test_biotype.py`:

```python
def test_fruitfly_canonical_biotype_column_exists():
    gt = _gene_table("fruit_fly")
    assert "canonical_biotype" in gt.columns


def test_fruitfly_canonical_biotype_no_nulls():
    gt = _gene_table("fruit_fly")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    assert gt["canonical_biotype"].isna().sum() == 0


def test_fruitfly_canonical_biotype_in_vocabulary():
    gt = _gene_table("fruit_fly")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    bad = set(gt["canonical_biotype"].unique()) - CANONICAL_BIOTYPES
    assert not bad, f"out-of-vocabulary: {bad}"


def test_fruitfly_has_protein_coding():
    gt = _gene_table("fruit_fly")
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    n = (gt["canonical_biotype"] == "protein_coding").sum()
    assert n > 8_000, f"unexpectedly few protein_coding genes: {n}"
```

- [ ] **Step 2: Verify tests fail**

```bash
python -m pytest tests/test_biotype.py -v -k "fruitfly"
```

Expected: FAIL — column absent.

- [ ] **Step 3: Update `_build_fruitfly_reference()` in `references.py`**

Find the `rows.append({...})` block. Add annotation_id extraction and store it as `gene_type`:

The loop iterates `map_df`. The annotation_id column is the 5th column in the FlyBase TSV (`annotation_ID`). Locate it with:

```python
# After the existing column detection lines in _build_fruitfly_reference(),
# add annotation_id column detection:
annotation_col = next((c for c in map_df.columns if "annotation_id" in c.lower()
                       and "secondary" not in c.lower()), None)
```

Then inside the loop, replace the existing `rows.append` with:

```python
annotation_id = ""
if annotation_col and pd.notna(m[annotation_col]):
    annotation_id = str(m[annotation_col]).strip()

rows.append({
    "ensembl_id": fbgn,
    "symbol": symbol,
    "alias_symbols": "|".join(alias_syms),
    "prev_symbols": "|".join(prev_fbgns),
    "gene_type": annotation_id,     # CG* = protein_coding, CR* = other_ncrna
    "status": "approved",
    "source": "FlyBase",
    "source_id": f"FlyBase:{fbgn}",
})
```

After `gene_table = pd.DataFrame(rows)`, add:

```python
gene_table["canonical_biotype"] = gene_table.apply(
    lambda r: normalize_biotype(r["gene_type"], source="FlyBase"), axis=1
)
```

- [ ] **Step 4: Rebuild fruit_fly parquet**

```bash
# From sh_dev interactive session:
python -c "
import stangene
stangene.build_reference('fruit_fly', force=True)
print('done')
"
```

- [ ] **Step 5: Run fruit_fly tests**

```bash
python -m pytest tests/test_biotype.py -v -k "fruitfly"
```

Expected: all PASS.

- [ ] **Step 6: Run full suite**

```bash
python -m pytest --tb=short -q
```

Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add src/stangene/references.py
git commit -m "feat(biotype): fruit_fly builder uses annotation_id prefix for canonical_biotype"
```

---

## Task 5: Update `harmonize.py` and `__init__.py`

**Files:**
- Modify: `src/stangene/harmonize.py` lines 153–155
- Modify: `src/stangene/__init__.py`

---

- [ ] **Step 1: Update `harmonize.py` lines 153–155**

Replace:
```python
gene_type = gene_info.get("gene_type", "")
if gene_type and "protein" not in str(gene_type).lower():
    notes.append(f"Non-protein-coding gene type: {gene_type}")
```

With:
```python
canonical = gene_info.get("canonical_biotype", "unknown")
if canonical not in ("protein_coding", "unknown"):
    notes.append(f"Non-protein-coding: {canonical}")
```

- [ ] **Step 2: Update `__init__.py`**

Add after the `from stangene.hb import hb_mask` line:

```python
from stangene.biotype import CANONICAL_BIOTYPES, normalize_biotype
```

- [ ] **Step 3: Verify the import works**

```bash
python -c "import stangene; print(stangene.CANONICAL_BIOTYPES)"
```

Expected: prints the frozenset of 13 values.

- [ ] **Step 4: Run harmonize tests**

```bash
python -m pytest tests/test_harmonize.py tests/test_harmonize_anndata.py -v --tb=short
```

Expected: all PASS (note text changed but no schema change).

- [ ] **Step 5: Run full suite**

```bash
python -m pytest --tb=short -q
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add src/stangene/harmonize.py src/stangene/__init__.py
git commit -m "feat(biotype): expose canonical_biotype API; use it in harmonize notes"
```

---

## Task 6: Commit regenerated parquet files + version bump

All 10 species' parquets now have the `canonical_biotype` column. Commit them and bump the version.

**Files:**
- Modify: `src/stangene/__init__.py` — `__version__`
- Modify: `src/stangene/data/refs/*/gene_table.parquet` (all 10 species)

---

- [ ] **Step 1: Confirm all 10 parquets have the column**

```bash
python - <<'EOF'
import pandas as pd, os
base = "src/stangene/data/refs"
for sp in os.listdir(base):
    p = f"{base}/{sp}/gene_table.parquet"
    if os.path.exists(p):
        gt = pd.read_parquet(p)
        has = "canonical_biotype" in gt.columns
        print(f"{sp:15s}  canonical_biotype={'YES' if has else 'NO ':3s}  rows={len(gt):,}")
EOF
```

Expected: all 10 species show `canonical_biotype=YES`.

- [ ] **Step 2: Bump version to 0.4.0**

In `src/stangene/__init__.py`, change:
```python
__version__ = "0.3.0"
```
to:
```python
__version__ = "0.4.0"
```

- [ ] **Step 3: Run full test suite one final time**

```bash
python -m pytest --tb=short -q
```

Expected: all PASS.

- [ ] **Step 4: Commit everything**

```bash
git add src/stangene/__init__.py src/stangene/data/refs/
git commit -m "feat(biotype): regenerate all 10 species parquets with canonical_biotype column, v0.4.0"
```
