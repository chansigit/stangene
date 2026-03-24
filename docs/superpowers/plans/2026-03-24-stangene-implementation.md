# Stangene Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `stangene` Python package — a gene identifier harmonization pipeline for single-cell transcriptomics datasets that preserves original information and tracks mapping provenance.

**Architecture:** Staged pipeline with composable modules: `io.py` (load/write) → `classify.py` (feature triage) → `references.py` (build/load annotation DBs) → `harmonize.py` (tiered matching cascade) → `merge.py` (conservative opt-in merge) → `report.py` (summaries/conflicts). A top-level `run()` function wires them together. CLI via `__main__.py`.

**Tech Stack:** Python 3.9+, pandas, anndata, pyarrow, argparse. No external HTTP library (stdlib `urllib`).

**Spec:** `docs/superpowers/specs/2026-03-24-gene-harmonization-design.md`

---

## File Map

| File | Responsibility |
|---|---|
| `pyproject.toml` | Package metadata, dependencies, CLI entry point |
| `src/stangene/__init__.py` | Public API: `run()`, `__version__`, re-exports |
| `src/stangene/species.py` | `SpeciesConfig` dataclass, human/mouse configs, classification patterns |
| `src/stangene/_logging.py` | Structured logger setup for the package |
| `src/stangene/io.py` | `load_features()`, `write_results()` — h5ad and TSV adapters |
| `src/stangene/classify.py` | `classify_features()` — feature type triage |
| `src/stangene/references.py` | `build_reference()`, `load_reference()`, `ReferenceNotFoundError` |
| `src/stangene/harmonize.py` | `harmonize()` — the 5-tier matching cascade, `HarmonizationResult` |
| `src/stangene/merge.py` | `merge_features()` — conservative merge, `MergeResult` |
| `src/stangene/report.py` | `summary()`, `conflict_report()`, `write_reports()` |
| `src/stangene/__main__.py` | CLI: `harmonize` and `build-refs` subcommands |
| `tests/conftest.py` | Shared pytest fixtures: mock reference data, sample DataFrames |
| `tests/test_species.py` | Tests for species configs |
| `tests/test_io.py` | Tests for load_features / write_results |
| `tests/test_classify.py` | Tests for classify_features |
| `tests/test_references.py` | Tests for build/load reference (with mocked downloads) |
| `tests/test_harmonize.py` | Tests for the matching cascade |
| `tests/test_merge.py` | Tests for merge logic |
| `tests/test_report.py` | Tests for summary/conflict/write_reports |
| `tests/test_run.py` | Integration test for run() end-to-end |
| `tests/fixtures/` | Small h5ad and TSV test files |

---

### Task 1: Project Scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `src/stangene/__init__.py`
- Create: `src/stangene/_logging.py`
- Create: `.gitignore`

- [ ] **Step 1: Create pyproject.toml**

```toml
[build-system]
requires = ["setuptools>=68.0", "setuptools-scm"]
build-backend = "setuptools.backends._legacy:_Backend"

[project]
name = "stangene"
version = "0.1.0"
description = "Gene identifier harmonization for single-cell transcriptomics"
requires-python = ">=3.9"
dependencies = [
    "pandas>=1.5",
    "anndata>=0.8",
    "pyarrow>=10.0",
]

[project.optional-dependencies]
dev = ["pytest>=7.0", "scanpy>=1.9"]

[project.scripts]
stangene = "stangene.__main__:main"

[tool.setuptools.packages.find]
where = ["src"]
```

- [ ] **Step 2: Create src/stangene/__init__.py**

```python
"""Stangene: gene identifier harmonization for single-cell transcriptomics."""

__version__ = "0.1.0"
```

The `run()` function and re-exports will be added in Task 10 after all modules exist.

- [ ] **Step 3: Create src/stangene/_logging.py**

```python
"""Structured logging for stangene."""

import logging


def get_logger(name: str) -> logging.Logger:
    """Get a logger under the stangene namespace."""
    logger = logging.getLogger(f"stangene.{name}")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(name)s] %(levelname)s: %(message)s")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
```

- [ ] **Step 4: Create .gitignore**

```
references/
*.egg-info/
__pycache__/
*.pyc
dist/
build/
.eggs/
*.parquet
```

- [ ] **Step 5: Install in dev mode and verify**

Run: `cd /scratch/users/chensj16/projects/stangene && pip install -e ".[dev]" 2>&1 | tail -3`
Expected: "Successfully installed stangene-0.1.0"

Run: `python -c "import stangene; print(stangene.__version__)"`
Expected: `0.1.0`

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml src/ .gitignore
git commit -m "feat: scaffold stangene package with pyproject.toml and logging"
```

---

### Task 2: Species Config (`species.py`)

**Files:**
- Create: `src/stangene/species.py`
- Create: `tests/test_species.py`

- [ ] **Step 1: Write failing test**

```python
# tests/test_species.py
from stangene.species import SpeciesConfig, get_species_config, CLASSIFICATION_PATTERNS


def test_get_human_config():
    cfg = get_species_config("human")
    assert cfg.name == "human"
    assert cfg.ensembl_prefix == "ENSG"
    assert cfg.transcript_prefix == "ENST"
    assert cfg.naming_convention == "uppercase"


def test_get_mouse_config():
    cfg = get_species_config("mouse")
    assert cfg.name == "mouse"
    assert cfg.ensembl_prefix == "ENSMUSG"
    assert cfg.transcript_prefix == "ENSMUST"
    assert cfg.naming_convention == "capitalized"


def test_unknown_species_raises():
    import pytest
    with pytest.raises(ValueError, match="Unknown species"):
        get_species_config("zebrafish")


def test_classification_patterns_detect_antibody():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "CD3_TotalSeqB"):
            assert ftype == "antibody_capture"
            return
    raise AssertionError("No pattern matched CD3_TotalSeqB")


def test_classification_patterns_detect_spike_in():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "ERCC-00002"):
            assert ftype == "spike_in"
            return
    raise AssertionError("No pattern matched ERCC-00002")


def test_classification_patterns_detect_peak():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "chr1:1000-2000"):
            assert ftype == "peak"
            return
    raise AssertionError("No pattern matched chr1:1000-2000")


def test_classification_patterns_detect_crispr():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "sg-TP53"):
            assert ftype == "crispr_guide"
            return
    raise AssertionError("No pattern matched sg-TP53")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_species.py -v 2>&1 | tail -5`
Expected: FAIL — `ModuleNotFoundError` or `ImportError`

- [ ] **Step 3: Implement species.py**

```python
# src/stangene/species.py
"""Species-specific configuration and classification patterns."""

import re
from dataclasses import dataclass, field


@dataclass
class SpeciesConfig:
    """Configuration for a species' gene naming and reference sources."""

    name: str
    ensembl_prefix: str
    transcript_prefix: str
    naming_convention: str  # "uppercase" (human) or "capitalized" (mouse)
    reference_sources: dict = field(default_factory=dict)


# Classification patterns: list of (compiled_regex_pattern, feature_type).
# Order matters — first match wins. These are checked against original_feature_name.
CLASSIFICATION_PATTERNS: list[tuple[re.Pattern, str]] = [
    # Transcript IDs (check before gene IDs since ENST is a subset prefix-wise)
    (re.compile(r"^ENST\d+", re.IGNORECASE), "transcript"),
    (re.compile(r"^ENSMUST\d+", re.IGNORECASE), "transcript"),
    # Ensembl gene IDs
    (re.compile(r"^ENSG\d+", re.IGNORECASE), "gene"),
    (re.compile(r"^ENSMUSG\d+", re.IGNORECASE), "gene"),
    # Antibody capture / protein tags
    (re.compile(r".*TotalSeq", re.IGNORECASE), "antibody_capture"),
    (re.compile(r".*_ADT$", re.IGNORECASE), "antibody_capture"),
    (re.compile(r".*_HTO$", re.IGNORECASE), "antibody_capture"),
    # CRISPR guides
    (re.compile(r"^sg-", re.IGNORECASE), "crispr_guide"),
    (re.compile(r"^gRNA-", re.IGNORECASE), "crispr_guide"),
    # Spike-ins
    (re.compile(r"^ERCC-", re.IGNORECASE), "spike_in"),
    # Genomic peaks (ATAC-seq style)
    (re.compile(r"^chr[\dXYMT]+:\d+-\d+$", re.IGNORECASE), "peak"),
]

# Excel-corruption date-like patterns that indicate corrupted gene names
EXCEL_DATE_PATTERN = re.compile(
    r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
    r"|^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{1,2}$"
    r"|^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{2,4}$",
    re.IGNORECASE,
)

# Gene symbols renamed by HGNC in 2020 due to Excel auto-conversion.
# Maps old symbol -> new symbol for flagging purposes.
EXCEL_RENAMED_GENES: dict[str, str] = {
    "MARCH1": "MARCHF1", "MARCH2": "MARCHF2", "MARCH3": "MARCHF3",
    "MARCH4": "MARCHF4", "MARCH5": "MARCHF5", "MARCH6": "MARCHF6",
    "MARCH7": "MARCHF7", "MARCH8": "MARCHF8", "MARCH9": "MARCHF9",
    "MARCH10": "MARCHF10", "MARCH11": "MARCHF11",
    "SEPT1": "SEPTIN1", "SEPT2": "SEPTIN2", "SEPT3": "SEPTIN3",
    "SEPT4": "SEPTIN4", "SEPT5": "SEPTIN5", "SEPT6": "SEPTIN6",
    "SEPT7": "SEPTIN7", "SEPT8": "SEPTIN8", "SEPT9": "SEPTIN9",
    "SEPT10": "SEPTIN10", "SEPT11": "SEPTIN11", "SEPT12": "SEPTIN12",
    "SEPT14": "SEPTIN14",
    "DEC1": "DELEC1", "DEC2": "BHLHE41",
}


_SPECIES_CONFIGS: dict[str, SpeciesConfig] = {
    "human": SpeciesConfig(
        name="human",
        ensembl_prefix="ENSG",
        transcript_prefix="ENST",
        naming_convention="uppercase",
        reference_sources={
            "hgnc": {
                "url": "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
                "description": "HGNC complete gene set with symbols, aliases, and Ensembl IDs",
            },
        },
    ),
    "mouse": SpeciesConfig(
        name="mouse",
        ensembl_prefix="ENSMUSG",
        transcript_prefix="ENSMUST",
        naming_convention="capitalized",
        reference_sources={
            "mgi_markers": {
                "url": "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt",
                "description": "MGI marker list with approved symbols and synonyms",
            },
            "mgi_ensembl": {
                "url": "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt",
                "description": "MGI to Ensembl ID mapping",
            },
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "description": "Ensembl BioMart mouse gene table (supplementary)",
            },
        },
    ),
}


def get_species_config(species: str) -> SpeciesConfig:
    """Get configuration for a species. Raises ValueError if unknown."""
    species_lower = species.lower()
    if species_lower not in _SPECIES_CONFIGS:
        raise ValueError(
            f"Unknown species: '{species}'. Supported: {list(_SPECIES_CONFIGS.keys())}"
        )
    return _SPECIES_CONFIGS[species_lower]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_species.py -v 2>&1 | tail -15`
Expected: All 7 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/species.py tests/test_species.py
git commit -m "feat: add species config with human/mouse and classification patterns"
```

---

### Task 3: Feature Classification (`classify.py`)

**Files:**
- Create: `src/stangene/classify.py`
- Create: `tests/test_classify.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_classify.py
import pandas as pd
from stangene.classify import classify_features


def _make_ft(names, feature_types=None, ids=None):
    """Helper to build a minimal FeatureTable."""
    data = {"original_feature_name": names, "species": "human", "dataset": "test"}
    if feature_types is not None:
        data["original_feature_type"] = feature_types
    if ids is not None:
        data["original_feature_id"] = ids
    return pd.DataFrame(data)


def test_explicit_labels_preserved():
    ft = _make_ft(
        ["TP53", "CD3_TotalSeqB"],
        feature_types=["Gene Expression", "Antibody Capture"],
    )
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "gene"
    assert result.loc[1, "original_feature_type"] == "antibody_capture"


def test_pattern_gene_ensembl():
    ft = _make_ft(["ENSG00000141510"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "gene"


def test_pattern_transcript():
    ft = _make_ft(["ENST00000269305"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "transcript"


def test_pattern_antibody():
    ft = _make_ft(["CD3_ADT"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "antibody_capture"


def test_pattern_spike_in():
    ft = _make_ft(["ERCC-00002"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "spike_in"


def test_pattern_peak():
    ft = _make_ft(["chr1:1000-2000"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "peak"


def test_pattern_crispr():
    ft = _make_ft(["sg-TP53"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "crispr_guide"


def test_default_gene():
    ft = _make_ft(["TP53"])
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "gene"


def test_non_gene_gets_mapping_status():
    ft = _make_ft(["ERCC-00002", "CD3_ADT", "chr1:100-200"])
    result = classify_features(ft)
    for i in range(3):
        assert result.loc[i, "mapping_status"] == "non_gene_feature"


def test_gene_features_no_mapping_status_yet():
    ft = _make_ft(["TP53"])
    result = classify_features(ft)
    assert pd.isna(result.loc[0, "mapping_status"])


def test_mixed_explicit_and_heuristic():
    ft = _make_ft(
        ["GeneA", "ERCC-00001", "GeneB"],
        feature_types=["Gene Expression", None, None],
    )
    result = classify_features(ft)
    assert result.loc[0, "original_feature_type"] == "gene"
    assert result.loc[1, "original_feature_type"] == "spike_in"
    assert result.loc[2, "original_feature_type"] == "gene"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_classify.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 3: Implement classify.py**

```python
# src/stangene/classify.py
"""Feature type classification for gene vs non-gene triage."""

import pandas as pd

from stangene._logging import get_logger
from stangene.species import CLASSIFICATION_PATTERNS

logger = get_logger("classify")

# Map 10x Cell Ranger feature_type labels to our internal types
_CELLRANGER_TYPE_MAP = {
    "Gene Expression": "gene",
    "Antibody Capture": "antibody_capture",
    "CRISPR Guide Capture": "crispr_guide",
    "Custom": "custom",
    "Peaks": "peak",
}

_NON_GENE_TYPES = frozenset([
    "transcript", "antibody_capture", "crispr_guide",
    "spike_in", "peak", "custom",
])


def classify_features(ft: pd.DataFrame) -> pd.DataFrame:
    """Classify features as gene or non-gene types.

    Adds/updates 'original_feature_type' and sets 'mapping_status' to
    'non_gene_feature' for non-gene rows. Returns a copy.
    """
    result = ft.copy()

    # Ensure columns exist
    if "original_feature_type" not in result.columns:
        result["original_feature_type"] = None
    if "mapping_status" not in result.columns:
        result["mapping_status"] = None
    if "mapping_notes" not in result.columns:
        result["mapping_notes"] = None

    for idx in result.index:
        existing_type = result.at[idx, "original_feature_type"]

        if pd.notna(existing_type) and existing_type in _CELLRANGER_TYPE_MAP:
            # Step 1: Trust explicit labels, normalize them
            result.at[idx, "original_feature_type"] = _CELLRANGER_TYPE_MAP[existing_type]
        elif pd.notna(existing_type) and existing_type in _NON_GENE_TYPES | {"gene"}:
            # Already in our internal format
            pass
        else:
            # Step 2: Pattern-based heuristic
            name = result.at[idx, "original_feature_name"]
            matched = False
            for pattern, ftype in CLASSIFICATION_PATTERNS:
                if pattern.match(str(name)):
                    result.at[idx, "original_feature_type"] = ftype
                    matched = True
                    break
            if not matched:
                # Default to gene, flag it
                result.at[idx, "original_feature_type"] = "gene"
                result.at[idx, "mapping_notes"] = "classified as gene by default (no pattern match)"

        # Step 3/4: Mark non-gene features
        if result.at[idx, "original_feature_type"] in _NON_GENE_TYPES:
            result.at[idx, "mapping_status"] = "non_gene_feature"

    classified_counts = result["original_feature_type"].value_counts().to_dict()
    logger.info("Feature classification: %s", classified_counts)
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_classify.py -v 2>&1 | tail -20`
Expected: All 12 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/classify.py tests/test_classify.py
git commit -m "feat: add feature classification with pattern-based triage"
```

---

### Task 4: Input Adapters (`io.py`)

**Files:**
- Create: `src/stangene/io.py`
- Create: `tests/test_io.py`
- Create: `tests/fixtures/` (test data)

- [ ] **Step 1: Create test fixtures**

Create `tests/fixtures/sample_features.tsv`:
```tsv
gene_name	gene_id
TP53	ENSG00000141510.18
BRCA1	ENSG00000012048
CD3_ADT
ERCC-00002
```

Create a small h5ad fixture programmatically in the test file (avoids committing binary).

- [ ] **Step 2: Write failing tests**

```python
# tests/test_io.py
import os
import tempfile

import anndata
import numpy as np
import pandas as pd
import pytest

from stangene.io import load_features


@pytest.fixture
def sample_h5ad(tmp_path):
    """Create a minimal h5ad file for testing."""
    n_obs, n_vars = 3, 5
    adata = anndata.AnnData(
        X=np.zeros((n_obs, n_vars)),
        var=pd.DataFrame(
            {
                "gene_ids": [
                    "ENSG00000141510.18",
                    "ENSG00000012048",
                    "",
                    "",
                    "ENSG00000139618.15",
                ],
                "feature_types": [
                    "Gene Expression",
                    "Gene Expression",
                    "Antibody Capture",
                    "Gene Expression",
                    "Gene Expression",
                ],
            },
            index=["TP53", "BRCA1", "CD3_ADT", "MYC", "BRCA2"],
        ),
    )
    path = str(tmp_path / "test.h5ad")
    adata.write_h5ad(path)
    return path


@pytest.fixture
def sample_tsv(tmp_path):
    """Create a minimal TSV file for testing."""
    path = str(tmp_path / "features.tsv")
    pd.DataFrame({
        "gene_name": ["TP53", "BRCA1", "CD3_ADT", "ERCC-00002"],
        "gene_id": ["ENSG00000141510.18", "ENSG00000012048", "", ""],
    }).to_csv(path, sep="\t", index=False)
    return path


def test_load_h5ad_basic(sample_h5ad):
    ft = load_features(sample_h5ad, species="human")
    assert len(ft) == 5
    assert "original_feature_name" in ft.columns
    assert "original_feature_id" in ft.columns
    assert "species" in ft.columns
    assert ft["species"].iloc[0] == "human"


def test_load_h5ad_preserves_names(sample_h5ad):
    ft = load_features(sample_h5ad, species="human")
    assert list(ft["original_feature_name"]) == ["TP53", "BRCA1", "CD3_ADT", "MYC", "BRCA2"]


def test_load_h5ad_extracts_ids(sample_h5ad):
    ft = load_features(sample_h5ad, species="human")
    assert ft["original_feature_id"].iloc[0] == "ENSG00000141510.18"


def test_load_h5ad_strips_version(sample_h5ad):
    ft = load_features(sample_h5ad, species="human")
    assert ft["feature_id_no_version"].iloc[0] == "ENSG00000141510"


def test_load_h5ad_dataset_name(sample_h5ad):
    ft = load_features(sample_h5ad, species="human", dataset_name="pbmc3k")
    assert ft["dataset"].iloc[0] == "pbmc3k"


def test_load_tsv_basic(sample_tsv):
    ft = load_features(sample_tsv, species="human")
    assert len(ft) == 4
    assert ft["original_feature_name"].iloc[0] == "TP53"


def test_load_tsv_with_column_map(sample_tsv):
    ft = load_features(
        sample_tsv,
        species="human",
        column_map={"gene_name": "original_feature_name", "gene_id": "original_feature_id"},
    )
    assert ft["original_feature_id"].iloc[0] == "ENSG00000141510.18"


def test_load_unsupported_format(tmp_path):
    path = str(tmp_path / "data.xyz")
    with open(path, "w") as f:
        f.write("junk")
    with pytest.raises(ValueError, match="Unsupported file format"):
        load_features(path, species="human")


def test_load_empty_ids_become_none(sample_h5ad):
    ft = load_features(sample_h5ad, species="human")
    # CD3_ADT has empty gene_id
    cd3_row = ft[ft["original_feature_name"] == "CD3_ADT"].iloc[0]
    assert pd.isna(cd3_row["original_feature_id"]) or cd3_row["original_feature_id"] == ""
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_io.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 4: Implement io.py**

```python
# src/stangene/io.py
"""Input/output adapters for loading feature metadata and writing results."""

import os
import re

import anndata
import pandas as pd

from stangene._logging import get_logger

logger = get_logger("io")

# Regex to strip Ensembl version suffix: ENSG00000141510.18 -> ENSG00000141510
_VERSION_SUFFIX = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+$")

# Common column names auto-detected for TSV/CSV files
_AUTO_COLUMN_MAP = {
    # source column name -> FeatureTable column name
    "gene": "original_feature_name",
    "gene_name": "original_feature_name",
    "feature_name": "original_feature_name",
    "gene_symbol": "original_feature_name",
    "symbol": "original_feature_name",
    "gene_id": "original_feature_id",
    "gene_ids": "original_feature_id",
    "ensembl_id": "original_feature_id",
    "ensembl_gene_id": "original_feature_id",
    "feature_id": "original_feature_id",
    "feature_types": "original_feature_type",
    "feature_type": "original_feature_type",
}


def _strip_version(eid: str) -> str:
    """Strip version suffix from an Ensembl ID. Returns empty string if no match."""
    if pd.isna(eid) or not eid:
        return ""
    m = _VERSION_SUFFIX.match(str(eid))
    return m.group(1) if m else ""


def _infer_reference_source(feature_ids: pd.Series) -> str:
    """Infer the reference source from ID patterns."""
    if feature_ids is None or feature_ids.isna().all():
        return ""
    sample = feature_ids.dropna().head(100)
    ensembl_count = sample.str.match(r"^ENS[A-Z]*G\d+").sum()
    if ensembl_count > len(sample) * 0.5:
        return "Ensembl/GENCODE"
    return ""


def load_features(
    path: str,
    species: str,
    dataset_name: str = None,
    column_map: dict = None,
) -> pd.DataFrame:
    """Load feature metadata from an h5ad or TSV/CSV file.

    Returns a standardized FeatureTable DataFrame. Does NOT load the expression matrix.
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in (".h5ad", ".h5"):
        ft = _load_h5ad(path)
    elif ext in (".tsv", ".csv", ".txt"):
        ft = _load_tabular(path, column_map)
    else:
        raise ValueError(
            f"Unsupported file format: '{ext}'. Supported: .h5ad, .tsv, .csv, .txt"
        )

    # Add metadata columns
    ft["species"] = species
    ft["dataset"] = dataset_name or os.path.splitext(os.path.basename(path))[0]

    # Derive version-stripped IDs
    if "original_feature_id" in ft.columns:
        ft["feature_id_no_version"] = ft["original_feature_id"].apply(_strip_version)
        ft.loc[ft["feature_id_no_version"] == "", "feature_id_no_version"] = None
    else:
        ft["original_feature_id"] = None
        ft["feature_id_no_version"] = None

    # Clean up empty strings to None
    ft["original_feature_id"] = ft["original_feature_id"].replace("", None)

    # Infer reference source
    ft["reference_source"] = _infer_reference_source(ft.get("original_feature_id"))
    ft["reference_release"] = None

    logger.info(
        "Loaded %d features from %s (species=%s, dataset=%s)",
        len(ft), path, species, ft["dataset"].iloc[0],
    )
    return ft


def _load_h5ad(path: str) -> pd.DataFrame:
    """Extract feature metadata from an h5ad file."""
    adata = anndata.read_h5ad(path, backed="r")
    var = adata.var.copy()

    ft = pd.DataFrame({"original_feature_name": var.index.tolist()})

    # Map known h5ad var columns
    if "gene_ids" in var.columns:
        ft["original_feature_id"] = var["gene_ids"].values
    if "feature_types" in var.columns:
        ft["original_feature_type"] = var["feature_types"].values

    if hasattr(adata, "file"):
        adata.file.close()

    return ft.reset_index(drop=True)


def _load_tabular(path: str, column_map: dict = None) -> pd.DataFrame:
    """Load feature metadata from a TSV/CSV file."""
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    raw = pd.read_csv(path, sep=sep)

    # Build column mapping
    if column_map is None:
        column_map = {}
        for col in raw.columns:
            col_lower = col.lower().strip()
            if col_lower in _AUTO_COLUMN_MAP:
                column_map[col] = _AUTO_COLUMN_MAP[col_lower]

    # Apply mapping
    ft = pd.DataFrame()
    for src_col, dst_col in column_map.items():
        if src_col in raw.columns:
            ft[dst_col] = raw[src_col].values

    # If no feature name was mapped, try using the first column or index
    if "original_feature_name" not in ft.columns:
        ft["original_feature_name"] = raw.iloc[:, 0].values

    return ft.reset_index(drop=True)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_io.py -v 2>&1 | tail -15`
Expected: All 10 tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/stangene/io.py tests/test_io.py
git commit -m "feat: add input adapters for h5ad and TSV with auto-detection"
```

---

### Task 5: Reference Build — Human (`references.py`, part 1)

**Files:**
- Create: `src/stangene/references.py`
- Create: `tests/test_references.py`

This task implements the reference building infrastructure and the human (HGNC) reference builder. Mouse is Task 6.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_references.py
import json
import os
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from stangene.references import (
    build_reference,
    load_reference,
    ReferenceNotFoundError,
    _build_human_reference,
    _default_reference_dir,
)


@pytest.fixture
def mock_hgnc_data():
    """Minimal HGNC-format TSV data."""
    return (
        "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tensembl_gene_id\tlocus_group\tstatus\n"
        'HGNC:11998\tTP53\tLFS1|p53\tLFS1\tENSG00000141510\tprotein-coding gene\tApproved\n'
        'HGNC:1100\tBRCA1\tRNF53|IRIS\tRNF53\tENSG00000012048\tprotein-coding gene\tApproved\n'
        'HGNC:1101\tBRCA2\tFACD|FANCD1\tFACD\tENSG00000139618\tprotein-coding gene\tApproved\n'
        'HGNC:999\tOLDGENE\t\t\t\tpseudogene\tEntry Withdrawn\n'
    )


@pytest.fixture
def ref_dir(tmp_path):
    return str(tmp_path / "references")


def test_reference_not_found_error(ref_dir):
    with pytest.raises(ReferenceNotFoundError):
        load_reference("human", reference_dir=ref_dir)


def test_build_human_creates_files(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    human_dir = os.path.join(ref_dir, "human")
    assert os.path.exists(os.path.join(human_dir, "gene_table.parquet"))
    assert os.path.exists(os.path.join(human_dir, "symbol_lookup.parquet"))
    assert os.path.exists(os.path.join(human_dir, "build_metadata.json"))


def test_build_human_gene_table(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    ref = load_reference("human", reference_dir=ref_dir)
    gt = ref["gene_table"]

    assert "ensembl_id" in gt.columns
    assert "symbol" in gt.columns
    assert "source_id" in gt.columns
    # TP53 should be present
    tp53 = gt[gt["symbol"] == "TP53"]
    assert len(tp53) == 1
    assert tp53.iloc[0]["ensembl_id"] == "ENSG00000141510"
    assert tp53.iloc[0]["source_id"] == "HGNC:11998"


def test_build_human_symbol_lookup(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    ref = load_reference("human", reference_dir=ref_dir)
    sl = ref["symbol_lookup"]

    # Check approved symbol lookup
    tp53_approved = sl[
        (sl["lookup_string"] == "TP53") & (sl["lookup_type"] == "approved_symbol")
    ]
    assert len(tp53_approved) == 1
    assert tp53_approved.iloc[0]["ensembl_id"] == "ENSG00000141510"

    # Check alias lookup
    lfs1_alias = sl[
        (sl["lookup_string"] == "LFS1") & (sl["lookup_type"] == "alias_symbol")
    ]
    assert len(lfs1_alias) == 1

    # Check uppercase column exists
    assert "lookup_string_upper" in sl.columns


def test_build_metadata(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    meta_path = os.path.join(ref_dir, "human", "build_metadata.json")
    with open(meta_path) as f:
        meta = json.load(f)
    assert "download_timestamp" in meta
    assert "sources" in meta


def test_build_skips_if_exists(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)
        build_reference("human", reference_dir=ref_dir)  # should skip
        assert mock_dl.call_count == 1  # only called once


def test_build_force_redownloads(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)
        build_reference("human", reference_dir=ref_dir, force=True)
        assert mock_dl.call_count >= 2
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_references.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 3: Implement references.py (infrastructure + human)**

```python
# src/stangene/references.py
"""Build and load species-specific gene annotation reference databases."""

import hashlib
import io
import json
import os
import urllib.request
from datetime import datetime, timezone

import pandas as pd

from stangene._logging import get_logger
from stangene.species import get_species_config

logger = get_logger("references")


class ReferenceNotFoundError(Exception):
    """Raised when reference data has not been built for a species."""

    pass


def _default_reference_dir() -> str:
    """Return the default reference directory (project-relative)."""
    return os.path.join(os.path.dirname(__file__), "..", "..", "references")


def _download_file(url: str) -> bytes:
    """Download a file from a URL and return its contents as bytes."""
    logger.info("Downloading %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "stangene/0.1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read()


def build_reference(
    species: str,
    reference_dir: str = None,
    force: bool = False,
) -> None:
    """Download and build reference tables for a species.

    Args:
        species: Species name (e.g., "human", "mouse").
        reference_dir: Directory to store references. Defaults to package-relative.
        force: If True, re-download and rebuild even if references exist.
    """
    config = get_species_config(species)
    ref_dir = os.path.join(reference_dir or _default_reference_dir(), config.name)

    # Check if already built
    gene_table_path = os.path.join(ref_dir, "gene_table.parquet")
    if os.path.exists(gene_table_path) and not force:
        logger.info("References for %s already exist at %s, skipping (use force=True to rebuild)", species, ref_dir)
        return

    os.makedirs(ref_dir, exist_ok=True)

    if config.name == "human":
        _build_human_reference(config, ref_dir)
    elif config.name == "mouse":
        _build_mouse_reference(config, ref_dir)
    else:
        raise ValueError(f"No reference builder for species: {species}")

    logger.info("Reference build complete for %s at %s", species, ref_dir)


def load_reference(
    species: str,
    reference_dir: str = None,
) -> dict:
    """Load built reference tables for a species.

    Returns:
        dict with keys: "gene_table" (DataFrame), "symbol_lookup" (DataFrame), "metadata" (dict).

    Raises:
        ReferenceNotFoundError: If references have not been built.
    """
    config = get_species_config(species)
    ref_dir = os.path.join(reference_dir or _default_reference_dir(), config.name)

    gene_table_path = os.path.join(ref_dir, "gene_table.parquet")
    lookup_path = os.path.join(ref_dir, "symbol_lookup.parquet")
    meta_path = os.path.join(ref_dir, "build_metadata.json")

    if not os.path.exists(gene_table_path):
        raise ReferenceNotFoundError(
            f"Reference data for '{species}' not found at {ref_dir}. "
            f"Run stangene.references.build_reference('{species}') first."
        )

    gene_table = pd.read_parquet(gene_table_path)
    symbol_lookup = pd.read_parquet(lookup_path)
    with open(meta_path) as f:
        metadata = json.load(f)

    logger.info("Loaded %s reference: %d genes, %d lookup entries", species, len(gene_table), len(symbol_lookup))
    return {
        "gene_table": gene_table,
        "symbol_lookup": symbol_lookup,
        "metadata": metadata,
    }


def _save_reference(ref_dir: str, gene_table: pd.DataFrame, symbol_lookup: pd.DataFrame, metadata: dict) -> None:
    """Save reference tables and metadata to disk."""
    gene_table.to_parquet(os.path.join(ref_dir, "gene_table.parquet"), index=False)
    symbol_lookup.to_parquet(os.path.join(ref_dir, "symbol_lookup.parquet"), index=False)
    with open(os.path.join(ref_dir, "build_metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)


def _build_symbol_lookup(gene_table: pd.DataFrame, source: str) -> pd.DataFrame:
    """Build the flattened symbol lookup index from a gene table."""
    rows = []

    for _, gene in gene_table.iterrows():
        eid = gene["ensembl_id"] if pd.notna(gene.get("ensembl_id")) else None
        sid = gene["source_id"]

        # Approved symbol
        if pd.notna(gene["symbol"]) and gene["symbol"]:
            rows.append({
                "lookup_string": gene["symbol"],
                "lookup_string_upper": gene["symbol"].upper(),
                "ensembl_id": eid,
                "source_id": sid,
                "lookup_type": "approved_symbol",
                "source": source,
            })

        # Alias symbols
        if pd.notna(gene.get("alias_symbols")) and gene["alias_symbols"]:
            for alias in str(gene["alias_symbols"]).split("|"):
                alias = alias.strip()
                if alias:
                    rows.append({
                        "lookup_string": alias,
                        "lookup_string_upper": alias.upper(),
                        "ensembl_id": eid,
                        "source_id": sid,
                        "lookup_type": "alias_symbol",
                        "source": source,
                    })

        # Previous symbols
        if pd.notna(gene.get("prev_symbols")) and gene["prev_symbols"]:
            for prev in str(gene["prev_symbols"]).split("|"):
                prev = prev.strip()
                if prev:
                    rows.append({
                        "lookup_string": prev,
                        "lookup_string_upper": prev.upper(),
                        "ensembl_id": eid,
                        "source_id": sid,
                        "lookup_type": "prev_symbol",
                        "source": source,
                    })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Human reference builder (HGNC)
# ---------------------------------------------------------------------------

def _build_human_reference(config, ref_dir: str) -> None:
    """Build human reference from HGNC complete gene set."""
    url = config.reference_sources["hgnc"]["url"]
    raw_data = _download_file(url)
    checksum = hashlib.sha256(raw_data).hexdigest()

    hgnc = pd.read_csv(io.BytesIO(raw_data), sep="\t", low_memory=False)

    # Build gene table
    gene_table = pd.DataFrame({
        "ensembl_id": hgnc["ensembl_gene_id"].where(hgnc["ensembl_gene_id"].notna(), None),
        "symbol": hgnc["symbol"],
        "alias_symbols": hgnc["alias_symbol"].fillna(""),
        "prev_symbols": hgnc["prev_symbol"].fillna(""),
        "gene_type": hgnc.get("locus_group", pd.Series(dtype=str)).fillna(""),
        "status": hgnc["status"].fillna(""),
        "source": "HGNC",
        "source_id": hgnc["hgnc_id"],
    })

    # Build symbol lookup
    symbol_lookup = _build_symbol_lookup(gene_table, source="HGNC")

    # Metadata
    metadata = {
        "species": "human",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "hgnc": {
                "url": url,
                "sha256": checksum,
                "rows": len(hgnc),
            }
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built human reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


# ---------------------------------------------------------------------------
# Mouse reference builder (MGI) — placeholder, implemented in Task 6
# ---------------------------------------------------------------------------

def _build_mouse_reference(config, ref_dir: str) -> None:
    """Build mouse reference from MGI marker files + Ensembl BioMart."""
    raise NotImplementedError("Mouse reference builder will be implemented in Task 6")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_references.py -v 2>&1 | tail -15`
Expected: All 7 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/references.py tests/test_references.py
git commit -m "feat: add reference build/load infrastructure with human HGNC builder"
```

---

### Task 6: Reference Build — Mouse (MGI)

**Files:**
- Modify: `src/stangene/references.py` (implement `_build_mouse_reference`)
- Modify: `tests/test_references.py` (add mouse tests)

- [ ] **Step 1: Write failing mouse tests**

Append to `tests/test_references.py`:

```python
@pytest.fixture
def mock_mgi_markers_data():
    """Minimal MGI MRK_List2.rpt format."""
    return (
        "MGI Accession ID\tChr\tcM Position\tgenome coordinate start\tgenome coordinate end\tstrand\tMarker Symbol\tStatus\tMarker Name\tMarker Type\tFeature Type\tMarker Synonyms (pipe-separated)\n"
        "MGI:87853\t11\t\t69580359\t69591872\t+\tTrp53\tO\ttransformation related protein 53\tGene\tprotein coding gene\tp53|Tp53\n"
        "MGI:104738\t11\t\t101453964\t101517817\t+\tBrca1\tO\tbreast cancer 1, early onset\tGene\tprotein coding gene\tBrca1\n"
        "MGI:12345\t1\t\t1000\t2000\t+\tFakeGene\tO\tfake gene\tGene\tprotein coding gene\t\n"
    )


@pytest.fixture
def mock_mgi_ensembl_data():
    """Minimal MGI MRK_ENSEMBL.rpt format."""
    return (
        "MGI Marker Accession ID\tMarker Symbol\tMarker Name\tcM Position\tChromosome\tEnsembl Gene ID\n"
        "MGI:87853\tTrp53\ttransformation related protein 53\t\t11\tENSMUSG00000059552\n"
        "MGI:104738\tBrca1\tbreast cancer 1, early onset\t\t11\tENSMUSG00000017146\n"
    )


def test_build_mouse_creates_files(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    mouse_dir = os.path.join(ref_dir, "mouse")
    assert os.path.exists(os.path.join(mouse_dir, "gene_table.parquet"))
    assert os.path.exists(os.path.join(mouse_dir, "symbol_lookup.parquet"))


def test_build_mouse_gene_table(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    ref = load_reference("mouse", reference_dir=ref_dir)
    gt = ref["gene_table"]

    # Trp53 should have Ensembl ID from MGI-Ensembl mapping
    trp53 = gt[gt["symbol"] == "Trp53"]
    assert len(trp53) == 1
    assert trp53.iloc[0]["ensembl_id"] == "ENSMUSG00000059552"
    assert trp53.iloc[0]["source_id"] == "MGI:87853"

    # FakeGene has no Ensembl ID — ensembl_id should be null
    fake = gt[gt["symbol"] == "FakeGene"]
    assert len(fake) == 1
    assert pd.isna(fake.iloc[0]["ensembl_id"])
    assert fake.iloc[0]["source_id"] == "MGI:12345"


def test_build_mouse_symbol_lookup(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    ref = load_reference("mouse", reference_dir=ref_dir)
    sl = ref["symbol_lookup"]

    # Approved symbol
    trp53_approved = sl[(sl["lookup_string"] == "Trp53") & (sl["lookup_type"] == "approved_symbol")]
    assert len(trp53_approved) == 1

    # Synonym (alias)
    p53_alias = sl[(sl["lookup_string"] == "p53") & (sl["lookup_type"] == "alias_symbol")]
    assert len(p53_alias) == 1
```

- [ ] **Step 2: Run new tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_references.py::test_build_mouse_creates_files -v 2>&1 | tail -5`
Expected: FAIL — `NotImplementedError`

- [ ] **Step 3: Implement _build_mouse_reference**

Replace the placeholder in `src/stangene/references.py`:

```python
def _build_mouse_reference(config, ref_dir: str) -> None:
    """Build mouse reference from MGI marker files + Ensembl mapping."""
    # Download MGI files
    markers_url = config.reference_sources["mgi_markers"]["url"]
    ensembl_url = config.reference_sources["mgi_ensembl"]["url"]

    markers_raw = _download_file(markers_url)
    ensembl_raw = _download_file(ensembl_url)

    markers_checksum = hashlib.sha256(markers_raw).hexdigest()
    ensembl_checksum = hashlib.sha256(ensembl_raw).hexdigest()

    # Parse MGI markers
    markers = pd.read_csv(io.BytesIO(markers_raw), sep="\t", low_memory=False)
    markers.columns = markers.columns.str.strip()

    # Parse MGI-Ensembl mapping
    ensembl_map = pd.read_csv(io.BytesIO(ensembl_raw), sep="\t", low_memory=False)
    ensembl_map.columns = ensembl_map.columns.str.strip()

    # Build MGI ID -> Ensembl ID lookup
    mgi_to_ensembl = {}
    eid_col = [c for c in ensembl_map.columns if "ensembl" in c.lower() and "id" in c.lower()]
    mid_col = [c for c in ensembl_map.columns if "mgi" in c.lower() and "accession" in c.lower()]
    if eid_col and mid_col:
        for _, row in ensembl_map.iterrows():
            mgi_id = row[mid_col[0]]
            ens_id = row[eid_col[0]]
            if pd.notna(mgi_id) and pd.notna(ens_id) and str(ens_id).startswith("ENSMUSG"):
                mgi_to_ensembl[str(mgi_id)] = str(ens_id)

    # Build gene table from markers
    sym_col = [c for c in markers.columns if "marker symbol" in c.lower()][0]
    status_col = [c for c in markers.columns if "status" in c.lower()][0]
    type_col = [c for c in markers.columns if "feature type" in c.lower()][0]
    mgi_col = [c for c in markers.columns if "mgi accession" in c.lower()][0]
    syn_col = [c for c in markers.columns if "synonym" in c.lower()][0]

    rows = []
    for _, m in markers.iterrows():
        mgi_id = str(m[mgi_col]) if pd.notna(m[mgi_col]) else ""
        symbol = str(m[sym_col]) if pd.notna(m[sym_col]) else ""
        ensembl_id = mgi_to_ensembl.get(mgi_id)
        synonyms = str(m[syn_col]) if pd.notna(m[syn_col]) else ""

        # MGI status: O = official, W = withdrawn
        status = "approved" if str(m[status_col]).strip() == "O" else "withdrawn"
        gene_type = str(m[type_col]) if pd.notna(m[type_col]) else ""

        rows.append({
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "alias_symbols": synonyms,
            "prev_symbols": "",
            "gene_type": gene_type,
            "status": status,
            "source": "MGI",
            "source_id": mgi_id,
        })

    gene_table = pd.DataFrame(rows)

    # Supplementary: fill Ensembl ID gaps using BioMart symbol-to-Ensembl mapping.
    # Build a symbol->ensembl_id lookup from BioMart for genes missing ensembl_id.
    biomart_url = config.reference_sources.get("ensembl_biomart", {}).get("url", "")
    if biomart_url:
        try:
            biomart_xml = (
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<!DOCTYPE Query>'
                '<Query virtualSchemaName="default" formatter="TSV" header="1">'
                '<Dataset name="mmusculus_gene_ensembl" interface="default">'
                '<Attribute name="ensembl_gene_id"/>'
                '<Attribute name="external_gene_name"/>'
                '<Attribute name="mgi_id"/>'
                '</Dataset></Query>'
            )
            import urllib.parse
            full_url = biomart_url + urllib.parse.quote(biomart_xml)
            biomart_raw = _download_file(full_url)
            biomart_df = pd.read_csv(io.BytesIO(biomart_raw), sep="\t", low_memory=False)
            biomart_df.columns = ["ensembl_id_bm", "symbol_bm", "mgi_id_bm"]

            # Fill gaps: for genes with no ensembl_id, try matching by MGI ID or symbol
            null_mask = gene_table["ensembl_id"].isna()
            if null_mask.any():
                # Match by MGI ID first
                bm_by_mgi = biomart_df.dropna(subset=["mgi_id_bm"]).drop_duplicates(subset=["mgi_id_bm"])
                bm_mgi_map = dict(zip(bm_by_mgi["mgi_id_bm"], bm_by_mgi["ensembl_id_bm"]))
                for idx in gene_table[null_mask].index:
                    sid = gene_table.at[idx, "source_id"]
                    if sid in bm_mgi_map:
                        gene_table.at[idx, "ensembl_id"] = bm_mgi_map[sid]

                # Match remaining gaps by symbol
                still_null = gene_table["ensembl_id"].isna()
                if still_null.any():
                    bm_by_sym = biomart_df.dropna(subset=["symbol_bm"]).drop_duplicates(subset=["symbol_bm"])
                    bm_sym_map = dict(zip(bm_by_sym["symbol_bm"], bm_by_sym["ensembl_id_bm"]))
                    for idx in gene_table[still_null].index:
                        sym = gene_table.at[idx, "symbol"]
                        if sym in bm_sym_map:
                            gene_table.at[idx, "ensembl_id"] = bm_sym_map[sym]

                filled = null_mask.sum() - gene_table["ensembl_id"].isna().sum()
                logger.info("BioMart supplementary fill: %d/%d Ensembl ID gaps filled", filled, null_mask.sum())
        except Exception as e:
            logger.warning("BioMart supplementary download failed (non-fatal): %s", e)

    # Build symbol lookup
    symbol_lookup = _build_symbol_lookup(gene_table, source="MGI")

    # Metadata
    metadata = {
        "species": "mouse",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "mgi_markers": {"url": markers_url, "sha256": markers_checksum, "rows": len(markers)},
            "mgi_ensembl": {"url": ensembl_url, "sha256": ensembl_checksum, "rows": len(ensembl_map)},
            "ensembl_biomart": {"url": biomart_url, "note": "supplementary Ensembl ID gap fill"},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built mouse reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))
```

- [ ] **Step 4: Run all reference tests**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_references.py -v 2>&1 | tail -20`
Expected: All 10 tests PASS (7 existing + 3 new mouse tests)

- [ ] **Step 5: Commit**

```bash
git add src/stangene/references.py tests/test_references.py
git commit -m "feat: add mouse reference builder from MGI markers + Ensembl mapping"
```

---

### Task 7: Harmonization Cascade (`harmonize.py`)

**Files:**
- Create: `src/stangene/harmonize.py`
- Create: `tests/test_harmonize.py`
- Create: `tests/conftest.py` (shared fixtures)

- [ ] **Step 1: Create shared test fixtures**

```python
# tests/conftest.py
"""Shared pytest fixtures for stangene tests."""

import pandas as pd
import pytest


@pytest.fixture
def mock_gene_table():
    """Minimal gene_table DataFrame mimicking HGNC reference."""
    return pd.DataFrame([
        {"ensembl_id": "ENSG00000141510", "symbol": "TP53", "alias_symbols": "LFS1|p53",
         "prev_symbols": "LFS1", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:11998"},
        {"ensembl_id": "ENSG00000012048", "symbol": "BRCA1", "alias_symbols": "RNF53|IRIS",
         "prev_symbols": "RNF53", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:1100"},
        {"ensembl_id": "ENSG00000139618", "symbol": "BRCA2", "alias_symbols": "FACD|FANCD1",
         "prev_symbols": "FACD", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:1101"},
        {"ensembl_id": "ENSG00000136997", "symbol": "MYC", "alias_symbols": "",
         "prev_symbols": "", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:7553"},
        {"ensembl_id": None, "symbol": "WITHDRAWN1", "alias_symbols": "",
         "prev_symbols": "", "gene_type": "pseudogene", "status": "Entry Withdrawn",
         "source": "HGNC", "source_id": "HGNC:99999"},
    ])


@pytest.fixture
def mock_symbol_lookup():
    """Minimal symbol_lookup DataFrame."""
    rows = [
        {"lookup_string": "TP53", "lookup_string_upper": "TP53", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "LFS1", "lookup_string_upper": "LFS1", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "p53", "lookup_string_upper": "P53", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "LFS1", "lookup_string_upper": "LFS1", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "prev_symbol", "source": "HGNC"},
        {"lookup_string": "BRCA1", "lookup_string_upper": "BRCA1", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "RNF53", "lookup_string_upper": "RNF53", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "RNF53", "lookup_string_upper": "RNF53", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "prev_symbol", "source": "HGNC"},
        {"lookup_string": "BRCA2", "lookup_string_upper": "BRCA2", "ensembl_id": "ENSG00000139618",
         "source_id": "HGNC:1101", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "FACD", "lookup_string_upper": "FACD", "ensembl_id": "ENSG00000139618",
         "source_id": "HGNC:1101", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "MYC", "lookup_string_upper": "MYC", "ensembl_id": "ENSG00000136997",
         "source_id": "HGNC:7553", "lookup_type": "approved_symbol", "source": "HGNC"},
        # Ambiguous alias: AMBIG maps to two different genes
        {"lookup_string": "AMBIG", "lookup_string_upper": "AMBIG", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "AMBIG", "lookup_string_upper": "AMBIG", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "alias_symbol", "source": "HGNC"},
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def mock_ref(mock_gene_table, mock_symbol_lookup):
    """Complete mock reference dict."""
    return {
        "gene_table": mock_gene_table,
        "symbol_lookup": mock_symbol_lookup,
        "metadata": {"species": "human", "download_timestamp": "2026-01-01T00:00:00Z"},
    }
```

- [ ] **Step 2: Write failing tests**

```python
# tests/test_harmonize.py
import pandas as pd
import pytest

from stangene.harmonize import harmonize, HarmonizationResult


def _make_ft(names, ids=None, types=None):
    """Helper to build a classified FeatureTable."""
    n = len(names)
    data = {
        "original_feature_name": names,
        "species": ["human"] * n,
        "dataset": ["test"] * n,
        "original_feature_type": types or (["gene"] * n),
        "mapping_status": [None] * n,
        "mapping_notes": [None] * n,
    }
    if ids is not None:
        data["original_feature_id"] = ids
        data["feature_id_no_version"] = [
            i.split(".")[0] if i and "." in i else (None if not i else "")
            for i in ids
        ]
    else:
        data["original_feature_id"] = [None] * n
        data["feature_id_no_version"] = [None] * n
    return pd.DataFrame(data)


class TestTier1ExactId:
    def test_exact_id_match(self, mock_ref):
        ft = _make_ft(["TP53"], ids=["ENSG00000141510"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_id"
        assert row["gene_id_harmonized"] == "ENSG00000141510"
        assert row["gene_symbol_harmonized"] == "TP53"
        assert row["mapping_confidence"] == "high"


class TestTier2VersionStripped:
    def test_version_stripped_match(self, mock_ref):
        ft = _make_ft(["TP53"], ids=["ENSG00000141510.18"])
        ft["feature_id_no_version"] = ["ENSG00000141510"]
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "id_no_version"
        assert row["gene_id_harmonized"] == "ENSG00000141510"
        assert row["mapping_confidence"] == "high"


class TestTier3ExactSymbol:
    def test_exact_symbol_match(self, mock_ref):
        ft = _make_ft(["MYC"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000136997"
        assert row["mapping_confidence"] == "high"


class TestTier4AliasAndPrevSymbol:
    def test_alias_match(self, mock_ref):
        ft = _make_ft(["p53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "alias_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000141510"
        assert row["mapping_confidence"] == "medium"

    def test_prev_symbol_match(self, mock_ref):
        ft = _make_ft(["RNF53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        # RNF53 exists as both alias and prev_symbol for BRCA1
        assert row["mapping_status"] in ("alias_symbol", "previous_symbol")
        assert row["gene_id_harmonized"] == "ENSG00000012048"

    def test_ambiguous_alias(self, mock_ref):
        ft = _make_ft(["AMBIG"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "ambiguous"
        assert row["mapping_confidence"] == "low"


class TestTier5Unmapped:
    def test_unmapped(self, mock_ref):
        ft = _make_ft(["NONEXISTENT_GENE"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "unmapped"
        assert pd.isna(row["gene_id_harmonized"]) or row["gene_id_harmonized"] is None

    def test_excel_date_unmapped(self, mock_ref):
        ft = _make_ft(["1-Mar"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "unmapped"
        assert "Excel" in str(row["mapping_notes"])

    def test_excel_date_format_sep(self, mock_ref):
        ft = _make_ft(["2-Sep"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "unmapped"
        assert "Excel" in str(row["mapping_notes"])


class TestWithdrawnGene:
    def test_withdrawn_gene_medium_confidence(self, mock_ref):
        """Matching a withdrawn gene should get mapping_confidence=medium."""
        # Add WITHDRAWN1 as an approved_symbol in the lookup so Tier 3 finds it
        import pandas as _pd
        extra_row = _pd.DataFrame([{
            "lookup_string": "WITHDRAWN1", "lookup_string_upper": "WITHDRAWN1",
            "ensembl_id": None, "source_id": "HGNC:99999",
            "lookup_type": "approved_symbol", "source": "HGNC",
        }])
        mock_ref["symbol_lookup"] = _pd.concat(
            [mock_ref["symbol_lookup"], extra_row], ignore_index=True,
        )
        ft = _make_ft(["WITHDRAWN1"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["mapping_confidence"] == "medium"
        assert "withdrawn" in str(row["mapping_notes"]).lower()


class TestNonGenePassthrough:
    def test_non_gene_skipped(self, mock_ref):
        ft = _make_ft(
            ["ERCC-00002"],
            types=["spike_in"],
        )
        ft["mapping_status"] = "non_gene_feature"
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "non_gene_feature"


class TestEarlyExit:
    def test_id_match_skips_symbol(self, mock_ref):
        """If ID matches at Tier 1, don't re-check symbol."""
        ft = _make_ft(["WRONG_NAME"], ids=["ENSG00000141510"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_id"
        assert row["gene_symbol_harmonized"] == "TP53"


class TestManyToOneDetection:
    def test_many_to_one_in_conflicts(self, mock_ref):
        ft = _make_ft(["TP53", "p53"], ids=["ENSG00000141510", None])
        ft["feature_id_no_version"] = [None, None]
        result = harmonize(ft, mock_ref)
        # Both map to ENSG00000141510 — should appear in conflicts
        assert len(result.conflicts) > 0


class TestResultStructure:
    def test_result_has_stats(self, mock_ref):
        ft = _make_ft(["TP53", "NONEXISTENT"], ids=["ENSG00000141510", None])
        result = harmonize(ft, mock_ref)
        assert isinstance(result, HarmonizationResult)
        assert "exact_id" in result.stats
        assert "unmapped" in result.stats

    def test_stangene_version_in_output(self, mock_ref):
        ft = _make_ft(["TP53"])
        result = harmonize(ft, mock_ref)
        assert "stangene_version" in result.mapping_table.columns
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_harmonize.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 4: Implement harmonize.py**

```python
# src/stangene/harmonize.py
"""Tiered harmonization cascade for gene identifier mapping."""

from dataclasses import dataclass

import pandas as pd

from stangene import __version__
from stangene._logging import get_logger
from stangene.species import EXCEL_DATE_PATTERN, EXCEL_RENAMED_GENES

logger = get_logger("harmonize")


@dataclass
class HarmonizationResult:
    """Result of the harmonization cascade."""

    mapping_table: pd.DataFrame
    conflicts: pd.DataFrame
    stats: dict


def harmonize(ft: pd.DataFrame, ref: dict) -> HarmonizationResult:
    """Run the tiered matching cascade on a classified FeatureTable.

    Args:
        ft: FeatureTable with classification columns (from classify_features).
        ref: Reference dict from load_reference (gene_table, symbol_lookup, metadata).

    Returns:
        HarmonizationResult with mapping_table, conflicts, and stats.
    """
    result = ft.copy()
    gene_table = ref["gene_table"]
    symbol_lookup = ref["symbol_lookup"]

    # Ensure output columns exist
    for col in [
        "gene_id_harmonized", "gene_symbol_harmonized",
        "mapping_status", "mapping_confidence", "mapping_source", "mapping_notes",
    ]:
        if col not in result.columns:
            result[col] = None

    # Build fast lookup sets
    ensembl_id_set = set(gene_table["ensembl_id"].dropna())
    ensembl_id_to_row = gene_table.dropna(subset=["ensembl_id"]).set_index("ensembl_id")

    # Pre-index symbol lookups by type for efficiency
    approved_lookup = symbol_lookup[symbol_lookup["lookup_type"] == "approved_symbol"]
    alias_prev_lookup = symbol_lookup[symbol_lookup["lookup_type"].isin(["alias_symbol", "prev_symbol"])]

    for idx in result.index:
        # Skip non-gene features (already have mapping_status set)
        if result.at[idx, "mapping_status"] == "non_gene_feature":
            continue

        feature_name = str(result.at[idx, "original_feature_name"]) if pd.notna(result.at[idx, "original_feature_name"]) else ""
        feature_id = result.at[idx, "original_feature_id"] if pd.notna(result.at[idx, "original_feature_id"]) else None
        feature_id_nv = result.at[idx, "feature_id_no_version"] if pd.notna(result.at[idx, "feature_id_no_version"]) else None

        notes = []

        # Check for Excel date corruption
        if EXCEL_DATE_PATTERN.match(feature_name):
            result.at[idx, "mapping_status"] = "unmapped"
            result.at[idx, "mapping_confidence"] = None
            result.at[idx, "mapping_notes"] = f"Likely Excel-corrupted date format: {feature_name}"
            continue

        # Note if this is a known Excel-renamed gene
        if feature_name.upper() in EXCEL_RENAMED_GENES:
            notes.append(f"Known Excel-renamed gene: {feature_name} -> {EXCEL_RENAMED_GENES[feature_name.upper()]}")

        # --- Tier 1: Exact stable ID match ---
        if feature_id and feature_id in ensembl_id_set:
            gene_row = ensembl_id_to_row.loc[feature_id]
            if isinstance(gene_row, pd.DataFrame):
                gene_row = gene_row.iloc[0]
            _apply_match(result, idx, gene_row, "exact_id", "high", f"HGNC:exact_id", notes)
            continue

        # --- Tier 2: Version-stripped ID match ---
        if feature_id_nv and feature_id_nv in ensembl_id_set:
            gene_row = ensembl_id_to_row.loc[feature_id_nv]
            if isinstance(gene_row, pd.DataFrame):
                gene_row = gene_row.iloc[0]
            notes.append(f"Matched via version-stripped ID: {feature_id} -> {feature_id_nv}")
            _apply_match(result, idx, gene_row, "id_no_version", "high", f"HGNC:id_no_version", notes)
            continue

        # --- Tier 3: Exact approved symbol match ---
        matches = approved_lookup[approved_lookup["lookup_string"] == feature_name]
        if len(matches) == 1:
            match = matches.iloc[0]
            eid = match["ensembl_id"]
            confidence = "high"
            # Check if matched gene is withdrawn or non-protein-coding
            gene_info = None
            if eid and eid in ensembl_id_set:
                gene_info = ensembl_id_to_row.loc[eid]
                if isinstance(gene_info, pd.DataFrame):
                    gene_info = gene_info.iloc[0]
            elif pd.notna(match.get("source_id")):
                sid_rows = gene_table[gene_table["source_id"] == match["source_id"]]
                if len(sid_rows) > 0:
                    gene_info = sid_rows.iloc[0]
            if gene_info is not None:
                if str(gene_info.get("status", "")).lower().startswith("entry withdrawn") or str(gene_info.get("status", "")).lower() == "withdrawn":
                    confidence = "medium"
                    notes.append("Matched withdrawn gene")
                gene_type = gene_info.get("gene_type", "")
                if gene_type and "protein" not in str(gene_type).lower():
                    notes.append(f"Non-protein-coding gene type: {gene_type}")
            _apply_match_from_lookup(result, idx, match, gene_table, "exact_symbol", confidence, notes)
            continue
        elif len(matches) > 1:
            candidates = matches["ensembl_id"].tolist()
            notes.append(f"Multiple approved symbol matches: {candidates}")
            result.at[idx, "mapping_status"] = "ambiguous"
            result.at[idx, "mapping_confidence"] = "low"
            result.at[idx, "mapping_notes"] = "; ".join(notes)
            continue

        # --- Tier 4: Alias / previous symbol match ---
        matches = alias_prev_lookup[alias_prev_lookup["lookup_string"] == feature_name]
        # Deduplicate by ensembl_id + source_id (same gene can appear as both alias and prev_symbol)
        unique_targets = matches.drop_duplicates(subset=["ensembl_id", "source_id"])
        if len(unique_targets) == 1:
            match = matches.iloc[0]
            lookup_type = match["lookup_type"]
            status = "alias_symbol" if lookup_type == "alias_symbol" else "previous_symbol"
            _apply_match_from_lookup(result, idx, match, gene_table, status, "medium", notes)
            continue
        elif len(unique_targets) > 1:
            candidates = unique_targets[["ensembl_id", "source_id", "lookup_type"]].to_dict("records")
            notes.append(f"Multiple alias/prev matches: {candidates}")
            result.at[idx, "mapping_status"] = "ambiguous"
            result.at[idx, "mapping_confidence"] = "low"
            result.at[idx, "mapping_notes"] = "; ".join(notes)
            continue

        # --- Tier 5: Unmapped ---
        result.at[idx, "mapping_status"] = "unmapped"
        result.at[idx, "mapping_confidence"] = None
        if notes:
            result.at[idx, "mapping_notes"] = "; ".join(notes)

    # Add stangene version and reference release
    result["stangene_version"] = __version__
    result["reference_release"] = ref.get("metadata", {}).get("download_timestamp", "")

    # Build stats
    stats = result["mapping_status"].value_counts().to_dict()

    # Detect conflicts: many-to-one (multiple original features -> same harmonized ID)
    harmonized_ids = result[result["gene_id_harmonized"].notna()]
    id_counts = harmonized_ids["gene_id_harmonized"].value_counts()
    duplicate_ids = id_counts[id_counts > 1].index.tolist()
    conflicts = result[result["gene_id_harmonized"].isin(duplicate_ids)].copy() if duplicate_ids else pd.DataFrame()

    logger.info("Harmonization complete: %s", stats)
    if len(conflicts) > 0:
        logger.info("Found %d features involved in %d many-to-one conflicts", len(conflicts), len(duplicate_ids))

    return HarmonizationResult(
        mapping_table=result,
        conflicts=conflicts,
        stats=stats,
    )


def _apply_match(result, idx, gene_row, status, confidence, source, notes):
    """Apply a match from gene_table directly."""
    result.at[idx, "gene_id_harmonized"] = gene_row.get("ensembl_id") or gene_row.get("source_id", "")
    result.at[idx, "gene_symbol_harmonized"] = gene_row.get("symbol", "")
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)


def _apply_match_from_lookup(result, idx, match, gene_table, status, confidence, notes):
    """Apply a match from symbol_lookup, looking up the full gene info."""
    eid = match.get("ensembl_id")
    sid = match.get("source_id")
    source_label = f"{match.get('source', '')}:{match.get('lookup_type', '')}"

    # Find the gene in gene_table for symbol
    symbol = ""
    if pd.notna(eid):
        gene_rows = gene_table[gene_table["ensembl_id"] == eid]
        if len(gene_rows) > 0:
            symbol = gene_rows.iloc[0]["symbol"]
    if not symbol and pd.notna(sid):
        gene_rows = gene_table[gene_table["source_id"] == sid]
        if len(gene_rows) > 0:
            symbol = gene_rows.iloc[0]["symbol"]

    harmonized_id = eid if pd.notna(eid) else sid
    result.at[idx, "gene_id_harmonized"] = harmonized_id
    result.at[idx, "gene_symbol_harmonized"] = symbol
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source_label
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_harmonize.py -v 2>&1 | tail -20`
Expected: All 12 tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/stangene/harmonize.py tests/test_harmonize.py tests/conftest.py
git commit -m "feat: add tiered harmonization cascade with conflict detection"
```

---

### Task 8: Reporting (`report.py`)

**Files:**
- Create: `src/stangene/report.py`
- Create: `tests/test_report.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_report.py
import json
import os

import pandas as pd
import pytest

from stangene.harmonize import HarmonizationResult
from stangene.report import summary, conflict_report, write_reports


@pytest.fixture
def sample_result():
    """A HarmonizationResult with diverse mapping statuses."""
    mapping_table = pd.DataFrame([
        {"original_feature_name": "TP53", "gene_id_harmonized": "ENSG00000141510",
         "gene_symbol_harmonized": "TP53", "mapping_status": "exact_id",
         "mapping_confidence": "high", "original_feature_type": "gene"},
        {"original_feature_name": "p53", "gene_id_harmonized": "ENSG00000141510",
         "gene_symbol_harmonized": "TP53", "mapping_status": "alias_symbol",
         "mapping_confidence": "medium", "original_feature_type": "gene"},
        {"original_feature_name": "MYC", "gene_id_harmonized": "ENSG00000136997",
         "gene_symbol_harmonized": "MYC", "mapping_status": "exact_symbol",
         "mapping_confidence": "high", "original_feature_type": "gene"},
        {"original_feature_name": "UNKNOWN", "gene_id_harmonized": None,
         "gene_symbol_harmonized": None, "mapping_status": "unmapped",
         "mapping_confidence": None, "original_feature_type": "gene"},
        {"original_feature_name": "ERCC-00002", "gene_id_harmonized": None,
         "gene_symbol_harmonized": None, "mapping_status": "non_gene_feature",
         "mapping_confidence": None, "original_feature_type": "spike_in"},
    ])
    conflicts = mapping_table[mapping_table["gene_id_harmonized"] == "ENSG00000141510"].copy()
    stats = mapping_table["mapping_status"].value_counts().to_dict()
    return HarmonizationResult(mapping_table=mapping_table, conflicts=conflicts, stats=stats)


def test_summary_keys(sample_result):
    s = summary(sample_result)
    assert "total_features" in s
    assert "gene_features" in s
    assert "non_gene_features" in s
    assert "status_counts" in s
    assert s["total_features"] == 5
    assert s["gene_features"] == 4
    assert s["non_gene_features"] == 1


def test_summary_duplicate_counts(sample_result):
    s = summary(sample_result)
    assert s["duplicate_harmonized_ids"] == 1  # ENSG00000141510 appears twice
    assert s["duplicate_harmonized_symbols"] == 1  # TP53 appears twice


def test_conflict_report_structure(sample_result):
    cr = conflict_report(sample_result)
    assert isinstance(cr, pd.DataFrame)
    assert "conflict_type" in cr.columns
    assert len(cr) > 0


def test_conflict_report_detects_many_to_one(sample_result):
    cr = conflict_report(sample_result)
    m2o = cr[cr["conflict_type"] == "many_to_one"]
    assert len(m2o) > 0


def test_conflict_report_includes_unmapped(sample_result):
    cr = conflict_report(sample_result)
    unmapped = cr[cr["conflict_type"] == "unmapped"]
    assert len(unmapped) == 1


def test_write_reports_creates_files(sample_result, tmp_path):
    output_dir = str(tmp_path / "output")
    write_reports(sample_result, output_dir)
    assert os.path.exists(os.path.join(output_dir, "harmonization_table.tsv"))
    assert os.path.exists(os.path.join(output_dir, "summary.json"))
    assert os.path.exists(os.path.join(output_dir, "conflicts.tsv"))
    assert os.path.exists(os.path.join(output_dir, "unmapped.tsv"))


def test_write_reports_summary_json(sample_result, tmp_path):
    output_dir = str(tmp_path / "output")
    write_reports(sample_result, output_dir)
    with open(os.path.join(output_dir, "summary.json")) as f:
        s = json.load(f)
    assert s["total_features"] == 5
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_report.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 3: Implement report.py**

```python
# src/stangene/report.py
"""Reporting: summaries, conflict detection, and output writing."""

import json
import os

import pandas as pd

from stangene._logging import get_logger
from stangene.species import EXCEL_DATE_PATTERN, EXCEL_RENAMED_GENES

logger = get_logger("report")


def summary(result) -> dict:
    """Generate dataset-level summary statistics.

    Args:
        result: HarmonizationResult from harmonize().

    Returns:
        dict with counts per category.
    """
    mt = result.mapping_table
    status_counts = mt["mapping_status"].value_counts().to_dict()

    gene_mask = mt["original_feature_type"].isin(["gene", "Gene Expression"])
    non_gene_mask = ~gene_mask

    # Count duplicates in harmonized IDs/symbols
    harmonized_ids = mt["gene_id_harmonized"].dropna()
    id_dupes = harmonized_ids.duplicated(keep=False).sum()
    harmonized_syms = mt["gene_symbol_harmonized"].dropna()
    sym_dupes = harmonized_syms.duplicated(keep=False).sum()

    # Count unique duplicate groups
    id_dupe_groups = harmonized_ids[harmonized_ids.duplicated(keep=False)].nunique()
    sym_dupe_groups = harmonized_syms[harmonized_syms.duplicated(keep=False)].nunique()

    return {
        "total_features": len(mt),
        "gene_features": int(gene_mask.sum()),
        "non_gene_features": int(non_gene_mask.sum()),
        "status_counts": status_counts,
        "duplicate_harmonized_ids": int(id_dupe_groups),
        "duplicate_harmonized_symbols": int(sym_dupe_groups),
    }


def conflict_report(result) -> pd.DataFrame:
    """Generate a conflict report table.

    Includes: many-to-one collisions, unmapped features, suspicious symbols.
    """
    mt = result.mapping_table
    rows = []

    # Many-to-one: multiple originals -> same harmonized ID
    harmonized_ids = mt[mt["gene_id_harmonized"].notna()]
    id_counts = harmonized_ids["gene_id_harmonized"].value_counts()
    for hid in id_counts[id_counts > 1].index:
        involved = mt[mt["gene_id_harmonized"] == hid]
        for _, row in involved.iterrows():
            rows.append({
                "conflict_type": "many_to_one",
                "original_feature_name": row["original_feature_name"],
                "gene_id_harmonized": hid,
                "gene_symbol_harmonized": row.get("gene_symbol_harmonized", ""),
                "mapping_status": row["mapping_status"],
                "details": f"{len(involved)} features map to {hid}",
            })

    # Unmapped features
    unmapped = mt[mt["mapping_status"] == "unmapped"]
    for _, row in unmapped.iterrows():
        details = ""
        name = str(row["original_feature_name"])
        if EXCEL_DATE_PATTERN.match(name):
            details = "Likely Excel date artifact"
        elif name.upper() in EXCEL_RENAMED_GENES:
            details = f"Known Excel-renamed gene (old: {name}, new: {EXCEL_RENAMED_GENES[name.upper()]})"
        rows.append({
            "conflict_type": "unmapped",
            "original_feature_name": name,
            "gene_id_harmonized": None,
            "gene_symbol_harmonized": None,
            "mapping_status": "unmapped",
            "details": details,
        })

    # Suspicious: previous_symbol matches (likely outdated names)
    prev_sym = mt[mt["mapping_status"] == "previous_symbol"]
    for _, row in prev_sym.iterrows():
        rows.append({
            "conflict_type": "outdated_name",
            "original_feature_name": row["original_feature_name"],
            "gene_id_harmonized": row.get("gene_id_harmonized"),
            "gene_symbol_harmonized": row.get("gene_symbol_harmonized"),
            "mapping_status": row["mapping_status"],
            "details": f"Mapped via previous symbol; current symbol is {row.get('gene_symbol_harmonized', '?')}",
        })

    # Ambiguous features
    ambiguous = mt[mt["mapping_status"] == "ambiguous"]
    for _, row in ambiguous.iterrows():
        rows.append({
            "conflict_type": "ambiguous",
            "original_feature_name": row["original_feature_name"],
            "gene_id_harmonized": None,
            "gene_symbol_harmonized": None,
            "mapping_status": "ambiguous",
            "details": row.get("mapping_notes", ""),
        })

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["conflict_type", "original_feature_name", "gene_id_harmonized",
                 "gene_symbol_harmonized", "mapping_status", "details"]
    )


def write_reports(result, output_dir: str, merge_result=None) -> None:
    """Write summary/conflict report files to output_dir.

    Note: harmonization_table.tsv and h5ad enrichment are handled by
    io.write_results(), not this function. This avoids duplicate writes.

    Args:
        result: HarmonizationResult.
        output_dir: Directory to write files.
        merge_result: Optional MergeResult to also write merge outputs.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Summary JSON
    s = summary(result)
    with open(os.path.join(output_dir, "summary.json"), "w") as f:
        json.dump(s, f, indent=2, default=str)

    # Conflict report
    cr = conflict_report(result)
    cr.to_csv(os.path.join(output_dir, "conflicts.tsv"), sep="\t", index=False)

    # Unmapped features
    unmapped = result.mapping_table[result.mapping_table["mapping_status"] == "unmapped"]
    unmapped.to_csv(os.path.join(output_dir, "unmapped.tsv"), sep="\t", index=False)

    # Merge outputs
    if merge_result is not None:
        merge_result.merged_table.to_csv(
            os.path.join(output_dir, "merged_table.tsv"), sep="\t", index=False,
        )
        merge_result.provenance.to_csv(
            os.path.join(output_dir, "merge_provenance.tsv"), sep="\t", index=False,
        )

    logger.info("Reports written to %s", output_dir)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_report.py -v 2>&1 | tail -15`
Expected: All 7 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/report.py tests/test_report.py
git commit -m "feat: add reporting with summary, conflict detection, and file output"
```

---

### Task 9: Conservative Merge (`merge.py`)

**Files:**
- Create: `src/stangene/merge.py`
- Create: `tests/test_merge.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_merge.py
import pandas as pd
import pytest

from stangene.harmonize import HarmonizationResult
from stangene.merge import merge_features, MergeResult


@pytest.fixture
def result_with_version_dupes():
    """Two features that map to the same gene via Tier 1 and Tier 2."""
    mt = pd.DataFrame([
        {"original_feature_name": "TP53_v1", "original_feature_id": "ENSG00000141510",
         "gene_id_harmonized": "ENSG00000141510", "gene_symbol_harmonized": "TP53",
         "mapping_status": "exact_id", "mapping_confidence": "high",
         "original_feature_type": "gene", "dataset": "ds1"},
        {"original_feature_name": "TP53_v2", "original_feature_id": "ENSG00000141510.18",
         "gene_id_harmonized": "ENSG00000141510", "gene_symbol_harmonized": "TP53",
         "mapping_status": "id_no_version", "mapping_confidence": "high",
         "original_feature_type": "gene", "dataset": "ds1"},
        {"original_feature_name": "MYC", "original_feature_id": "ENSG00000136997",
         "gene_id_harmonized": "ENSG00000136997", "gene_symbol_harmonized": "MYC",
         "mapping_status": "exact_id", "mapping_confidence": "high",
         "original_feature_type": "gene", "dataset": "ds1"},
    ])
    conflicts = mt[mt["gene_id_harmonized"] == "ENSG00000141510"].copy()
    return HarmonizationResult(mt, conflicts, mt["mapping_status"].value_counts().to_dict())


@pytest.fixture
def result_with_alias_dupes():
    """Two features mapping to the same gene, one via alias."""
    mt = pd.DataFrame([
        {"original_feature_name": "TP53", "gene_id_harmonized": "ENSG00000141510",
         "gene_symbol_harmonized": "TP53", "mapping_status": "exact_symbol",
         "mapping_confidence": "high", "original_feature_type": "gene", "dataset": "ds1"},
        {"original_feature_name": "p53", "gene_id_harmonized": "ENSG00000141510",
         "gene_symbol_harmonized": "TP53", "mapping_status": "alias_symbol",
         "mapping_confidence": "medium", "original_feature_type": "gene", "dataset": "ds1"},
    ])
    conflicts = mt.copy()
    return HarmonizationResult(mt, conflicts, mt["mapping_status"].value_counts().to_dict())


def test_strict_merges_id_based(result_with_version_dupes):
    mr = merge_features(result_with_version_dupes, policy="strict")
    assert isinstance(mr, MergeResult)
    # TP53 rows should be merged, MYC stays
    assert len(mr.merged_table) == 2
    assert "ENSG00000141510" in mr.merged_table["gene_id_harmonized"].values


def test_strict_does_not_merge_alias(result_with_alias_dupes):
    mr = merge_features(result_with_alias_dupes, policy="strict")
    # alias_symbol is not Tier 1/2, so strict should not merge
    assert len(mr.merged_table) == 2


def test_symbol_policy_merges_exact_symbol(result_with_alias_dupes):
    # Modify: make both exact_symbol (Tier 3)
    result_with_alias_dupes.mapping_table.at[1, "mapping_status"] = "exact_symbol"
    mr = merge_features(result_with_alias_dupes, policy="symbol")
    assert len(mr.merged_table) == 1


def test_provenance_tracks_originals(result_with_version_dupes):
    mr = merge_features(result_with_version_dupes, policy="strict")
    # Provenance should have 2 entries for the TP53 merge
    tp53_prov = mr.provenance[mr.provenance["gene_id_harmonized"] == "ENSG00000141510"]
    assert len(tp53_prov) == 2


def test_merge_log_not_empty(result_with_version_dupes):
    mr = merge_features(result_with_version_dupes, policy="strict")
    assert len(mr.merge_log) > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_merge.py -v 2>&1 | tail -5`
Expected: FAIL — `ImportError`

- [ ] **Step 3: Implement merge.py**

```python
# src/stangene/merge.py
"""Conservative merge logic for harmonized features."""

from dataclasses import dataclass, field

import pandas as pd

from stangene._logging import get_logger

logger = get_logger("merge")

# Mapping statuses eligible for merge under each policy
_STRICT_ELIGIBLE = frozenset(["exact_id", "id_no_version"])
_SYMBOL_ELIGIBLE = frozenset(["exact_id", "id_no_version", "exact_symbol"])


@dataclass
class MergeResult:
    """Result of conservative merge."""

    merged_table: pd.DataFrame
    provenance: pd.DataFrame
    merge_log: list = field(default_factory=list)


def merge_features(result, policy: str = "strict") -> MergeResult:
    """Merge features sharing the same gene_id_harmonized under a conservative policy.

    Args:
        result: HarmonizationResult from harmonize().
        policy: "strict" (Tier 1-2 only) or "symbol" (Tier 1-3).

    Returns:
        MergeResult with merged_table, provenance, and merge_log.
    """
    if policy == "strict":
        eligible = _STRICT_ELIGIBLE
    elif policy == "symbol":
        eligible = _SYMBOL_ELIGIBLE
    else:
        raise ValueError(f"Unknown merge policy: {policy}. Use 'strict' or 'symbol'.")

    mt = result.mapping_table.copy()
    merge_log = []
    provenance_rows = []

    # Identify merge candidates: rows with non-null gene_id_harmonized and eligible status
    eligible_mask = mt["mapping_status"].isin(eligible) & mt["gene_id_harmonized"].notna()
    eligible_df = mt[eligible_mask]
    ineligible_df = mt[~eligible_mask]

    # Group eligible rows by gene_id_harmonized
    merged_rows = []
    for hid, group in eligible_df.groupby("gene_id_harmonized"):
        if len(group) == 1:
            # No merge needed
            merged_rows.append(group.iloc[0].to_dict())
            provenance_rows.append({
                "gene_id_harmonized": hid,
                "original_feature_name": group.iloc[0]["original_feature_name"],
                "dataset": group.iloc[0].get("dataset", ""),
                "mapping_status": group.iloc[0]["mapping_status"],
                "merge_action": "kept_single",
            })
        else:
            # Merge: keep first row as representative, log all contributors
            rep = group.iloc[0].to_dict()
            originals = group["original_feature_name"].tolist()
            merge_log.append(
                f"Merged {len(group)} features into {hid} ({rep.get('gene_symbol_harmonized', '')}): {originals}"
            )
            merged_rows.append(rep)
            for _, row in group.iterrows():
                provenance_rows.append({
                    "gene_id_harmonized": hid,
                    "original_feature_name": row["original_feature_name"],
                    "dataset": row.get("dataset", ""),
                    "mapping_status": row["mapping_status"],
                    "merge_action": "merged",
                })

    # Ineligible rows pass through as-is
    for _, row in ineligible_df.iterrows():
        merged_rows.append(row.to_dict())
        provenance_rows.append({
            "gene_id_harmonized": row.get("gene_id_harmonized"),
            "original_feature_name": row["original_feature_name"],
            "dataset": row.get("dataset", ""),
            "mapping_status": row["mapping_status"],
            "merge_action": "not_eligible",
        })

    merged_table = pd.DataFrame(merged_rows)
    provenance = pd.DataFrame(provenance_rows)

    logger.info("Merge complete (policy=%s): %d -> %d rows, %d merges", policy, len(mt), len(merged_table), len(merge_log))
    return MergeResult(merged_table=merged_table, provenance=provenance, merge_log=merge_log)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_merge.py -v 2>&1 | tail -15`
Expected: All 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/merge.py tests/test_merge.py
git commit -m "feat: add conservative merge with strict and symbol policies"
```

---

### Task 10: Top-level API and CLI

**Files:**
- Modify: `src/stangene/__init__.py` (add `run()` and re-exports)
- Create: `src/stangene/__main__.py`
- Create: `tests/test_run.py`

- [ ] **Step 1: Write failing integration test**

```python
# tests/test_run.py
import os
from unittest.mock import patch

import anndata
import numpy as np
import pandas as pd
import pytest

import stangene
from stangene.harmonize import HarmonizationResult
from stangene.references import ReferenceNotFoundError


@pytest.fixture
def sample_h5ad(tmp_path):
    """Create a minimal h5ad file for integration testing."""
    adata = anndata.AnnData(
        X=np.zeros((3, 4)),
        var=pd.DataFrame(
            {"gene_ids": ["ENSG00000141510", "ENSG00000012048", "", ""]},
            index=["TP53", "BRCA1", "ERCC-00002", "NONEXISTENT"],
        ),
    )
    path = str(tmp_path / "test.h5ad")
    adata.write_h5ad(path)
    return path


@pytest.fixture
def ref_dir_with_human(tmp_path, mock_ref):
    """Create a reference directory pre-populated with mock human data."""
    ref_dir = str(tmp_path / "references")
    human_dir = os.path.join(ref_dir, "human")
    os.makedirs(human_dir)
    mock_ref["gene_table"].to_parquet(os.path.join(human_dir, "gene_table.parquet"))
    mock_ref["symbol_lookup"].to_parquet(os.path.join(human_dir, "symbol_lookup.parquet"))
    import json
    with open(os.path.join(human_dir, "build_metadata.json"), "w") as f:
        json.dump(mock_ref["metadata"], f)
    return ref_dir


def test_run_returns_result(sample_h5ad, ref_dir_with_human):
    result = stangene.run(sample_h5ad, species="human", reference_dir=ref_dir_with_human)
    assert isinstance(result, HarmonizationResult)
    assert len(result.mapping_table) == 4


def test_run_writes_output(sample_h5ad, ref_dir_with_human, tmp_path):
    output_dir = str(tmp_path / "output")
    result = stangene.run(
        sample_h5ad, species="human",
        output_dir=output_dir, reference_dir=ref_dir_with_human,
    )
    assert os.path.exists(os.path.join(output_dir, "harmonization_table.tsv"))
    assert os.path.exists(os.path.join(output_dir, "summary.json"))


def test_run_raises_without_refs(sample_h5ad, tmp_path):
    with pytest.raises(ReferenceNotFoundError):
        stangene.run(sample_h5ad, species="human", reference_dir=str(tmp_path / "empty"))


def test_run_mapping_statuses(sample_h5ad, ref_dir_with_human):
    result = stangene.run(sample_h5ad, species="human", reference_dir=ref_dir_with_human)
    statuses = set(result.mapping_table["mapping_status"])
    # Should have at least exact_id (TP53, BRCA1), non_gene_feature (ERCC), and unmapped
    assert "exact_id" in statuses
    assert "non_gene_feature" in statuses
    assert "unmapped" in statuses
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_run.py -v 2>&1 | tail -5`
Expected: FAIL — `stangene.run` not defined

- [ ] **Step 3: Update __init__.py with run() and re-exports**

```python
# src/stangene/__init__.py
"""Stangene: gene identifier harmonization for single-cell transcriptomics."""

__version__ = "0.1.0"

from stangene.classify import classify_features
from stangene.harmonize import HarmonizationResult, harmonize
from stangene.io import load_features
from stangene.merge import MergeResult, merge_features
from stangene.references import (
    ReferenceNotFoundError,
    build_reference,
    load_reference,
)
from stangene.report import conflict_report, summary, write_reports
from stangene._logging import get_logger

_logger = get_logger("run")


def run(
    path: str,
    species: str,
    output_dir: str = None,
    dataset_name: str = None,
    reference_dir: str = None,
) -> HarmonizationResult:
    """Run the full harmonization pipeline on a single dataset.

    Args:
        path: Path to input file (.h5ad, .tsv, .csv).
        species: Species name (e.g., "human", "mouse").
        output_dir: If provided, write report files here.
        dataset_name: Optional name for the dataset.
        reference_dir: Optional custom reference directory.

    Returns:
        HarmonizationResult with mapping_table, conflicts, and stats.

    Raises:
        ReferenceNotFoundError: If references have not been built for the species.
    """
    _logger.info("Starting harmonization: path=%s, species=%s", path, species)

    # Step 1: Load features
    ft = load_features(path, species=species, dataset_name=dataset_name)

    # Step 2: Classify features
    ft = classify_features(ft)

    # Step 3: Load reference (raises if not built)
    ref = load_reference(species, reference_dir=reference_dir)

    # Step 4: Harmonize
    result = harmonize(ft, ref)

    # Step 5: Write outputs if output_dir provided
    if output_dir:
        from stangene.io import write_results
        write_results(result, output_dir, input_path=path)
        write_reports(result, output_dir)

    _logger.info("Harmonization complete: %s", result.stats)
    return result
```

- [ ] **Step 4: Implement __main__.py**

```python
# src/stangene/__main__.py
"""CLI entry point for stangene."""

import argparse
import sys

from stangene._logging import get_logger

logger = get_logger("cli")


def main():
    parser = argparse.ArgumentParser(
        prog="stangene",
        description="Gene identifier harmonization for single-cell transcriptomics",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # harmonize subcommand
    harm_parser = subparsers.add_parser("harmonize", help="Harmonize gene identifiers in a dataset")
    harm_parser.add_argument("--input", required=True, help="Path to input file (.h5ad, .tsv, .csv)")
    harm_parser.add_argument("--species", required=True, help="Species name (e.g., human, mouse)")
    harm_parser.add_argument("--output-dir", required=True, help="Directory to write output files")
    harm_parser.add_argument("--dataset-name", default=None, help="Optional dataset name")
    harm_parser.add_argument("--reference-dir", default=None, help="Custom reference directory")

    # build-refs subcommand
    build_parser = subparsers.add_parser("build-refs", help="Build reference databases")
    build_parser.add_argument("--species", required=True, help="Species name (e.g., human, mouse)")
    build_parser.add_argument("--reference-dir", default=None, help="Custom reference directory")
    build_parser.add_argument("--force", action="store_true", help="Force re-download and rebuild")

    args = parser.parse_args()

    if args.command == "harmonize":
        from stangene import run
        try:
            result = run(
                path=args.input,
                species=args.species,
                output_dir=args.output_dir,
                dataset_name=args.dataset_name,
                reference_dir=args.reference_dir,
            )
            print(f"Harmonization complete. {len(result.mapping_table)} features processed.")
            print(f"Status counts: {result.stats}")
            print(f"Reports written to: {args.output_dir}")
        except Exception as e:
            logger.error("Harmonization failed: %s", e)
            sys.exit(1)

    elif args.command == "build-refs":
        from stangene.references import build_reference
        try:
            build_reference(
                species=args.species,
                reference_dir=args.reference_dir,
                force=args.force,
            )
            print(f"Reference build complete for {args.species}.")
        except Exception as e:
            logger.error("Reference build failed: %s", e)
            sys.exit(1)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_run.py -v 2>&1 | tail -15`
Expected: All 4 tests PASS

- [ ] **Step 6: Run all tests**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/ -v 2>&1 | tail -30`
Expected: All tests PASS

- [ ] **Step 7: Commit**

```bash
git add src/stangene/__init__.py src/stangene/__main__.py tests/test_run.py
git commit -m "feat: add top-level run() API and CLI entry point"
```

---

### Task 11: Write Results to h5ad (`io.py` — `write_results`)

**Files:**
- Modify: `src/stangene/io.py` (implement `write_results`)
- Modify: `tests/test_io.py` (add write tests)

- [ ] **Step 1: Write failing tests**

Append to `tests/test_io.py`:

```python
from stangene.io import write_results
from stangene.harmonize import HarmonizationResult


@pytest.fixture
def sample_result():
    mt = pd.DataFrame([
        {"original_feature_name": "TP53", "gene_id_harmonized": "ENSG00000141510",
         "gene_symbol_harmonized": "TP53", "mapping_status": "exact_id",
         "mapping_confidence": "high", "original_feature_type": "gene",
         "species": "human", "dataset": "test"},
    ])
    return HarmonizationResult(mt, pd.DataFrame(), {"exact_id": 1})


def test_write_results_creates_tsv(sample_result, tmp_path):
    output_dir = str(tmp_path / "out")
    write_results(sample_result, output_dir)
    assert os.path.exists(os.path.join(output_dir, "harmonization_table.tsv"))


def test_write_results_enriches_h5ad(sample_result, sample_h5ad, tmp_path):
    output_dir = str(tmp_path / "out")
    write_results(sample_result, output_dir, input_path=sample_h5ad)
    # Check that enriched h5ad was written
    enriched_path = os.path.join(output_dir, "test_harmonized.h5ad")
    assert os.path.exists(enriched_path)
    adata = anndata.read_h5ad(enriched_path)
    assert "gene_id_harmonized" in adata.var.columns
    # Original var_names should be preserved
    assert "TP53" in adata.var_names
```

- [ ] **Step 2: Run new tests to verify they fail**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_io.py::test_write_results_creates_tsv tests/test_io.py::test_write_results_enriches_h5ad -v 2>&1 | tail -5`
Expected: FAIL — `write_results` not implemented or import error

- [ ] **Step 3: Implement write_results in io.py**

Add to `src/stangene/io.py`:

```python
def write_results(
    result,
    output_dir: str,
    input_path: str = None,
    overwrite_h5ad: bool = False,
) -> None:
    """Write harmonization results to disk.

    Always writes harmonization_table.tsv. If input_path is an h5ad file,
    also writes an enriched copy with harmonization columns in adata.var.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Write TSV mapping table
    result.mapping_table.to_csv(
        os.path.join(output_dir, "harmonization_table.tsv"), sep="\t", index=False,
    )
    logger.info("Wrote harmonization_table.tsv to %s", output_dir)

    # Optionally enrich h5ad
    if input_path and os.path.splitext(input_path)[1].lower() in (".h5ad", ".h5"):
        _write_enriched_h5ad(result, input_path, output_dir, overwrite_h5ad)


def _write_enriched_h5ad(result, input_path: str, output_dir: str, overwrite: bool) -> None:
    """Add harmonization columns to adata.var and save."""
    adata = anndata.read_h5ad(input_path)

    # Build a mapping from original_feature_name to harmonization columns
    mt = result.mapping_table.set_index("original_feature_name")
    harm_cols = [
        "gene_id_harmonized", "gene_symbol_harmonized",
        "mapping_status", "mapping_confidence", "mapping_source", "mapping_notes",
        "original_feature_type", "feature_id_no_version",
    ]

    for col in harm_cols:
        if col in mt.columns:
            # Map by var_names
            adata.var[col] = adata.var_names.map(
                mt[col].to_dict() if not mt.index.duplicated().any()
                else mt[~mt.index.duplicated(keep="first")][col].to_dict()
            )

    if overwrite:
        adata.write_h5ad(input_path)
        logger.info("Overwrote %s with harmonization columns", input_path)
    else:
        basename = os.path.splitext(os.path.basename(input_path))[0]
        out_path = os.path.join(output_dir, f"{basename}_harmonized.h5ad")
        adata.write_h5ad(out_path)
        logger.info("Wrote enriched h5ad to %s", out_path)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/test_io.py -v 2>&1 | tail -20`
Expected: All 12 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/stangene/io.py tests/test_io.py
git commit -m "feat: add write_results with h5ad enrichment support"
```

---

### Task 12: Claude Code Skill File

**Files:**
- Create: `skill.md` (at project root, the Claude Code skill definition)

- [ ] **Step 1: Write the skill markdown**

```markdown
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
```

- [ ] **Step 2: Commit**

```bash
git add skill.md
git commit -m "feat: add Claude Code skill definition for gene harmonization"
```

---

### Task 13: Final Integration Test and README

**Files:**
- Create: `README.md`
- Run full test suite

- [ ] **Step 1: Run full test suite**

Run: `cd /scratch/users/chensj16/projects/stangene && python -m pytest tests/ -v --tb=short 2>&1`
Expected: All tests PASS

- [ ] **Step 2: Write README.md**

```markdown
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
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: add README with usage instructions and design overview"
```

- [ ] **Step 4: Final commit log review**

Run: `cd /scratch/users/chensj16/projects/stangene && git log --oneline`

Verify the commit history is clean and incremental.
