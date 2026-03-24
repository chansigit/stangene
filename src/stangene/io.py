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

    ft["species"] = species
    ft["dataset"] = dataset_name or os.path.splitext(os.path.basename(path))[0]

    if "original_feature_id" in ft.columns:
        ft["feature_id_no_version"] = ft["original_feature_id"].apply(_strip_version)
        ft.loc[ft["feature_id_no_version"] == "", "feature_id_no_version"] = None
    else:
        ft["original_feature_id"] = None
        ft["feature_id_no_version"] = None

    ft["original_feature_id"] = ft["original_feature_id"].replace("", None)
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

    if "gene_ids" in var.columns:
        ft["original_feature_id"] = var["gene_ids"].values
    if "feature_types" in var.columns:
        ft["original_feature_type"] = var["feature_types"].values

    if hasattr(adata, "file"):
        adata.file.close()

    return ft.reset_index(drop=True)


def _load_tabular(path: str, column_map: dict = None) -> pd.DataFrame:
    """Load feature metadata from a TSV/CSV file."""
    ext = os.path.splitext(path)[1].lower()
    sep = "\t" if ext in (".tsv", ".txt") else ","
    raw = pd.read_csv(path, sep=sep)

    if column_map is None:
        column_map = {}
        for col in raw.columns:
            col_lower = col.lower().strip()
            if col_lower in _AUTO_COLUMN_MAP:
                column_map[col] = _AUTO_COLUMN_MAP[col_lower]

    ft = pd.DataFrame()
    for src_col, dst_col in column_map.items():
        if src_col in raw.columns:
            ft[dst_col] = raw[src_col].values

    if "original_feature_name" not in ft.columns:
        ft["original_feature_name"] = raw.iloc[:, 0].values

    return ft.reset_index(drop=True)
