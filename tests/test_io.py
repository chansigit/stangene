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
    cd3_row = ft[ft["original_feature_name"] == "CD3_ADT"].iloc[0]
    assert pd.isna(cd3_row["original_feature_id"]) or cd3_row["original_feature_id"] == ""


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
    enriched_path = os.path.join(output_dir, "test_harmonized.h5ad")
    assert os.path.exists(enriched_path)
    adata = anndata.read_h5ad(enriched_path)
    assert "gene_id_harmonized" in adata.var.columns
    assert "TP53" in adata.var_names
