"""Edge case tests for robustness: empty inputs, case sensitivity, format detection."""

import os

import anndata
import numpy as np
import pandas as pd
import pytest

from stangene.classify import classify_features
from stangene.harmonize import harmonize, HarmonizationResult
from stangene.io import load_features
from stangene.report import summary, conflict_report


# ---------------------------------------------------------------------------
# Empty input tests (#12)
# ---------------------------------------------------------------------------

class TestEmptyInputs:
    def test_load_empty_h5ad(self, tmp_path):
        """h5ad with zero features should load without crashing."""
        adata = anndata.AnnData(
            X=np.zeros((0, 0)),
            var=pd.DataFrame(index=[]),
        )
        path = str(tmp_path / "empty.h5ad")
        adata.write_h5ad(path)
        ft = load_features(path, species="human", dataset_name="empty")
        assert len(ft) == 0
        assert "original_feature_name" in ft.columns

    def test_load_empty_tsv(self, tmp_path):
        """TSV with header only should load without crashing."""
        path = str(tmp_path / "empty.tsv")
        pd.DataFrame({"gene_name": [], "gene_id": []}).to_csv(path, sep="\t", index=False)
        ft = load_features(path, species="human", dataset_name="empty")
        assert len(ft) == 0

    def test_classify_empty_df(self):
        """classify_features on empty DataFrame should return empty."""
        ft = pd.DataFrame({
            "original_feature_name": [],
            "species": [],
            "dataset": [],
        })
        result = classify_features(ft)
        assert len(result) == 0

    def test_harmonize_empty_df(self, mock_ref):
        """harmonize on empty DataFrame should return empty result."""
        ft = pd.DataFrame({
            "original_feature_name": [],
            "original_feature_id": [],
            "feature_id_no_version": [],
            "original_feature_type": [],
            "mapping_status": [],
            "mapping_notes": [],
            "species": [],
            "dataset": [],
        })
        result = harmonize(ft, mock_ref)
        assert isinstance(result, HarmonizationResult)
        assert len(result.mapping_table) == 0
        assert len(result.conflicts) == 0

    def test_summary_empty_result(self, mock_ref):
        """summary on empty result should not crash."""
        ft = pd.DataFrame({
            "original_feature_name": [],
            "original_feature_id": [],
            "feature_id_no_version": [],
            "original_feature_type": [],
            "mapping_status": [],
            "mapping_notes": [],
            "species": [],
            "dataset": [],
        })
        result = harmonize(ft, mock_ref)
        s = summary(result)
        assert s["total_features"] == 0


# ---------------------------------------------------------------------------
# Case-sensitivity tests (#13)
# ---------------------------------------------------------------------------

def _make_ft(names, ids=None, types=None):
    n = len(names)
    data = {
        "original_feature_name": names,
        "species": ["human"] * n,
        "dataset": ["test"] * n,
        "original_feature_type": types or (["gene"] * n),
        "mapping_status": [None] * n,
        "mapping_notes": [None] * n,
        "original_feature_id": ids or ([None] * n),
        "feature_id_no_version": [None] * n,
    }
    return pd.DataFrame(data)


class TestCaseSensitivity:
    def test_exact_case_matches_first(self, mock_ref):
        """TP53 (exact case) should match via exact-case lookup."""
        ft = _make_ft(["TP53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000141510"

    def test_lowercase_matches_via_fallback(self, mock_ref):
        """tp53 (lowercase) should match TP53 via uppercase fallback."""
        ft = _make_ft(["tp53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000141510"

    def test_mixed_case_matches_via_fallback(self, mock_ref):
        """Tp53 (mixed case) should match TP53 via uppercase fallback."""
        ft = _make_ft(["Tp53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000141510"

    def test_lowercase_alias_matches_via_fallback(self, mock_ref):
        """P53 (uppercase alias) should match via alias lookup."""
        ft = _make_ft(["P53"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        # p53 is the alias in mock data; P53 should match via uppercase fallback
        assert row["mapping_status"] == "alias_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000141510"


# ---------------------------------------------------------------------------
# .txt extension test (#14)
# ---------------------------------------------------------------------------

class TestTxtExtension:
    def test_load_txt_as_tsv(self, tmp_path):
        """A .txt file should be parsed as tab-separated."""
        path = str(tmp_path / "features.txt")
        pd.DataFrame({
            "gene_name": ["TP53", "BRCA1"],
            "gene_id": ["ENSG00000141510", "ENSG00000012048"],
        }).to_csv(path, sep="\t", index=False)
        ft = load_features(path, species="human")
        assert len(ft) == 2
        assert ft["original_feature_name"].iloc[0] == "TP53"


# ---------------------------------------------------------------------------
# Auto-detection fallback test (#15)
# ---------------------------------------------------------------------------

class TestAutoDetectionFallback:
    def test_unrecognized_columns_use_first(self, tmp_path):
        """When no columns match auto-detection, first column is used as feature name."""
        path = str(tmp_path / "weird.tsv")
        pd.DataFrame({
            "MyCustomGeneCol": ["TP53", "BRCA1"],
            "SomeOtherCol": ["x", "y"],
        }).to_csv(path, sep="\t", index=False)
        ft = load_features(path, species="human")
        assert len(ft) == 2
        assert ft["original_feature_name"].iloc[0] == "TP53"

    def test_column_map_overrides_auto(self, tmp_path):
        """Explicit column_map should override auto-detection."""
        path = str(tmp_path / "custom.csv")
        pd.DataFrame({
            "Symbol": ["TP53", "BRCA1"],
            "EnsemblID": ["ENSG00000141510", "ENSG00000012048"],
        }).to_csv(path, index=False)
        ft = load_features(
            path, species="human",
            column_map={"Symbol": "original_feature_name", "EnsemblID": "original_feature_id"},
        )
        assert ft["original_feature_name"].iloc[0] == "TP53"
        assert ft["original_feature_id"].iloc[0] == "ENSG00000141510"
