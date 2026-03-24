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
    assert len(mr.merged_table) == 2
    assert "ENSG00000141510" in mr.merged_table["gene_id_harmonized"].values


def test_strict_does_not_merge_alias(result_with_alias_dupes):
    mr = merge_features(result_with_alias_dupes, policy="strict")
    assert len(mr.merged_table) == 2


def test_symbol_policy_merges_exact_symbol(result_with_alias_dupes):
    result_with_alias_dupes.mapping_table.at[1, "mapping_status"] = "exact_symbol"
    mr = merge_features(result_with_alias_dupes, policy="symbol")
    assert len(mr.merged_table) == 1


def test_provenance_tracks_originals(result_with_version_dupes):
    mr = merge_features(result_with_version_dupes, policy="strict")
    tp53_prov = mr.provenance[mr.provenance["gene_id_harmonized"] == "ENSG00000141510"]
    assert len(tp53_prov) == 2


def test_merge_log_not_empty(result_with_version_dupes):
    mr = merge_features(result_with_version_dupes, policy="strict")
    assert len(mr.merge_log) > 0
