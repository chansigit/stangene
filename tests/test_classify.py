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
