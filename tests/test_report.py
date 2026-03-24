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
    assert s["duplicate_harmonized_ids"] == 1
    assert s["duplicate_harmonized_symbols"] == 1


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
