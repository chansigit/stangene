"""Tests for markdown report generation."""

import os

import pandas as pd
import pytest

from stangene.harmonize import HarmonizationResult
from stangene.report import generate_markdown_report


@pytest.fixture
def rich_result():
    """A HarmonizationResult with all mapping statuses represented."""
    mt = pd.DataFrame([
        # Tier 1: exact ID
        {"original_feature_name": "TP53", "original_feature_id": "ENSG00000141510",
         "gene_id_harmonized": "ENSG00000141510", "gene_symbol_harmonized": "TP53",
         "mapping_status": "exact_id", "mapping_confidence": "high",
         "mapping_source": "HGNC:exact_id", "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Tier 2: version-stripped
        {"original_feature_name": "BRCA1", "original_feature_id": "ENSG00000012048.23",
         "gene_id_harmonized": "ENSG00000012048", "gene_symbol_harmonized": "BRCA1",
         "mapping_status": "id_no_version", "mapping_confidence": "high",
         "mapping_source": "HGNC:id_no_version",
         "mapping_notes": "Matched via version-stripped ID: ENSG00000012048.23 -> ENSG00000012048",
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Tier 3: exact symbol
        {"original_feature_name": "MYC", "original_feature_id": None,
         "gene_id_harmonized": "ENSG00000136997", "gene_symbol_harmonized": "MYC",
         "mapping_status": "exact_symbol", "mapping_confidence": "high",
         "mapping_source": "HGNC:approved_symbol", "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Tier 3 duplicate: another feature -> same MYC
        {"original_feature_name": "MYC_dup", "original_feature_id": None,
         "gene_id_harmonized": "ENSG00000136997", "gene_symbol_harmonized": "MYC",
         "mapping_status": "exact_symbol", "mapping_confidence": "high",
         "mapping_source": "HGNC:approved_symbol", "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Tier 4: alias
        {"original_feature_name": "p53", "original_feature_id": None,
         "gene_id_harmonized": "ENSG00000141510", "gene_symbol_harmonized": "TP53",
         "mapping_status": "alias_symbol", "mapping_confidence": "medium",
         "mapping_source": "HGNC:alias_symbol", "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Tier 4: previous symbol
        {"original_feature_name": "RNF53", "original_feature_id": None,
         "gene_id_harmonized": "ENSG00000012048", "gene_symbol_harmonized": "BRCA1",
         "mapping_status": "previous_symbol", "mapping_confidence": "medium",
         "mapping_source": "HGNC:prev_symbol",
         "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Ambiguous
        {"original_feature_name": "AMBIG_GENE", "original_feature_id": None,
         "gene_id_harmonized": None, "gene_symbol_harmonized": None,
         "mapping_status": "ambiguous", "mapping_confidence": "low",
         "mapping_source": None,
         "mapping_notes": "Multiple alias/prev matches: [ENSG001, ENSG002]",
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Unmapped
        {"original_feature_name": "FAKE_GENE1", "original_feature_id": None,
         "gene_id_harmonized": None, "gene_symbol_harmonized": None,
         "mapping_status": "unmapped", "mapping_confidence": None,
         "mapping_source": None, "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        {"original_feature_name": "FAKE_GENE2", "original_feature_id": None,
         "gene_id_harmonized": None, "gene_symbol_harmonized": None,
         "mapping_status": "unmapped", "mapping_confidence": None,
         "mapping_source": None, "mapping_notes": None,
         "original_feature_type": "gene", "species": "human", "dataset": "test_ds"},
        # Non-gene feature
        {"original_feature_name": "ERCC-00002", "original_feature_id": None,
         "gene_id_harmonized": None, "gene_symbol_harmonized": None,
         "mapping_status": "non_gene_feature", "mapping_confidence": None,
         "mapping_source": None, "mapping_notes": None,
         "original_feature_type": "spike_in", "species": "human", "dataset": "test_ds"},
        {"original_feature_name": "CD3_ADT", "original_feature_id": None,
         "gene_id_harmonized": None, "gene_symbol_harmonized": None,
         "mapping_status": "non_gene_feature", "mapping_confidence": None,
         "mapping_source": None, "mapping_notes": None,
         "original_feature_type": "antibody_capture", "species": "human", "dataset": "test_ds"},
    ])
    conflicts = mt[mt["gene_id_harmonized"].isin(["ENSG00000141510", "ENSG00000136997"])].copy()
    stats = mt["mapping_status"].value_counts().to_dict()
    return HarmonizationResult(mapping_table=mt, conflicts=conflicts, stats=stats)


class TestMarkdownReportContent:
    def test_returns_string(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert isinstance(md, str)

    def test_contains_title(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "# Gene Harmonization Report" in md

    def test_contains_dataset_name(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "test_ds" in md

    def test_contains_species(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "human" in md

    def test_contains_summary_table(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "Total features" in md
        assert "| 11 |" in md  # total features

    def test_contains_tier_breakdown(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "exact_id" in md
        assert "exact_symbol" in md
        assert "alias_symbol" in md

    def test_contains_confidence_breakdown(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "high" in md.lower()
        assert "medium" in md.lower()

    def test_contains_conflict_section(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "Many-to-One" in md or "many-to-one" in md.lower()

    def test_contains_unmapped_section(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "FAKE_GENE1" in md
        assert "FAKE_GENE2" in md

    def test_contains_non_gene_section(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "ERCC-00002" in md or "Non-Gene" in md or "non_gene" in md

    def test_contains_outdated_names(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "RNF53" in md

    def test_contains_ambiguous_section(self, rich_result):
        md = generate_markdown_report(rich_result)
        assert "AMBIG_GENE" in md


class TestMarkdownReportOptions:
    def test_max_rows_limits_unmapped(self, rich_result):
        md = generate_markdown_report(rich_result, max_unmapped_rows=1)
        # Should show 1 unmapped + "and N more"
        assert "and 1 more" in md

    def test_custom_title(self, rich_result):
        md = generate_markdown_report(rich_result, title="Custom Report Title")
        assert "# Custom Report Title" in md


class TestMarkdownReportWriteToFile:
    def test_write_to_file(self, rich_result, tmp_path):
        out = str(tmp_path / "report.md")
        generate_markdown_report(rich_result, output_path=out)
        assert os.path.exists(out)
        with open(out) as f:
            content = f.read()
        assert "# Gene Harmonization Report" in content

    def test_write_reports_includes_markdown(self, rich_result, tmp_path):
        """write_reports() should also produce report.md."""
        from stangene.report import write_reports
        output_dir = str(tmp_path / "output")
        write_reports(rich_result, output_dir)
        assert os.path.exists(os.path.join(output_dir, "report.md"))
