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
        ft = _make_ft(["ERCC-00002"], types=["spike_in"])
        ft["mapping_status"] = "non_gene_feature"
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "non_gene_feature"


class TestEarlyExit:
    def test_id_match_skips_symbol(self, mock_ref):
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
