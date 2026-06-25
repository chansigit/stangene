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


# ---------------------------------------------------------------------------
# Integration tests for non-Ensembl species (FlyBase, Ensembl BioMart)
# ---------------------------------------------------------------------------

def _make_flybase_mock_ref():
    gene_table = pd.DataFrame([
        {"ensembl_id": "FBgn0003996", "symbol": "w", "alias_symbols": "white|CG2759",
         "prev_symbols": "FBgn0000058|FBgn0000083", "gene_type": "",
         "status": "approved", "source": "FlyBase", "source_id": "FlyBase:FBgn0003996"},
        {"ensembl_id": "FBgn0011676", "symbol": "Nos", "alias_symbols": "dNOS|NOS1",
         "prev_symbols": "", "gene_type": "", "status": "approved",
         "source": "FlyBase", "source_id": "FlyBase:FBgn0011676"},
    ])
    rows = [
        {"lookup_string": "w", "lookup_string_upper": "W", "ensembl_id": "FBgn0003996",
         "source_id": "FlyBase:FBgn0003996", "lookup_type": "approved_symbol", "source": "FlyBase"},
        {"lookup_string": "white", "lookup_string_upper": "WHITE", "ensembl_id": "FBgn0003996",
         "source_id": "FlyBase:FBgn0003996", "lookup_type": "alias_symbol", "source": "FlyBase"},
        {"lookup_string": "FBgn0000058", "lookup_string_upper": "FBGN0000058", "ensembl_id": "FBgn0003996",
         "source_id": "FlyBase:FBgn0003996", "lookup_type": "prev_symbol", "source": "FlyBase"},
        {"lookup_string": "Nos", "lookup_string_upper": "NOS", "ensembl_id": "FBgn0011676",
         "source_id": "FlyBase:FBgn0011676", "lookup_type": "approved_symbol", "source": "FlyBase"},
        {"lookup_string": "dNOS", "lookup_string_upper": "DNOS", "ensembl_id": "FBgn0011676",
         "source_id": "FlyBase:FBgn0011676", "lookup_type": "alias_symbol", "source": "FlyBase"},
    ]
    return {
        "gene_table": gene_table,
        "symbol_lookup": pd.DataFrame(rows),
        "metadata": {"species": "fruit_fly", "download_timestamp": "2026-01-01T00:00:00Z"},
    }


def _make_biomart_mock_ref():
    gene_table = pd.DataFrame([
        {"ensembl_id": "ENSMFAG00000001234", "symbol": "TP53", "alias_symbols": "LFS1|p53",
         "prev_symbols": "", "gene_type": "protein_coding", "status": "approved",
         "source": "Ensembl", "source_id": "Ensembl:ENSMFAG00000001234"},
        {"ensembl_id": "ENSMFAG00000009999", "symbol": "", "alias_symbols": "",
         "prev_symbols": "", "gene_type": "protein_coding", "status": "approved",
         "source": "Ensembl", "source_id": "Ensembl:ENSMFAG00000009999"},
    ])
    rows = [
        {"lookup_string": "TP53", "lookup_string_upper": "TP53", "ensembl_id": "ENSMFAG00000001234",
         "source_id": "Ensembl:ENSMFAG00000001234", "lookup_type": "approved_symbol", "source": "Ensembl"},
        {"lookup_string": "LFS1", "lookup_string_upper": "LFS1", "ensembl_id": "ENSMFAG00000001234",
         "source_id": "Ensembl:ENSMFAG00000001234", "lookup_type": "alias_symbol", "source": "Ensembl"},
    ]
    return {
        "gene_table": gene_table,
        "symbol_lookup": pd.DataFrame(rows),
        "metadata": {"species": "cynomolgus", "download_timestamp": "2026-01-01T00:00:00Z"},
    }


def _make_ft_species(names, ids=None, species="fruit_fly"):
    n = len(names)
    return pd.DataFrame({
        "original_feature_name": names,
        "species": [species] * n,
        "dataset": ["test"] * n,
        "original_feature_type": ["gene"] * n,
        "mapping_status": [None] * n,
        "mapping_notes": [None] * n,
        "original_feature_id": ids if ids else [None] * n,
        "feature_id_no_version": [None] * n,
    })


class TestFlyBaseHarmonization:
    def test_exact_fbgn_id_matches_tier1(self):
        """FBgn ID as original_feature_id should match Tier 1 (exact_id)."""
        ref = _make_flybase_mock_ref()
        ft = _make_ft_species(["w"], ids=["FBgn0003996"])
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_id"
        assert row["gene_id_harmonized"] == "FBgn0003996"
        assert row["gene_symbol_harmonized"] == "w"

    def test_fbgn_symbol_matches_tier3(self):
        """FlyBase symbol 'w' should match via exact_symbol."""
        ref = _make_flybase_mock_ref()
        ft = _make_ft_species(["w"])
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "FBgn0003996"

    def test_fbgn_alias_matches_tier4(self):
        """FlyBase alias 'white' should match via alias_symbol."""
        ref = _make_flybase_mock_ref()
        ft = _make_ft_species(["white"])
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "alias_symbol"
        assert row["gene_id_harmonized"] == "FBgn0003996"

    def test_versioned_fbgn_matches_tier2(self):
        """FBgn0003996.2 (versioned) should match via id_no_version."""
        from stangene.io import _strip_version
        # First verify _strip_version handles FBgn:
        assert _strip_version("FBgn0003996.2") == "FBgn0003996"

        ref = _make_flybase_mock_ref()
        ft = pd.DataFrame({
            "original_feature_name": ["w_versioned"],
            "species": ["fruit_fly"],
            "dataset": ["test"],
            "original_feature_type": ["gene"],
            "mapping_status": [None],
            "mapping_notes": [None],
            "original_feature_id": ["FBgn0003996.2"],
            "feature_id_no_version": ["FBgn0003996"],
        })
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "id_no_version"
        assert row["gene_id_harmonized"] == "FBgn0003996"


class TestNonProteinCodingNote:
    def test_non_protein_coding_note_in_mapping_notes(self, mock_ref):
        """Matching a non-protein-coding gene should add 'Non-protein-coding: <biotype>' to mapping_notes."""
        import pandas as _pd

        # Add a non-withdrawn lncRNA gene to gene_table and symbol_lookup
        extra_gene = _pd.DataFrame([{
            "ensembl_id": "ENSG00000260612", "symbol": "LINC00261",
            "alias_symbols": "", "prev_symbols": "",
            "gene_type": "lncRNA gene", "status": "Approved",
            "source": "HGNC", "source_id": "HGNC:27518",
            "canonical_biotype": "lncRNA",
        }])
        extra_lookup = _pd.DataFrame([{
            "lookup_string": "LINC00261", "lookup_string_upper": "LINC00261",
            "ensembl_id": "ENSG00000260612", "source_id": "HGNC:27518",
            "lookup_type": "approved_symbol", "source": "HGNC",
        }])
        mock_ref["gene_table"] = _pd.concat(
            [mock_ref["gene_table"], extra_gene], ignore_index=True,
        )
        mock_ref["symbol_lookup"] = _pd.concat(
            [mock_ref["symbol_lookup"], extra_lookup], ignore_index=True,
        )

        ft = _make_ft(["LINC00261"])
        result = harmonize(ft, mock_ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_symbol"
        assert row["gene_id_harmonized"] == "ENSG00000260612"
        assert "Non-protein-coding: lncRNA" in str(row["mapping_notes"])


class TestBioMartHarmonization:
    def test_exact_ensmfag_id_matches_tier1(self):
        """ENSMFAG (cynomolgus) ID should match Tier 1."""
        ref = _make_biomart_mock_ref()
        ft = _make_ft_species(["TP53"], ids=["ENSMFAG00000001234"], species="cynomolgus")
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_id"
        assert row["gene_id_harmonized"] == "ENSMFAG00000001234"

    def test_biomart_symbolless_gene_reachable_by_id(self):
        """Gene with empty symbol in BioMart can still be Tier-1 matched by ID."""
        ref = _make_biomart_mock_ref()
        ft = _make_ft_species(["orphan"], ids=["ENSMFAG00000009999"], species="cynomolgus")
        result = harmonize(ft, ref)
        row = result.mapping_table.iloc[0]
        assert row["mapping_status"] == "exact_id"
        assert row["gene_id_harmonized"] == "ENSMFAG00000009999"
