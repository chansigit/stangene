"""Tests for stangene.harmonize_anndata (in-memory AnnData harmonization)."""

import anndata
import numpy as np

from stangene import harmonize_anndata


def _adata(var_names, gene_ids=None):
    a = anndata.AnnData(X=np.zeros((2, len(var_names)), dtype="float32"))
    a.var_names = list(var_names)
    if gene_ids is not None:
        a.var["gene_ids"] = list(gene_ids)
    return a


def test_rows_aligned_and_symbols_resolved():
    # TP53 (exact symbol), p53 (alias), ENSG00000141510 (TP53's Ensembl id),
    # and a bogus feature that must stay unmapped.
    names = ["TP53", "p53", "ENSG00000141510", "NOT_A_GENE_XYZ"]
    res = harmonize_anndata(_adata(names), "human")
    mt = res.mapping_table

    # row order preserved -> positional alignment with adata.var_names
    assert list(mt["original_feature_name"]) == names

    by_name = mt.set_index("original_feature_name")
    assert by_name.loc["TP53", "gene_symbol_harmonized"] == "TP53"
    assert by_name.loc["ENSG00000141510", "gene_symbol_harmonized"] == "TP53"
    assert by_name.loc["NOT_A_GENE_XYZ", "mapping_status"] == "unmapped"


def test_uses_gene_ids_column_when_present():
    # var_names are junk, but gene_ids carry the real Ensembl id -> still maps.
    a = _adata(["junk1"], gene_ids=["ENSG00000141510"])
    res = harmonize_anndata(a, "human")
    row = res.mapping_table.iloc[0]
    assert row["gene_symbol_harmonized"] == "TP53"


def test_does_not_mutate_adata():
    a = _adata(["TP53", "p53"])
    before = list(a.var_names)
    harmonize_anndata(a, "human")
    assert list(a.var_names) == before  # var_names untouched
