import os

import pandas as pd
import pytest

from stangene.biotype import CANONICAL_BIOTYPES, normalize_biotype

_REFS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "src", "stangene", "data", "refs"
)


def _gene_table(species: str) -> pd.DataFrame:
    path = os.path.join(_REFS_DIR, species, "gene_table.parquet")
    if not os.path.exists(path):
        pytest.skip(f"Reference not built for {species}")
    return pd.read_parquet(path)


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_column_exists(species):
    gt = _gene_table(species)
    assert "canonical_biotype" in gt.columns, f"{species}: missing canonical_biotype"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_no_nulls(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    assert gt["canonical_biotype"].isna().sum() == 0, f"{species}: NaN in canonical_biotype"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
])
def test_canonical_biotype_in_vocabulary(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    bad = set(gt["canonical_biotype"].unique()) - CANONICAL_BIOTYPES
    assert not bad, f"{species}: out-of-vocabulary values: {bad}"


@pytest.mark.parametrize("species", [
    "human", "mouse", "rat",
])
def test_canonical_biotype_has_protein_coding(species):
    gt = _gene_table(species)
    if "canonical_biotype" not in gt.columns:
        pytest.skip("column absent")
    n = (gt["canonical_biotype"] == "protein_coding").sum()
    assert n > 10_000, f"{species}: unexpectedly few protein_coding genes ({n})"


# --- vocabulary ---

def test_canonical_biotypes_is_frozenset():
    assert isinstance(CANONICAL_BIOTYPES, frozenset)


def test_canonical_biotypes_contains_expected():
    expected = {
        "protein_coding", "lncRNA", "pseudogene",
        "miRNA", "snoRNA", "snRNA", "rRNA", "tRNA",
        "IG_gene", "TR_gene", "other_ncrna", "other", "unknown",
    }
    assert expected == CANONICAL_BIOTYPES


# --- normalize_biotype: HGNC (locus_type values) ---

@pytest.mark.parametrize("raw,expected", [
    ("gene with protein product", "protein_coding"),
    ("RNA, micro",                "miRNA"),
    ("RNA, small nucleolar",      "snoRNA"),
    ("RNA, small nuclear",        "snRNA"),
    ("RNA, ribosomal",            "rRNA"),
    ("RNA, transfer",             "tRNA"),
    ("RNA, long non-coding",      "lncRNA"),
    ("immunoglobulin gene",       "IG_gene"),
    ("T cell receptor gene",      "TR_gene"),
    ("pseudogene",                "pseudogene"),
    ("endogenous retrovirus",     "other"),
    ("unknown",                   "unknown"),
])
def test_hgnc_mappings(raw, expected):
    assert normalize_biotype(raw, "HGNC") == expected


def test_hgnc_case_insensitive():
    assert normalize_biotype("Gene With Protein Product", "HGNC") == "protein_coding"
    assert normalize_biotype("  RNA, Micro  ", "HGNC") == "miRNA"


# --- normalize_biotype: MGI ---

@pytest.mark.parametrize("raw,expected", [
    ("protein coding gene",              "protein_coding"),
    ("lncRNA gene",                      "lncRNA"),
    ("pseudogene",                       "pseudogene"),
    ("miRNA gene",                       "miRNA"),
    ("snoRNA gene",                      "snoRNA"),
    ("snRNA gene",                       "snRNA"),
    ("rRNA gene",                        "rRNA"),
    ("tRNA gene",                        "tRNA"),
    ("unclassified non-coding RNA gene", "other_ncrna"),
    ("transgene",                        "other"),
    ("unclassified gene",                "other"),
    ("gene segment",                     "other"),
])
def test_mgi_mappings(raw, expected):
    assert normalize_biotype(raw, "MGI") == expected


# --- normalize_biotype: RGD ---

@pytest.mark.parametrize("raw,expected", [
    ("protein-coding", "protein_coding"),
    ("ncrna",          "other_ncrna"),
    ("lincrna",        "lncRNA"),
    ("lncrna",         "lncRNA"),
    ("mirna",          "miRNA"),
    ("snorna",         "snoRNA"),
    ("snrna",          "snRNA"),
    ("rrna",           "rRNA"),
    ("trna",           "tRNA"),
    ("pseudo",         "pseudogene"),
    ("pseudogene",     "pseudogene"),
])
def test_rgd_mappings(raw, expected):
    assert normalize_biotype(raw, "RGD") == expected


# --- normalize_biotype: Ensembl (BioMart gene_biotype) ---

@pytest.mark.parametrize("raw,expected", [
    ("protein_coding",                          "protein_coding"),
    ("lincRNA",                                 "lncRNA"),
    ("lncRNA",                                  "lncRNA"),
    ("miRNA",                                   "miRNA"),
    ("snoRNA",                                  "snoRNA"),
    ("snRNA",                                   "snRNA"),
    ("rRNA",                                    "rRNA"),
    ("Mt_rRNA",                                 "rRNA"),
    ("tRNA",                                    "tRNA"),
    ("Mt_tRNA",                                 "tRNA"),
    ("pseudogene",                              "pseudogene"),
    ("processed_pseudogene",                    "pseudogene"),
    ("transcribed_unprocessed_pseudogene",      "pseudogene"),
    ("IG_C_gene",                               "IG_gene"),
    ("IG_V_gene",                               "IG_gene"),
    ("TR_C_gene",                               "TR_gene"),
    ("TR_J_gene",                               "TR_gene"),
    ("misc_RNA",                                "other_ncrna"),
    ("ribozyme",                                "other_ncrna"),
    ("vault_RNA",                               "other_ncrna"),
    ("TEC",                                     "other"),
])
def test_ensembl_mappings(raw, expected):
    assert normalize_biotype(raw, "Ensembl") == expected


# --- normalize_biotype: WormBase ---

@pytest.mark.parametrize("raw,expected", [
    ("protein_coding_gene", "protein_coding"),
    ("pseudogene",          "pseudogene"),
    ("lncRNA_gene",         "lncRNA"),
    ("miRNA_gene",          "miRNA"),
    ("snoRNA_gene",         "snoRNA"),
    ("snRNA_gene",          "snRNA"),
    ("rRNA_gene",           "rRNA"),
    ("tRNA_gene",           "tRNA"),
    ("ncRNA_gene",          "other_ncrna"),
    ("transposable_element_gene", "other"),
])
def test_wormbase_mappings(raw, expected):
    assert normalize_biotype(raw, "WormBase") == expected


# --- normalize_biotype: ZFIN (SO IDs) ---

@pytest.mark.parametrize("raw,expected", [
    ("SO:0001217", "protein_coding"),
    ("SO:0000336", "pseudogene"),
    ("SO:0001263", "miRNA"),
    ("SO:0001267", "snoRNA"),
    ("SO:0001268", "snRNA"),
    ("SO:0001637", "rRNA"),
    ("SO:0001272", "tRNA"),
    ("SO:0001877", "lncRNA"),
])
def test_zfin_mappings(raw, expected):
    assert normalize_biotype(raw, "ZFIN") == expected


# --- normalize_biotype: FlyBase (annotation_id prefix) ---

@pytest.mark.parametrize("raw,expected", [
    ("CG12345",  "protein_coding"),
    ("CR12345",  "other_ncrna"),
    ("FBtr0000", "other"),
    ("",         "unknown"),
])
def test_flybase_mappings(raw, expected):
    assert normalize_biotype(raw, "FlyBase") == expected


# --- unknown source / unknown value ---

def test_unknown_source_returns_unknown():
    assert normalize_biotype("protein-coding gene", "UnknownDB") == "unknown"


def test_empty_raw_returns_unknown():
    assert normalize_biotype("", "HGNC") == "unknown"


def test_whitespace_raw_returns_unknown():
    assert normalize_biotype("   ", "MGI") == "unknown"


# --- all return values are in vocabulary ---

@pytest.mark.parametrize("source,raw", [
    ("HGNC",      "gene with protein product"),
    ("MGI",       "miRNA gene"),
    ("RGD",       "lincrna"),
    ("Ensembl",   "processed_pseudogene"),
    ("WormBase",  "rRNA_gene"),
    ("ZFIN",      "SO:0001217"),
    ("FlyBase",   "CG9999"),
    ("HGNC",      "something_totally_new"),
])
def test_return_value_always_in_vocabulary(source, raw):
    result = normalize_biotype(raw, source)
    assert result in CANONICAL_BIOTYPES
