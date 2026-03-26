import json
import os
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from stangene.references import (
    build_reference,
    load_reference,
    ReferenceNotFoundError,
    _build_human_reference,
    _default_reference_dir,
)


@pytest.fixture
def mock_hgnc_data():
    """Minimal HGNC-format TSV data."""
    return (
        "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tensembl_gene_id\tlocus_group\tstatus\n"
        'HGNC:11998\tTP53\tLFS1|p53\tLFS1\tENSG00000141510\tprotein-coding gene\tApproved\n'
        'HGNC:1100\tBRCA1\tRNF53|IRIS\tRNF53\tENSG00000012048\tprotein-coding gene\tApproved\n'
        'HGNC:1101\tBRCA2\tFACD|FANCD1\tFACD\tENSG00000139618\tprotein-coding gene\tApproved\n'
        'HGNC:999\tOLDGENE\t\t\t\tpseudogene\tEntry Withdrawn\n'
    )


@pytest.fixture
def ref_dir(tmp_path):
    return str(tmp_path / "references")


def test_reference_not_found_error(ref_dir):
    with pytest.raises(ReferenceNotFoundError):
        load_reference("human", reference_dir=ref_dir)


def test_build_human_creates_files(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    human_dir = os.path.join(ref_dir, "human")
    assert os.path.exists(os.path.join(human_dir, "gene_table.parquet"))
    assert os.path.exists(os.path.join(human_dir, "symbol_lookup.parquet"))
    assert os.path.exists(os.path.join(human_dir, "build_metadata.json"))


def test_build_human_gene_table(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    ref = load_reference("human", reference_dir=ref_dir)
    gt = ref["gene_table"]

    assert "ensembl_id" in gt.columns
    assert "symbol" in gt.columns
    assert "source_id" in gt.columns
    tp53 = gt[gt["symbol"] == "TP53"]
    assert len(tp53) == 1
    assert tp53.iloc[0]["ensembl_id"] == "ENSG00000141510"
    assert tp53.iloc[0]["source_id"] == "HGNC:11998"


def test_build_human_symbol_lookup(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    ref = load_reference("human", reference_dir=ref_dir)
    sl = ref["symbol_lookup"]

    tp53_approved = sl[
        (sl["lookup_string"] == "TP53") & (sl["lookup_type"] == "approved_symbol")
    ]
    assert len(tp53_approved) == 1
    assert tp53_approved.iloc[0]["ensembl_id"] == "ENSG00000141510"

    lfs1_alias = sl[
        (sl["lookup_string"] == "LFS1") & (sl["lookup_type"] == "alias_symbol")
    ]
    assert len(lfs1_alias) == 1

    assert "lookup_string_upper" in sl.columns


def test_build_metadata(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)

    meta_path = os.path.join(ref_dir, "human", "build_metadata.json")
    with open(meta_path) as f:
        meta = json.load(f)
    assert "download_timestamp" in meta
    assert "sources" in meta


def test_build_skips_if_exists(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)
        build_reference("human", reference_dir=ref_dir)
        assert mock_dl.call_count == 1


def test_build_force_redownloads(ref_dir, mock_hgnc_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_hgnc_data.encode("utf-8")
        build_reference("human", reference_dir=ref_dir)
        build_reference("human", reference_dir=ref_dir, force=True)
        assert mock_dl.call_count >= 2


@pytest.fixture
def mock_mgi_markers_data():
    """Minimal MGI MRK_List2.rpt format."""
    return (
        "MGI Accession ID\tChr\tcM Position\tgenome coordinate start\tgenome coordinate end\tstrand\tMarker Symbol\tStatus\tMarker Name\tMarker Type\tFeature Type\tMarker Synonyms (pipe-separated)\n"
        "MGI:87853\t11\t\t69580359\t69591872\t+\tTrp53\tO\ttransformation related protein 53\tGene\tprotein coding gene\tp53|Tp53\n"
        "MGI:104738\t11\t\t101453964\t101517817\t+\tBrca1\tO\tbreast cancer 1, early onset\tGene\tprotein coding gene\tBrca1\n"
        "MGI:12345\t1\t\t1000\t2000\t+\tFakeGene\tO\tfake gene\tGene\tprotein coding gene\t\n"
    )


@pytest.fixture
def mock_mgi_ensembl_data():
    """Minimal MGI MRK_ENSEMBL.rpt format."""
    return (
        "MGI Marker Accession ID\tMarker Symbol\tMarker Name\tcM Position\tChromosome\tEnsembl Gene ID\n"
        "MGI:87853\tTrp53\ttransformation related protein 53\t\t11\tENSMUSG00000059552\n"
        "MGI:104738\tBrca1\tbreast cancer 1, early onset\t\t11\tENSMUSG00000017146\n"
    )


def test_build_mouse_creates_files(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    mouse_dir = os.path.join(ref_dir, "mouse")
    assert os.path.exists(os.path.join(mouse_dir, "gene_table.parquet"))
    assert os.path.exists(os.path.join(mouse_dir, "symbol_lookup.parquet"))


def test_build_mouse_gene_table(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    ref = load_reference("mouse", reference_dir=ref_dir)
    gt = ref["gene_table"]

    trp53 = gt[gt["symbol"] == "Trp53"]
    assert len(trp53) == 1
    assert trp53.iloc[0]["ensembl_id"] == "ENSMUSG00000059552"
    assert trp53.iloc[0]["source_id"] == "MGI:87853"

    fake = gt[gt["symbol"] == "FakeGene"]
    assert len(fake) == 1
    assert pd.isna(fake.iloc[0]["ensembl_id"])
    assert fake.iloc[0]["source_id"] == "MGI:12345"


def test_build_mouse_symbol_lookup(ref_dir, mock_mgi_markers_data, mock_mgi_ensembl_data):
    def mock_download(url):
        if "MRK_List2" in url:
            return mock_mgi_markers_data.encode("utf-8")
        elif "MRK_ENSEMBL" in url:
            return mock_mgi_ensembl_data.encode("utf-8")
        return b""

    with patch("stangene.references._download_file", side_effect=mock_download):
        build_reference("mouse", reference_dir=ref_dir)

    ref = load_reference("mouse", reference_dir=ref_dir)
    sl = ref["symbol_lookup"]

    trp53_approved = sl[(sl["lookup_string"] == "Trp53") & (sl["lookup_type"] == "approved_symbol")]
    assert len(trp53_approved) == 1

    p53_alias = sl[(sl["lookup_string"] == "p53") & (sl["lookup_type"] == "alias_symbol")]
    assert len(p53_alias) == 1


# ---------------------------------------------------------------------------
# Rat (RGD) tests
# ---------------------------------------------------------------------------

@pytest.fixture
def mock_rgd_data():
    """Minimal RGD GENES_RAT.txt format with comment header."""
    return (
        "# RGD-PIPELINE: ftp-file-extracts\n"
        "# MODULE: genes  build 2024-06-24\n"
        "#\n"
        "GENE_RGD_ID\tSYMBOL\tNAME\tGENE_DESC\tCHROMOSOME_CELERA\tCHROMOSOME_mRatBN7.2\tCHROMOSOME_RGSC_v3.4\tFISH_BAND\tSTART_POS_CELERA\tSTOP_POS_CELERA\tSTRAND_CELERA\tSTART_POS_mRatBN7.2\tSTOP_POS_mRatBN7.2\tSTRAND_mRatBN7.2\tSTART_POS_RGSC_v3.4\tSTOP_POS_RGSC_v3.4\tSTRAND_RGSC_v3.4\tCURATED_REF_RGD_ID\tCURATED_REF_PUBMED_ID\tUNCURATED_PUBMED_ID\tNCBI_GENE_ID\tUNIPROT_ID\tGENE_REFSEQ_STATUS\tGENBANK_NUCLEOTIDE\tTIGR_ID\tGENBANK_PROTEIN\tCANONICAL_PROTEIN\tMARKER_RGD_ID\tMARKER_SYMBOL\tOLD_SYMBOL\tOLD_NAME\tQTL_RGD_ID\tQTL_SYMBOL\tNOMENCLATURE_STATUS\tSPLICE_RGD_ID\tSPLICE_SYMBOL\tGENE_TYPE\tENSEMBL_ID\n"
        "2003\tAsip\tagouti signaling protein\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tA;ASP\t\t\t\tAPPROVED\t\t\tprotein-coding\tENSRNOG00000017701\n"
        "2004\tA2m\talpha-2-macroglobulin\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tA2MAC1;Mam\t\t\t\tPROVISIONAL\t\t\tprotein-coding\tENSRNOG00000028896;ENSRNOG00000045772\n"
        "9999\tNoEnsGene\tgene without ensembl\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAPPROVED\t\t\tprotein-coding\t\n"
    )


def test_build_rat_creates_files(ref_dir, mock_rgd_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_rgd_data.encode("utf-8")
        build_reference("rat", reference_dir=ref_dir)

    rat_dir = os.path.join(ref_dir, "rat")
    assert os.path.exists(os.path.join(rat_dir, "gene_table.parquet"))
    assert os.path.exists(os.path.join(rat_dir, "symbol_lookup.parquet"))
    assert os.path.exists(os.path.join(rat_dir, "build_metadata.json"))


def test_build_rat_gene_table(ref_dir, mock_rgd_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_rgd_data.encode("utf-8")
        build_reference("rat", reference_dir=ref_dir)

    ref = load_reference("rat", reference_dir=ref_dir)
    gt = ref["gene_table"]

    # Asip should have Ensembl ID (first ENSRNOG from the field)
    asip = gt[gt["symbol"] == "Asip"]
    assert len(asip) == 1
    assert asip.iloc[0]["ensembl_id"] == "ENSRNOG00000017701"
    assert asip.iloc[0]["source_id"] == "RGD:2003"
    assert asip.iloc[0]["status"] == "approved"

    # A2m has multiple Ensembl IDs — should pick the first ENSRNOG one
    a2m = gt[gt["symbol"] == "A2m"]
    assert len(a2m) == 1
    assert a2m.iloc[0]["ensembl_id"] == "ENSRNOG00000028896"
    assert a2m.iloc[0]["status"] == "provisional"

    # NoEnsGene has no Ensembl ID — ensembl_id should be null
    no_ens = gt[gt["symbol"] == "NoEnsGene"]
    assert len(no_ens) == 1
    assert pd.isna(no_ens.iloc[0]["ensembl_id"])
    assert no_ens.iloc[0]["source_id"] == "RGD:9999"


def test_build_rat_symbol_lookup(ref_dir, mock_rgd_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_rgd_data.encode("utf-8")
        build_reference("rat", reference_dir=ref_dir)

    ref = load_reference("rat", reference_dir=ref_dir)
    sl = ref["symbol_lookup"]

    # Approved symbol
    asip_approved = sl[(sl["lookup_string"] == "Asip") & (sl["lookup_type"] == "approved_symbol")]
    assert len(asip_approved) == 1

    # OLD_SYMBOL should be in prev_symbol (not alias_symbol)
    asp_prev = sl[(sl["lookup_string"] == "ASP") & (sl["lookup_type"] == "prev_symbol")]
    assert len(asp_prev) == 1
    assert asp_prev.iloc[0]["source_id"] == "RGD:2003"

    # A2m's old symbols should also be prev_symbol
    mam_prev = sl[(sl["lookup_string"] == "Mam") & (sl["lookup_type"] == "prev_symbol")]
    assert len(mam_prev) == 1


def test_build_rat_metadata(ref_dir, mock_rgd_data):
    with patch("stangene.references._download_file") as mock_dl:
        mock_dl.return_value = mock_rgd_data.encode("utf-8")
        build_reference("rat", reference_dir=ref_dir)

    meta_path = os.path.join(ref_dir, "rat", "build_metadata.json")
    with open(meta_path) as f:
        meta = json.load(f)
    assert meta["species"] == "rat"
    assert "rgd_genes" in meta["sources"]
