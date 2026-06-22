"""Tests for stangene.hb_mask (species-aware hemoglobin detection)."""

from stangene import hb_mask


def _masked(symbols, species):
    m = hb_mask(symbols, species)
    return [s for s, keep in zip(symbols, m) if keep]


def test_human_curated_set_excludes_lookalikes():
    syms = ["HBA1", "HBB", "HBE1", "HBZ", "HBEGF", "HBP1", "HBS1L", "TP53"]
    # real hemoglobins match; HBEGF / HBP1 / HBS1L (non-hemoglobin) are excluded
    assert _masked(syms, "hs") == ["HBA1", "HBB", "HBE1", "HBZ"]


def test_primates_share_human_set():
    syms = ["HBA", "HBA1", "HBE1", "HBM", "PTGS1"]
    assert _masked(syms, "rhesus") == ["HBA", "HBA1", "HBE1", "HBM"]
    assert _masked(["HBZ", "HBE1", "ACTB"], "marmoset") == ["HBZ", "HBE1"]


def test_mouse_dashed_cluster_and_theta():
    syms = ["Hba-a1", "Hbb-bs", "Hbb-y", "Hbq1a", "Hba-ps3", "Actb"]
    # dashed alpha/beta cluster + theta; pseudogene Hba-ps3 excluded
    assert _masked(syms, "mm") == ["Hba-a1", "Hbb-bs", "Hbb-y", "Hbq1a"]


def test_rat_set():
    syms = ["Hba-a1", "Hbb-b1", "Hbe1", "Hbz", "Actb"]
    assert _masked(syms, "rn") == ["Hba-a1", "Hbb-b1", "Hbe1", "Hbz"]


def test_zebrafish_isoforms():
    syms = ["hbaa1", "hbba1", "hbae1.1", "hbbe2", "actb1"]
    assert _masked(syms, "zebrafish") == ["hbaa1", "hbba1", "hbae1.1", "hbbe2"]


def test_invertebrates_empty():
    # insects/nematodes have no RBC hemoglobin; fly "hb" is hunchback, not hemoglobin
    assert _masked(["hb", "glob1", "glob2"], "dm") == []
    assert _masked(["glb-1", "glb-2"], "ce") == []
    # human curated names must NOT leak into invertebrates
    assert _masked(["HBA1", "HBB"], "dm") == []


def test_empty_and_dtype():
    m = hb_mask([], "human")
    assert m.dtype == bool and len(m) == 0
    m2 = hb_mask(["HBB", "TP53"], "human")
    assert m2.dtype == bool and list(m2) == [True, False]
