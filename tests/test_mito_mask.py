"""Tests for stangene.mito_mask (species-aware mitochondrial detection)."""

import numpy as np

from stangene import mito_mask


def _masked(symbols, species):
    m = mito_mask(symbols, species)
    return [s for s, keep in zip(symbols, m) if keep]


def test_human_mt_dash_prefix_only():
    syms = ["MT-CO1", "MT-CYB", "TP53", "MT1A", "COX1"]
    # MT- prefix only; MT1A (metallothionein) and bare COX1 (=PTGS1 alias) excluded
    assert _masked(syms, "hs") == ["MT-CO1", "MT-CYB"]


def test_mouse_and_rat_uppercased_prefix():
    assert _masked(["mt-Co1", "mt-Nd1", "Actb"], "mm") == ["mt-Co1", "mt-Nd1"]
    assert _masked(["Mt-co1", "Mt-cyb", "Actb"], "rn") == ["Mt-co1", "Mt-cyb"]


def test_primate_bare_names():
    syms = ["ND1", "COX1", "CYTB", "ATP6", "PTGS1", "TP53"]
    # bare mito names match; PTGS1 (separate nuclear gene) does NOT
    assert _masked(syms, "rhesus") == ["ND1", "COX1", "CYTB", "ATP6"]
    # primates also accept an MT- prefix (Ensembl-style data)
    assert _masked(["MT-CO1", "ND1"], "cynomolgus") == ["MT-CO1", "ND1"]


def test_fruit_fly_colon_prefix():
    syms = ["mt:CoI", "mt:Cyt-b", "mt:ND1", "Act5C"]
    assert _masked(syms, "dm") == ["mt:CoI", "mt:Cyt-b", "mt:ND1"]
    # human MT- rule must NOT apply to fly
    assert _masked(["MT-CO1"], "dm") == []


def test_c_elegans_named_set():
    syms = ["nduo-1", "ctc-1", "ctb-1", "atp-6", "unc-22"]
    assert _masked(syms, "ce") == ["nduo-1", "ctc-1", "ctb-1", "atp-6"]


def test_empty_and_dtype():
    m = mito_mask([], "human")
    assert m.dtype == bool and len(m) == 0
    m2 = mito_mask(["MT-CO1", "TP53"], "human")
    assert m2.dtype == bool and list(m2) == [True, False]
