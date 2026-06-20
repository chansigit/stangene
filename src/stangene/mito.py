"""Species-aware mitochondrial-gene detection on (harmonized) gene symbols.

Why per-species: the human ``MT-`` prefix does NOT generalize. Verified against
the bundled references:
  - human / mouse / rat : approved symbols are MT-/mt-/Mt- prefixed (MT- after upper)
  - cynomolgus / rhesus / marmoset / mouse_lemur : bare names (ND1, COX1, CYTB, ...)
  - fruit_fly (FlyBase) : ``mt:`` prefix (e.g. mt:CoI, mt:Cyt-b)
  - c_elegans (WormBase) : idiosyncratic (nduo-*, ctc-*, ctb-1, atp-6)
  - zebrafish : ZFIN uses an ``mt-`` prefix in data (MT- after upper)

Bare primate names and the worm set were checked to be collision-free with
nuclear genes in their references (e.g. bare ``COX1`` is an *alias* of PTGS1 in
human, so it is deliberately NOT used for mammals — only as an approved symbol
in primates, where PTGS1 is a separate entry).
"""

from __future__ import annotations

import numpy as np

from stangene.species import resolve_species

# Mito-specific "MT-" prefix holds across mammals, primates, and zebrafish
# (nuclear metallothioneins are MT1A/MT2A — no hyphen — so they don't match).
_MT_DASH_SPECIES = frozenset({
    "human", "mouse", "rat", "zebrafish",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
})
# Primates additionally encode mito genes as bare approved symbols.
_BARE_SPECIES = frozenset({"cynomolgus", "rhesus", "marmoset", "mouse_lemur"})
_BARE_MITO = frozenset({
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8",
})
_FLY_SPECIES = frozenset({"fruit_fly"})
_WORM_SPECIES = frozenset({"c_elegans"})
_WORM_MITO = frozenset({
    "ATP-6", "CTB-1", "CTC-1", "CTC-2", "CTC-3",
    "NDUO-1", "NDUO-2", "NDUO-3", "NDUO-4", "NDUO-5", "NDUO-6",
})


def mito_mask(symbols, species: str) -> np.ndarray:
    """Boolean mask over ``symbols`` marking mitochondrial genes for a species.

    ``species`` accepts a short code or full name (resolved via
    :func:`resolve_species`). ``symbols`` is any iterable of gene symbols
    (ideally already harmonized to canonical symbols). Returns a numpy bool
    array aligned to ``symbols``.
    """
    canon = resolve_species(species)
    up = [str(s).upper() for s in symbols]
    mask = np.zeros(len(up), dtype=bool)
    if not up:
        return mask
    if canon in _MT_DASH_SPECIES:
        mask |= np.array([s.startswith("MT-") for s in up], dtype=bool)
    if canon in _BARE_SPECIES:
        mask |= np.array([s in _BARE_MITO for s in up], dtype=bool)
    if canon in _FLY_SPECIES:
        mask |= np.array([s.startswith("MT:") for s in up], dtype=bool)
    if canon in _WORM_SPECIES:
        mask |= np.array([s in _WORM_MITO for s in up], dtype=bool)
    return mask
