"""Species-aware hemoglobin-gene detection on (harmonized) gene symbols.

Why per-species (and why explicit sets, not a prefix): hemoglobin symbols differ
sharply across species, and a naive ``HB`` prefix over-captures non-hemoglobin
genes. So each species uses an explicit, curated set of approved hemoglobin
symbols, verified against the bundled references:
  - human + primates (cyno/rhesus/marmoset/mouse_lemur) : HGNC-style symbols
    (HBA1, HBB, HBE1, ...). The reference also carries HBEGF / HBP1 / HBS1L as
    separate, non-hemoglobin genes — deliberately NOT in the set.
  - mouse (MGI) : dashed alpha/beta cluster (Hba-a1, Hbb-bs, ...) + theta Hbq1a/b
  - rat (RGD)   : dashed cluster + Hbe1/Hbe2/Hbz/Hbq1
  - zebrafish (ZFIN) : hbaa*/hbba*/hbae*/hbbe* isoforms
  - fruit_fly / c_elegans : NOT APPLICABLE — insects and nematodes have no
    red-blood-cell hemoglobin, so the mask is always empty. (Drosophila ``hb``
    is hunchback, ``glob1-3`` are intracellular globins — neither is a blood
    contamination marker, so they must never match.)

Hemoglobin fraction is a red-blood-cell / lysis contamination QC signal, only
meaningful for vertebrates. Pseudogenes (e.g. mouse ``Hba-ps*``) are omitted.
"""

from __future__ import annotations

import numpy as np

from stangene.species import resolve_species

# HGNC-style symbols shared by human and the curated primate references. Bare
# "HBA" appears as an approved symbol in the rhesus reference. HBEGF/HBP1/HBS1L
# are separate non-hemoglobin genes and are deliberately absent.
_HB_PRIMATE = frozenset({
    "HBA", "HBA1", "HBA2", "HBB", "HBD",
    "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ",
})
# Mouse: dashed alpha/beta cluster (incl. embryonic Hba-x / Hbb-bh* / Hbb-y) +
# theta Hbq1a/b. Pseudogenes (Hba-ps*) excluded.
_HB_MOUSE = frozenset({
    "HBA-A1", "HBA-A2", "HBA-X",
    "HBB-B1", "HBB-B2", "HBB-BH0", "HBB-BH1", "HBB-BH2", "HBB-BH3",
    "HBB-BS", "HBB-BT", "HBB-Y", "HBQ1A", "HBQ1B",
})
_HB_RAT = frozenset({
    "HBA-A1", "HBA-A2", "HBA-A3",
    "HBB-B1", "HBB-B2", "HBB-BS", "HBB-BT",
    "HBE1", "HBE2", "HBQ1", "HBQ1B", "HBZ",
})
_HB_ZEBRAFISH = frozenset({
    "HBAA1", "HBAA2", "HBAE1.1", "HBAE1.2", "HBAE1.3", "HBAE3", "HBAE4", "HBAE5",
    "HBBA1", "HBBA2", "HBBE1.1", "HBBE1.2", "HBBE1.3", "HBBE2", "HBBE3",
})

_PRIMATE_SPECIES = frozenset({
    "human", "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
})

# Species -> curated set. Species absent from this map (fruit_fly, c_elegans)
# have no RBC hemoglobin and get an all-False mask.
_HB_BY_SPECIES = {
    **{sp: _HB_PRIMATE for sp in _PRIMATE_SPECIES},
    "mouse": _HB_MOUSE,
    "rat": _HB_RAT,
    "zebrafish": _HB_ZEBRAFISH,
}


def hb_mask(symbols, species: str) -> np.ndarray:
    """Boolean mask over ``symbols`` marking hemoglobin genes for a species.

    ``species`` accepts a short code or full name (resolved via
    :func:`resolve_species`). ``symbols`` is any iterable of gene symbols
    (ideally already harmonized to canonical symbols). Returns a numpy bool
    array aligned to ``symbols``. For species with no red-blood-cell hemoglobin
    (fruit fly, C. elegans) the mask is all-False.
    """
    canon = resolve_species(species)
    up = [str(s).upper() for s in symbols]
    mask = np.zeros(len(up), dtype=bool)
    if not up:
        return mask
    hb_set = _HB_BY_SPECIES.get(canon)
    if hb_set is None:
        return mask
    return np.array([s in hb_set for s in up], dtype=bool)
