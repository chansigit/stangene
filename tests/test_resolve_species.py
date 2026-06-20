"""Tests for stangene.resolve_species (short code / full name -> canonical)."""

import pytest

from stangene import resolve_species

CANONICAL = {
    "human", "mouse", "rat", "zebrafish", "fruit_fly", "c_elegans",
    "cynomolgus", "rhesus", "marmoset", "mouse_lemur",
}


@pytest.mark.parametrize("code,expected", [
    ("hs", "human"), ("human", "human"), ("HS", "human"),
    ("mm", "mouse"), ("mouse", "mouse"),
    ("rn", "rat"), ("rat", "rat"),
    ("dr", "zebrafish"), ("zebrafish", "zebrafish"),
    ("dm", "fruit_fly"), ("fruit_fly", "fruit_fly"),
    ("ce", "c_elegans"), ("c_elegans", "c_elegans"),
    ("cyno", "cynomolgus"), ("cynomolgus", "cynomolgus"),
    ("rhesus", "rhesus"),
    ("marmoset", "marmoset"),
    ("lemur", "mouse_lemur"), ("mouse_lemur", "mouse_lemur"),
    ("  Mouse_Lemur  ", "mouse_lemur"),  # whitespace + case
])
def test_resolve_known_codes(code, expected):
    assert resolve_species(code) == expected


def test_every_canonical_species_resolvable():
    # All 10 stangene species are reachable by their full name.
    for sp in CANONICAL:
        assert resolve_species(sp) == sp


def test_unknown_code_raises():
    with pytest.raises(ValueError) as e:
        resolve_species("xx")
    assert "Supported codes" in str(e.value)
