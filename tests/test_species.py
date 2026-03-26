from stangene.species import SpeciesConfig, get_species_config, CLASSIFICATION_PATTERNS


def test_get_human_config():
    cfg = get_species_config("human")
    assert cfg.name == "human"
    assert cfg.ensembl_prefix == "ENSG"
    assert cfg.transcript_prefix == "ENST"
    assert cfg.naming_convention == "uppercase"


def test_get_mouse_config():
    cfg = get_species_config("mouse")
    assert cfg.name == "mouse"
    assert cfg.ensembl_prefix == "ENSMUSG"
    assert cfg.transcript_prefix == "ENSMUST"
    assert cfg.naming_convention == "capitalized"


def test_get_rat_config():
    cfg = get_species_config("rat")
    assert cfg.name == "rat"
    assert cfg.ensembl_prefix == "ENSRNOG"
    assert cfg.transcript_prefix == "ENSRNOT"
    assert cfg.naming_convention == "capitalized"


def test_unknown_species_raises():
    import pytest
    with pytest.raises(ValueError, match="Unknown species"):
        get_species_config("zebrafish")


def test_classification_patterns_detect_antibody():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "CD3_TotalSeqB"):
            assert ftype == "antibody_capture"
            return
    raise AssertionError("No pattern matched CD3_TotalSeqB")


def test_classification_patterns_detect_spike_in():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "ERCC-00002"):
            assert ftype == "spike_in"
            return
    raise AssertionError("No pattern matched ERCC-00002")


def test_classification_patterns_detect_peak():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "chr1:1000-2000"):
            assert ftype == "peak"
            return
    raise AssertionError("No pattern matched chr1:1000-2000")


def test_classification_patterns_detect_crispr():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "sg-TP53"):
            assert ftype == "crispr_guide"
            return
    raise AssertionError("No pattern matched sg-TP53")
