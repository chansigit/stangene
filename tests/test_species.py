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
        get_species_config("platypus")


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


def test_get_zebrafish_config():
    cfg = get_species_config("zebrafish")
    assert cfg.name == "zebrafish"
    assert cfg.ensembl_prefix == "ENSDARG"
    assert cfg.transcript_prefix == "ENSDART"
    assert "zfin_genes" in cfg.reference_sources


def test_get_fruitfly_config():
    cfg = get_species_config("fruit_fly")
    assert cfg.name == "fruit_fly"
    assert cfg.ensembl_prefix == "FBgn"
    assert cfg.transcript_prefix == "FBtr"
    assert "flybase_gene_map" in cfg.reference_sources


def test_get_celegans_config():
    cfg = get_species_config("c_elegans")
    assert cfg.name == "c_elegans"
    assert cfg.ensembl_prefix == "WBGene"
    assert "wormbase_gene_ids" in cfg.reference_sources


def test_get_cynomolgus_config():
    cfg = get_species_config("cynomolgus")
    assert cfg.name == "cynomolgus"
    assert cfg.ensembl_prefix == "ENSMFAG"
    assert cfg.transcript_prefix == "ENSMFAT"
    assert cfg.reference_sources["ensembl_biomart"]["dataset"] == "mfascicularis_gene_ensembl"


def test_get_rhesus_config():
    cfg = get_species_config("rhesus")
    assert cfg.name == "rhesus"
    assert cfg.ensembl_prefix == "ENSMMUG"
    assert cfg.reference_sources["ensembl_biomart"]["dataset"] == "mmulatta_gene_ensembl"


def test_get_marmoset_config():
    cfg = get_species_config("marmoset")
    assert cfg.name == "marmoset"
    assert cfg.ensembl_prefix == "ENSCJAG"
    assert cfg.reference_sources["ensembl_biomart"]["dataset"] == "cjacchus_gene_ensembl"


def test_get_mouse_lemur_config():
    cfg = get_species_config("mouse_lemur")
    assert cfg.name == "mouse_lemur"
    assert cfg.ensembl_prefix == "ENSMICG"
    assert cfg.reference_sources["ensembl_biomart"]["dataset"] == "mmurinus_gene_ensembl"


def test_classification_patterns_detect_zebrafish_gene():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "ENSDARG00000012345"):
            assert ftype == "gene"
            return
    raise AssertionError("No pattern matched ENSDARG00000012345")


def test_classification_patterns_detect_macaque_gene():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "ENSMFAG00000056789"):
            assert ftype == "gene"
            return
    raise AssertionError("No pattern matched ENSMFAG00000056789")


def test_classification_patterns_detect_flybase_gene():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "FBgn0000008"):
            assert ftype == "gene"
            return
    raise AssertionError("No pattern matched FBgn0000008")


def test_classification_patterns_detect_flybase_transcript():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "FBtr0078163"):
            assert ftype == "transcript"
            return
    raise AssertionError("No pattern matched FBtr0078163")


def test_classification_patterns_detect_wormbase_gene():
    import re
    for pattern, ftype in CLASSIFICATION_PATTERNS:
        if re.match(pattern, "WBGene00000001"):
            assert ftype == "gene"
            return
    raise AssertionError("No pattern matched WBGene00000001")


def test_classification_patterns_transcripts_not_genes():
    """ENST/ENSMUST/ENSDART should match transcript, not gene."""
    import re
    for feature_id in ["ENST00000269305", "ENSMUST00000017146", "ENSDART00000123456"]:
        for pattern, ftype in CLASSIFICATION_PATTERNS:
            if re.match(pattern, feature_id):
                assert ftype == "transcript", f"{feature_id} matched {ftype}, expected transcript"
                break
