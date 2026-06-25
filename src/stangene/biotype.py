"""Canonical gene biotype vocabulary and per-source normalization."""

from __future__ import annotations

CANONICAL_BIOTYPES: frozenset[str] = frozenset({
    "protein_coding",
    "lncRNA",
    "pseudogene",
    "miRNA",
    "snoRNA",
    "snRNA",
    "rRNA",
    "tRNA",
    "IG_gene",
    "TR_gene",
    "other_ncrna",
    "other",
    "unknown",
})

# HGNC: locus_type field (more granular than locus_group)
_HGNC_MAP: dict[str, str] = {
    "gene with protein product": "protein_coding",
    "rna, micro":                "miRNA",
    "rna, small nucleolar":      "snoRNA",
    "rna, small nuclear":        "snRNA",
    "rna, ribosomal":            "rRNA",
    "rna, transfer":             "tRNA",
    "rna, long non-coding":      "lncRNA",
    "rna, y":                    "other_ncrna",
    "rna, vault":                "other_ncrna",
    "rna, small cytoplasmic":    "other_ncrna",
    "immunoglobulin gene":       "IG_gene",
    "t cell receptor gene":      "TR_gene",
    "pseudogene":                "pseudogene",
    "endogenous retrovirus":     "other",
    "complex locus constituent": "other",
    "fragile site":              "other",
    "readthrough":               "other",
    "region":                    "other",
    "unknown":                   "unknown",
}

# MGI: feature_type field
_MGI_MAP: dict[str, str] = {
    "protein coding gene":               "protein_coding",
    "lncrna gene":                       "lncRNA",
    "pseudogene":                        "pseudogene",
    "mirna gene":                        "miRNA",
    "snorna gene":                       "snoRNA",
    "snrna gene":                        "snRNA",
    "rrna gene":                         "rRNA",
    "trna gene":                         "tRNA",
    "ig gene":                           "IG_gene",
    "tr gene":                           "TR_gene",
    "unclassified non-coding rna gene":  "other_ncrna",
    "pirna gene":                        "other_ncrna",
    "scrna gene":                        "other_ncrna",
    "rnase mrp rna gene":                "other_ncrna",
    "rnase p rna gene":                  "other_ncrna",
    "srp rna gene":                      "other_ncrna",
    "telomerase rna gene":               "other_ncrna",
    "antisense lncrna gene":             "lncRNA",
    "bidirectional promoter lncrna":     "lncRNA",
    "transgene":                         "other",
    "unclassified gene":                 "other",
    "gene segment":                      "other",
    "other feature types":               "other",
}

# RGD: GENE_TYPE field (lowercase in data)
_RGD_MAP: dict[str, str] = {
    "protein-coding":  "protein_coding",
    "ncrna":           "other_ncrna",
    "lincrna":         "lncRNA",
    "lncrna":          "lncRNA",
    "mirna":           "miRNA",
    "snorna":          "snoRNA",
    "snrna":           "snRNA",
    "rrna":            "rRNA",
    "trna":            "tRNA",
    "pseudo":          "pseudogene",
    "pseudogene":      "pseudogene",
    "gene":            "other",
}

# Ensembl BioMart: gene_biotype field
# Pseudogene subtypes all end in "pseudogene" — handled by suffix fallback below.
# IG/TR subtypes all match IG_*_gene / TR_*_gene — handled by prefix fallback.
_ENSEMBL_MAP: dict[str, str] = {
    "protein_coding":   "protein_coding",
    "lincrna":          "lncRNA",
    "lncrna":           "lncRNA",
    "mirna":            "miRNA",
    "snorna":           "snoRNA",
    "snrna":            "snRNA",
    "rrna":             "rRNA",
    "mt_rrna":          "rRNA",
    "trna":             "tRNA",
    "mt_trna":          "tRNA",
    "pseudogene":       "pseudogene",
    "misc_rna":         "other_ncrna",
    "ribozyme":         "other_ncrna",
    "vault_rna":        "other_ncrna",
    "scrna":            "other_ncrna",
    "pirna":            "other_ncrna",
    "scarna":           "other_ncrna",
    "known_ncrna":      "other_ncrna",
    "processed_transcript": "other_ncrna",
    "tec":              "other",
    "artifact":         "other",
    "disrupted_domain": "other",
}

# WormBase: biotype field
_WORMBASE_MAP: dict[str, str] = {
    "protein_coding_gene":       "protein_coding",
    "pseudogene":                "pseudogene",
    "lncrna_gene":               "lncRNA",
    "mirna_gene":                "miRNA",
    "snorna_gene":               "snoRNA",
    "snrna_gene":                "snRNA",
    "rrna_gene":                 "rRNA",
    "trna_gene":                 "tRNA",
    "ncrna_gene":                "other_ncrna",
    "pirna_gene":                "other_ncrna",
    "scrna_gene":                "other_ncrna",
    "transposable_element_gene": "other",
    "operon":                    "other",
    "gene":                      "other",        # generic unclassified WormBase gene
    "lincrna_gene":              "lncRNA",       # long intergenic ncRNA
    "antisense_lncrna_gene":     "lncRNA",       # antisense lncRNA
    "transcript":                "other",        # data artifact
}

# ZFIN: SO (Sequence Ontology) ID, e.g. "SO:0001217"
_ZFIN_SO_MAP: dict[str, str] = {
    "so:0001217": "protein_coding",   # protein_coding_gene
    "so:0000336": "pseudogene",       # pseudogene
    "so:0001263": "miRNA",            # miRNA_gene
    "so:0001267": "snoRNA",           # snoRNA_gene
    "so:0001268": "snRNA",            # snRNA_gene
    "so:0001637": "rRNA",             # rRNA_gene
    "so:0001272": "tRNA",             # tRNA_gene
    "so:0001877": "lncRNA",           # lncRNA_gene
    "so:0001900": "other_ncrna",      # ncRNA_gene
    "so:0001265": "other_ncrna",      # miRNA_primary_transcript (pre-miRNA)
    "so:0001032": "other_ncrna",      # other ncRNA
    "so:0000704": "other",            # gene (unspecified)
    "so:0000233": "other",            # processed_transcript
}

_SOURCE_MAPS: dict[str, dict[str, str]] = {
    "HGNC":      _HGNC_MAP,
    "MGI":       _MGI_MAP,
    "RGD":       _RGD_MAP,
    "Ensembl":   _ENSEMBL_MAP,
    "WormBase":  _WORMBASE_MAP,
    "ZFIN":      _ZFIN_SO_MAP,
}


def normalize_biotype(raw: str, source: str) -> str:
    """Map a raw gene_type string to a canonical biotype.

    Args:
        raw: The raw gene_type string from the upstream database.
        source: One of "HGNC", "MGI", "RGD", "Ensembl", "WormBase",
                "ZFIN", or "FlyBase".

    Returns:
        A member of CANONICAL_BIOTYPES. Unmapped values return "unknown".
    """
    key = raw.strip().lower()
    if not key:
        return "unknown"

    # FlyBase uses annotation_id prefix convention (no dict lookup needed)
    if source == "FlyBase":
        if key.startswith("cg"):
            return "protein_coding"
        if key.startswith("cr"):
            return "other_ncrna"
        return "other"

    source_map = _SOURCE_MAPS.get(source)
    if source_map is None:
        return "unknown"

    result = source_map.get(key)
    if result is not None:
        return result

    # Ensembl fallbacks: many pseudogene subtypes end in "pseudogene";
    # IG/TR gene subtypes start with ig_ or tr_ and end with _gene.
    if source == "Ensembl":
        if key.endswith("pseudogene"):
            return "pseudogene"
        if key.startswith("ig_") and key.endswith("_gene"):
            return "IG_gene"
        if key.startswith("tr_") and key.endswith("_gene"):
            return "TR_gene"

    return "unknown"
