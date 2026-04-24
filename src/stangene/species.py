"""Species-specific configuration and classification patterns."""

import re
from dataclasses import dataclass, field


@dataclass
class SpeciesConfig:
    """Configuration for a species' gene naming and reference sources."""

    name: str
    ensembl_prefix: str
    transcript_prefix: str
    naming_convention: str  # "uppercase", "capitalized", or "lowercase"
    reference_sources: dict = field(default_factory=dict)


# Classification patterns: list of (compiled_regex_pattern, feature_type).
# Order matters — first match wins. These are checked against original_feature_name.
CLASSIFICATION_PATTERNS: list[tuple[re.Pattern, str]] = [
    # Ensembl transcript IDs (check before gene IDs since both start with ENS)
    # Matches ENST, ENSMUST, ENSRNOT, ENSDART, ENSMFAT, ENSMMUT, ENSCJAT, ENSMICT, etc.
    (re.compile(r"^ENS[A-Z]*T\d+", re.IGNORECASE), "transcript"),
    # Ensembl gene IDs
    # Matches ENSG, ENSMUSG, ENSRNOG, ENSDARG, ENSMFAG, ENSMMUG, ENSCJAG, ENSMICG, etc.
    (re.compile(r"^ENS[A-Z]*G\d+", re.IGNORECASE), "gene"),
    # FlyBase IDs (transcripts first)
    (re.compile(r"^FBtr\d+", re.IGNORECASE), "transcript"),
    (re.compile(r"^FBgn\d+", re.IGNORECASE), "gene"),
    # WormBase IDs
    (re.compile(r"^WBGene\d+", re.IGNORECASE), "gene"),
    # Antibody capture / protein tags
    (re.compile(r".*TotalSeq", re.IGNORECASE), "antibody_capture"),
    (re.compile(r".*_ADT$", re.IGNORECASE), "antibody_capture"),
    (re.compile(r".*_HTO$", re.IGNORECASE), "antibody_capture"),
    # CRISPR guides
    (re.compile(r"^sg-", re.IGNORECASE), "crispr_guide"),
    (re.compile(r"^gRNA-", re.IGNORECASE), "crispr_guide"),
    # Spike-ins
    (re.compile(r"^ERCC-", re.IGNORECASE), "spike_in"),
    # Genomic peaks (ATAC-seq style)
    (re.compile(r"^chr[\dXYMT]+:\d+-\d+$", re.IGNORECASE), "peak"),
]

# Excel-corruption date-like patterns that indicate corrupted gene names
EXCEL_DATE_PATTERN = re.compile(
    r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
    r"|^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{1,2}$"
    r"|^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d{2,4}$",
    re.IGNORECASE,
)

# Gene symbols renamed by HGNC in 2020 due to Excel auto-conversion.
# Maps old symbol -> new symbol for flagging purposes.
EXCEL_RENAMED_GENES: dict[str, str] = {
    "MARCH1": "MARCHF1", "MARCH2": "MARCHF2", "MARCH3": "MARCHF3",
    "MARCH4": "MARCHF4", "MARCH5": "MARCHF5", "MARCH6": "MARCHF6",
    "MARCH7": "MARCHF7", "MARCH8": "MARCHF8", "MARCH9": "MARCHF9",
    "MARCH10": "MARCHF10", "MARCH11": "MARCHF11",
    "SEPT1": "SEPTIN1", "SEPT2": "SEPTIN2", "SEPT3": "SEPTIN3",
    "SEPT4": "SEPTIN4", "SEPT5": "SEPTIN5", "SEPT6": "SEPTIN6",
    "SEPT7": "SEPTIN7", "SEPT8": "SEPTIN8", "SEPT9": "SEPTIN9",
    "SEPT10": "SEPTIN10", "SEPT11": "SEPTIN11", "SEPT12": "SEPTIN12",
    "SEPT14": "SEPTIN14",
    "DEC1": "DELEC1", "DEC2": "BHLHE41",
}


_SPECIES_CONFIGS: dict[str, SpeciesConfig] = {
    "human": SpeciesConfig(
        name="human",
        ensembl_prefix="ENSG",
        transcript_prefix="ENST",
        naming_convention="uppercase",
        reference_sources={
            "hgnc": {
                "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                "description": "HGNC complete gene set with symbols, aliases, and Ensembl IDs",
            },
        },
    ),
    "mouse": SpeciesConfig(
        name="mouse",
        ensembl_prefix="ENSMUSG",
        transcript_prefix="ENSMUST",
        naming_convention="capitalized",
        reference_sources={
            "mgi_markers": {
                "url": "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt",
                "description": "MGI marker list with approved symbols and synonyms",
            },
            "mgi_ensembl": {
                "url": "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt",
                "description": "MGI to Ensembl ID mapping",
            },
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "description": "Ensembl BioMart mouse gene table (supplementary)",
            },
        },
    ),
    "rat": SpeciesConfig(
        name="rat",
        ensembl_prefix="ENSRNOG",
        transcript_prefix="ENSRNOT",
        naming_convention="capitalized",
        reference_sources={
            "rgd_genes": {
                "url": "https://download.rgd.mcw.edu/data_release/RAT/GENES_RAT.txt",
                "description": "RGD rat gene file with approved symbols, synonyms, and Ensembl IDs",
            },
        },
    ),
    "zebrafish": SpeciesConfig(
        name="zebrafish",
        ensembl_prefix="ENSDARG",
        transcript_prefix="ENSDART",
        naming_convention="capitalized",
        reference_sources={
            "zfin_genes": {
                "url": "https://zfin.org/downloads/gene.txt",
                "description": "ZFIN zebrafish gene list with Ensembl IDs",
            },
            "zfin_aliases": {
                "url": "https://zfin.org/downloads/aliases.txt",
                "description": "ZFIN gene aliases and previous symbols",
            },
            "zfin_ensembl": {
                "url": "https://zfin.org/downloads/ensembl_1_to_1.txt",
                "description": "ZFIN to Ensembl 1:1 mapping",
            },
        },
    ),
    "fruit_fly": SpeciesConfig(
        name="fruit_fly",
        ensembl_prefix="FBgn",
        transcript_prefix="FBtr",
        naming_convention="capitalized",
        reference_sources={
            "flybase_gene_map": {
                "url": "https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID.tsv.gz",
                "description": "FlyBase FBgn to annotation ID mapping with current symbols",
            },
            "flybase_synonyms": {
                # FlyBase synonyms file has no "current"-named symlink; only a
                # versioned filename exists under /releases/current/. Users will
                # need to bump the release suffix when FlyBase publishes a new
                # release (e.g. fb_2026_02, fb_2026_03, ...).
                "url": "https://s3ftp.flybase.org/releases/current/precomputed_files/synonyms/fb_synonym_fb_2026_01.tsv.gz",
                "description": "FlyBase synonyms (aliases + secondary FBgns = previous IDs)",
            },
        },
    ),
    "c_elegans": SpeciesConfig(
        name="c_elegans",
        ensembl_prefix="WBGene",
        transcript_prefix="",  # C. elegans transcripts have varied naming; no single prefix
        naming_convention="lowercase",
        reference_sources={
            "wormbase_gene_ids": {
                "url": "https://downloads.wormbase.org/releases/current-production-release/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.current.geneIDs.txt.gz",
                "description": "WormBase C. elegans gene IDs with public names and biotypes",
            },
            "wormbase_other_ids": {
                "url": "https://downloads.wormbase.org/releases/current-production-release/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.current.geneOtherIDs.txt.gz",
                "description": "WormBase other IDs (aliases, previous symbols)",
            },
        },
    ),
    "cynomolgus": SpeciesConfig(
        name="cynomolgus",
        ensembl_prefix="ENSMFAG",
        transcript_prefix="ENSMFAT",
        naming_convention="uppercase",
        reference_sources={
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "dataset": "mfascicularis_gene_ensembl",
                "description": "Ensembl BioMart Macaca fascicularis gene table",
            },
        },
    ),
    "rhesus": SpeciesConfig(
        name="rhesus",
        ensembl_prefix="ENSMMUG",
        transcript_prefix="ENSMMUT",
        naming_convention="uppercase",
        reference_sources={
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "dataset": "mmulatta_gene_ensembl",
                "description": "Ensembl BioMart Macaca mulatta gene table",
            },
        },
    ),
    "marmoset": SpeciesConfig(
        name="marmoset",
        ensembl_prefix="ENSCJAG",
        transcript_prefix="ENSCJAT",
        naming_convention="uppercase",
        reference_sources={
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "dataset": "cjacchus_gene_ensembl",
                "description": "Ensembl BioMart Callithrix jacchus gene table",
            },
        },
    ),
    "mouse_lemur": SpeciesConfig(
        name="mouse_lemur",
        ensembl_prefix="ENSMICG",
        transcript_prefix="ENSMICT",
        naming_convention="uppercase",
        reference_sources={
            "ensembl_biomart": {
                "url": "https://www.ensembl.org/biomart/martservice?query=",
                "dataset": "mmurinus_gene_ensembl",
                "description": "Ensembl BioMart Microcebus murinus gene table",
            },
        },
    ),
}


def get_species_config(species: str) -> SpeciesConfig:
    """Get configuration for a species. Raises ValueError if unknown."""
    species_lower = species.lower()
    if species_lower not in _SPECIES_CONFIGS:
        raise ValueError(
            f"Unknown species: '{species}'. Supported: {list(_SPECIES_CONFIGS.keys())}"
        )
    return _SPECIES_CONFIGS[species_lower]
