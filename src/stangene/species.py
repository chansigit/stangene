"""Species-specific configuration and classification patterns."""

import re
from dataclasses import dataclass, field


@dataclass
class SpeciesConfig:
    """Configuration for a species' gene naming and reference sources."""

    name: str
    ensembl_prefix: str
    transcript_prefix: str
    naming_convention: str  # "uppercase" (human) or "capitalized" (mouse)
    reference_sources: dict = field(default_factory=dict)


# Classification patterns: list of (compiled_regex_pattern, feature_type).
# Order matters — first match wins. These are checked against original_feature_name.
CLASSIFICATION_PATTERNS: list[tuple[re.Pattern, str]] = [
    # Transcript IDs (check before gene IDs since ENST is a subset prefix-wise)
    (re.compile(r"^ENST\d+", re.IGNORECASE), "transcript"),
    (re.compile(r"^ENSMUST\d+", re.IGNORECASE), "transcript"),
    (re.compile(r"^ENSRNOT\d+", re.IGNORECASE), "transcript"),
    # Ensembl gene IDs
    (re.compile(r"^ENSG\d+", re.IGNORECASE), "gene"),
    (re.compile(r"^ENSMUSG\d+", re.IGNORECASE), "gene"),
    (re.compile(r"^ENSRNOG\d+", re.IGNORECASE), "gene"),
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
}


def get_species_config(species: str) -> SpeciesConfig:
    """Get configuration for a species. Raises ValueError if unknown."""
    species_lower = species.lower()
    if species_lower not in _SPECIES_CONFIGS:
        raise ValueError(
            f"Unknown species: '{species}'. Supported: {list(_SPECIES_CONFIGS.keys())}"
        )
    return _SPECIES_CONFIGS[species_lower]
