"""Feature type classification for gene vs non-gene triage."""

import pandas as pd

from stangene._logging import get_logger
from stangene.species import CLASSIFICATION_PATTERNS

logger = get_logger("classify")

# Map 10x Cell Ranger feature_type labels to our internal types
_CELLRANGER_TYPE_MAP = {
    "Gene Expression": "gene",
    "Antibody Capture": "antibody_capture",
    "CRISPR Guide Capture": "crispr_guide",
    "Custom": "custom",
    "Peaks": "peak",
}

_NON_GENE_TYPES = frozenset([
    "transcript", "antibody_capture", "crispr_guide",
    "spike_in", "peak", "custom",
])


def classify_features(ft: pd.DataFrame) -> pd.DataFrame:
    """Classify features as gene or non-gene types.

    Adds/updates 'original_feature_type' and sets 'mapping_status' to
    'non_gene_feature' for non-gene rows. Returns a copy.
    """
    result = ft.copy()

    # Ensure columns exist
    if "original_feature_type" not in result.columns:
        result["original_feature_type"] = None
    if "mapping_status" not in result.columns:
        result["mapping_status"] = None
    if "mapping_notes" not in result.columns:
        result["mapping_notes"] = None

    for idx in result.index:
        existing_type = result.at[idx, "original_feature_type"]

        if pd.notna(existing_type) and existing_type in _CELLRANGER_TYPE_MAP:
            result.at[idx, "original_feature_type"] = _CELLRANGER_TYPE_MAP[existing_type]
        elif pd.notna(existing_type) and existing_type in _NON_GENE_TYPES | {"gene"}:
            pass
        else:
            name = result.at[idx, "original_feature_name"]
            matched = False
            for pattern, ftype in CLASSIFICATION_PATTERNS:
                if pattern.match(str(name)):
                    result.at[idx, "original_feature_type"] = ftype
                    matched = True
                    break
            if not matched:
                result.at[idx, "original_feature_type"] = "gene"
                result.at[idx, "mapping_notes"] = "classified as gene by default (no pattern match)"

        if result.at[idx, "original_feature_type"] in _NON_GENE_TYPES:
            result.at[idx, "mapping_status"] = "non_gene_feature"

    classified_counts = result["original_feature_type"].value_counts().to_dict()
    logger.info("Feature classification: %s", classified_counts)
    return result
