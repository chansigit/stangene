"""Stangene: gene identifier harmonization for single-cell transcriptomics."""

__version__ = "0.1.0"

from stangene.classify import classify_features
from stangene.harmonize import HarmonizationResult, harmonize
from stangene.io import load_features
from stangene.merge import MergeResult, merge_features
from stangene.references import (
    ReferenceNotFoundError,
    build_reference,
    load_reference,
)
from stangene.report import conflict_report, generate_markdown_report, summary, write_reports
from stangene._logging import get_logger

_logger = get_logger("run")


def run(
    path: str,
    species: str,
    output_dir: str = None,
    dataset_name: str = None,
    reference_dir: str = None,
) -> HarmonizationResult:
    """Run the full harmonization pipeline on a single dataset."""
    _logger.info("Starting harmonization: path=%s, species=%s", path, species)

    ft = load_features(path, species=species, dataset_name=dataset_name)
    ft = classify_features(ft)
    ref = load_reference(species, reference_dir=reference_dir)
    result = harmonize(ft, ref)

    if output_dir:
        from stangene.io import write_results
        write_results(result, output_dir, input_path=path)
        write_reports(result, output_dir)

    _logger.info("Harmonization complete: %s", result.stats)
    return result
