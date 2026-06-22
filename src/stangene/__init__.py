"""Stangene: gene identifier harmonization for single-cell transcriptomics."""

__version__ = "0.3.0"

import pandas as pd

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
from stangene.species import resolve_species
from stangene.mito import mito_mask
from stangene.hb import hb_mask
from stangene._logging import get_logger

_logger = get_logger("run")


def harmonize_anndata(adata, species: str, *, reference_dir: str = None) -> HarmonizationResult:
    """Harmonize an in-memory AnnData's features without touching the matrix.

    Builds a FeatureTable from ``adata.var`` (mirroring the h5ad loader:
    ``original_feature_name`` from var_names, plus ``gene_ids``/``feature_types``
    columns if present), then runs classify -> load_reference -> harmonize. The
    returned ``HarmonizationResult.mapping_table`` rows are positionally aligned
    with ``adata.var_names``. Does NOT modify ``adata`` or write any file —
    applying the mapping (e.g. renaming var_names) is the caller's policy.
    """
    from stangene.io import _infer_reference_source, _strip_version

    var = adata.var
    ft = pd.DataFrame({"original_feature_name": list(var.index)})
    if "gene_ids" in var.columns:
        ft["original_feature_id"] = var["gene_ids"].values
    else:
        # No dedicated ID column: var_names may themselves be Ensembl/stable IDs.
        # Feed them as the ID too so ID tiers can match (a real symbol is not in
        # the Ensembl set, so it simply falls through to the symbol tier — no
        # false ID match), while original_feature_name still drives symbol matching.
        ft["original_feature_id"] = list(var.index)
    if "feature_types" in var.columns:
        ft["original_feature_type"] = var["feature_types"].values
    ft = ft.reset_index(drop=True)

    # Mirror load_features' column normalization (harmonize() requires these).
    ft["species"] = species
    ft["dataset"] = "anndata"
    if "original_feature_id" in ft.columns:
        ft["feature_id_no_version"] = ft["original_feature_id"].apply(_strip_version)
        ft.loc[ft["feature_id_no_version"] == "", "feature_id_no_version"] = None
    else:
        ft["original_feature_id"] = None
        ft["feature_id_no_version"] = None
    ft["original_feature_id"] = ft["original_feature_id"].replace("", None)
    ft["reference_source"] = _infer_reference_source(ft.get("original_feature_id"))
    ft["reference_release"] = None

    ft = classify_features(ft)
    ref = load_reference(species, reference_dir=reference_dir)
    return harmonize(ft, ref)


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
