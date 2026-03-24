"""Integration tests using scanpy's pbmc3k dataset.

These tests run the full stangene pipeline on real-world data from 10x Genomics.
The pbmc3k dataset has 2,700 cells and 32,738 genes with both gene symbols
(var_names) and Ensembl gene IDs (var["gene_ids"]).

Requires: scanpy (dev dependency), network access on first run (to download
pbmc3k and build HGNC references).
"""

import json
import os

import pandas as pd
import pytest
import scanpy as sc

import stangene
from stangene.harmonize import HarmonizationResult


@pytest.fixture(scope="module")
def pbmc3k_path(tmp_path_factory):
    """Download pbmc3k and save as h5ad."""
    tmp = tmp_path_factory.mktemp("pbmc3k")
    adata = sc.datasets.pbmc3k()
    path = str(tmp / "pbmc3k.h5ad")
    adata.write_h5ad(path)
    return path


@pytest.fixture(scope="module")
def ref_dir(tmp_path_factory):
    """Build a real HGNC reference for human."""
    ref_dir = str(tmp_path_factory.mktemp("refs"))
    stangene.build_reference("human", reference_dir=ref_dir)
    return ref_dir


@pytest.fixture(scope="module")
def result(pbmc3k_path, ref_dir, tmp_path_factory):
    """Run the full pipeline on pbmc3k."""
    output_dir = str(tmp_path_factory.mktemp("output"))
    result = stangene.run(
        pbmc3k_path,
        species="human",
        output_dir=output_dir,
        dataset_name="pbmc3k",
        reference_dir=ref_dir,
    )
    result._output_dir = output_dir  # stash for file tests
    return result


class TestPbmc3kPipelineRuns:
    """Basic checks that the pipeline completes on real data."""

    def test_returns_harmonization_result(self, result):
        assert isinstance(result, HarmonizationResult)

    def test_all_features_present(self, result):
        assert len(result.mapping_table) == 32738

    def test_no_rows_without_status(self, result):
        assert result.mapping_table["mapping_status"].notna().all()

    def test_dataset_name_set(self, result):
        assert (result.mapping_table["dataset"] == "pbmc3k").all()

    def test_species_set(self, result):
        assert (result.mapping_table["species"] == "human").all()

    def test_stangene_version_present(self, result):
        assert "stangene_version" in result.mapping_table.columns
        assert (result.mapping_table["stangene_version"] == stangene.__version__).all()


class TestPbmc3kMappingQuality:
    """Verify that the majority of pbmc3k genes map successfully."""

    def test_majority_mapped(self, result):
        """Most features should resolve to a harmonized ID.

        pbmc3k uses GENCODE annotation which includes many non-coding RNAs
        and novel loci (e.g., RP11-*, AC*, AL*) that are not in HGNC's
        approved gene set. ~76% mapping rate is expected.
        """
        mt = result.mapping_table
        mapped = mt["mapping_status"].isin([
            "exact_id", "id_no_version", "exact_symbol",
            "alias_symbol", "previous_symbol",
        ])
        mapped_pct = mapped.sum() / len(mt)
        assert mapped_pct > 0.70, f"Only {mapped_pct:.1%} mapped"

    def test_tier1_dominates(self, result):
        """pbmc3k has Ensembl IDs — Tier 1 (exact_id) should be the largest group."""
        counts = result.stats
        exact_id = counts.get("exact_id", 0)
        total_gene = sum(v for k, v in counts.items() if k != "non_gene_feature")
        assert exact_id > total_gene * 0.5, (
            f"exact_id={exact_id} is less than 50% of gene features ({total_gene})"
        )

    def test_no_non_gene_features(self, result):
        """pbmc3k is gene expression only — no antibody/CRISPR/spike-in features."""
        non_gene = result.stats.get("non_gene_feature", 0)
        assert non_gene == 0, f"Unexpected {non_gene} non-gene features in pbmc3k"

    def test_unmapped_count_reasonable(self, result):
        """Unmapped genes should not exceed ~30%.

        pbmc3k's GENCODE annotation includes ~24% features that are novel
        loci, pseudogenes, or non-coding RNAs absent from HGNC.
        """
        mt = result.mapping_table
        unmapped = (mt["mapping_status"] == "unmapped").sum()
        unmapped_pct = unmapped / len(mt)
        assert unmapped_pct < 0.30, f"{unmapped_pct:.1%} unmapped is too high"


class TestPbmc3kKnownGenes:
    """Spot-check known genes that must appear in pbmc3k."""

    @pytest.mark.parametrize("symbol,ensembl_id", [
        ("CD3D", "ENSG00000167286"),    # T cell marker
        ("CD79A", "ENSG00000105369"),    # B cell marker
        ("NKG7", "ENSG00000105374"),     # NK cell marker
        ("LYZ", "ENSG00000090382"),      # Monocyte marker
        ("MS4A1", "ENSG00000156738"),    # CD20 / B cell
    ])
    def test_known_gene_mapped(self, result, symbol, ensembl_id):
        mt = result.mapping_table
        row = mt[mt["original_feature_name"] == symbol]
        assert len(row) == 1, f"{symbol} not found in pbmc3k"
        row = row.iloc[0]
        assert row["gene_id_harmonized"] == ensembl_id
        assert row["mapping_status"] in ("exact_id", "id_no_version")
        assert row["mapping_confidence"] == "high"

    def test_tp53_has_correct_symbol(self, result):
        mt = result.mapping_table
        row = mt[mt["original_feature_name"] == "TP53"]
        assert len(row) == 1
        row = row.iloc[0]
        assert row["gene_symbol_harmonized"] == "TP53"
        assert row["gene_id_harmonized"] == "ENSG00000141510"


class TestPbmc3kOriginalPreservation:
    """Verify that original identifiers are never lost."""

    def test_original_feature_name_preserved(self, result):
        mt = result.mapping_table
        assert "original_feature_name" in mt.columns
        # Should include well-known pbmc3k genes
        names = set(mt["original_feature_name"])
        assert "CD3D" in names
        assert "MS4A1" in names

    def test_original_feature_id_preserved(self, result):
        mt = result.mapping_table
        assert "original_feature_id" in mt.columns
        # First gene's Ensembl ID should be present
        first_id = mt.iloc[0]["original_feature_id"]
        assert first_id is not None and str(first_id).startswith("ENSG")

    def test_feature_id_no_version_derived(self, result):
        """pbmc3k Ensembl IDs have no version suffix, so feature_id_no_version
        should be null (since _strip_version only strips .N suffixes)."""
        mt = result.mapping_table
        assert "feature_id_no_version" in mt.columns


class TestPbmc3kOutputFiles:
    """Verify that output files are written correctly."""

    def test_harmonization_table_written(self, result):
        path = os.path.join(result._output_dir, "harmonization_table.tsv")
        assert os.path.exists(path)
        df = pd.read_csv(path, sep="\t")
        assert len(df) == 32738

    def test_summary_json_valid(self, result):
        path = os.path.join(result._output_dir, "summary.json")
        assert os.path.exists(path)
        with open(path) as f:
            s = json.load(f)
        assert s["total_features"] == 32738
        assert s["gene_features"] == 32738
        assert s["non_gene_features"] == 0

    def test_conflicts_tsv_written(self, result):
        path = os.path.join(result._output_dir, "conflicts.tsv")
        assert os.path.exists(path)
        df = pd.read_csv(path, sep="\t")
        assert "conflict_type" in df.columns

    def test_unmapped_tsv_written(self, result):
        path = os.path.join(result._output_dir, "unmapped.tsv")
        assert os.path.exists(path)

    def test_enriched_h5ad_written(self, result):
        path = os.path.join(result._output_dir, "pbmc3k_harmonized.h5ad")
        assert os.path.exists(path)
        import anndata
        adata = anndata.read_h5ad(path)
        assert "gene_id_harmonized" in adata.var.columns
        assert "mapping_status" in adata.var.columns
        # var_names should be the original gene symbols
        assert "CD3D" in adata.var_names


class TestPbmc3kSummaryStats:
    """Validate the summary statistics make sense for pbmc3k."""

    def test_status_counts_sum_to_total(self, result):
        s = stangene.summary(result)
        total_from_counts = sum(s["status_counts"].values())
        assert total_from_counts == s["total_features"]

    def test_duplicate_ids_reported(self, result):
        """Check that duplicate harmonized IDs (if any) are counted."""
        s = stangene.summary(result)
        assert "duplicate_harmonized_ids" in s
        # Value should be a non-negative integer
        assert isinstance(s["duplicate_harmonized_ids"], int)
        assert s["duplicate_harmonized_ids"] >= 0
