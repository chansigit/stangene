import os
from unittest.mock import patch

import anndata
import numpy as np
import pandas as pd
import pytest

import stangene
from stangene.harmonize import HarmonizationResult
from stangene.references import ReferenceNotFoundError


@pytest.fixture
def sample_h5ad(tmp_path):
    adata = anndata.AnnData(
        X=np.zeros((3, 4)),
        var=pd.DataFrame(
            {"gene_ids": ["ENSG00000141510", "ENSG00000012048", "", ""]},
            index=["TP53", "BRCA1", "ERCC-00002", "NONEXISTENT"],
        ),
    )
    path = str(tmp_path / "test.h5ad")
    adata.write_h5ad(path)
    return path


@pytest.fixture
def ref_dir_with_human(tmp_path, mock_ref):
    ref_dir = str(tmp_path / "references")
    human_dir = os.path.join(ref_dir, "human")
    os.makedirs(human_dir)
    mock_ref["gene_table"].to_parquet(os.path.join(human_dir, "gene_table.parquet"))
    mock_ref["symbol_lookup"].to_parquet(os.path.join(human_dir, "symbol_lookup.parquet"))
    import json
    with open(os.path.join(human_dir, "build_metadata.json"), "w") as f:
        json.dump(mock_ref["metadata"], f)
    return ref_dir


def test_run_returns_result(sample_h5ad, ref_dir_with_human):
    result = stangene.run(sample_h5ad, species="human", reference_dir=ref_dir_with_human)
    assert isinstance(result, HarmonizationResult)
    assert len(result.mapping_table) == 4


def test_run_writes_output(sample_h5ad, ref_dir_with_human, tmp_path):
    output_dir = str(tmp_path / "output")
    result = stangene.run(
        sample_h5ad, species="human",
        output_dir=output_dir, reference_dir=ref_dir_with_human,
    )
    assert os.path.exists(os.path.join(output_dir, "harmonization_table.tsv"))
    assert os.path.exists(os.path.join(output_dir, "summary.json"))


def test_run_raises_without_refs(sample_h5ad, tmp_path):
    with pytest.raises(ReferenceNotFoundError):
        stangene.run(sample_h5ad, species="human", reference_dir=str(tmp_path / "empty"))


def test_run_mapping_statuses(sample_h5ad, ref_dir_with_human):
    result = stangene.run(sample_h5ad, species="human", reference_dir=ref_dir_with_human)
    statuses = set(result.mapping_table["mapping_status"])
    assert "exact_id" in statuses
    assert "non_gene_feature" in statuses
    assert "unmapped" in statuses
