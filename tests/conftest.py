"""Shared pytest fixtures for stangene tests."""

import pandas as pd
import pytest


@pytest.fixture
def mock_gene_table():
    """Minimal gene_table DataFrame mimicking HGNC reference."""
    return pd.DataFrame([
        {"ensembl_id": "ENSG00000141510", "symbol": "TP53", "alias_symbols": "LFS1|p53",
         "prev_symbols": "LFS1", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:11998", "canonical_biotype": "protein_coding"},
        {"ensembl_id": "ENSG00000012048", "symbol": "BRCA1", "alias_symbols": "RNF53|IRIS",
         "prev_symbols": "RNF53", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:1100", "canonical_biotype": "protein_coding"},
        {"ensembl_id": "ENSG00000139618", "symbol": "BRCA2", "alias_symbols": "FACD|FANCD1",
         "prev_symbols": "FACD", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:1101", "canonical_biotype": "protein_coding"},
        {"ensembl_id": "ENSG00000136997", "symbol": "MYC", "alias_symbols": "",
         "prev_symbols": "", "gene_type": "protein-coding gene", "status": "Approved",
         "source": "HGNC", "source_id": "HGNC:7553", "canonical_biotype": "protein_coding"},
        {"ensembl_id": None, "symbol": "WITHDRAWN1", "alias_symbols": "",
         "prev_symbols": "", "gene_type": "pseudogene", "status": "Entry Withdrawn",
         "source": "HGNC", "source_id": "HGNC:99999", "canonical_biotype": "pseudogene"},
    ])


@pytest.fixture
def mock_symbol_lookup():
    """Minimal symbol_lookup DataFrame."""
    rows = [
        {"lookup_string": "TP53", "lookup_string_upper": "TP53", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "LFS1", "lookup_string_upper": "LFS1", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "p53", "lookup_string_upper": "P53", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "LFS1", "lookup_string_upper": "LFS1", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "prev_symbol", "source": "HGNC"},
        {"lookup_string": "BRCA1", "lookup_string_upper": "BRCA1", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "RNF53", "lookup_string_upper": "RNF53", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "RNF53", "lookup_string_upper": "RNF53", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "prev_symbol", "source": "HGNC"},
        {"lookup_string": "BRCA2", "lookup_string_upper": "BRCA2", "ensembl_id": "ENSG00000139618",
         "source_id": "HGNC:1101", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "FACD", "lookup_string_upper": "FACD", "ensembl_id": "ENSG00000139618",
         "source_id": "HGNC:1101", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "MYC", "lookup_string_upper": "MYC", "ensembl_id": "ENSG00000136997",
         "source_id": "HGNC:7553", "lookup_type": "approved_symbol", "source": "HGNC"},
        {"lookup_string": "AMBIG", "lookup_string_upper": "AMBIG", "ensembl_id": "ENSG00000141510",
         "source_id": "HGNC:11998", "lookup_type": "alias_symbol", "source": "HGNC"},
        {"lookup_string": "AMBIG", "lookup_string_upper": "AMBIG", "ensembl_id": "ENSG00000012048",
         "source_id": "HGNC:1100", "lookup_type": "alias_symbol", "source": "HGNC"},
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def mock_ref(mock_gene_table, mock_symbol_lookup):
    """Complete mock reference dict."""
    return {
        "gene_table": mock_gene_table,
        "symbol_lookup": mock_symbol_lookup,
        "metadata": {"species": "human", "download_timestamp": "2026-01-01T00:00:00Z"},
    }
