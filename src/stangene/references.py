"""Build and load species-specific gene annotation reference databases."""

import hashlib
import io
import json
import os
import urllib.request
from datetime import datetime, timezone

import pandas as pd

from stangene._logging import get_logger
from stangene.species import get_species_config

logger = get_logger("references")


class ReferenceNotFoundError(Exception):
    """Raised when reference data has not been built for a species."""
    pass


def _default_reference_dir() -> str:
    """Return the default reference directory (project-relative)."""
    return os.path.join(os.path.dirname(__file__), "..", "..", "references")


def _download_file(url: str) -> bytes:
    """Download a file from a URL and return its contents as bytes."""
    logger.info("Downloading %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "stangene/0.1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read()


def build_reference(
    species: str,
    reference_dir: str = None,
    force: bool = False,
) -> None:
    config = get_species_config(species)
    ref_dir = os.path.join(reference_dir or _default_reference_dir(), config.name)

    gene_table_path = os.path.join(ref_dir, "gene_table.parquet")
    if os.path.exists(gene_table_path) and not force:
        logger.info("References for %s already exist, skipping", species)
        return

    os.makedirs(ref_dir, exist_ok=True)

    if config.name == "human":
        _build_human_reference(config, ref_dir)
    elif config.name == "mouse":
        _build_mouse_reference(config, ref_dir)
    else:
        raise ValueError(f"No reference builder for species: {species}")

    logger.info("Reference build complete for %s at %s", species, ref_dir)


def load_reference(
    species: str,
    reference_dir: str = None,
) -> dict:
    config = get_species_config(species)
    ref_dir = os.path.join(reference_dir or _default_reference_dir(), config.name)

    gene_table_path = os.path.join(ref_dir, "gene_table.parquet")
    lookup_path = os.path.join(ref_dir, "symbol_lookup.parquet")
    meta_path = os.path.join(ref_dir, "build_metadata.json")

    if not os.path.exists(gene_table_path):
        raise ReferenceNotFoundError(
            f"Reference data for '{species}' not found at {ref_dir}. "
            f"Run stangene.references.build_reference('{species}') first."
        )

    gene_table = pd.read_parquet(gene_table_path)
    symbol_lookup = pd.read_parquet(lookup_path)
    with open(meta_path) as f:
        metadata = json.load(f)

    logger.info("Loaded %s reference: %d genes, %d lookup entries", species, len(gene_table), len(symbol_lookup))
    return {"gene_table": gene_table, "symbol_lookup": symbol_lookup, "metadata": metadata}


def _save_reference(ref_dir, gene_table, symbol_lookup, metadata):
    gene_table.to_parquet(os.path.join(ref_dir, "gene_table.parquet"), index=False)
    symbol_lookup.to_parquet(os.path.join(ref_dir, "symbol_lookup.parquet"), index=False)
    with open(os.path.join(ref_dir, "build_metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)


def _build_symbol_lookup(gene_table: pd.DataFrame, source: str) -> pd.DataFrame:
    rows = []
    for _, gene in gene_table.iterrows():
        eid = gene["ensembl_id"] if pd.notna(gene.get("ensembl_id")) else None
        sid = gene["source_id"]

        if pd.notna(gene["symbol"]) and gene["symbol"]:
            rows.append({
                "lookup_string": gene["symbol"],
                "lookup_string_upper": gene["symbol"].upper(),
                "ensembl_id": eid,
                "source_id": sid,
                "lookup_type": "approved_symbol",
                "source": source,
            })

        if pd.notna(gene.get("alias_symbols")) and gene["alias_symbols"]:
            for alias in str(gene["alias_symbols"]).split("|"):
                alias = alias.strip()
                if alias:
                    rows.append({
                        "lookup_string": alias,
                        "lookup_string_upper": alias.upper(),
                        "ensembl_id": eid,
                        "source_id": sid,
                        "lookup_type": "alias_symbol",
                        "source": source,
                    })

        if pd.notna(gene.get("prev_symbols")) and gene["prev_symbols"]:
            for prev in str(gene["prev_symbols"]).split("|"):
                prev = prev.strip()
                if prev:
                    rows.append({
                        "lookup_string": prev,
                        "lookup_string_upper": prev.upper(),
                        "ensembl_id": eid,
                        "source_id": sid,
                        "lookup_type": "prev_symbol",
                        "source": source,
                    })

    return pd.DataFrame(rows)


def _build_human_reference(config, ref_dir: str) -> None:
    url = config.reference_sources["hgnc"]["url"]
    raw_data = _download_file(url)
    checksum = hashlib.sha256(raw_data).hexdigest()

    hgnc = pd.read_csv(io.BytesIO(raw_data), sep="\t", low_memory=False)

    gene_table = pd.DataFrame({
        "ensembl_id": hgnc["ensembl_gene_id"].where(hgnc["ensembl_gene_id"].notna(), None),
        "symbol": hgnc["symbol"],
        "alias_symbols": hgnc["alias_symbol"].fillna(""),
        "prev_symbols": hgnc["prev_symbol"].fillna(""),
        "gene_type": hgnc.get("locus_group", pd.Series(dtype=str)).fillna(""),
        "status": hgnc["status"].fillna(""),
        "source": "HGNC",
        "source_id": hgnc["hgnc_id"],
    })

    symbol_lookup = _build_symbol_lookup(gene_table, source="HGNC")

    metadata = {
        "species": "human",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "hgnc": {"url": url, "sha256": checksum, "rows": len(hgnc)}
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built human reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


def _build_mouse_reference(config, ref_dir: str) -> None:
    """Placeholder — implemented in Task 6."""
    raise NotImplementedError("Mouse reference builder will be implemented in Task 6")
