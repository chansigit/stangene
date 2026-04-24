"""Build and load species-specific gene annotation reference databases."""

import hashlib
import io
import json
import os
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone

import pandas as pd

from stangene._logging import get_logger
from stangene.species import get_species_config

logger = get_logger("references")


class ReferenceNotFoundError(Exception):
    """Raised when reference data has not been built for a species."""
    pass


def _default_reference_dir() -> str:
    """Return the default reference directory.

    Uses project-relative 'references/' if it exists (dev mode),
    otherwise falls back to ~/.cache/stangene/references.
    """
    project_dir = os.path.join(os.path.dirname(__file__), "..", "..", "references")
    if os.path.isdir(project_dir):
        return project_dir
    cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "stangene", "references")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


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
    elif config.name == "rat":
        _build_rat_reference(config, ref_dir)
    elif config.name == "zebrafish":
        _build_zebrafish_reference(config, ref_dir)
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

    missing = [p for p in [gene_table_path, lookup_path, meta_path] if not os.path.exists(p)]
    if missing:
        raise ReferenceNotFoundError(
            f"Reference data for '{species}' incomplete at {ref_dir}. "
            f"Missing: {[os.path.basename(p) for p in missing]}. "
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
    """Build mouse reference from MGI marker files + Ensembl mapping."""
    markers_url = config.reference_sources["mgi_markers"]["url"]
    ensembl_url = config.reference_sources["mgi_ensembl"]["url"]

    markers_raw = _download_file(markers_url)
    ensembl_raw = _download_file(ensembl_url)

    markers_checksum = hashlib.sha256(markers_raw).hexdigest()
    ensembl_checksum = hashlib.sha256(ensembl_raw).hexdigest()

    markers = pd.read_csv(io.BytesIO(markers_raw), sep="\t", low_memory=False)
    markers.columns = markers.columns.str.strip()

    ensembl_map = pd.read_csv(io.BytesIO(ensembl_raw), sep="\t", low_memory=False)
    ensembl_map.columns = ensembl_map.columns.str.strip()

    # Build MGI ID -> Ensembl ID lookup
    mgi_to_ensembl = {}
    eid_col = [c for c in ensembl_map.columns if "ensembl" in c.lower() and "id" in c.lower()]
    mid_col = [c for c in ensembl_map.columns if "mgi" in c.lower() and "accession" in c.lower()]
    if eid_col and mid_col:
        for _, row in ensembl_map.iterrows():
            mgi_id = row[mid_col[0]]
            ens_id = row[eid_col[0]]
            if pd.notna(mgi_id) and pd.notna(ens_id) and str(ens_id).startswith("ENSMUSG"):
                mgi_to_ensembl[str(mgi_id)] = str(ens_id)
    else:
        logger.warning("MGI-Ensembl mapping columns not found. Expected 'Ensembl Gene ID' and 'MGI Marker Accession ID'. Found: %s", list(ensembl_map.columns))

    def _find_col(df, pattern, label):
        matches = [c for c in df.columns if pattern in c.lower()]
        if not matches:
            raise ValueError(f"Expected column containing '{label}' not found in {list(df.columns)}")
        return matches[0]

    # Build gene table
    sym_col = _find_col(markers, "marker symbol", "Marker Symbol")
    status_col = _find_col(markers, "status", "Status")
    type_col = _find_col(markers, "feature type", "Feature Type")
    mgi_col = _find_col(markers, "mgi accession", "MGI Accession")
    syn_col = [c for c in markers.columns if "synonym" in c.lower()][0]

    rows = []
    for _, m in markers.iterrows():
        mgi_id = str(m[mgi_col]) if pd.notna(m[mgi_col]) else ""
        symbol = str(m[sym_col]) if pd.notna(m[sym_col]) else ""
        ensembl_id = mgi_to_ensembl.get(mgi_id)
        synonyms = str(m[syn_col]) if pd.notna(m[syn_col]) else ""

        status = "approved" if str(m[status_col]).strip() == "O" else "withdrawn"
        gene_type = str(m[type_col]) if pd.notna(m[type_col]) else ""

        rows.append({
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "alias_symbols": synonyms,
            "prev_symbols": "",
            "gene_type": gene_type,
            "status": status,
            "source": "MGI",
            "source_id": mgi_id,
        })

    gene_table = pd.DataFrame(rows)

    # Supplementary BioMart fill (non-fatal if it fails)
    biomart_url = config.reference_sources.get("ensembl_biomart", {}).get("url", "")
    if biomart_url:
        try:
            biomart_xml = (
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<!DOCTYPE Query>'
                '<Query virtualSchemaName="default" formatter="TSV" header="1">'
                '<Dataset name="mmusculus_gene_ensembl" interface="default">'
                '<Attribute name="ensembl_gene_id"/>'
                '<Attribute name="external_gene_name"/>'
                '<Attribute name="mgi_id"/>'
                '</Dataset></Query>'
            )
            import urllib.parse
            full_url = biomart_url + urllib.parse.quote(biomart_xml)
            biomart_raw = _download_file(full_url)
            biomart_df = pd.read_csv(io.BytesIO(biomart_raw), sep="\t", low_memory=False)
            biomart_df.columns = ["ensembl_id_bm", "symbol_bm", "mgi_id_bm"]

            null_mask = gene_table["ensembl_id"].isna()
            if null_mask.any():
                bm_by_mgi = biomart_df.dropna(subset=["mgi_id_bm"]).drop_duplicates(subset=["mgi_id_bm"])
                bm_mgi_map = dict(zip(bm_by_mgi["mgi_id_bm"], bm_by_mgi["ensembl_id_bm"]))
                for idx in gene_table[null_mask].index:
                    sid = gene_table.at[idx, "source_id"]
                    if sid in bm_mgi_map:
                        gene_table.at[idx, "ensembl_id"] = bm_mgi_map[sid]

                still_null = gene_table["ensembl_id"].isna()
                if still_null.any():
                    bm_by_sym = biomart_df.dropna(subset=["symbol_bm"]).drop_duplicates(subset=["symbol_bm"])
                    bm_sym_map = dict(zip(bm_by_sym["symbol_bm"], bm_by_sym["ensembl_id_bm"]))
                    for idx in gene_table[still_null].index:
                        sym = gene_table.at[idx, "symbol"]
                        if sym in bm_sym_map:
                            gene_table.at[idx, "ensembl_id"] = bm_sym_map[sym]

                filled = null_mask.sum() - gene_table["ensembl_id"].isna().sum()
                logger.info("BioMart supplementary fill: %d/%d gaps filled", filled, null_mask.sum())
        except Exception as e:
            logger.warning("BioMart supplementary download failed (non-fatal): %s", e)

    symbol_lookup = _build_symbol_lookup(gene_table, source="MGI")

    metadata = {
        "species": "mouse",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "mgi_markers": {"url": markers_url, "sha256": markers_checksum, "rows": len(markers)},
            "mgi_ensembl": {"url": ensembl_url, "sha256": ensembl_checksum, "rows": len(ensembl_map)},
            "ensembl_biomart": {"url": biomart_url, "note": "supplementary Ensembl ID gap fill"},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built mouse reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


def _build_rat_reference(config, ref_dir: str) -> None:
    """Build rat reference from RGD (Rat Genome Database) gene file."""
    url = config.reference_sources["rgd_genes"]["url"]
    raw_data = _download_file(url)
    checksum = hashlib.sha256(raw_data).hexdigest()

    # RGD file has comment lines starting with '#'; actual header is the first non-comment line
    lines = raw_data.decode("utf-8", errors="replace").split("\n")
    data_start = next(i for i, l in enumerate(lines) if not l.startswith("#") and l.strip())
    rgd = pd.read_csv(
        io.StringIO("\n".join(lines[data_start:])),
        sep="\t",
        low_memory=False,
    )
    rgd.columns = rgd.columns.str.strip()

    # Map RGD nomenclature status to internal status
    _RGD_STATUS_MAP = {
        "APPROVED": "approved",
        "PROVISIONAL": "provisional",
        "INTERIM": "provisional",
    }

    rows = []
    for _, gene in rgd.iterrows():
        symbol = str(gene["SYMBOL"]).strip() if pd.notna(gene.get("SYMBOL")) else ""
        if not symbol:
            continue

        # RGD ENSEMBL_ID column may contain multiple IDs separated by ';'
        # Use the first canonical ENSRNOG ID
        ensembl_id = None
        raw_ens = str(gene.get("ENSEMBL_ID", "")).strip() if pd.notna(gene.get("ENSEMBL_ID")) else ""
        if raw_ens:
            candidates = [e.strip() for e in raw_ens.split(";") if e.strip().startswith("ENSRNOG")]
            if candidates:
                ensembl_id = candidates[0]

        gene_type = str(gene.get("GENE_TYPE", "")).strip() if pd.notna(gene.get("GENE_TYPE")) else ""

        # OLD_SYMBOL contains previous names (';'-separated) — map to prev_symbols
        old_symbols = str(gene.get("OLD_SYMBOL", "")).strip() if pd.notna(gene.get("OLD_SYMBOL")) else ""
        prev_symbols = "|".join(
            s.strip() for s in old_symbols.split(";") if s.strip()
        ) if old_symbols else ""

        # Map nomenclature status
        raw_status = str(gene.get("NOMENCLATURE_STATUS", "")).strip().upper() if pd.notna(gene.get("NOMENCLATURE_STATUS")) else ""
        status = _RGD_STATUS_MAP.get(raw_status, "provisional")

        # Build source_id with RGD: prefix
        rgd_id = ""
        if pd.notna(gene.get("GENE_RGD_ID")):
            try:
                rgd_id = f"RGD:{int(gene['GENE_RGD_ID'])}"
            except (ValueError, TypeError):
                rgd_id = f"RGD:{gene['GENE_RGD_ID']}"

        rows.append({
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "alias_symbols": "",
            "prev_symbols": prev_symbols,
            "gene_type": gene_type,
            "status": status,
            "source": "RGD",
            "source_id": rgd_id,
        })

    gene_table = pd.DataFrame(rows)
    symbol_lookup = _build_symbol_lookup(gene_table, source="RGD")

    metadata = {
        "species": "rat",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "rgd_genes": {"url": url, "sha256": checksum, "rows": len(rgd)},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built rat reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


def _build_zebrafish_reference(config, ref_dir: str) -> None:
    """Build zebrafish reference from ZFIN (Zebrafish Information Network)."""
    genes_url = config.reference_sources["zfin_genes"]["url"]
    aliases_url = config.reference_sources["zfin_aliases"]["url"]
    ensembl_url = config.reference_sources["zfin_ensembl"]["url"]

    genes_raw = _download_file(genes_url)
    aliases_raw = _download_file(aliases_url)
    ensembl_raw = _download_file(ensembl_url)

    genes_checksum = hashlib.sha256(genes_raw).hexdigest()
    aliases_checksum = hashlib.sha256(aliases_raw).hexdigest()
    ensembl_checksum = hashlib.sha256(ensembl_raw).hexdigest()

    # ZFIN files are headerless TSVs
    genes_df = pd.read_csv(
        io.BytesIO(genes_raw), sep="\t", header=None,
        names=["zfin_id", "so_id", "symbol", "ensembl_id"],
        low_memory=False, dtype=str,
    )

    # Supplementary Ensembl mapping (used only if genes file lacks it)
    ensembl_df = pd.read_csv(
        io.BytesIO(ensembl_raw), sep="\t", header=None,
        names=["zfin_id", "symbol_em", "ensembl_id_em"],
        low_memory=False, dtype=str,
    )
    ensembl_map = dict(zip(ensembl_df["zfin_id"], ensembl_df["ensembl_id_em"]))

    # Aliases file: zfin_id, current_symbol, alias_string, alias_type
    aliases_df = pd.read_csv(
        io.BytesIO(aliases_raw), sep="\t", header=None,
        names=["zfin_id", "current_symbol", "alias_string", "alias_type"],
        low_memory=False, dtype=str,
    )

    # Group aliases by zfin_id — separate PREVIOUS NAME vs ALIAS
    prev_by_id = defaultdict(list)
    alias_by_id = defaultdict(list)
    for _, a in aliases_df.iterrows():
        zid = a["zfin_id"]
        alias_str = str(a["alias_string"]).strip()
        atype = str(a["alias_type"]).strip().upper()
        if not zid or not alias_str:
            continue
        if "PREVIOUS" in atype:
            prev_by_id[zid].append(alias_str)
        else:
            alias_by_id[zid].append(alias_str)

    rows = []
    for _, g in genes_df.iterrows():
        zid = str(g["zfin_id"]).strip()
        symbol = str(g["symbol"]).strip() if pd.notna(g["symbol"]) else ""
        if not symbol or not zid:
            continue

        ensembl_id = str(g["ensembl_id"]).strip() if pd.notna(g["ensembl_id"]) and str(g["ensembl_id"]).strip() else None
        if not ensembl_id:
            ensembl_id = ensembl_map.get(zid)
        if ensembl_id and not ensembl_id.startswith("ENSDARG"):
            ensembl_id = None

        alias_symbols = "|".join(alias_by_id.get(zid, []))
        prev_symbols = "|".join(prev_by_id.get(zid, []))

        rows.append({
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "alias_symbols": alias_symbols,
            "prev_symbols": prev_symbols,
            "gene_type": "",
            "status": "approved",
            "source": "ZFIN",
            "source_id": f"ZFIN:{zid}",
        })

    gene_table = pd.DataFrame(rows)
    symbol_lookup = _build_symbol_lookup(gene_table, source="ZFIN")

    metadata = {
        "species": "zebrafish",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "zfin_genes": {"url": genes_url, "sha256": genes_checksum, "rows": len(genes_df)},
            "zfin_aliases": {"url": aliases_url, "sha256": aliases_checksum, "rows": len(aliases_df)},
            "zfin_ensembl": {"url": ensembl_url, "sha256": ensembl_checksum, "rows": len(ensembl_df)},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built zebrafish reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))
