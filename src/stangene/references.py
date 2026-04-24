"""Build and load species-specific gene annotation reference databases."""

import gzip
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
    # BioMart queries for large datasets can take 3-5 minutes
    with urllib.request.urlopen(req, timeout=300) as resp:
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
    elif config.name == "fruit_fly":
        _build_fruitfly_reference(config, ref_dir)
    elif config.name == "c_elegans":
        _build_celegans_reference(config, ref_dir)
    elif config.name in ("cynomolgus", "rhesus", "marmoset", "mouse_lemur"):
        _build_ensembl_biomart_reference(config, ref_dir)
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
    ensembl_map = {
        k: v for k, v in zip(ensembl_df["zfin_id"], ensembl_df["ensembl_id_em"])
        if pd.notna(k) and pd.notna(v)
    }

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
        zid = str(a["zfin_id"]).strip() if pd.notna(a["zfin_id"]) else ""
        alias_str = str(a["alias_string"]).strip() if pd.notna(a["alias_string"]) else ""
        atype = str(a["alias_type"]).strip().upper() if pd.notna(a["alias_type"]) else ""
        if not zid or not alias_str:
            continue
        # ZFIN alias_type enum is small and stable; use exact match on known values.
        if atype == "PREVIOUS NAME":
            prev_by_id[zid].append(alias_str)
        else:
            alias_by_id[zid].append(alias_str)

    rows = []
    for _, g in genes_df.iterrows():
        zid = str(g["zfin_id"]).strip() if pd.notna(g["zfin_id"]) else ""
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


def _build_fruitfly_reference(config, ref_dir: str) -> None:
    """Build fruit fly reference from FlyBase gene annotation + synonyms files."""
    map_url = config.reference_sources["flybase_gene_map"]["url"]
    syn_url = config.reference_sources["flybase_synonyms"]["url"]

    map_raw = _download_file(map_url)
    syn_raw = _download_file(syn_url)

    # Decompress gzip
    map_data = gzip.decompress(map_raw)
    syn_data = gzip.decompress(syn_raw)

    map_checksum = hashlib.sha256(map_raw).hexdigest()
    syn_checksum = hashlib.sha256(syn_raw).hexdigest()

    # FlyBase files have ## comment lines at top, then a ##-prefixed header line.
    # Anchor on column-token content to tolerate FlyBase metadata lines that
    # also lack a trailing space (e.g., "##date\t2024-01-01").
    _HEADER_TOKENS = ("gene_symbol", "primary_FBgn", "primary_FBid", "symbol_synonym")

    def _read_flybase_tsv(data: bytes) -> pd.DataFrame:
        lines = data.decode("utf-8", errors="replace").split("\n")
        header_idx = None
        for i, line in enumerate(lines):
            if (line.startswith("##") and "\t" in line
                    and any(tok in line for tok in _HEADER_TOKENS)):
                header_idx = i
                break
        if header_idx is None:
            raise ValueError("FlyBase TSV header line with expected column tokens not found")
        header = lines[header_idx].lstrip("#").strip().split("\t")
        data_lines = [l for l in lines[header_idx + 1:] if l.strip() and not l.startswith("#")]
        df = pd.DataFrame(
            [line.split("\t") for line in data_lines],
            columns=header,
        )
        return df

    def _require_col(df: pd.DataFrame, needle: str, label: str) -> str:
        match = next((c for c in df.columns if needle in c), None)
        if match is None:
            raise ValueError(f"FlyBase: expected column containing '{label}' not found in {list(df.columns)}")
        return match

    map_df = _read_flybase_tsv(map_data)
    syn_df = _read_flybase_tsv(syn_data)

    sym_col = _require_col(map_df, "gene_symbol", "gene_symbol")
    primary_col = _require_col(map_df, "primary_FBgn", "primary_FBgn")
    secondary_col = next((c for c in map_df.columns if "secondary_FBgn" in c), None)

    # Filter to Dmel only (some files include other Drosophila species)
    if "organism_abbreviation" in map_df.columns:
        map_df = map_df[map_df["organism_abbreviation"] == "Dmel"]

    syn_fbgn_col = next((c for c in syn_df.columns if "primary_FBid" in c or "primary_FBgn" in c), None)
    if syn_fbgn_col is None:
        raise ValueError(f"FlyBase: expected 'primary_FBid' or 'primary_FBgn' column in synonyms file, got {list(syn_df.columns)}")
    syn_symbol_col = _require_col(syn_df, "symbol_synonym", "symbol_synonym")
    if "organism_abbreviation" in syn_df.columns:
        syn_df = syn_df[syn_df["organism_abbreviation"] == "Dmel"]

    # Build FBgn -> list of symbol synonyms
    syn_by_fbgn = {}
    for _, r in syn_df.iterrows():
        fbgn = str(r[syn_fbgn_col]).strip() if pd.notna(r[syn_fbgn_col]) else ""
        raw = str(r[syn_symbol_col]).strip() if pd.notna(r[syn_symbol_col]) else ""
        syns = [s.strip() for s in raw.split(",") if s.strip()] if raw else []
        if fbgn:
            syn_by_fbgn[fbgn] = syns

    rows = []
    for _, m in map_df.iterrows():
        symbol = str(m[sym_col]).strip() if pd.notna(m[sym_col]) else ""
        fbgn = str(m[primary_col]).strip() if pd.notna(m[primary_col]) else ""
        if not symbol or not fbgn.startswith("FBgn"):
            continue

        # secondary FBgn# are "previous IDs" for the same gene
        prev_fbgns = []
        if secondary_col and pd.notna(m[secondary_col]):
            raw_sec = str(m[secondary_col]).strip()
            prev_fbgns = [s.strip() for s in raw_sec.split(",") if s.strip().startswith("FBgn")]

        alias_syms = syn_by_fbgn.get(fbgn, [])

        rows.append({
            "ensembl_id": fbgn,  # FlyBase FBgn is the primary gene ID
            "symbol": symbol,
            "alias_symbols": "|".join(alias_syms),
            "prev_symbols": "|".join(prev_fbgns),
            "gene_type": "",
            "status": "approved",
            "source": "FlyBase",
            "source_id": f"FlyBase:{fbgn}",
        })

    gene_table = pd.DataFrame(rows)
    symbol_lookup = _build_symbol_lookup(gene_table, source="FlyBase")

    metadata = {
        "species": "fruit_fly",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "flybase_gene_map": {"url": map_url, "sha256": map_checksum, "rows": len(map_df)},
            "flybase_synonyms": {"url": syn_url, "sha256": syn_checksum, "rows": len(syn_df)},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built fruit fly reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


def _build_celegans_reference(config, ref_dir: str) -> None:
    """Build C. elegans reference from WormBase geneIDs + geneOtherIDs files."""
    ids_url = config.reference_sources["wormbase_gene_ids"]["url"]
    other_url = config.reference_sources["wormbase_other_ids"]["url"]

    ids_raw = _download_file(ids_url)
    other_raw = _download_file(other_url)

    ids_data = gzip.decompress(ids_raw)
    other_data = gzip.decompress(other_raw)

    ids_checksum = hashlib.sha256(ids_raw).hexdigest()
    other_checksum = hashlib.sha256(other_raw).hexdigest()

    # geneIDs.txt is CSV (not TSV): taxon_id, WBGene, public_name, sequence_name, status, biotype
    ids_df = pd.read_csv(
        io.BytesIO(ids_data), sep=",", header=None,
        names=["taxon_id", "wbgene", "public_name", "sequence_name", "status", "biotype"],
        low_memory=False, dtype=str,
    )

    # geneOtherIDs.txt is TSV: WBGene, public_name, sequence_name, other_names, other_ids
    other_df = pd.read_csv(
        io.BytesIO(other_data), sep="\t", header=None,
        names=["wbgene", "public_name", "sequence_name", "other_names", "other_ids"],
        low_memory=False, dtype=str,
    )

    other_by_wb = {}
    for _, r in other_df.iterrows():
        wb = str(r["wbgene"]).strip() if pd.notna(r["wbgene"]) else ""
        raw = str(r["other_names"]).strip() if pd.notna(r["other_names"]) else ""
        others = [s.strip() for s in raw.split(" ") if s.strip()] if raw else []
        if wb:
            other_by_wb[wb] = others

    rows = []
    for _, g in ids_df.iterrows():
        wb = str(g["wbgene"]).strip() if pd.notna(g["wbgene"]) else ""
        public = str(g["public_name"]).strip() if pd.notna(g["public_name"]) else ""
        seq = str(g["sequence_name"]).strip() if pd.notna(g["sequence_name"]) else ""
        raw_status = str(g["status"]).strip() if pd.notna(g["status"]) else ""
        biotype = str(g["biotype"]).strip() if pd.notna(g["biotype"]) else ""

        if not wb or not wb.startswith("WBGene"):
            continue

        symbol = public if public else seq
        if not symbol:
            continue

        status = "approved" if raw_status.lower() == "live" else "withdrawn"

        # Aliases: sequence_name (if differs from symbol) + other_names
        aliases = []
        if seq and seq != symbol:
            aliases.append(seq)
        aliases.extend(other_by_wb.get(wb, []))

        rows.append({
            "ensembl_id": wb,  # WBGene is the primary gene ID
            "symbol": symbol,
            "alias_symbols": "|".join(aliases),
            "prev_symbols": "",
            "gene_type": biotype,
            "status": status,
            "source": "WormBase",
            "source_id": f"WormBase:{wb}",
        })

    gene_table = pd.DataFrame(rows)
    symbol_lookup = _build_symbol_lookup(gene_table, source="WormBase")

    metadata = {
        "species": "c_elegans",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "wormbase_gene_ids": {"url": ids_url, "sha256": ids_checksum, "rows": len(ids_df)},
            "wormbase_other_ids": {"url": other_url, "sha256": other_checksum, "rows": len(other_df)},
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built C. elegans reference: %d genes, %d lookup entries", len(gene_table), len(symbol_lookup))


def _build_ensembl_biomart_reference(config, ref_dir: str) -> None:
    """Build reference from Ensembl BioMart (Tier 2 species with no dedicated authority).

    Used for cynomolgus macaque, rhesus macaque, common marmoset, and mouse lemur.
    Symbol lookups come from BioMart's external_gene_name + external_synonym attributes.
    No previous_symbols (BioMart doesn't track those). Status hardcoded to "approved".
    """
    import urllib.parse

    src = config.reference_sources["ensembl_biomart"]
    base_url = src["url"]
    dataset = src["dataset"]

    biomart_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1">'
        f'<Dataset name="{dataset}" interface="default">'
        '<Attribute name="ensembl_gene_id"/>'
        '<Attribute name="external_gene_name"/>'
        '<Attribute name="external_synonym"/>'
        '<Attribute name="gene_biotype"/>'
        '</Dataset></Query>'
    )
    full_url = base_url + urllib.parse.quote(biomart_xml)

    raw_data = _download_file(full_url)
    checksum = hashlib.sha256(raw_data).hexdigest()

    bm_df = pd.read_csv(io.BytesIO(raw_data), sep="\t", low_memory=False, dtype=str)
    bm_df.columns = bm_df.columns.str.strip()

    if len(bm_df) == 0:
        raise RuntimeError(
            f"Ensembl BioMart returned empty result for dataset '{dataset}'. "
            "Check network connectivity and dataset name."
        )

    # BioMart may return different column names depending on Ensembl version
    eid_col = next(c for c in bm_df.columns if "gene stable id" in c.lower() or c.lower() == "ensembl_gene_id")
    sym_col = next(c for c in bm_df.columns if "gene name" in c.lower() or c.lower() == "external_gene_name")
    syn_col = next((c for c in bm_df.columns if "synonym" in c.lower()), None)
    type_col = next((c for c in bm_df.columns if "gene type" in c.lower() or c.lower() == "gene_biotype"), None)

    # Group by ensembl_gene_id (BioMart returns one row per synonym)
    by_eid = {}
    for _, r in bm_df.iterrows():
        eid = str(r[eid_col]).strip() if pd.notna(r[eid_col]) else ""
        if not eid:
            continue
        sym = str(r[sym_col]).strip() if pd.notna(r[sym_col]) else ""
        syn = str(r[syn_col]).strip() if syn_col and pd.notna(r[syn_col]) else ""
        gtype = str(r[type_col]).strip() if type_col and pd.notna(r[type_col]) else ""

        if eid not in by_eid:
            by_eid[eid] = {"symbol": sym, "synonyms": set(), "gene_type": gtype}
        if syn:
            by_eid[eid]["synonyms"].add(syn)
        # Prefer non-empty symbol/gene_type if first row had them empty
        if sym and not by_eid[eid]["symbol"]:
            by_eid[eid]["symbol"] = sym
        if gtype and not by_eid[eid]["gene_type"]:
            by_eid[eid]["gene_type"] = gtype

    rows = []
    for eid, info in by_eid.items():
        # A synonym that matches the symbol is not a synonym
        synonyms = sorted(info["synonyms"] - {info["symbol"]})
        rows.append({
            "ensembl_id": eid,
            "symbol": info["symbol"],
            "alias_symbols": "|".join(synonyms),
            "prev_symbols": "",
            "gene_type": info["gene_type"],
            "status": "approved",
            "source": "Ensembl",
            "source_id": f"Ensembl:{eid}",
        })

    gene_table = pd.DataFrame(rows)
    symbol_lookup = _build_symbol_lookup(gene_table, source="Ensembl")

    metadata = {
        "species": config.name,
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "ensembl_biomart": {
                "url": full_url,
                "dataset": dataset,
                "sha256": checksum,
                "rows": len(bm_df),
            },
        },
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
    logger.info("Built %s reference from BioMart: %d genes, %d lookup entries", config.name, len(gene_table), len(symbol_lookup))
