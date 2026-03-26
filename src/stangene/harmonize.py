"""Tiered harmonization cascade for gene identifier mapping."""

from collections import defaultdict
from dataclasses import dataclass

import pandas as pd

from stangene import __version__
from stangene._logging import get_logger
from stangene.species import EXCEL_DATE_PATTERN, EXCEL_RENAMED_GENES

logger = get_logger("harmonize")


@dataclass
class HarmonizationResult:
    """Result of the harmonization cascade."""
    mapping_table: pd.DataFrame
    conflicts: pd.DataFrame
    stats: dict


def _build_lookup_dicts(symbol_lookup: pd.DataFrame):
    """Pre-build dictionary lookups for O(1) access per feature.

    Returns:
        approved_dict: {lookup_string: [list of matching rows as dicts]}
        alias_prev_dict: {lookup_string: [list of matching rows as dicts]}
    """
    approved_dict = defaultdict(list)       # exact case
    approved_dict_upper = defaultdict(list)  # uppercase fallback
    alias_prev_dict = defaultdict(list)
    alias_prev_dict_upper = defaultdict(list)

    for row in symbol_lookup.itertuples(index=False):
        entry = {
            "lookup_string": row.lookup_string,
            "ensembl_id": row.ensembl_id if pd.notna(row.ensembl_id) else None,
            "source_id": row.source_id,
            "lookup_type": row.lookup_type,
            "source": row.source,
        }
        upper_key = row.lookup_string_upper if pd.notna(row.lookup_string_upper) else row.lookup_string.upper()
        if row.lookup_type == "approved_symbol":
            approved_dict[row.lookup_string].append(entry)
            approved_dict_upper[upper_key].append(entry)
        else:
            alias_prev_dict[row.lookup_string].append(entry)
            alias_prev_dict_upper[upper_key].append(entry)

    return (dict(approved_dict), dict(approved_dict_upper),
            dict(alias_prev_dict), dict(alias_prev_dict_upper))


def _build_gene_info_dicts(gene_table: pd.DataFrame):
    """Pre-build gene info lookups by ensembl_id and source_id.

    Returns:
        eid_to_gene: {ensembl_id: {symbol, status, gene_type, source_id, ...}}
        sid_to_gene: {source_id: {symbol, status, gene_type, ensembl_id, ...}}
    """
    eid_to_gene = {}
    sid_to_gene = {}

    for row in gene_table.itertuples(index=False):
        gene_info = {
            "ensembl_id": row.ensembl_id if pd.notna(row.ensembl_id) else None,
            "symbol": row.symbol,
            "status": row.status if pd.notna(row.status) else "",
            "gene_type": row.gene_type if pd.notna(row.gene_type) else "",
            "source_id": row.source_id,
        }
        if pd.notna(row.ensembl_id):
            eid_to_gene[row.ensembl_id] = gene_info
        if pd.notna(row.source_id):
            sid_to_gene[row.source_id] = gene_info

    return eid_to_gene, sid_to_gene


def harmonize(ft: pd.DataFrame, ref: dict) -> HarmonizationResult:
    result = ft.copy()
    gene_table = ref["gene_table"]
    symbol_lookup = ref["symbol_lookup"]

    for col in ["gene_id_harmonized", "gene_symbol_harmonized",
                "mapping_status", "mapping_confidence", "mapping_source", "mapping_notes"]:
        if col not in result.columns:
            result[col] = None

    # Build O(1) lookup structures
    ensembl_id_set = set(gene_table["ensembl_id"].dropna())
    eid_to_gene, sid_to_gene = _build_gene_info_dicts(gene_table)
    approved_dict, approved_dict_upper, alias_prev_dict, alias_prev_dict_upper = _build_lookup_dicts(symbol_lookup)

    logger.info("Built lookup dicts: %d ensembl IDs, %d approved symbols, %d alias/prev entries",
                len(ensembl_id_set), len(approved_dict), len(alias_prev_dict))

    for idx in result.index:
        if result.at[idx, "mapping_status"] == "non_gene_feature":
            continue

        feature_name = str(result.at[idx, "original_feature_name"]) if pd.notna(result.at[idx, "original_feature_name"]) else ""
        feature_id = result.at[idx, "original_feature_id"] if pd.notna(result.at[idx, "original_feature_id"]) else None
        feature_id_nv = result.at[idx, "feature_id_no_version"] if pd.notna(result.at[idx, "feature_id_no_version"]) else None

        notes = []

        # Check for Excel date corruption
        if EXCEL_DATE_PATTERN.match(feature_name):
            result.at[idx, "mapping_status"] = "unmapped"
            result.at[idx, "mapping_confidence"] = None
            result.at[idx, "mapping_notes"] = f"Likely Excel-corrupted date format: {feature_name}"
            continue

        if feature_name.upper() in EXCEL_RENAMED_GENES:
            notes.append(f"Known Excel-renamed gene: {feature_name} -> {EXCEL_RENAMED_GENES[feature_name.upper()]}")

        # Tier 1: Exact stable ID match
        if feature_id and feature_id in ensembl_id_set:
            gene_info = eid_to_gene[feature_id]
            _apply_gene_match(result, idx, gene_info, "exact_id", "high",
                              f"{gene_info.get('source_id', 'HGNC')}:exact_id", notes)
            continue

        # Tier 2: Version-stripped ID match
        if feature_id_nv and feature_id_nv in ensembl_id_set:
            gene_info = eid_to_gene[feature_id_nv]
            notes.append(f"Matched via version-stripped ID: {feature_id} -> {feature_id_nv}")
            _apply_gene_match(result, idx, gene_info, "id_no_version", "high",
                              f"{gene_info.get('source_id', 'HGNC')}:id_no_version", notes)
            continue

        # Tier 3: Exact approved symbol match (case-sensitive first, then uppercase fallback)
        matches = approved_dict.get(feature_name, [])
        if not matches:
            matches = approved_dict_upper.get(feature_name.upper(), [])
        if len(matches) == 1:
            match = matches[0]
            eid = match["ensembl_id"]
            confidence = "high"
            # Check if matched gene is withdrawn or non-protein-coding
            gene_info = None
            if eid and eid in eid_to_gene:
                gene_info = eid_to_gene[eid]
            elif match.get("source_id") and match["source_id"] in sid_to_gene:
                gene_info = sid_to_gene[match["source_id"]]
            if gene_info is not None:
                status_lower = str(gene_info.get("status", "")).lower()
                if status_lower.startswith("entry withdrawn") or status_lower == "withdrawn":
                    confidence = "medium"
                    notes.append("Matched withdrawn gene")
                gene_type = gene_info.get("gene_type", "")
                if gene_type and "protein" not in str(gene_type).lower():
                    notes.append(f"Non-protein-coding gene type: {gene_type}")
            _apply_lookup_match(result, idx, match, eid_to_gene, sid_to_gene,
                                "exact_symbol", confidence, notes)
            continue
        elif len(matches) > 1:
            candidates = [m["ensembl_id"] for m in matches]
            notes.append(f"Multiple approved symbol matches: {candidates}")
            result.at[idx, "mapping_status"] = "ambiguous"
            result.at[idx, "mapping_confidence"] = "low"
            result.at[idx, "mapping_notes"] = "; ".join(notes)
            continue

        # Tier 4: Alias / previous symbol match (case-sensitive first, then uppercase fallback)
        matches = alias_prev_dict.get(feature_name, [])
        if not matches:
            matches = alias_prev_dict_upper.get(feature_name.upper(), [])
        # Deduplicate by (ensembl_id, source_id) — same gene can appear as both alias and prev
        seen_targets = {}
        for m in matches:
            key = (m["ensembl_id"], m["source_id"])
            if key not in seen_targets:
                seen_targets[key] = m

        unique_targets = list(seen_targets.values())
        if len(unique_targets) == 1:
            # Prefer alias_symbol over previous_symbol deterministically
            has_alias = any(m["lookup_type"] == "alias_symbol" for m in matches)
            status = "alias_symbol" if has_alias else "previous_symbol"
            match = matches[0]
            _apply_lookup_match(result, idx, match, eid_to_gene, sid_to_gene,
                                status, "medium", notes)
            continue
        elif len(unique_targets) > 1:
            candidates = [{"ensembl_id": m["ensembl_id"], "source_id": m["source_id"],
                           "lookup_type": m["lookup_type"]} for m in unique_targets]
            notes.append(f"Multiple alias/prev matches: {candidates}")
            result.at[idx, "mapping_status"] = "ambiguous"
            result.at[idx, "mapping_confidence"] = "low"
            result.at[idx, "mapping_notes"] = "; ".join(notes)
            continue

        # Tier 5: Unmapped
        result.at[idx, "mapping_status"] = "unmapped"
        result.at[idx, "mapping_confidence"] = None
        if notes:
            result.at[idx, "mapping_notes"] = "; ".join(notes)

    result["stangene_version"] = __version__
    result["reference_release"] = ref.get("metadata", {}).get("download_timestamp", "")

    stats = result["mapping_status"].value_counts().to_dict()

    harmonized_ids = result[result["gene_id_harmonized"].notna()]
    id_counts = harmonized_ids["gene_id_harmonized"].value_counts()
    duplicate_ids = id_counts[id_counts > 1].index.tolist()
    conflicts = result[result["gene_id_harmonized"].isin(duplicate_ids)].copy() if duplicate_ids else pd.DataFrame()

    logger.info("Harmonization complete: %s", stats)
    if len(conflicts) > 0:
        logger.info("Found %d features in %d many-to-one conflicts", len(conflicts), len(duplicate_ids))

    return HarmonizationResult(mapping_table=result, conflicts=conflicts, stats=stats)


def _apply_gene_match(result, idx, gene_info, status, confidence, source, notes):
    """Apply a match from pre-built gene_info dict."""
    eid = gene_info.get("ensembl_id")
    result.at[idx, "gene_id_harmonized"] = eid or gene_info.get("source_id", "")
    result.at[idx, "gene_symbol_harmonized"] = gene_info.get("symbol", "")
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)


def _apply_lookup_match(result, idx, match, eid_to_gene, sid_to_gene, status, confidence, notes):
    """Apply a match from symbol lookup, resolving symbol via gene info dicts."""
    eid = match.get("ensembl_id")
    sid = match.get("source_id")
    source_label = f"{match.get('source', '')}:{match.get('lookup_type', '')}"

    # Resolve symbol
    symbol = ""
    if eid and eid in eid_to_gene:
        symbol = eid_to_gene[eid].get("symbol", "")
    if not symbol and sid and sid in sid_to_gene:
        symbol = sid_to_gene[sid].get("symbol", "")

    harmonized_id = eid if eid else sid
    result.at[idx, "gene_id_harmonized"] = harmonized_id
    result.at[idx, "gene_symbol_harmonized"] = symbol
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source_label
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)
