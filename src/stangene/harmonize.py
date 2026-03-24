"""Tiered harmonization cascade for gene identifier mapping."""

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


def harmonize(ft: pd.DataFrame, ref: dict) -> HarmonizationResult:
    result = ft.copy()
    gene_table = ref["gene_table"]
    symbol_lookup = ref["symbol_lookup"]

    for col in ["gene_id_harmonized", "gene_symbol_harmonized",
                "mapping_status", "mapping_confidence", "mapping_source", "mapping_notes"]:
        if col not in result.columns:
            result[col] = None

    ensembl_id_set = set(gene_table["ensembl_id"].dropna())
    ensembl_id_to_row = gene_table.dropna(subset=["ensembl_id"]).set_index("ensembl_id")

    approved_lookup = symbol_lookup[symbol_lookup["lookup_type"] == "approved_symbol"]
    alias_prev_lookup = symbol_lookup[symbol_lookup["lookup_type"].isin(["alias_symbol", "prev_symbol"])]

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
            gene_row = ensembl_id_to_row.loc[feature_id]
            if isinstance(gene_row, pd.DataFrame):
                gene_row = gene_row.iloc[0]
            _apply_match(result, idx, gene_row, "exact_id", "high", "HGNC:exact_id", notes)
            continue

        # Tier 2: Version-stripped ID match
        if feature_id_nv and feature_id_nv in ensembl_id_set:
            gene_row = ensembl_id_to_row.loc[feature_id_nv]
            if isinstance(gene_row, pd.DataFrame):
                gene_row = gene_row.iloc[0]
            notes.append(f"Matched via version-stripped ID: {feature_id} -> {feature_id_nv}")
            _apply_match(result, idx, gene_row, "id_no_version", "high", "HGNC:id_no_version", notes)
            continue

        # Tier 3: Exact approved symbol match
        matches = approved_lookup[approved_lookup["lookup_string"] == feature_name]
        if len(matches) == 1:
            match = matches.iloc[0]
            eid = match["ensembl_id"]
            confidence = "high"
            # Check if matched gene is withdrawn or non-protein-coding
            gene_info = None
            if eid and eid in ensembl_id_set:
                gene_info = ensembl_id_to_row.loc[eid]
                if isinstance(gene_info, pd.DataFrame):
                    gene_info = gene_info.iloc[0]
            elif pd.notna(match.get("source_id")):
                sid_rows = gene_table[gene_table["source_id"] == match["source_id"]]
                if len(sid_rows) > 0:
                    gene_info = sid_rows.iloc[0]
            if gene_info is not None:
                if str(gene_info.get("status", "")).lower().startswith("entry withdrawn") or str(gene_info.get("status", "")).lower() == "withdrawn":
                    confidence = "medium"
                    notes.append("Matched withdrawn gene")
                gene_type = gene_info.get("gene_type", "")
                if gene_type and "protein" not in str(gene_type).lower():
                    notes.append(f"Non-protein-coding gene type: {gene_type}")
            _apply_match_from_lookup(result, idx, match, gene_table, "exact_symbol", confidence, notes)
            continue
        elif len(matches) > 1:
            candidates = matches["ensembl_id"].tolist()
            notes.append(f"Multiple approved symbol matches: {candidates}")
            result.at[idx, "mapping_status"] = "ambiguous"
            result.at[idx, "mapping_confidence"] = "low"
            result.at[idx, "mapping_notes"] = "; ".join(notes)
            continue

        # Tier 4: Alias / previous symbol match
        matches = alias_prev_lookup[alias_prev_lookup["lookup_string"] == feature_name]
        unique_targets = matches.drop_duplicates(subset=["ensembl_id", "source_id"])
        if len(unique_targets) == 1:
            match = matches.iloc[0]
            lookup_type = match["lookup_type"]
            status = "alias_symbol" if lookup_type == "alias_symbol" else "previous_symbol"
            _apply_match_from_lookup(result, idx, match, gene_table, status, "medium", notes)
            continue
        elif len(unique_targets) > 1:
            candidates = unique_targets[["ensembl_id", "source_id", "lookup_type"]].to_dict("records")
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


def _apply_match(result, idx, gene_row, status, confidence, source, notes):
    # gene_row is a pandas Series from set_index("ensembl_id"), so the ensembl_id
    # is the Series name (index label), not a column. Fall back to source_id if needed.
    eid = gene_row.name if gene_row.name and pd.notna(gene_row.name) else gene_row.get("ensembl_id")
    result.at[idx, "gene_id_harmonized"] = eid or gene_row.get("source_id", "")
    result.at[idx, "gene_symbol_harmonized"] = gene_row.get("symbol", "")
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)


def _apply_match_from_lookup(result, idx, match, gene_table, status, confidence, notes):
    eid = match.get("ensembl_id")
    sid = match.get("source_id")
    source_label = f"{match.get('source', '')}:{match.get('lookup_type', '')}"

    symbol = ""
    if pd.notna(eid):
        gene_rows = gene_table[gene_table["ensembl_id"] == eid]
        if len(gene_rows) > 0:
            symbol = gene_rows.iloc[0]["symbol"]
    if not symbol and pd.notna(sid):
        gene_rows = gene_table[gene_table["source_id"] == sid]
        if len(gene_rows) > 0:
            symbol = gene_rows.iloc[0]["symbol"]

    harmonized_id = eid if pd.notna(eid) else sid
    result.at[idx, "gene_id_harmonized"] = harmonized_id
    result.at[idx, "gene_symbol_harmonized"] = symbol
    result.at[idx, "mapping_status"] = status
    result.at[idx, "mapping_confidence"] = confidence
    result.at[idx, "mapping_source"] = source_label
    if notes:
        result.at[idx, "mapping_notes"] = "; ".join(notes)
