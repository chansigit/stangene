"""Reporting: summaries, conflict detection, and output writing."""

import json
import os

import pandas as pd

from stangene._logging import get_logger
from stangene.species import EXCEL_DATE_PATTERN, EXCEL_RENAMED_GENES

logger = get_logger("report")


def summary(result) -> dict:
    mt = result.mapping_table
    status_counts = mt["mapping_status"].value_counts().to_dict()

    gene_mask = mt["original_feature_type"].isin(["gene", "Gene Expression"])
    non_gene_mask = ~gene_mask

    harmonized_ids = mt["gene_id_harmonized"].dropna()
    harmonized_syms = mt["gene_symbol_harmonized"].dropna()

    id_dupe_groups = harmonized_ids[harmonized_ids.duplicated(keep=False)].nunique()
    sym_dupe_groups = harmonized_syms[harmonized_syms.duplicated(keep=False)].nunique()

    return {
        "total_features": len(mt),
        "gene_features": int(gene_mask.sum()),
        "non_gene_features": int(non_gene_mask.sum()),
        "status_counts": status_counts,
        "duplicate_harmonized_ids": int(id_dupe_groups),
        "duplicate_harmonized_symbols": int(sym_dupe_groups),
    }


def conflict_report(result) -> pd.DataFrame:
    mt = result.mapping_table
    rows = []

    harmonized_ids = mt[mt["gene_id_harmonized"].notna()]
    id_counts = harmonized_ids["gene_id_harmonized"].value_counts()
    for hid in id_counts[id_counts > 1].index:
        involved = mt[mt["gene_id_harmonized"] == hid]
        for _, row in involved.iterrows():
            rows.append({
                "conflict_type": "many_to_one",
                "original_feature_name": row["original_feature_name"],
                "gene_id_harmonized": hid,
                "gene_symbol_harmonized": row.get("gene_symbol_harmonized", ""),
                "mapping_status": row["mapping_status"],
                "details": f"{len(involved)} features map to {hid}",
            })

    unmapped = mt[mt["mapping_status"] == "unmapped"]
    for _, row in unmapped.iterrows():
        details = ""
        name = str(row["original_feature_name"])
        if EXCEL_DATE_PATTERN.match(name):
            details = "Likely Excel date artifact"
        elif name.upper() in EXCEL_RENAMED_GENES:
            details = f"Known Excel-renamed gene (old: {name}, new: {EXCEL_RENAMED_GENES[name.upper()]})"
        rows.append({
            "conflict_type": "unmapped",
            "original_feature_name": name,
            "gene_id_harmonized": None,
            "gene_symbol_harmonized": None,
            "mapping_status": "unmapped",
            "details": details,
        })

    prev_sym = mt[mt["mapping_status"] == "previous_symbol"]
    for _, row in prev_sym.iterrows():
        rows.append({
            "conflict_type": "outdated_name",
            "original_feature_name": row["original_feature_name"],
            "gene_id_harmonized": row.get("gene_id_harmonized"),
            "gene_symbol_harmonized": row.get("gene_symbol_harmonized"),
            "mapping_status": row["mapping_status"],
            "details": f"Mapped via previous symbol; current symbol is {row.get('gene_symbol_harmonized', '?')}",
        })

    ambiguous = mt[mt["mapping_status"] == "ambiguous"]
    for _, row in ambiguous.iterrows():
        rows.append({
            "conflict_type": "ambiguous",
            "original_feature_name": row["original_feature_name"],
            "gene_id_harmonized": None,
            "gene_symbol_harmonized": None,
            "mapping_status": "ambiguous",
            "details": row.get("mapping_notes", ""),
        })

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["conflict_type", "original_feature_name", "gene_id_harmonized",
                 "gene_symbol_harmonized", "mapping_status", "details"]
    )


def write_reports(result, output_dir: str, merge_result=None) -> None:
    os.makedirs(output_dir, exist_ok=True)

    # Harmonization table
    result.mapping_table.to_csv(
        os.path.join(output_dir, "harmonization_table.tsv"), sep="\t", index=False,
    )

    # Summary JSON
    s = summary(result)
    with open(os.path.join(output_dir, "summary.json"), "w") as f:
        json.dump(s, f, indent=2, default=str)

    # Conflict report
    cr = conflict_report(result)
    cr.to_csv(os.path.join(output_dir, "conflicts.tsv"), sep="\t", index=False)

    # Unmapped features
    unmapped = result.mapping_table[result.mapping_table["mapping_status"] == "unmapped"]
    unmapped.to_csv(os.path.join(output_dir, "unmapped.tsv"), sep="\t", index=False)

    # Merge outputs
    if merge_result is not None:
        merge_result.merged_table.to_csv(
            os.path.join(output_dir, "merged_table.tsv"), sep="\t", index=False,
        )
        merge_result.provenance.to_csv(
            os.path.join(output_dir, "merge_provenance.tsv"), sep="\t", index=False,
        )

    logger.info("Reports written to %s", output_dir)
