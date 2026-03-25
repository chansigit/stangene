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


def generate_markdown_report(
    result,
    title: str = "Gene Harmonization Report",
    max_unmapped_rows: int = 50,
    max_conflict_rows: int = 50,
    max_outdated_rows: int = 50,
    max_ambiguous_rows: int = 50,
    output_path: str = None,
) -> str:
    """Generate a comprehensive markdown report of harmonization results.

    Args:
        result: HarmonizationResult from harmonize().
        title: Report title.
        max_unmapped_rows: Max unmapped features to list individually.
        max_conflict_rows: Max many-to-one collision groups to list.
        max_outdated_rows: Max outdated name rows to list.
        max_ambiguous_rows: Max ambiguous features to list.
        output_path: If provided, write the report to this file.

    Returns:
        The markdown report as a string.
    """
    mt = result.mapping_table
    s = summary(result)
    cr = conflict_report(result)

    lines = []

    def _add(text=""):
        lines.append(text)

    # --- Header ---
    dataset_name = mt["dataset"].iloc[0] if "dataset" in mt.columns else "unknown"
    species = mt["species"].iloc[0] if "species" in mt.columns else "unknown"
    version = mt["stangene_version"].iloc[0] if "stangene_version" in mt.columns else "?"
    ref_release = mt["reference_release"].iloc[0] if "reference_release" in mt.columns and pd.notna(mt["reference_release"].iloc[0]) else "N/A"

    _add(f"# {title}")
    _add()
    _add("| Field | Value |")
    _add("|---|---|")
    _add(f"| Dataset | `{dataset_name}` |")
    _add(f"| Species | {species} |")
    _add(f"| stangene version | {version} |")
    _add(f"| Reference date | {ref_release} |")
    _add()

    # --- Summary ---
    _add("## Summary")
    _add()
    total = s["total_features"]
    gene = s["gene_features"]
    non_gene = s["non_gene_features"]
    mapped = sum(v for k, v in s["status_counts"].items()
                 if k in ("exact_id", "id_no_version", "exact_symbol", "alias_symbol", "previous_symbol"))
    mapped_pct = f"{mapped / total * 100:.1f}" if total > 0 else "0.0"

    _add("| Metric | Count | % |")
    _add("|---|---|---|")
    _add(f"| Total features | {total} | 100% |")
    _add(f"| Gene features | {gene} | {gene / total * 100:.1f}% |" if total > 0 else f"| Gene features | {gene} | - |")
    _add(f"| Non-gene features | {non_gene} | {non_gene / total * 100:.1f}% |" if total > 0 else f"| Non-gene features | {non_gene} | - |")
    _add(f"| **Successfully mapped** | **{mapped}** | **{mapped_pct}%** |")
    _add(f"| Unmapped | {s['status_counts'].get('unmapped', 0)} | {s['status_counts'].get('unmapped', 0) / total * 100:.1f}% |" if total > 0 else "| Unmapped | 0 | - |")
    _add(f"| Ambiguous | {s['status_counts'].get('ambiguous', 0)} | {s['status_counts'].get('ambiguous', 0) / total * 100:.1f}% |" if total > 0 else "| Ambiguous | 0 | - |")
    _add(f"| Duplicate harmonized IDs | {s['duplicate_harmonized_ids']} | |")
    _add(f"| Duplicate harmonized symbols | {s['duplicate_harmonized_symbols']} | |")
    _add()

    # --- Tier Breakdown ---
    _add("## Mapping Tier Breakdown")
    _add()

    tier_order = [
        ("exact_id", "Tier 1: Exact Ensembl ID", "high"),
        ("id_no_version", "Tier 2: Version-stripped ID", "high"),
        ("exact_symbol", "Tier 3: Exact approved symbol", "high"),
        ("alias_symbol", "Tier 4a: Alias symbol", "medium"),
        ("previous_symbol", "Tier 4b: Previous symbol", "medium"),
        ("ambiguous", "Ambiguous (multiple candidates)", "low"),
        ("unmapped", "Unmapped", "-"),
        ("non_gene_feature", "Non-gene feature (skipped)", "-"),
    ]

    _add("| Tier | Status | Count | % | Confidence |")
    _add("|---|---|---|---|---|")
    for status_key, label, conf in tier_order:
        count = s["status_counts"].get(status_key, 0)
        if count > 0:
            pct = f"{count / total * 100:.1f}%" if total > 0 else "-"
            _add(f"| {label} | `{status_key}` | {count} | {pct} | {conf} |")
    _add()

    # --- Confidence Distribution ---
    _add("## Confidence Distribution")
    _add()
    if "mapping_confidence" in mt.columns:
        conf_counts = mt["mapping_confidence"].value_counts().to_dict()
        null_count = int(mt["mapping_confidence"].isna().sum())
        _add("| Confidence | Count |")
        _add("|---|---|")
        for conf_level in ["high", "medium", "low"]:
            if conf_level in conf_counts:
                _add(f"| {conf_level} | {conf_counts[conf_level]} |")
        if null_count > 0:
            _add(f"| N/A (unmapped/non-gene) | {null_count} |")
        _add()

    # --- Many-to-One Conflicts ---
    m2o = cr[cr["conflict_type"] == "many_to_one"] if len(cr) > 0 else pd.DataFrame()
    if len(m2o) > 0:
        # Group by harmonized ID
        m2o_groups = m2o.groupby("gene_id_harmonized")
        _add("## Many-to-One Conflicts")
        _add()
        _add(f"**{len(m2o_groups)} canonical genes** have multiple original features mapping to them.")
        _add()

        shown = 0
        for hid, group in m2o_groups:
            if shown >= max_conflict_rows:
                remaining = len(m2o_groups) - shown
                _add(f"*... and {remaining} more collision groups.*")
                break
            symbol = group["gene_symbol_harmonized"].iloc[0]
            _add(f"**{symbol}** (`{hid}`) -- {len(group)} features:")
            _add()
            for _, row in group.iterrows():
                status_badge = f"`{row['mapping_status']}`"
                _add(f"- `{row['original_feature_name']}` ({status_badge})")
            _add()
            shown += 1

    # --- Outdated Names (previous_symbol) ---
    outdated = cr[cr["conflict_type"] == "outdated_name"] if len(cr) > 0 else pd.DataFrame()
    if len(outdated) > 0:
        _add("## Outdated Gene Names")
        _add()
        _add(f"**{len(outdated)} features** were mapped via a previous (outdated) symbol.")
        _add("These genes have been renamed in the current reference.")
        _add()
        _add("| Original Name | Current Symbol | Harmonized ID |")
        _add("|---|---|---|")
        shown = 0
        for _, row in outdated.iterrows():
            if shown >= max_outdated_rows:
                remaining = len(outdated) - shown
                _add(f"| *... and {remaining} more* | | |")
                break
            _add(f"| `{row['original_feature_name']}` | {row['gene_symbol_harmonized']} | `{row['gene_id_harmonized']}` |")
            shown += 1
        _add()

    # --- Ambiguous Features ---
    ambiguous = cr[cr["conflict_type"] == "ambiguous"] if len(cr) > 0 else pd.DataFrame()
    if len(ambiguous) > 0:
        _add("## Ambiguous Features")
        _add()
        _add(f"**{len(ambiguous)} features** could not be uniquely resolved (multiple candidate genes).")
        _add()
        _add("| Feature | Details |")
        _add("|---|---|")
        shown = 0
        for _, row in ambiguous.iterrows():
            if shown >= max_ambiguous_rows:
                remaining = len(ambiguous) - shown
                _add(f"| *... and {remaining} more* | |")
                break
            details = str(row.get("details", "")) if pd.notna(row.get("details")) else ""
            _add(f"| `{row['original_feature_name']}` | {details} |")
            shown += 1
        _add()

    # --- Unmapped Features ---
    unmapped_df = mt[mt["mapping_status"] == "unmapped"]
    if len(unmapped_df) > 0:
        _add("## Unmapped Features")
        _add()
        _add(f"**{len(unmapped_df)} features** could not be matched to any reference gene.")
        _add()

        # Check for common patterns in unmapped names
        unmapped_names = unmapped_df["original_feature_name"].tolist()
        excel_dates = [n for n in unmapped_names if EXCEL_DATE_PATTERN.match(str(n))]
        excel_renamed = [n for n in unmapped_names if str(n).upper() in EXCEL_RENAMED_GENES]

        if excel_dates:
            _add(f"**Warning:** {len(excel_dates)} features appear to be Excel-corrupted dates: {', '.join(f'`{n}`' for n in excel_dates[:10])}")
            _add()
        if excel_renamed:
            _add(f"**Warning:** {len(excel_renamed)} features use gene names that were renamed by HGNC due to Excel auto-conversion:")
            for n in excel_renamed[:10]:
                _add(f"- `{n}` (now `{EXCEL_RENAMED_GENES[str(n).upper()]}`)")
            _add()

        # Sample unmapped
        _add("| # | Feature Name | Notes |")
        _add("|---|---|---|")
        shown = 0
        for _, row in unmapped_df.iterrows():
            if shown >= max_unmapped_rows:
                remaining = len(unmapped_df) - shown
                _add(f"| | *... and {remaining} more* | |")
                break
            notes = str(row.get("mapping_notes", "")) if pd.notna(row.get("mapping_notes")) else ""
            _add(f"| {shown + 1} | `{row['original_feature_name']}` | {notes} |")
            shown += 1
        _add()

    # --- Non-Gene Features ---
    non_gene_df = mt[mt["mapping_status"] == "non_gene_feature"]
    if len(non_gene_df) > 0:
        _add("## Non-Gene Features")
        _add()
        type_counts = non_gene_df["original_feature_type"].value_counts().to_dict()
        _add(f"**{len(non_gene_df)} non-gene features** were excluded from harmonization:")
        _add()
        _add("| Feature Type | Count |")
        _add("|---|---|")
        for ftype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
            _add(f"| {ftype} | {count} |")
        _add()

        # Show a few examples per type
        for ftype in type_counts:
            examples = non_gene_df[non_gene_df["original_feature_type"] == ftype]["original_feature_name"].head(5).tolist()
            _add(f"**{ftype}** examples: {', '.join(f'`{e}`' for e in examples)}")
            _add()

    # --- Footer ---
    _add("---")
    _add()
    _add("*Report generated by [stangene](https://github.com/chansigit/stangene). "
         "Ambiguous and unmapped features are preserved for manual review. "
         "No original identifiers were modified.*")

    report = "\n".join(lines)

    if output_path:
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info("Markdown report written to %s", output_path)

    return report


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

    # Markdown report
    generate_markdown_report(result, output_path=os.path.join(output_dir, "report.md"))

    # Merge outputs
    if merge_result is not None:
        merge_result.merged_table.to_csv(
            os.path.join(output_dir, "merged_table.tsv"), sep="\t", index=False,
        )
        merge_result.provenance.to_csv(
            os.path.join(output_dir, "merge_provenance.tsv"), sep="\t", index=False,
        )

    logger.info("Reports written to %s", output_dir)
