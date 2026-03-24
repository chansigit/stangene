"""Conservative merge logic for harmonized features."""

from dataclasses import dataclass, field

import pandas as pd

from stangene._logging import get_logger

logger = get_logger("merge")

_STRICT_ELIGIBLE = frozenset(["exact_id", "id_no_version"])
_SYMBOL_ELIGIBLE = frozenset(["exact_id", "id_no_version", "exact_symbol"])


@dataclass
class MergeResult:
    """Result of conservative merge."""
    merged_table: pd.DataFrame
    provenance: pd.DataFrame
    merge_log: list = field(default_factory=list)


def merge_features(result, policy: str = "strict") -> MergeResult:
    if policy == "strict":
        eligible = _STRICT_ELIGIBLE
    elif policy == "symbol":
        eligible = _SYMBOL_ELIGIBLE
    else:
        raise ValueError(f"Unknown merge policy: {policy}. Use 'strict' or 'symbol'.")

    mt = result.mapping_table.copy()
    merge_log = []
    provenance_rows = []

    eligible_mask = mt["mapping_status"].isin(eligible) & mt["gene_id_harmonized"].notna()
    eligible_df = mt[eligible_mask]
    ineligible_df = mt[~eligible_mask]

    merged_rows = []
    for hid, group in eligible_df.groupby("gene_id_harmonized"):
        if len(group) == 1:
            merged_rows.append(group.iloc[0].to_dict())
            provenance_rows.append({
                "gene_id_harmonized": hid,
                "original_feature_name": group.iloc[0]["original_feature_name"],
                "dataset": group.iloc[0].get("dataset", ""),
                "mapping_status": group.iloc[0]["mapping_status"],
                "merge_action": "kept_single",
            })
        else:
            rep = group.iloc[0].to_dict()
            originals = group["original_feature_name"].tolist()
            merge_log.append(
                f"Merged {len(group)} features into {hid} ({rep.get('gene_symbol_harmonized', '')}): {originals}"
            )
            merged_rows.append(rep)
            for _, row in group.iterrows():
                provenance_rows.append({
                    "gene_id_harmonized": hid,
                    "original_feature_name": row["original_feature_name"],
                    "dataset": row.get("dataset", ""),
                    "mapping_status": row["mapping_status"],
                    "merge_action": "merged",
                })

    for _, row in ineligible_df.iterrows():
        merged_rows.append(row.to_dict())
        provenance_rows.append({
            "gene_id_harmonized": row.get("gene_id_harmonized"),
            "original_feature_name": row["original_feature_name"],
            "dataset": row.get("dataset", ""),
            "mapping_status": row["mapping_status"],
            "merge_action": "not_eligible",
        })

    merged_table = pd.DataFrame(merged_rows)
    provenance = pd.DataFrame(provenance_rows)

    logger.info("Merge complete (policy=%s): %d -> %d rows, %d merges", policy, len(mt), len(merged_table), len(merge_log))
    return MergeResult(merged_table=merged_table, provenance=provenance, merge_log=merge_log)
