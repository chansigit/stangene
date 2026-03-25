# How It Works

## Pipeline Overview

```
Input Dataset (.h5ad / .tsv / .csv)
        |
        v
   load_features()          Extract feature metadata (not expression matrix)
        |
        v
   classify_features()      Triage: gene vs non-gene features
        |
    +---+---+
    |       |
    v       v
  Genes   Non-gene features → labeled and passed through
    |
    v
   load_reference()          Load species-specific annotation tables
        |
        v
   harmonize()               5-tier matching cascade
        |
        v
   HarmonizationResult       Mapping table + conflicts + stats
        |
    +---+---+
    |       |
    v       v
  write_results()     write_reports()
  (TSV + h5ad)        (summary + conflicts + report.md)
```

## Feature Classification

Before harmonization, every feature is classified to separate true gene features from non-gene features that should not go through gene-name matching.

| Pattern | Type | Action |
|---|---|---|
| `ENSG*` / `ENSMUSG*` | `gene` | Harmonize |
| `ENST*` / `ENSMUST*` | `transcript` | Skip (v1) |
| `*_ADT`, `*_HTO`, `*TotalSeq*` | `antibody_capture` | Skip |
| `sg-*`, `gRNA-*` | `crispr_guide` | Skip |
| `ERCC-*` | `spike_in` | Skip |
| `chrN:start-end` | `peak` | Skip |
| 10x Cell Ranger labels | Mapped directly | Depends on type |
| Anything else | `gene` (default) | Harmonize |

Non-gene features receive `mapping_status = non_gene_feature` and are preserved in the output without modification.

## Matching Cascade

Gene features are matched against the reference database in strict priority order. Once a feature is resolved at a tier, lower tiers are skipped (**early exit**).

### Tier 1: Exact Ensembl ID

If the feature has an Ensembl gene ID (e.g., `ENSG00000141510`) and it exactly matches an ID in the reference, the mapping is immediate.

- **Confidence:** high
- **Example:** `ENSG00000141510` → TP53

### Tier 2: Version-stripped ID

If the Ensembl ID includes a version suffix (e.g., `ENSG00000141510.18`), the suffix is stripped and the base ID is matched.

- **Confidence:** high
- **Example:** `ENSG00000141510.18` → strip to `ENSG00000141510` → TP53

### Tier 3: Exact approved symbol

If no stable ID is available, the feature name is matched against official approved gene symbols from HGNC (human) or MGI (mouse).

- **Confidence:** high (downgraded to medium if the matched gene is withdrawn)
- **Example:** `TP53` matches the HGNC approved symbol for ENSG00000141510

### Tier 4: Alias / previous symbol

If the approved symbol doesn't match, the feature name is checked against:
- **Alias symbols** (alternative names) → `mapping_status = alias_symbol`
- **Previous symbols** (old names) → `mapping_status = previous_symbol`

If the alias maps to exactly one gene, it's used. If multiple genes share the alias, the feature is marked ambiguous.

- **Confidence:** medium
- **Example:** `p53` is an alias for TP53; `RNF53` is a previous name for BRCA1

### Tier 5: Unmapped

If no match is found at any tier, the feature is marked `unmapped`. It is preserved in the output for manual review.

- **Confidence:** N/A

### Additional checks

- **Excel date detection:** Features like `1-Mar`, `2-Sep` are detected as likely spreadsheet artifacts and marked unmapped with a warning.
- **Known renamed genes:** Symbols that HGNC renamed in 2020 due to Excel auto-conversion (e.g., MARCH1→MARCHF1, SEPT2→SEPTIN2) are flagged in mapping notes.
- **Withdrawn genes:** Matches against withdrawn genes are downgraded to medium confidence.

## Design Principles

1. **Never destroy original information.** All original identifiers are preserved in dedicated columns. Harmonized identifiers live in separate columns.

2. **Stable IDs over symbols.** Ensembl gene IDs are the preferred canonical key. Gene symbols are a display layer, not the identity key.

3. **Separate identity from display.** `gene_id_harmonized` (canonical identity) and `gene_symbol_harmonized` (display name) are distinct outputs.

4. **Within-species only.** Each species is harmonized independently. Cross-species linkage via orthology is a separate concern.

5. **Conservative by default.** An explicit `ambiguous` or `unmapped` is always better than an incorrect forced mapping.

6. **Full traceability.** Every mapping records its tier, confidence, source, and notes.
