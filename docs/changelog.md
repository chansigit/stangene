# Changelog

## v0.4.0 (2026-06-24)

### Features

- `canonical_biotype` column added to all 10 `gene_table.parquet` files
- New `biotype.py` module with `CANONICAL_BIOTYPES` frozenset (13 categories: `protein_coding`, `lncRNA`, `pseudogene`, `miRNA`, `snoRNA`, `snRNA`, `rRNA`, `tRNA`, `IG_gene`, `TR_gene`, `other_ncrna`, `other`, `unknown`)
- `normalize_biotype(raw, source)` function maps raw `gene_type` strings to the canonical vocabulary for all 7 upstream sources (HGNC, MGI, RGD, Ensembl, WormBase, ZFIN, FlyBase)
- Harmonization notes now include `canonical_biotype` for all mapped genes

## v0.3.0 (2026-06-20)

### Features

- `mito_mask()`: species-aware mitochondrial gene detection across all 10 species; handles per-species prefix conventions (`MT-` for mammals/zebrafish, `mt:` for fruit fly, curated symbol sets for C. elegans and bare-name primates)
- `hb_mask()`: species-aware hemoglobin gene detection across all 10 species; uses curated explicit symbol sets per species rather than a naive prefix to avoid false positives

## v0.2.0 (2026-04-23)

### Features

- 7-species expansion: rat, zebrafish, fruit_fly, c_elegans, cynomolgus, rhesus, marmoset, mouse_lemur (10 species total)
- `resolve_species()`: resolves short codes and full names to canonical species identifiers (e.g. `"hs"` → `"human"`, `"cyno"` → `"cynomolgus"`)
- `harmonize_anndata()`: in-memory / pipeline harmonization that modifies an existing AnnData object rather than writing to disk

## v0.1.0 (2026-03-24)

Initial release.

### Features

- 5-tier harmonization cascade: exact Ensembl ID, version-stripped ID, exact approved symbol, alias/previous symbol, unmapped
- Human gene harmonization via HGNC
- Mouse gene harmonization via MGI + Ensembl BioMart supplementary
- Feature classification: auto-detects gene, transcript, antibody capture, CRISPR guide, spike-in, and peak features
- h5ad and TSV/CSV input support
- Conservative merge with strict and symbol policies
- Comprehensive reporting: summary JSON, conflict TSV, unmapped TSV, markdown report
- Enriched h5ad output with harmonization columns in `adata.var`
- Excel corruption detection: date-format artifacts and HGNC-renamed gene symbols
- Withdrawn gene flagging
- CLI: `stangene harmonize` and `stangene build-refs`
- Claude Code skill integration
