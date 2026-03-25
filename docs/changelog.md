# Changelog

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
