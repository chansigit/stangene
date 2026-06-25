# Reference Data

stangene harmonizes gene names against official annotation databases. References must be built once per species before harmonization.

## Building references

```bash
# CLI
stangene build-refs --species human
stangene build-refs --species mouse
stangene build-refs --species human --force  # re-download and rebuild

# Python
from stangene import build_reference
build_reference("human")
build_reference("mouse", force=True)
```

## Human (HGNC)

**Source:** [HGNC complete gene set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) (~15 MB)

Provides:
- Approved gene symbols
- Alias (alternative) symbols
- Previous (old) symbols
- Ensembl gene IDs
- HGNC IDs
- Gene types (protein-coding, lncRNA, pseudogene, etc.)
- Approval status (Approved, Entry Withdrawn)

## Mouse (MGI + Ensembl BioMart)

**Sources:**
- [MGI marker list](https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt) (~7 MB) — approved symbols, synonyms, feature types
- [MGI-to-Ensembl mapping](https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt) — MGI ID to Ensembl ID links
- Ensembl BioMart (supplementary, non-fatal if unavailable) — fills Ensembl ID gaps for mouse genes not covered by MGI mapping

## Rat (RGD)

**Source:** [RGD rat gene file](https://download.rgd.mcw.edu/data_release/RAT/GENES_RAT.txt) (~5 MB)

Provides:
- Approved gene symbols (SYMBOL)
- Previous symbols (OLD_SYMBOL, semicolon-separated)
- Ensembl gene IDs (ENSRNOG prefix; multi-valued field, first canonical ID used)
- RGD gene IDs
- Gene types (protein-coding, ncRNA, lncRNA, etc.)
- Nomenclature status (APPROVED, PROVISIONAL, INTERIM)

## Zebrafish (ZFIN)

**Sources:**
- [ZFIN gene list](https://zfin.org/downloads/gene.txt) — primary gene identifiers with Ensembl IDs (ENSDARG prefix)
- [ZFIN aliases](https://zfin.org/downloads/aliases.txt) — alias and previous symbols, typed by `PREVIOUS NAME` vs other alias types
- [ZFIN Ensembl 1:1 mapping](https://zfin.org/downloads/ensembl_1_to_1.txt) — supplementary ZFIN-to-Ensembl mapping for gap filling

## Fruit Fly (FlyBase)

**Sources:**
- [FlyBase FBgn annotation ID mapping](https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID.tsv.gz) — primary FBgn IDs with current gene symbols and annotation IDs; Drosophila melanogaster only
- [FlyBase synonyms](https://s3ftp.flybase.org/releases/current/precomputed_files/synonyms/fb_synonym_fb_2026_01.tsv.gz) — symbol synonyms and secondary FBgn IDs (previous IDs); versioned filename, bump the release suffix when FlyBase publishes a new release

Note: FBgn IDs serve as the primary Ensembl-equivalent identifiers. Biotype is inferred from annotation ID prefix: `CG*` → protein_coding, `CR*` → other_ncrna.

## C. elegans (WormBase)

**Sources:**
- [WormBase geneIDs](https://downloads.wormbase.org/releases/current-production-release/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS298.geneIDs.txt.gz) — WBGene IDs with public names, sequence names, biotypes, and live/dead status
- [WormBase geneOtherIDs](https://downloads.wormbase.org/releases/current-production-release/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS298.geneOtherIDs.txt.gz) — other names (aliases) per WBGene

Note: WBGene IDs serve as the primary identifiers (Ensembl-equivalent). File URLs carry a WormBase release number (WS298); bump when a new release ships. As of 2026-06, `downloads.wormbase.org` returns HTTP 403 — update the URL when this is resolved.

## Cynomolgus macaque (Ensembl BioMart)

**Source:** [Ensembl BioMart](https://www.ensembl.org/biomart/martservice?query=) — dataset `mfascicularis_gene_ensembl` (*Macaca fascicularis*)

Provides Ensembl gene IDs (ENSMFAG prefix), gene symbols, synonyms, and gene biotypes via BioMart's `external_gene_name`, `external_synonym`, and `gene_biotype` attributes. No dedicated nomenclature authority exists for this species. Previous symbols are not tracked by BioMart.

## Rhesus macaque (Ensembl BioMart)

**Source:** [Ensembl BioMart](https://www.ensembl.org/biomart/martservice?query=) — dataset `mmulatta_gene_ensembl` (*Macaca mulatta*)

Provides Ensembl gene IDs (ENSMMUG prefix), gene symbols, synonyms, and gene biotypes. No dedicated nomenclature authority exists for this species.

## Common marmoset (Ensembl BioMart)

**Source:** [Ensembl BioMart](https://www.ensembl.org/biomart/martservice?query=) — dataset `cjacchus_gene_ensembl` (*Callithrix jacchus*)

Provides Ensembl gene IDs (ENSCJAG prefix), gene symbols, synonyms, and gene biotypes. No dedicated nomenclature authority exists for this species.

## Mouse lemur (Ensembl BioMart)

**Source:** [Ensembl BioMart](https://www.ensembl.org/biomart/martservice?query=) — dataset `mmurinus_gene_ensembl` (*Microcebus murinus*)

Provides Ensembl gene IDs (ENSMICG prefix), gene symbols, synonyms, and gene biotypes. No dedicated nomenclature authority exists for this species.

## Internal format

Built references are stored as parquet files:

```
references/<species>/
├── gene_table.parquet       # one row per gene
├── symbol_lookup.parquet    # flattened symbol → gene index
└── build_metadata.json      # source URLs, timestamps, checksums
```

### gene_table columns

| Column | Description |
|---|---|
| `ensembl_id` | Ensembl gene ID (nullable for some mouse genes) |
| `symbol` | Approved gene symbol |
| `alias_symbols` | Pipe-delimited alias symbols |
| `prev_symbols` | Pipe-delimited previous symbols |
| `gene_type` | Gene biotype (protein-coding, lncRNA, etc.) |
| `canonical_biotype` | Normalised biotype using a 13-category vocabulary (protein_coding, lncRNA, pseudogene, miRNA, snoRNA, snRNA, rRNA, tRNA, IG_gene, TR_gene, other_ncrna, other, unknown) |
| `status` | Approval status (Approved / Entry Withdrawn) |
| `source` | Reference authority (HGNC / MGI) |
| `source_id` | Authority-specific ID (HGNC:12345 / MGI:12345) |

### symbol_lookup columns

| Column | Description |
|---|---|
| `lookup_string` | The symbol/alias/prev string (original case) |
| `lookup_string_upper` | Uppercased for case-insensitive matching |
| `ensembl_id` | Target Ensembl gene ID (nullable) |
| `source_id` | Target authority ID (always present) |
| `lookup_type` | `approved_symbol`, `alias_symbol`, or `prev_symbol` |
| `source` | Reference authority |

### build_metadata.json

Records exactly what was downloaded and when, for reproducibility:
- Source URLs
- SHA-256 checksums of downloaded files
- Download timestamps
- Row counts

## Custom reference directory

By default, references are stored in a `references/` directory relative to the package. To use a custom location:

```python
build_reference("human", reference_dir="/path/to/my/refs")
result = stangene.run("data.h5ad", species="human", reference_dir="/path/to/my/refs")
```

## Versioning references

The `references/` directory is gitignored by default. If you want to version-control your references (recommended for reproducibility), you can:

1. Commit the parquet files to a separate git repo or GitHub release
2. Or remove `references/` from `.gitignore` and commit them directly
