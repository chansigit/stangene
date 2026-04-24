# Adding a New Species

stangene is designed to be extensible. Adding support for a new model organism requires three steps.

## Step 1: Add a SpeciesConfig

In `src/stangene/species.py`, add an entry to `_SPECIES_CONFIGS`:

```python
_SPECIES_CONFIGS["rat"] = SpeciesConfig(
    name="rat",
    ensembl_prefix="ENSRNOG",
    transcript_prefix="ENSRNOT",
    naming_convention="capitalized",  # or "uppercase"
    reference_sources={
        "rgd": {
            "url": "https://download.rgd.mcw.edu/data_release/GENES_RAT.txt",
            "description": "RGD rat gene annotations",
        },
    },
)
```

If the new species has unique feature classification patterns (e.g., a new Ensembl prefix), add them to `CLASSIFICATION_PATTERNS`:

```python
CLASSIFICATION_PATTERNS.insert(0, (re.compile(r"^ENSRNOG\d+"), "gene"))
CLASSIFICATION_PATTERNS.insert(0, (re.compile(r"^ENSRNOT\d+"), "transcript"))
```

## Step 2: Implement a reference builder

In `src/stangene/references.py`, add a `_build_rat_reference()` function:

```python
def _build_rat_reference(config, ref_dir: str) -> None:
    """Build rat reference from RGD gene data."""
    url = config.reference_sources["rgd"]["url"]
    raw_data = _download_file(url)
    checksum = hashlib.sha256(raw_data).hexdigest()

    # Parse the source file into a gene_table DataFrame with columns:
    # ensembl_id, symbol, alias_symbols, prev_symbols, gene_type, status, source, source_id
    gene_table = ...  # your parsing logic here

    # Build the flattened symbol lookup index
    symbol_lookup = _build_symbol_lookup(gene_table, source="RGD")

    # Save metadata
    metadata = {
        "species": "rat",
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "sources": {"rgd": {"url": url, "sha256": checksum, "rows": len(gene_table)}},
        "gene_count": len(gene_table),
        "lookup_count": len(symbol_lookup),
    }

    _save_reference(ref_dir, gene_table, symbol_lookup, metadata)
```

## Step 3: Register the builder

In the `build_reference()` function, add the dispatch:

```python
elif config.name == "rat":
    _build_rat_reference(config, ref_dir)
```

## Key requirements

- The `gene_table` must have columns: `ensembl_id`, `symbol`, `alias_symbols`, `prev_symbols`, `gene_type`, `status`, `source`, `source_id`
- `ensembl_id` can be null for genes without Ensembl annotation; `source_id` serves as fallback
- `alias_symbols` and `prev_symbols` are pipe-delimited strings
- `_build_symbol_lookup()` handles the flattening automatically

## Testing

Add test fixtures with mock data for your species in `tests/test_references.py`, following the pattern of the existing human (HGNC) and mouse (MGI) tests.

## Tier 2: Ensembl BioMart (no dedicated authority)

For species without a curated nomenclature authority (e.g., most non-human primates), use the shared `_build_ensembl_biomart_reference()` instead of writing a custom parser.

```python
"my_species": SpeciesConfig(
    name="my_species",
    ensembl_prefix="ENS...G",
    transcript_prefix="ENS...T",
    naming_convention="uppercase",
    reference_sources={
        "ensembl_biomart": {
            "url": "https://www.ensembl.org/biomart/martservice?query=",
            "dataset": "<species>_gene_ensembl",
            "description": "Ensembl BioMart gene table for <species>",
        },
    },
),
```

Then in `build_reference()`, add your species to the tuple:

```python
elif config.name in ("cynomolgus", "rhesus", "marmoset", "mouse_lemur", "my_species"):
    _build_ensembl_biomart_reference(config, ref_dir)
```

BioMart returns no `previous_symbols` and no `status` — the builder hardcodes `status="approved"` and `prev_symbols=""`. Genes with no symbol are still included in the gene_table (so Tier 1 exact_id matching works for raw Ensembl IDs), but skipped from `symbol_lookup`.

See `_build_ensembl_biomart_reference` in `src/stangene/references.py` for the full implementation.
