# Command Line Interface

stangene provides a CLI via `python -m stangene` or the `stangene` command (if installed via pip).

## harmonize

Harmonize gene identifiers in a dataset.

```bash
stangene harmonize --input <file> --species <species> --output-dir <dir> [options]
```

| Argument | Required | Description |
|---|---|---|
| `--input` | Yes | Path to input file (`.h5ad`, `.tsv`, `.csv`) |
| `--species` | Yes | Species name (`human` or `mouse`) |
| `--output-dir` | Yes | Directory to write output files |
| `--dataset-name` | No | Optional dataset label (defaults to filename) |
| `--reference-dir` | No | Custom reference directory |

**Example:**

```bash
stangene harmonize \
    --input data/pbmc3k.h5ad \
    --species human \
    --output-dir results/pbmc3k/ \
    --dataset-name pbmc3k
```

## build-refs

Build reference annotation databases for a species.

```bash
stangene build-refs --species <species> [options]
```

| Argument | Required | Description |
|---|---|---|
| `--species` | Yes | Species name (`human` or `mouse`) |
| `--reference-dir` | No | Custom reference directory |
| `--force` | No | Force re-download and rebuild |

**Examples:**

```bash
# First-time build
stangene build-refs --species human
stangene build-refs --species mouse

# Update to latest references
stangene build-refs --species human --force
```
