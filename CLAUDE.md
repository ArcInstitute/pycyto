# CLAUDE.md

This file provides guidance to Claude Code when working with pycyto code.

## Project Overview

pycyto is a Python companion tool for cyto, providing utilities to work with cyto outputs:

- **convert**: Transform Matrix Market (MTX) format to h5ad (AnnData) for downstream analysis
- **aggregate**: Combine multi-probe cyto outputs into unified sample-level datasets

pycyto is designed specifically for 10x Genomics Flex workflows where a single biological sample is split across multiple flex probe barcodes (BC001-BC016 for GEX, CR001-CR016 for CRISPR). The main challenge is intelligently aggregating data across these technical replicates while handling multi-modal experiments (GEX + CRISPR).

## Installation

```bash
uv tool install pycyto
```

Or from source:

```bash
cd pycyto
uv pip install -e .
```

## Architecture

### Module Organization

- `__main__.py` - CLI entry point using Typer (convert, aggregate commands)
- `config.py` - Configuration parsing and barcode expansion DSL
- `aggregate.py` - Multi-modal sample aggregation logic
- `convert.py` - Simple MTX to h5ad conversion utilities

### Key Dependencies

- `polars` - Fast dataframe operations for configuration handling
- `anndata` - Single-cell data structures
- `anndata.experimental` - Lazy loading for memory efficiency
- `pandas` - Legacy support for anndata `.obs` manipulation

## The `config` Module

**Purpose**: Parse JSON configuration files that specify how to map flex probe barcodes to biological samples.

**Key challenge**: Users need a concise way to specify barcode ranges and pairings without listing hundreds of individual barcodes.

### Supported Barcode Formats

Pycyto supports two barcode naming conventions:

**Flex-V1 (16-plex)**:

- **Format**: `BC001`-`BC016`, `CR001`-`CR016`, `AB001`-`AB016`
- **Use case**: Standard 16-plex multiplexing for GEX, CRISPR, or antibody capture

**Flex-V2 (384-plex)**:

- **Format**: `[ABCD]-[ABCDEFGH][01-12]` or `[ABCD]_[ABCDEFGH][01-12]`
  - Examples: `A-A01`, `B-C05`, `D-H12` or `A_A01`, `B_C05`, `D_H12`
  - Both hyphen and underscore separators are supported
- **Structure**: 4 sets × 8 rows × 12 columns = 384 unique barcodes
- **Sets**: A, B, C, D
- **Rows**: A-H (8 rows per set, like a 96-well plate)
- **Columns**: 01-12 (12 columns per row)
- **Use case**: High-throughput 384-plex multiplexing
- **Note**: No BC/CR dual naming for Flex-V2 (single naming scheme per barcode set)

The barcode format is automatically detected from the naming convention - no config flag needed.

### Barcode Expansion DSL

The config parser implements a mini-language for specifying barcodes:

**Flex-V1 Range syntax**:

```json
"BC1..8"  → expands to BC001, BC002, BC003, ..., BC008
```

**Flex-V2 Range syntax**:

```json
"A-A01..A-A12"  → expands to A-A01, A-A02, ..., A-A12
"A-A01..A-H12"  → expands to all 96 barcodes in set A (full plate)
```

**Selection syntax** (both V1 and V2):

```json
"BC1|BC3|BC5"  → expands to BC001, BC003, BC005
"A-A01|A-A05|A-C03"  → expands to A-A01, A-A05, A-C03
```

**Combined syntax**:

```json
"BC1..4|BC7|BC9..12"  → expands to BC001, BC002, BC003, BC004, BC007, BC009, BC010, BC011, BC012
"A-A01..A-A12|B-B01..B-B12"  → expands to ranges from sets A and B
```

**Pairing syntax** (for multi-modal experiments, Flex-V1 only):

```json
"BC1..8+CR1..8"  → pairs BC001+CR001, BC002+CR002, ..., BC008+CR008
```

**Multiple independent pairings**:

```json
"BC001+CR001|BC002+CR002"  → two separate pairings
```

**Flex-V2 Multi-Modal**: Since Flex-V2 uses a single naming scheme (no BC/CR dual naming), no `+` pairing is needed:

```json
"A-A01..A-H12"  → uses same barcodes for both GEX and CRISPR
```

### Ambiguity Resolution

The parser must handle ambiguity with the `|` operator:

- Does `BC1..4|BC5+CR1..4|CR5` mean two pairings OR ranges within each component?
- Solution: Check if each `|`-separated part has the expected number of `+` separators
  - If all parts are complete pairings (correct number of `+`), treat as multiple combinations
  - Otherwise, treat `|` as selection within components

### Important Functions

- `parse_config(config_path)` - Main entry point, returns polars DataFrame with expanded barcodes
- `_parse_barcodes(entry, nlib)` - Handles the barcode expansion DSL
- `_expand_barcode_component(component)` - Expands a single component like `BC1..8`
- `_expand_range(range_str)` - Expands `5..7` to `[5, 6, 7]`
- `_expand_selection(selection)` - Expands `1|3|5..7` to `[1, 3, 5, 6, 7]`

### Output Schema

Returns a polars DataFrame where each row represents a single flex barcode assignment:

```
experiment | sample | mode | bc_component | bc_idx | features | probe_set | feature_path | expected_prefix
```

**Critical columns**:

- `bc_component`: The specific flex barcode (BC001, CR001, etc.)
- `bc_idx`: Index within the pairing (for tracking which barcodes go together)
- `expected_prefix`: Used for directory matching (e.g., "experiment_GEX_Lane")

## The `aggregate` Module

**Purpose**: Aggregate cyto outputs across multiple flex barcodes and lanes into unified sample-level datasets.

**Key challenge**: Cyto outputs are organized by `experiment_MODE_Lane#` directories, with one file per flex barcode. For multi-modal experiments (GEX + CRISPR), data must be merged at the cell level.

### High-Level Workflow

1. **Directory Discovery**: For each sample, find all matching cyto output directories
   - Uses regex to match `{experiment}_{MODE}_Lane*` patterns
   - Aggregates across all lanes for an experiment+mode combination

2. **Data Loading**: Load data for all barcodes assigned to a sample
   - GEX: Filtered h5ad files (`BC001.filt.h5ad`)
   - CRISPR: Unfiltered h5ad files (`CR001.h5ad`) + assignment files (`CR001.assignments.tsv`)
   - Reads: Per-barcode read statistics (`BC001.reads.tsv.zst`)

3. **Multi-Modal Merging** (GEX + CRISPR case):
   - Concatenate all GEX data across barcodes
   - Concatenate all CRISPR data and assignments
   - **Convert CRISPR barcodes**: CR → BC for cell matching
   - Merge guide assignments into GEX `.obs` metadata
   - Filter CRISPR data to only cells present in filtered GEX data

4. **Output Generation**: Write sample-level files
   - `{sample}_gex.h5ad` - Gene expression with guide annotations
   - `{sample}_crispr.h5ad` - Guide counts (GEX-filtered)
   - `{sample}_assignments.parquet` - Guide assignments per cell
   - `{sample}_reads.parquet` - Read/UMI statistics

### Barcode Format Conversion

**The problem (Flex-V1 only)**: CRISPR libraries can use either BC or CR flex barcode prefixes depending on the experimental design. When pairing GEX (which always uses BC) with CRISPR data that uses CR prefixes, the cell barcodes need to be matched correctly.

**Why separate prefixes (Flex-V1)**: CRISPR and GEX can be run with different flex barcode sets:

- Same barcodes: `BC1..8+BC9..16` (CRISPR uses BC)
- Different barcodes: `BC1..8+CR1..8` (CRISPR uses CR)

**The solution**:

1. Detect barcode format (Flex-V1 vs Flex-V2) from GEX data
2. **For Flex-V1**:
   - Detect if CRISPR uses CR prefixes by checking assignment data
   - If CR detected: cell barcodes like `ACGTACGT-CR001-1` are converted to `ACGTACGT-BC001-1`
   - If BC detected: no conversion needed, barcodes already match
3. **For Flex-V2**:
   - No conversion needed - single naming scheme per barcode set
   - Both GEX and CRISPR use the same barcode identifiers (e.g., `A-A01`)
4. Matching is always done on `cell_barcode + flex_barcode + lane_id`

See `_process_gex_crispr_set()` around line 140-160 for the detection and conversion logic.

### Memory Management

The aggregation process is designed to minimize memory usage:

- Uses `anndata.experimental.read_lazy()` to avoid loading full matrices into memory
- Explicit `del` statements after concatenation to free memory immediately
- Parallel processing at the sample level (not barcode level) to avoid memory explosion

### Important Functions

- `aggregate_data(config, cyto_outdir, outdir, ...)` - Main entry point, orchestrates parallel processing
- `process_sample(sample, config, ...)` - Processes a single sample (runs in parallel)
- `_process_gex_crispr_set(...)` - Handles GEX + CRISPR merging and filtering
- `_filter_crispr_adata_to_gex_barcodes(...)` - Filters CRISPR to GEX cells
- `_load_*_for_experiment_sample(...)` - Load specific data types from cyto outputs

### Parallel Processing

Samples are processed in parallel using multiprocessing with `spawn` context:

- One process per sample (controlled by `--threads`)
- Each worker initializes its own logger
- No shared state between workers

## Configuration Format

Configurations are JSON files with two main sections:

### Libraries Section

Maps logical names to feature file paths:

```json
{
  "libraries": {
    "GEX_PROBES": "./gex_probes.tsv",
    "GUIDE_LIBRARY": "./guides.tsv"
  }
}
```

### Samples Section

Specifies how to aggregate barcodes into samples:

```json
{
  "samples": [
    {
      "experiment": "exp1",
      "sample": "sample_name",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDE_LIBRARY",
      "barcodes": "BC1..8+CR1..8"
    }
  ]
}
```

**Key fields**:

- `experiment`: Must match cyto output directory prefix
- `mode`: `gex`, `crispr`, or `gex+crispr`
- `features`: Library names from `libraries` section (use `+` for multi-modal)
- `barcodes`: Barcode specification using the DSL (use `+` for pairing)

## Expected Cyto Output Structure

The aggregate command expects cyto outputs organized as:

```
cyto_outdir/
├── {experiment}_GEX_Lane1/
│   ├── counts/
│   │   ├── BC001.filt.h5ad
│   │   └── BC002.filt.h5ad
│   └── stats/reads/
│       ├── BC001.reads.tsv.zst
│       └── BC002.reads.tsv.zst
└── {experiment}_CRISPR_Lane1/
    ├── counts/
    │   ├── CR001.h5ad
    │   └── CR002.h5ad
    └── assignments/
        ├── CR001.assignments.tsv
        └── CR002.assignments.tsv
```

## Common Patterns

### Single-Modal GEX Experiment (Flex-V1)

```json
{
  "libraries": { "GEX": "./gex_probes.tsv" },
  "samples": [
    {
      "experiment": "exp1",
      "mode": "gex",
      "features": "GEX",
      "sample": "control",
      "barcodes": "BC1..4"
    }
  ]
}
```

### Multi-Modal Perturb-seq (Flex-V1)

```json
{
  "libraries": {
    "GEX": "./gex_probes.tsv",
    "GUIDES": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "perturbseq",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "screen_rep1",
      "barcodes": "BC1..8+CR1..8"
    }
  ]
}
```

### Single-Modal GEX Experiment (Flex-V2, 384-plex)

```json
{
  "libraries": { "GEX": "./gex_probes.tsv" },
  "samples": [
    {
      "experiment": "exp1",
      "mode": "gex",
      "features": "GEX",
      "sample": "poolA_sample1",
      "barcodes": "A-A01..A-H12"
    },
    {
      "experiment": "exp1",
      "mode": "gex",
      "features": "GEX",
      "sample": "poolB_sample1",
      "barcodes": "B-A01..B-H12"
    }
  ]
}
```

### Multi-Modal Perturb-seq (Flex-V2)

```json
{
  "libraries": {
    "GEX": "./gex_probes.tsv",
    "GUIDES": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "perturbseq",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "screen_poolA",
      "barcodes": "A-A01..A-H12"
    }
  ]
}
```

**Note**: Flex-V2 uses a single naming scheme, so no `+` pairing is needed for multi-modal experiments.

## Testing

Run tests with:

```bash
pytest tests/
```

Example configurations are available in `examples/`:

- `aggregation.json` - Simple multi-modal example with explicit pairings
- `ngn2_agg.json` - Developmental timecourse with range syntax

## See Also

- **cyto**: Main processing pipeline - https://github.com/arcinstitute/cyto
- **anndata**: Data structures - https://anndata.readthedocs.io
