# pycyto

Python utilities for cyto - format conversion and sample aggregation.

## Overview

`pycyto` provides command-line tools for working with cyto outputs:

- **convert**: Transform MTX format to h5ad (AnnData)
- **aggregate**: Combine multi-probe cyto outputs into unified sample-level datasets

## Installation

```bash
uv tool install pycyto
pycyto --help
```

## Commands

### convert

Convert Matrix Market (MTX) format from cyto to h5ad for downstream analysis.

```bash
pycyto convert <mtx_directory> <output.h5ad>
```

**Arguments**:
- `mtx_directory`: Path to MTX directory (containing `matrix.mtx`, `features.tsv`, `barcodes.tsv`) or direct path to `.mtx` file
- `output.h5ad`: Output h5ad file path

**Options**:
- `--compress / --no-compress`: Enable gzip compression in h5ad (default: enabled)
- `--integer / --no-integer`: Store counts as int32 instead of float32 (default: float32)

**Examples**:

```bash
# Convert MTX directory to h5ad
pycyto convert cyto_out/counts/BC001.counts.mtx sample.h5ad

# Convert without compression
pycyto convert cyto_out/counts/BC001.counts.mtx sample.h5ad --no-compress

# Store as integers for memory efficiency
pycyto convert cyto_out/counts/BC001.counts.mtx sample.h5ad --integer
```

**Input structure** (from `cyto ibu count --format mtx`):
```
mtx_directory/
├── matrix.mtx      # Sparse count matrix (gene × cell)
├── features.tsv    # Gene/feature names
└── barcodes.tsv    # Cell barcodes
```

**Output**: AnnData object (cell × gene) in CSR sparse format, ready for scanpy/seurat workflows

### aggregate

Aggregate multi-probe cyto outputs across flex barcodes into unified sample-level datasets. Designed specifically for single-modal (GEX) and multi-modal Flex experiments (GEX + CRISPR).

```bash
pycyto aggregate <config.json> <cyto_outdir> <output_directory>
```

**Arguments**:
- `config.json`: JSON configuration specifying sample structure and barcode assignments
- `cyto_outdir`: Directory containing cyto workflow outputs
- `output_directory`: Where to write aggregated files

**Options**:
- `--compress / --no-compress`: Compress output h5ad files (default: no compression)
- `--threads INT`: Number of parallel sample processing threads (default: -1 for all cores)
- `--verbose`: Enable detailed logging

**What it does**:
- Concatenates data across multiple flex probe barcodes per sample
- Merges GEX and CRISPR modalities when both present
- Adds guide assignments to GEX cell metadata
- Filters CRISPR data to match filtered GEX cells (useful if using alternative guide-assignment algorithms)
- Preserves per-cell read/UMI statistics

**Output structure** (per sample):
```
output_directory/
└── sample_name/
    ├── sample_name_gex.h5ad              # Gene expression data
    ├── sample_name_crispr.h5ad           # Guide RNA counts (GEX-filtered)
    ├── sample_name_assignments.parquet   # Guide assignments per cell
    └── sample_name_reads.parquet         # Read/UMI statistics per barcode
```

## Configuration Format

### Basic Structure

Configuration files require:
- `libraries`: Named paths to feature files (probe lists, guide libraries)
- `samples`: Array of sample specifications

### Sample Specification

Each sample requires:
- `experiment`: Experiment identifier (must match cyto output directory names)
- `sample`: Unique sample name (used for output files)
- `mode`: Processing mode (`gex`, `crispr`, or `gex+crispr`)
- `features`: Which library/libraries to use (must match mode with `+`)
- `barcodes`: Flex probe barcode assignments

Matches an output directory of the following path structure:

```text
cyto_output_directory/
└── [experiment_name]_[mode]_Lane*/
    └── ...
```

> Note: All `Lane`s for an `Experiment` + `Mode` will be concatenated.
> If you have differing barcode poolings based on `Lane` you will need to adjust your `Experiment` name to reflect that.


### Barcode Syntax

**Single barcode**:
```json
"barcodes": "BC001"
```

**Barcode range** (expands to BC001, BC002, BC003):
```json
"barcodes": "BC1..3"
```

**Non-contiguous selection**:
```json
"barcodes": "BC001|BC003|BC005"
```

**Combined range and selection**:
```json
"barcodes": "BC1..3|BC005|BC7..9"
```

**Multi-modal pairing** (GEX + CRISPR on same cells):
```json
"mode": "gex+crispr",
"features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
"barcodes": "BC1..8+CR1..8"
```
Pairs BC001+CR001, BC002+CR002, ..., BC008+CR008

**Multiple independent combinations**:
```json
"barcodes": "BC001+CR001|BC002+CR002"
```
Processes BC001+CR001 as one pair, BC002+CR002 as another

### Configuration Examples

#### Example 1: Simple GEX experiment

```json
{
  "libraries": {
    "GEX_PROBES": "./gex_probes.tsv"
  },
  "samples": [
    {
      "experiment": "exp1",
      "mode": "gex",
      "features": "GEX_PROBES",
      "sample": "control",
      "barcodes": "BC1..4"
    },
    {
      "experiment": "exp1",
      "mode": "gex",
      "features": "GEX_PROBES",
      "sample": "treatment",
      "barcodes": "BC5..8"
    }
  ]
}
```

#### Example 2: Multi-modal Perturb-seq

```json
{
  "libraries": {
    "GEX_PROBES": "./gex_probes.tsv",
    "GUIDE_LIBRARY": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "perturbseq_screen",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDE_LIBRARY",
      "sample": "screen_replicate1",
      "barcodes": "BC1..8+CR1..8"
    },
    {
      "experiment": "perturbseq_screen",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDE_LIBRARY",
      "sample": "screen_replicate2",
      "barcodes": "BC9..16+CR9..16"
    }
  ]
}
```

#### Example 3: Developmental timecourse

```json
{
  "libraries": {
    "GEX_PROBES": "./gex_probes.tsv",
    "GUIDES": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "timecourse_20250101",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDES",
      "sample": "day0",
      "barcodes": "BC1..4+CR1..4"
    },
    {
      "experiment": "timecourse_20250101",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDES",
      "sample": "day7",
      "barcodes": "BC5..7+CR5..7"
    },
    {
      "experiment": "timecourse_20250101",
      "mode": "gex+crispr",
      "features": "GEX_PROBES+GUIDES",
      "sample": "day14",
      "barcodes": "BC8..12+CR8..12"
    }
  ]
}
```

## Typical Workflows

### Single Sample Conversion

```bash
# 1. Run cyto workflow
cyto workflow gex -c probes.tsv -w whitelist.txt -o cyto_out sample.vbq

# 2. Convert probe barcode to h5ad
pycyto convert cyto_out/counts/BC001.counts.mtx sample_BC001.h5ad
```

### Multi-Modal Perturb-seq Aggregation

```bash
# 1. Run cyto for GEX and CRISPR
cyto workflow gex -c gex_probes.tsv -w whitelist.txt -p probes.txt -o cyto_out/perturbseq_GEX_Lane1 sample.vbq
cyto workflow crispr -c guides.tsv -w whitelist.txt -p probes.txt -o cyto_out/perturbseq_CRISPR_Lane1 sample.vbq

# 2. Create aggregation config
cat > config.json << 'EOF'
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
      "sample": "mysample",
      "barcodes": "BC1..16+CR1..16"
    }
  ]
}
EOF

# 3. Aggregate (merges GEX + CRISPR modalities)
pycyto aggregate config.json ./cyto_out ./aggr

# Output:
# ./aggr/mysample/
#   ├── mysample_gex.h5ad            # Gene expression with guide annotations
#   ├── mysample_crispr.h5ad         # Guide counts (filtered to GEX cells)
#   ├── mysample_assignments.parquet # Guide assignments
#   └── mysample_reads.parquet       # QC statistics
```

### Multiple Samples from Same Experiment

```bash
# 1. Run cyto workflows
cyto workflow gex -c probes.tsv -w whitelist.txt -p probes.txt -o cyto_out/exp_GEX_Lane1 samples.vbq
cyto workflow crispr -c guides.tsv -w whitelist.txt -p probes.txt -o cyto_out/exp_CRISPR_Lane1 samples.vbq

# 2. Create config assigning barcodes to biological samples
cat > config.json << 'EOF'
{
  "libraries": {
    "GEX": "./gex_probes.tsv",
    "GUIDES": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "exp",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "control_rep1",
      "barcodes": "BC1..2+CR1..2"
    },
    {
      "experiment": "exp",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "control_rep2",
      "barcodes": "BC3..4+CR3..4"
    },
    {
      "experiment": "exp",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "treatment_rep1",
      "barcodes": "BC5..6+CR5..6"
    },
    {
      "experiment": "exp",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "treatment_rep2",
      "barcodes": "BC7..8+CR7..8"
    }
  ]
}
EOF

# 3. Aggregate all samples in parallel
pycyto aggregate config.json ./cyto_out ./aggr
```

## Aggregation Details

### What Gets Aggregated

For each sample, `aggregate` combines data across all specified flex barcodes:

**GEX mode**: Concatenates gene expression matrices  
**CRISPR mode**: Concatenates guide count matrices  
**GEX + CRISPR mode**: 
- Merges guide assignments into GEX `.obs` metadata
- Filters CRISPR data to cells present in filtered GEX data
- Adds read/UMI statistics for both modalities

### Cell Metadata in Aggregated h5ad

The aggregated GEX h5ad includes:
- `experiment`: Experiment identifier
- `sample`: Sample name
- `flex_barcode`: Original flex probe barcode (e.g., "BC001")
- `lane_id`: Sequencing lane identifier
- `assignment`: Assigned guide(s) from CRISPR data
- `moi`: Multiplicity of infection (number of guides per cell)
- `umis`: Guide UMI counts
- `n_reads_gex`: Total GEX reads per cell
- `n_umis_gex`: Total GEX UMIs per cell
- `n_reads_crispr`: Total CRISPR reads per cell
- `n_umis_crispr`: Total CRISPR UMIs per cell

### Barcode Matching

When pairing GEX and CRISPR data:
- CRISPR barcodes (CR) are automatically converted to match GEX format (BC) for cell matching
- Cells are matched on: `cell_barcode + flex_barcode + lane_id`
- Only cells present in filtered GEX data are retained in CRISPR output

### Cyto Output Directory Structure

`aggregate` expects cyto outputs organized as:

```
cyto_outdir/
├── {experiment}_GEX_Lane*/
│   └── counts/
│       ├── BC001.h5ad
│       ├── BC002.h5ad
│       └── ...
└── {experiment}_CRISPR_Lane*/
    ├── counts/
    │   ├── CR001.h5ad
    │   └── ...
    └── assignments/
        ├── CR001.assignments.tsv
        └── ...
```

Where `{experiment}` matches the `experiment` field in your config.

## Advanced Usage

### Processing Subset of Barcodes

Process only specific barcode combinations by modifying the config:

```json
{
  "samples": [
    {
      "experiment": "exp1",
      "mode": "gex+crispr",
      "features": "GEX+GUIDES",
      "sample": "high_quality_subset",
      "barcodes": "BC001+CR001|BC003+CR003|BC007+CR007"
    }
  ]
}
```

### Different Feature Libraries per Sample

Reference different probe/guide libraries:

```json
{
  "libraries": {
    "TISSUE_PANEL": "./tissue_probes.tsv",
    "IMMUNE_PANEL": "./immune_probes.tsv",
    "SCREEN_GUIDES": "./guides.tsv"
  },
  "samples": [
    {
      "experiment": "exp1",
      "mode": "gex+crispr",
      "features": "TISSUE_PANEL+SCREEN_GUIDES",
      "sample": "tissue_sample",
      "barcodes": "BC1..8+CR1..8"
    },
    {
      "experiment": "exp1",
      "mode": "gex+crispr",
      "features": "IMMUNE_PANEL+SCREEN_GUIDES",
      "sample": "pbmc_sample",
      "barcodes": "BC9..16+CR9..16"
    }
  ]
}
```

### Parallel Processing Control

```bash
# Use all available cores
pycyto aggregate config.json cyto_out output --threads -1

# Limit to 8 samples processed simultaneously
pycyto aggregate config.json cyto_out output --threads 8

# Single-threaded (minimal memory)
pycyto aggregate config.json cyto_out output --threads 1
```

## Troubleshooting

### Missing Files Error

**Problem**: "Expected Feature file does not exist"  
**Solution**: Ensure MTX directory contains `matrix.mtx`, `features.tsv`, and `barcodes.tsv`

### Barcode Format Errors

**Problem**: "Invalid barcode format"  
**Solution**: Barcodes must be BC/CR/AB followed by numbers. Use range syntax: `BC1..8` not `BC1-8`

### No Data Found

**Problem**: "No data found to process for sample"  
**Solution**: 
- Verify experiment name in config matches cyto output directory prefix
- Check that barcode specifications match available probe barcodes in cyto output
- Ensure cyto outputs are in expected directory structure

### Memory Issues

**Problem**: Out of memory during aggregation  
**Solution**: 
- Reduce `--threads` to process fewer samples concurrently
- Use `--compress` to reduce output file sizes
- Process samples in smaller batches with separate configs

## Performance Notes

- **Parallel aggregation**: Samples are processed independently in parallel (one per thread)
- **Lazy loading**: Uses anndata experimental lazy loading to minimize memory overhead
- **Compression**: Optional compression reduces h5ad file sizes ~40-50%
- **Integer storage**: `--integer` flag in `convert` reduces memory vs float32

## See Also

- **cyto**: Main processing pipeline - https://github.com/arcinstitute/cyto
- **anndata**: Annotated data structures - https://anndata.readthedocs.io

## Citation

If you use pycyto in your research, please cite:

```
Teyssier, N. and Dobin, A. (2025). cyto: ultra-high throughput processing 
of 10x-flex single cell sequencing. bioRxiv.
```

## Support

- **Issues**: https://github.com/arcinstitute/pycyto/issues
- **Documentation**: Run `pycyto --help` or `pycyto <command> --help`
