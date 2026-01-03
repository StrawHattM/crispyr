
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crispyr: CRISPR Screen FastQ Processing

<!-- badges: start -->
<!-- badges: end -->

## Overview

`crispyr` is an R package for processing multiplexed FASTQ files from pooled CRISPR/Cas9 screens into count matrices suitable for downstream statistical analysis. It implements barcode extraction, construct library matching, gene mapping, and quality control metrics similar to the Broad Institute's PoolQ tool.

### Features

- **Flexible barcode extraction** with configurable extraction policies (default: Broad Institute format)
- **Dual barcode demultiplexing**: Construct barcodes (20bp guides) + Sample barcodes (8bp)
- **Library matching** with exact or approximate matching (allows 0-2 mismatches)
- **Gene annotation** using chip files with 1-to-many barcode-to-gene mappings
- **Count matrices** at construct and gene levels
- **Log-normalization** using PoolQ method (depth-corrected)
- **Comprehensive QC metrics** including coverage, depth, and match rates
- **Complete output files** in TSV format for downstream analysis

## Installation

Install the development version from GitHub:

```r
# Using pak (recommended)
# install.packages("pak")
pak::pak("StrawHattM/crispyr")

# Or using devtools
# install.packages("devtools")
devtools::install_github("StrawHattM/crispyr")
```

## Quick Start

### Basic Workflow

```r
library(crispyr)

# Process a pooled CRISPR screen
result <- process_pooled_screen(
  fastq_files = "sample.fastq.gz",
  reference_library_csv = "CP0045_reference.csv",
  chip_file = "CP0045_GRCm38_strict.chip",
  sample_manifest_csv = "samples.csv",
  output_dir = "./results"
)

# Access results
construct_counts <- result$construct_counts$matrix
gene_counts <- result$gene_counts$matrix
qc_metrics <- result$qc_metrics
```

### Step-by-Step Pipeline

```r
library(crispyr)

# 1. Load reference data
ref_lib <- load_reference_library("reference.csv")
chip <- load_chip_file("annotation.chip")
manifest <- load_sample_manifest("samples.csv")

# 2. Extract barcodes from FASTQ
extracted <- extract_barcodes_from_fastq(
  fastq_files = "sample.fastq.gz",
  barcode_policy = "PREFIX:CACCG@6"
)

# 3. Match to reference library
matched <- match_constructs(extracted, ref_lib)

# 4. Map to genes
with_genes <- map_constructs_to_genes(matched, chip)

# 5. Aggregate by sample
aggregated <- aggregate_by_sample_barcode(with_genes, manifest)

# 6. Create count matrices
construct_counts <- create_count_matrix(aggregated, level = "construct")
gene_counts <- create_count_matrix(aggregated, level = "gene")

# 7. Normalize
lognorm <- lognormalize_counts(construct_counts$matrix)

# 8. Quality control
qc <- calculate_qc_metrics(extracted, matched, aggregated)
```

## Key Functions

| Function | Description |
|----------|-------------|
| `process_pooled_screen()` | Main wrapper - complete pipeline in one call |
| `extract_barcodes_from_fastq()` | Extract construct and sample barcodes from reads |
| `match_constructs()` | Map barcodes to reference library (exact or approximate) |
| `map_constructs_to_genes()` | Assign barcodes to target genes |
| `create_count_matrix()` | Generate count matrices at construct or gene level |
| `lognormalize_counts()` | Apply depth-normalization using PoolQ method |
| `calculate_qc_metrics()` | Compute comprehensive quality control statistics |

## Example Data

The package includes example data from the Broad Institute's CP0045 BRIE AiO library:

```r
# These files are available in the example_data folder:
# - CP0045_reference_20160120.csv       # Barcode reference
# - CP0045_GRCm38_*_strict.chip         # Gene annotations
# - GP P-6350_conditions.csv            # Sample manifest
# - example.fastq.gz                    # Sample FASTQ data
```

## Documentation

Complete documentation is available through R's help system:

```r
# Get help on specific functions
?process_pooled_screen
?extract_barcodes_from_fastq
?match_constructs

# Read the vignette
vignette("crispyr-workflow")
```

## Requirements

### Dependencies

- `ShortRead` - FASTQ file I/O
- `Biostrings` - DNA sequence handling  
- `data.table` - Fast data manipulation
- `dplyr` / `tidyr` - Data wrangling
- `stringr` - String processing
- `stringdist` - Approximate string matching (for mismatches)

### System Requirements

- R ≥ 4.0
- ~2GB RAM for typical CRISPR screens (500M reads)
- Gzip-compressed FASTQ files recommended

## Output Files

When saving results with `output_dir` parameter, the following files are created:

- `*_construct_counts.tsv` - Raw construct barcode counts
- `*_gene_counts.tsv` - Gene-aggregated counts
- `*_construct_lognorm.tsv` - Log-normalized construct counts
- `*_gene_lognorm.tsv` - Log-normalized gene counts
- `*_qc_overall.txt` - Overall quality metrics
- `*_qc_per_sample.tsv` - Per-sample quality metrics

## Advanced Usage

### Custom Barcode Policies

Define custom barcode extraction for different library designs:

```r
# Broad Institute default (most common)
result <- process_pooled_screen(
  ...,
  barcode_policy = "PREFIX:CACCG@6"
)

# Custom flanking sequence
result <- process_pooled_screen(
  ...,
  barcode_policy = "PREFIX:GCCTG@5"
)

# Policy format: "PREFIX:SEQUENCE@POSITION"
# - Find SEQUENCE in read
# - Extract 20bp starting at POSITION after the sequence
```

### Approximate Matching

Allow sequencing errors in barcode matching:

```r
result <- process_pooled_screen(
  ...,
  max_mismatches = 1  # Allow 1bp mismatch
)
```

## Validation

The implementation is validated against the Broad Institute's PoolQ tool. Output matrices should be highly similar to PoolQ results (within expected rounding differences).

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use crispyr in your research, please cite:

```
González Fernández, M. (2024). crispyr: CRISPR Screen FastQ Processing in R. 
GitHub: https://github.com/StrawHattM/crispyr
```

## License

This package is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

- Broad GPP PoolQ: https://portals.broadinstitute.org/gpp/
- Hart et al. (2015) "High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities." *Nature Methods* 12(5):494-498
- Original PoolQ publication details available on Broad website

## Contact

For questions, issues, or suggestions, please open an issue on GitHub:
https://github.com/StrawHattM/crispyr/issues

