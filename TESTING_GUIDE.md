# Testing Guide for crispyr Package

## Quick Start - Test Scripts Ready to Use

Two example scripts are available in the repository root:

### 1. Simple Pipeline (`01_pipeline_wrapper.R`)
- **Size:** 1.7 KB
- **Purpose:** Run the complete pipeline with a single function call
- **Best for:** Quick testing and validation
- **Time complexity:** Linear - same for all dataset sizes
- **Usage:**
  ```r
  source("01_pipeline_wrapper.R")
  ```

### 2. Step-by-Step Pipeline (`02_pipeline_stepwise.R`)
- **Size:** 5.0 KB
- **Purpose:** Run each pipeline step individually with detailed timing
- **Best for:** Understanding the pipeline and debugging
- **Time complexity:** Same total time, but broken down by step
- **Usage:**
  ```r
  source("02_pipeline_stepwise.R")
  ```

## Performance Expectations

All timing estimates are for **standard consumer hardware:**
- SSD storage (500 MB/s read/write)
- 4-8 core CPU
- 8+ GB RAM

### Dataset Size vs Runtime

| Dataset Size | Expected Runtime | Step Duration |
|--------------|-----------------|---------|
| 1-5M reads | 30-60 sec | ~15 sec/1M reads |
| 10-20M reads | 2-4 min | ~12 sec/1M reads |
| 50M reads | 8-12 min | ~10 sec/1M reads |
| 100M+ reads | 20-40 min | ~10 sec/1M reads |

### Hardware Variations

**SSD vs HDD:** ±50% variation (SSD strongly recommended)
- SSD: baseline performance
- HDD: 1.5-2x slower

**CPU Cores:** Minor impact (<10% variation)
- Current implementation is single-threaded
- Multi-threading could be added in future versions

**RAM:** Minimal impact
- Typically uses <2GB even for 100M reads
- Larger screens benefit from available cache

## Expected Timing Breakdown (per 1 million reads)

```
Barcode extraction:     15-20 seconds
Library matching:        5-10 seconds
Gene mapping:            2-3 seconds
Count matrix creation:   2-3 seconds
Normalization:           <1 second
QC metrics:              2-3 seconds
File I/O:                5-10 seconds
                        ───────────
TOTAL:                  32-50 seconds per 1M reads
```

## How to Test

### Test 1: Load and Verify Installation

```r
# Load the package
devtools::load_all()

# Check that all functions are available
ls("package:crispyr")

# Should show 14+ exported functions
```

**Expected output:** List of all exported functions

---

### Test 2: Run Simple Example (30-60 seconds)

```r
# Edit paths in 01_pipeline_wrapper.R
source("01_pipeline_wrapper.R")

# This will:
# - Load reference data
# - Extract barcodes from FASTQ files
# - Match constructs
# - Create count matrices
# - Save results to output directory
```

**Expected output:**
- Progress messages showing each pipeline stage
- Timing information at the end
- Count matrices printed to console
- Files saved to `output_dir`

---

### Test 3: Run Step-by-Step Example (30-60 seconds)

```r
# Edit paths in 02_pipeline_stepwise.R
source("02_pipeline_stepwise.R")

# This will:
# - Run each pipeline step individually
# - Print detailed timing for each step
# - Generate detailed timing report
```

**Expected output:**
- Detailed progress for each step
- Timing table showing:
  - Step name
  - Duration in seconds
  - Cumulative time
- Summary statistics
- All results saved to `output_dir`

---

### Test 4: Validate Results (5 minutes)

After running either script, validate the outputs:

```r
# Load results
construct_counts <- read.csv("construct_counts_raw.csv", row.names = 1)
gene_counts <- read.csv("gene_counts_raw.csv", row.names = 1)
qc_metrics <- read.csv("qc_metrics_overall.csv")

# Check dimensions
dim(construct_counts)  # Should be constructs x samples
dim(gene_counts)       # Should be genes x samples

# Check for NAs
sum(is.na(construct_counts))  # Should be 0
sum(is.na(gene_counts))       # Should be 0

# Check count distributions
summary(construct_counts[, 1])
hist(log2(construct_counts[, 1] + 1))

# Review QC metrics
print(qc_metrics)
```

---

### Test 5: Benchmark Your Hardware (Optional)

To measure actual performance on your system:

```r
# Modify 02_pipeline_stepwise.R to print timing at each step
# Run the script and note the total time

# Compare to expected:
expected_rate <- 10  # seconds per 1M reads
actual_reads <- 2000000  # your dataset size
expected_time <- (actual_reads / 1000000) * expected_rate

# Your system multiplier:
actual_multiplier <- actual_time / expected_time
# < 1.0 = faster than baseline
# > 1.0 = slower than baseline
```

---

## Troubleshooting

### Issue: "Some FASTQ files not found"
**Solution:**
- Check file paths are correct
- Ensure files exist and are readable
- Verify file extensions (.fastq.gz or .fastq)

### Issue: Low extraction rates (<80%)
**Possible causes:**
- Barcode policy doesn't match your library
- FASTQ reads don't contain flanking sequence
- Quality issues in the sequencing

**Solutions:**
- Check that `barcode_policy` matches your construct layout
- Verify flanking sequence is present in reads:
  ```r
  # Check first few reads
  fq <- ShortRead::readFastq("sample.fastq.gz")
  head(as.character(Biostrings::sread(fq)))
  ```

### Issue: Low construct matching rates (<90%)
**Possible causes:**
- Reference library incomplete
- Systematic sequencing errors
- Barcode degradation

**Solutions:**
- Try increasing `max_mismatches` (0, 1, or 2)
- Verify reference library matches sequenced barcodes
- Check for systematic errors in your FASTQ data

### Issue: Memory errors with large datasets
**Solutions:**
- Process files in batches
- Close other applications to free RAM
- Consider downsampling FASTQ files for testing

---

## Comparing to PoolQ

If you have PoolQ results to compare:

```r
# Load both crispyr and PoolQ results
crispyr_genes <- read.csv("gene_counts_raw.csv", row.names = 1)
poolq_genes <- read.csv("poolq_results.csv", row.names = 1)

# Compare dimensions
dim(crispyr_genes)
dim(poolq_genes)

# Correlation analysis
common_genes <- intersect(rownames(crispyr_genes), rownames(poolq_genes))
common_samples <- intersect(colnames(crispyr_genes), colnames(poolq_genes))

if(length(common_genes) > 0 & length(common_samples) > 0) {
  # Subset to common features
  crispyr_sub <- crispyr_genes[common_genes, common_samples]
  poolq_sub <- poolq_genes[common_genes, common_samples]
  
  # Calculate correlation
  cor_matrix <- cor(as.vector(crispyr_sub), as.vector(poolq_sub))
  print(paste("Pearson correlation:", round(cor_matrix, 3)))
  
  # Plot comparison
  plot(as.vector(crispyr_sub) + 1, as.vector(poolq_sub) + 1, 
       log = "xy", xlab = "crispyr", ylab = "PoolQ")
  abline(0, 1, col = "red")
}
```

---

## Success Criteria

Your testing is successful if:

✓ Both example scripts run without errors
✓ Output files are created
✓ Count matrices have expected dimensions (constructs/genes × samples)
✓ QC metrics are generated
✓ Runtime matches expected estimates (±50%)
✓ Results are biologically reasonable
✓ Correlation with PoolQ >0.95 (if comparing)

---

## Next Steps

1. **Test with your CP0045 data**
   - Edit file paths in example scripts
   - Run 01_pipeline_wrapper.R first
   - Then try 02_pipeline_stepwise.R

2. **Customize for your screen**
   - Adjust `barcode_policy` if needed
   - Modify `max_mismatches` for your quality level
   - Set `min_quality_score` if filtering low-quality reads

3. **Validate biological results**
   - Compare with PoolQ or other tools
   - Check expected genes are present
   - Verify condition effects are reasonable

4. **Scale up to full dataset**
   - Test with subsampled FASTQ files first
   - Monitor timing as you increase dataset size
   - Adjust parameters based on real performance

