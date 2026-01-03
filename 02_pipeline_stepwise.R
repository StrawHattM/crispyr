# ====================================================================
# CRISPYR: Complete CRISPR Screen Processing Pipeline (Step-by-Step)
# ====================================================================
# This script runs each pipeline step individually with timing.
#
# TIMING BREAKDOWN (for ~2M reads):
# - Load reference: 5-10 sec
# - Extract barcodes: 20-30 sec
# - Match barcodes: 10-15 sec
# - Aggregate samples: 2-3 sec
# - Count matrices: 2-3 sec
# - Normalize: <1 sec
# - QC metrics: 2-3 sec
# - Save: 1-2 sec
# TOTAL: ~45-70 seconds
#
# For 50M reads: expect 8-12 minutes
# For 100M+ reads: expect 20-40 minutes
# ====================================================================

library(crispyr)
library(data.table)

# UPDATE THESE PATHS
data_dir <- "path/to/your/CP0045/data"
fastq_files <- file.path(data_dir, c("sample1.fastq.gz", "sample2.fastq.gz"))
reference_csv <- file.path(data_dir, "CP0045_reference.csv")
chip_file <- file.path(data_dir, "CP0045_GRCm38_NCBI_CRISPRko_strict_gene.chip")
sample_manifest_csv <- file.path(data_dir, "samples.csv")
output_dir <- file.path(data_dir, "crispyr_results_stepwise")
dir.create(output_dir, showWarnings = FALSE)

cat("\n========== CRISPYR Step-by-Step Pipeline ==========\n\n")

# STEP 1: Load reference data
cat("Step 1: Loading reference data...\n")
t1 <- Sys.time()
ref_library <- load_reference_library(reference_csv, quiet = TRUE)
chip <- load_chip_file(chip_file, quiet = TRUE)
sample_manifest <- load_sample_manifest(sample_manifest_csv, quiet = TRUE)
cat(sprintf("  Loaded %d constructs, %d gene mappings, %d samples\n", 
            nrow(ref_library), nrow(chip), nrow(sample_manifest)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t1, units="secs")))

# STEP 2: Extract barcodes
cat("Step 2: Extracting barcodes...\n")
t2 <- Sys.time()
extracted <- extract_barcodes_from_fastq(fastq_files, quiet = TRUE)
cat(sprintf("  Extracted from %d reads\n", nrow(extracted)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t2, units="secs")))

# STEP 3: Match barcodes
cat("Step 3: Matching to library...\n")
t3 <- Sys.time()
matched <- match_constructs(extracted, ref_library, max_mismatches = 1, quiet = TRUE)
mapped <- map_constructs_to_genes(matched, chip, quiet = TRUE)
cat(sprintf("  Matched %.1f%% of reads\n", 
            100 * sum(!is.na(mapped$construct_id)) / nrow(mapped)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t3, units="secs")))

# STEP 4: Aggregate by sample
cat("Step 4: Aggregating by sample...\n")
t4 <- Sys.time()
aggregated <- aggregate_by_sample_barcode(mapped, sample_manifest, quiet = TRUE)
cat(sprintf("  Conditions: %d\n", data.table::uniqueN(aggregated$condition, na.rm=TRUE)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t4, units="secs")))

# STEP 5: Create count matrices
cat("Step 5: Creating count matrices...\n")
t5 <- Sys.time()
construct_counts <- create_count_matrix(aggregated, level = "construct", quiet = TRUE)
gene_counts <- create_count_matrix(aggregated, level = "gene", quiet = TRUE)
cat(sprintf("  Constructs: %d x %d\n", nrow(construct_counts$matrix), ncol(construct_counts$matrix)))
cat(sprintf("  Genes: %d x %d\n", nrow(gene_counts$matrix), ncol(gene_counts$matrix)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t5, units="secs")))

# STEP 6: Normalize
cat("Step 6: Normalizing counts...\n")
t6 <- Sys.time()
norm_construct <- lognormalize_counts(construct_counts$matrix, quiet = TRUE)
norm_gene <- lognormalize_counts(gene_counts$matrix, quiet = TRUE)
cat(sprintf("  Normalization complete\n"))
cat(sprintf("  Time: %.2f sec\n\n", as.numeric(difftime(Sys.time(), t6, units="secs")))

# STEP 7: QC metrics
cat("Step 7: Calculating QC metrics...\n")
t7 <- Sys.time()
matched_dt <- extracted[!is.na(construct_barcode)]
qc <- calculate_qc_metrics(extracted, matched_dt, aggregated, quiet = TRUE)
cat(sprintf("  Overall metrics calculated\n"))
cat(sprintf("  Per-sample metrics: %d samples\n", nrow(qc$per_sample_metrics)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t7, units="secs")))

# STEP 8: Save results
cat("Step 8: Saving results...\n")
t8 <- Sys.time()
write.csv(construct_counts$matrix, file.path(output_dir, "construct_counts_raw.csv"))
write.csv(norm_construct, file.path(output_dir, "construct_counts_normalized.csv"))
write.csv(gene_counts$matrix, file.path(output_dir, "gene_counts_raw.csv"))
write.csv(norm_gene, file.path(output_dir, "gene_counts_normalized.csv"))
write.csv(qc$overall_metrics, file.path(output_dir, "qc_metrics_overall.csv"), row.names = FALSE)
write.csv(qc$per_sample_metrics, file.path(output_dir, "qc_metrics_per_sample.csv"), row.names = FALSE)
cat(sprintf("  Results saved to: %s\n", output_dir))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t8, units="secs")))

# Summary
total_time <- Sys.time() - t1
cat("========== PIPELINE COMPLETE ==========\n")
cat(sprintf("Total time: %.1f seconds (%.2f minutes)\n\n", as.numeric(total_time, units="secs"), as.numeric(total_time, units="mins")))
