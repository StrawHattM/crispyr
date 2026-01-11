# ====================================================================
# CRISPYR: Complete CRISPR Screen Processing Pipeline (Step-by-Step)
# ====================================================================
# This script runs each pipeline step individually with timing.
#
# TIMING BREAKDOWN (for ~100k reads):
# - Load reference: 1-2 sec
# - Extract barcodes: 5-10 sec
# - Match barcodes: 2-5 sec
# - Aggregate samples: 1 sec
# - Count matrices: 1 sec
# - Normalize: <1 sec
# - QC metrics: 1 sec
# - Save: 1 sec
# TOTAL: ~15-25 seconds
# ====================================================================

devtools::load_all()
library(data.table)

# Data directory with subset example
data_dir <- "D:/Bibliotecas/Biologia/crispyr_fastq/example_data/subset_example"

# File paths
fastq_files <- paste0(data_dir, "/subset_sample.fastq.gz")
reference_csv <- paste0(data_dir, "/CP0045_reference_subset.csv")
chip_file <- paste0(data_dir, "/CP0045_GRCm38_NCBI_CRISPRko_strict_gene_subset.chip")
sample_manifest_csv <- paste0(data_dir, "/sample_manifest_subset.csv")
output_dir <- paste0(data_dir, "/crispyr_results_stepwise")
dir.create(output_dir, showWarnings = FALSE)

cat("\n========== CRISPYR Step-by-Step Pipeline ==========\n\n")

# STEP 1: Load reference data
cat("Step 1: Loading reference data...\n")
t1 <- Sys.time()
ref_library <- load_reference_library(reference_csv, quiet = FALSE)
chip <- load_chip_file(chip_file, quiet = TRUE)
sample_manifest <- load_sample_manifest(sample_manifest_csv, quiet = FALSE)
cat(sprintf("  Loaded %d constructs, %d gene mappings, %d samples\n",
            nrow(ref_library), nrow(chip), nrow(sample_manifest)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t1, units = "secs"))))

# STEP 2: Extract barcodes
cat("Step 2: Extracting barcodes...\n")
t2 <- Sys.time()
extracted <- extract_barcodes_from_fastq(fastq_files, quiet = TRUE)
cat(sprintf("  Extracted from %d reads\n", nrow(extracted)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t2, units = "secs"))))

# STEP 3: Match barcodes
cat("Step 3: Matching to library...\n")
t3 <- Sys.time()
matched <- match_constructs(extracted, ref_library, max_mismatches = 1, quiet = FALSE)
mapped <- map_constructs_to_genes(matched, chip, quiet = FALSE)
cat(sprintf("  Matched %.1f%% of reads\n",
            100 * sum(!is.na(mapped$construct_id)) / nrow(mapped)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t3, units = "secs"))))

# STEP 4: Aggregate by sample
cat("Step 4: Aggregating by sample...\n")
t4 <- Sys.time()
aggregated <- aggregate_by_sample_barcode(mapped, sample_manifest, quiet = FALSE, keep_unmatched_samples = TRUE)
cat(sprintf("  Conditions: %d\n", data.table::uniqueN(aggregated$condition, na.rm = TRUE)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t4, units = "secs"))))

# STEP 5: Create count matrices
cat("Step 5: Creating count matrices...\n")
t5 <- Sys.time()
construct_counts <- create_count_matrix(aggregated, level = "construct", quiet = FALSE)
gene_counts <- create_count_matrix(aggregated, level = "gene", quiet = FALSE)
cat(sprintf("  Constructs: %d x %d\n", nrow(construct_counts$matrix), ncol(construct_counts$matrix)))
cat(sprintf("  Genes: %d x %d\n", nrow(gene_counts$matrix), ncol(gene_counts$matrix)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t5, units = "secs"))))

# STEP 6: Normalize
cat("Step 6: Normalizing counts...\n")
t6 <- Sys.time()
norm_construct <- lognormalize_counts(construct_counts$matrix, quiet = FALSE)
norm_gene <- lognormalize_counts(gene_counts$matrix, quiet = FALSE)
cat(sprintf("  Normalization complete\n"))
cat(sprintf("  Time: %.2f sec\n\n", as.numeric(difftime(Sys.time(), t6, units = "secs"))))

# STEP 7: QC metrics
cat("Step 7: Calculating QC metrics...\n")
t7 <- Sys.time()
matched_dt <- extracted[!is.na(construct_barcode)]
qc <- calculate_qc_metrics(extracted, matched_dt, aggregated, quiet = FALSE)
cat(sprintf("  Overall metrics calculated\n"))
cat(sprintf("  Per-sample metrics: %d samples\n", nrow(qc$per_sample_metrics)))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t7, units = "secs"))))

# STEP 8: Save results
cat("Step 8: Saving results...\n")
t8 <- Sys.time()
write.csv(construct_counts$matrix, paste0(output_dir, "/construct_counts_raw.csv"))
write.csv(norm_construct, paste0(output_dir, "/construct_counts_normalized.csv"))
write.csv(gene_counts$matrix, paste0(output_dir, "/gene_counts_raw.csv"))
write.csv(norm_gene, paste0(output_dir, "/gene_counts_normalized.csv"))
write.csv(qc$overall_metrics, paste0(output_dir, "/qc_metrics_overall.csv"), row.names = FALSE)
write.csv(qc$per_sample_metrics, paste0(output_dir, "/qc_metrics_per_sample.csv"), row.names = FALSE)
cat(sprintf("  Results saved to: %s\n", output_dir))
cat(sprintf("  Time: %.1f sec\n\n", as.numeric(difftime(Sys.time(), t8, units = "secs"))))

# Summary
total_time <- Sys.time() - t1
cat("========== PIPELINE COMPLETE ==========\n")
cat(sprintf("Total time: %.1f seconds (%.2f minutes)\n\n", as.numeric(total_time, units = "secs"), as.numeric(total_time, units = "mins")))
