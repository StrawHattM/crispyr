# ====================================================================
# CRISPYR: Complete CRISPR Screen Processing Pipeline (Wrapper)
# ====================================================================
# This script runs the entire pipeline with a single function call.
#
# SYSTEM PERFORMANCE ESTIMATE:
# - Small screens (1-5M reads): 30-60 seconds
# - Medium screens (10-20M reads): 2-4 minutes
# - Large screens (50M reads): 8-12 minutes
# - Very large screens (100M+ reads): 20-40 minutes
# ====================================================================

library(crispyr)

# Data directory with subset example
data_dir <- "/d/Bibliotecas/Biologia/crispyr_fastq/example_data/subset_example"

# File paths
fastq_files <- paste0(data_dir, "/subset_sample.fastq.gz")
reference_csv <- paste0(data_dir, "/CP0045_reference_subset.csv")
chip_file <- paste0(data_dir, "/CP0045_GRCm38_NCBI_CRISPRko_strict_gene_subset.chip")
sample_manifest <- paste0(data_dir, "/sample_manifest_subset.csv")
output_dir <- paste0(data_dir, "/crispyr_results")

dir.create(output_dir, showWarnings = FALSE)

cat("\n========== CRISPYR Pipeline (Wrapper) ==========\n\n")

start_time <- Sys.time()

result <- process_pooled_screen(
  fastq_files = fastq_files,
  reference_library_csv = reference_csv,
  chip_file = chip_file,
  sample_manifest_csv = sample_manifest,
  barcode_policy = "PREFIX:CACCG@6",
  max_mismatches = 1,
  normalize = TRUE,
  output_dir = output_dir,
  quiet = FALSE
)

elapsed <- difftime(Sys.time(), start_time, units = "mins")

cat("\n========== PIPELINE COMPLETE ==========\n")
cat(sprintf("Total time: %.2f minutes\n\n", as.numeric(elapsed)))

# View results
cat("Construct-level counts (first 5 rows):\n")
print(head(result$construct_counts$matrix, 5))
cat("\nGene-level counts (first 5 rows):\n")
print(head(result$gene_counts$matrix, 5))
cat("\nQC Metrics:\n")
print(result$qc_metrics)
