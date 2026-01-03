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

# UPDATE THESE PATHS TO YOUR DATA DIRECTORY
data_dir <- "path/to/your/CP0045/data"
fastq_files <- file.path(data_dir, c("sample1.fastq.gz", "sample2.fastq.gz"))
reference_csv <- file.path(data_dir, "CP0045_reference.csv")
chip_file <- file.path(data_dir, "CP0045_GRCm38_NCBI_CRISPRko_strict_gene.chip")
sample_manifest <- file.path(data_dir, "samples.csv")
output_dir <- file.path(data_dir, "crispyr_results")

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
cat("Construct-level counts:\n")
print(head(result$construct_counts$matrix))
cat("\nGene-level counts:\n")
print(head(result$gene_counts$matrix))
