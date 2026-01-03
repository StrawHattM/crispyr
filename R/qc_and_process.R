#' Calculate Quality Control Metrics
#'
#' Computes QC metrics for the pooled CRISPR screen analysis,
#' including match percentages, coverage depth, and library complexity.
#'
#' @param extracted_barcodes A data.table from [extract_barcodes_from_fastq()]
#'   containing extraction status for each read.
#' @param matched_reads A data.table from [match_constructs()] with barcode
#'   matching information.
#' @param aggregated_reads A data.table from [aggregate_by_sample_barcode()]
#'   with sample condition assignments.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A list with elements:
#'   - `overall_metrics`: Data.frame with overall statistics
#'     * total_reads
#'     * successful_extractions (%)
#'     * matched_constructs (%)
#'     * unique_constructs
#'     * median_depth
#'   - `per_sample_metrics`: Data.frame with per-sample statistics
#'     * condition
#'     * total_reads
#'     * matched_reads (%)
#'     * unique_constructs
#'     * median_depth
#'     * coverage (% of reference constructs detected)
#'
#' @examples
#' \dontrun{
#'   # Calculate QC metrics
#'   qc <- calculate_qc_metrics(
#'     extracted_barcodes,
#'     matched_reads,
#'     aggregated_reads
#'   )
#'
#'   # View overall metrics
#'   qc$overall_metrics
#'
#'   # View per-sample metrics
#'   qc$per_sample_metrics
#' }
#'
#' @importFrom data.table as.data.table
#' @importFrom dplyr summarise group_by n_distinct
#' @export
calculate_qc_metrics <- function(extracted_barcodes,
                                  matched_reads,
                                  aggregated_reads,
                                  quiet = FALSE) {
  # Convert to data.table if needed
  if (!data.table::is.data.table(extracted_barcodes)) {
    extracted_barcodes <- data.table::as.data.table(extracted_barcodes)
  }
  if (!data.table::is.data.table(matched_reads)) {
    matched_reads <- data.table::as.data.table(matched_reads)
  }
  if (!data.table::is.data.table(aggregated_reads)) {
    aggregated_reads <- data.table::as.data.table(aggregated_reads)
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Calculating quality control metrics"
    ))
  }

  # ===== OVERALL METRICS =====
  total_reads <- nrow(extracted_barcodes)

  successful_extractions <- extracted_barcodes[
    extraction_status %in% c("success", "no_sample_barcode"),
    .N
  ]
  successful_pct <- 100 * successful_extractions / total_reads

  matched_constructs_n <- matched_reads[
    !is.na(construct_id),
    .N
  ]
  matched_pct <- 100 * matched_constructs_n / successful_extractions

  unique_constructs <- matched_reads[
    !is.na(construct_id),
    uniqueN(construct_id)
  ]

  # Calculate depth statistics from construct-level counts
  construct_counts <- matched_reads[
    !is.na(construct_id),
    .(read_count = .N),
    by = "construct_id"
  ]

  median_depth <- median(construct_counts$read_count)
  mean_depth <- mean(construct_counts$read_count)
  min_depth <- min(construct_counts$read_count)
  max_depth <- max(construct_counts$read_count)

  overall_metrics <- data.frame(
    metric = c(
      "Total Reads", "Successful Extractions (%)", "Matched Constructs (%)",
      "Unique Constructs", "Median Depth", "Mean Depth", "Min Depth", "Max Depth"
    ),
    value = c(
      total_reads, round(successful_pct, 2), round(matched_pct, 2),
      unique_constructs, median_depth, round(mean_depth, 2),
      min_depth, max_depth
    )
  )

  # ===== PER-SAMPLE METRICS =====
  per_sample <- aggregated_reads[
    !is.na(condition),
    .(
      total_reads = .N,
      matched_constructs = sum(!is.na(construct_id)),
      unique_constructs = n_distinct(construct_id, na.rm = TRUE),
      unique_genes = n_distinct(gene_id, na.rm = TRUE)
    ),
    by = "condition"
  ] %>%
    dplyr::mutate(
      matched_pct = round(100 * matched_constructs / total_reads, 2),
      condition = as.character(condition)
    ) %>%
    dplyr::select(
      condition, total_reads, matched_constructs, matched_pct,
      unique_constructs, unique_genes
    ) %>%
    as.data.frame()

  # Calculate per-sample depth
  sample_depths <- aggregated_reads[
    !is.na(condition) & !is.na(construct_id),
    .(median_depth = median(.N), mean_depth = mean(.N)),
    by = c("condition", "construct_id")
  ][
    ,
    .(
      sample_median_depth = median(median_depth),
      sample_mean_depth = mean(mean_depth)
    ),
    by = "condition"
  ]

  per_sample <- merge(
    per_sample,
    sample_depths,
    by = "condition",
    all.x = TRUE
  )

  per_sample <- per_sample %>%
    dplyr::mutate(
      sample_median_depth = round(sample_median_depth, 2),
      sample_mean_depth = round(sample_mean_depth, 2)
    ) %>%
    as.data.frame()

  if (!quiet) {
    cli::cli_inform(c(
      "✓" = "QC metrics calculated",
      "i" = "Overall: {total_reads} reads, {unique_constructs} unique constructs",
      "i" = "{nrow(per_sample)} samples analyzed"
    ))
  }

  list(
    overall_metrics = overall_metrics,
    per_sample_metrics = per_sample
  )
}


#' Process Pooled CRISPR Screen FASTQ Files
#'
#' Main wrapper function that orchestrates the complete CRISPR screen analysis pipeline:
#' FASTQ processing → barcode extraction → library matching → gene mapping →
#' count matrix generation → normalization → QC reporting.
#'
#' @param fastq_files Character vector of FASTQ file paths. Can be gzip-compressed.
#' @param reference_library_csv Path to the reference library CSV file
#'   (see [load_reference_library()]).
#' @param chip_file Path to the gene annotation chip file (see [load_chip_file()]).
#' @param sample_manifest_csv Path to the sample manifest CSV file
#'   (see [load_sample_manifest()]).
#' @param barcode_policy Character. Barcode extraction policy string.
#'   Format: "PREFIX:SEQUENCE@POSITION". Default: "PREFIX:CACCG@6" (Broad Institute).
#' @param max_mismatches Integer. Maximum mismatches for barcode matching. Default: 0 (exact).
#' @param normalize Logical. If `TRUE`, output log-normalized counts. Default: `TRUE`.
#' @param output_dir Character. Directory path for saving output files.
#'   If NULL (default), files are not saved.
#' @param output_prefix Character. Prefix for output file names. Default: "crispyr_output".
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A list with elements:
#'   - `construct_counts`: Construct-level count matrix
#'   - `gene_counts`: Gene-level count matrix
#'   - `construct_lognorm`: Log-normalized construct counts (if `normalize = TRUE`)
#'   - `gene_lognorm`: Log-normalized gene counts (if `normalize = TRUE`)
#'   - `qc_metrics`: Quality control metrics
#'   - `feature_metadata`: Metadata for features (constructs/genes)
#'   - `sample_metadata`: Metadata for samples
#'
#' @details
#' ## Pipeline Steps
#'
#' 1. Extract construct (20bp) and sample (8bp) barcodes from FASTQ reads
#' 2. Match construct barcodes to reference library
#' 3. Map construct barcodes to target genes using chip file
#' 4. Group reads by sample barcode and condition
#' 5. Create construct-level and gene-level count matrices
#' 6. Log-normalize counts (optional)
#' 7. Calculate quality control metrics
#'
#' ## Output Files
#'
#' If `output_dir` is specified, the following files are saved:
#' - `{prefix}_construct_counts.tsv`: Raw construct counts
#' - `{prefix}_gene_counts.tsv`: Raw gene counts
#' - `{prefix}_construct_lognorm.tsv`: Log-normalized construct counts
#' - `{prefix}_gene_lognorm.tsv`: Log-normalized gene counts
#' - `{prefix}_qc_overall.txt`: Overall QC metrics
#' - `{prefix}_qc_per_sample.tsv`: Per-sample QC metrics
#'
#' @examples
#' \dontrun{
#'   # Process a pooled CRISPR screen
#'   result <- process_pooled_screen(
#'     fastq_files = c("sample1.fastq.gz", "sample2.fastq.gz"),
#'     reference_library_csv = "CP0045_reference.csv",
#'     chip_file = "CP0045_GRCm38_strict.chip",
#'     sample_manifest_csv = "samples.csv",
#'     output_dir = "./results"
#'   )
#'
#'   # Access results
#'   construct_counts <- result$construct_counts$matrix
#'   gene_counts <- result$gene_counts$matrix
#' }
#'
#' @export
process_pooled_screen <- function(fastq_files,
                                   reference_library_csv,
                                   chip_file,
                                   sample_manifest_csv,
                                   barcode_policy = "PREFIX:CACCG@6",
                                   max_mismatches = 0,
                                   normalize = TRUE,
                                   output_dir = NULL,
                                   output_prefix = "crispyr_output",
                                   quiet = FALSE) {
  if (!quiet) {
    cli::cli_h1("CRISPR Screen Processing Pipeline")
  }

  # Step 1: Load reference data
  if (!quiet) {
    cli::cli_h2("Loading reference files")
  }

  ref_library <- load_reference_library(reference_library_csv, quiet = quiet)
  chip <- load_chip_file(chip_file, quiet = quiet)
  manifest <- load_sample_manifest(sample_manifest_csv, quiet = quiet)

  # Step 2: Extract barcodes
  if (!quiet) {
    cli::cli_h2("Extracting barcodes from FASTQ files")
  }

  extracted <- extract_barcodes_from_fastq(
    fastq_files = fastq_files,
    barcode_policy = barcode_policy,
    quiet = quiet
  )

  # Step 3: Match constructs
  if (!quiet) {
    cli::cli_h2("Matching barcodes to reference library")
  }

  matched <- match_constructs(
    extracted,
    ref_library,
    max_mismatches = max_mismatches,
    quiet = quiet
  )

  # Step 4: Map genes
  if (!quiet) {
    cli::cli_h2("Mapping constructs to genes")
  }

  with_genes <- map_constructs_to_genes(matched, chip, quiet = quiet)

  # Step 5: Aggregate by sample
  if (!quiet) {
    cli::cli_h2("Aggregating reads by sample")
  }

  aggregated <- aggregate_by_sample_barcode(with_genes, manifest, quiet = quiet)

  # Step 6: Create count matrices
  if (!quiet) {
    cli::cli_h2("Creating count matrices")
  }

  construct_result <- create_count_matrix(
    aggregated,
    level = "construct",
    quiet = quiet
  )

  gene_result <- create_count_matrix(
    aggregated,
    level = "gene",
    quiet = quiet
  )

  # Step 7: Normalize (optional)
  results <- list(
    construct_counts = construct_result,
    gene_counts = gene_result
  )

  if (normalize) {
    if (!quiet) {
      cli::cli_h2("Normalizing counts")
    }

    construct_lognorm <- lognormalize_counts(
      construct_result$matrix,
      method = "poolq",
      quiet = quiet
    )

    gene_lognorm <- lognormalize_counts(
      gene_result$matrix,
      method = "poolq",
      quiet = quiet
    )

    results$construct_lognorm <- construct_lognorm
    results$gene_lognorm <- gene_lognorm
  }

  # Step 8: QC metrics
  if (!quiet) {
    cli::cli_h2("Calculating quality control metrics")
  }

  qc <- calculate_qc_metrics(
    extracted,
    matched,
    aggregated,
    quiet = quiet
  )

  results$qc_metrics <- qc
  results$feature_metadata <- construct_result$feature_metadata
  results$sample_metadata <- construct_result$sample_metadata

  # Step 9: Save results (optional)
  if (!is.null(output_dir)) {
    if (!quiet) {
      cli::cli_h2("Saving results")
    }

    save_results(
      results,
      output_dir = output_dir,
      output_prefix = output_prefix,
      quiet = quiet
    )
  }

  if (!quiet) {
    cli::cli_h1("Processing complete!")
  }

  results
}


#' Save Processing Results
#'
#' Saves the results from [process_pooled_screen()] to TSV and text files.
#'
#' @param results List returned from [process_pooled_screen()].
#' @param output_dir Character path to output directory.
#' @param output_prefix Character prefix for output files.
#' @param quiet Logical. If `TRUE`, suppress messages. Default: `FALSE`.
#'
#' @return Invisibly returns the results list (called for side effects).
#'
#' @keywords internal
save_results <- function(results,
                          output_dir,
                          output_prefix = "crispyr_output",
                          quiet = FALSE) {
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save count matrices
  if (!is.null(results$construct_counts)) {
    write.table(
      results$construct_counts$matrix,
      file = file.path(output_dir, paste0(output_prefix, "_construct_counts.tsv")),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  if (!is.null(results$gene_counts)) {
    write.table(
      results$gene_counts$matrix,
      file = file.path(output_dir, paste0(output_prefix, "_gene_counts.tsv")),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  # Save normalized matrices
  if (!is.null(results$construct_lognorm)) {
    write.table(
      results$construct_lognorm,
      file = file.path(output_dir, paste0(output_prefix, "_construct_lognorm.tsv")),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  if (!is.null(results$gene_lognorm)) {
    write.table(
      results$gene_lognorm,
      file = file.path(output_dir, paste0(output_prefix, "_gene_lognorm.tsv")),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  # Save QC metrics
  if (!is.null(results$qc_metrics)) {
    write.table(
      results$qc_metrics$overall_metrics,
      file = file.path(output_dir, paste0(output_prefix, "_qc_overall.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    write.table(
      results$qc_metrics$per_sample_metrics,
      file = file.path(output_dir, paste0(output_prefix, "_qc_per_sample.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (!quiet) {
    cli::cli_inform(c(
      "✓" = "Results saved to {output_dir}"
    ))
  }

  invisible(results)
}
