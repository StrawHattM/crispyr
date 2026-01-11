#' Aggregate Reads by Sample Barcode
#'
#' Groups extracted barcodes by their 8bp sample barcode and maps them to
#' experimental conditions/samples using the sample manifest.
#'
#' @param mapped_reads A data.table from [map_constructs_to_genes()] containing
#'   barcode mapping results with a `sample_barcode` column.
#' @param sample_manifest A tibble/data.frame from [load_sample_manifest()] mapping
#'   8bp sample barcodes to condition names.
#' @param keep_unmatched_samples Logical. If `TRUE`, keep reads with sample barcodes
#'   not in the manifest. Default: `FALSE`.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A [data.table::data.table] with input columns plus:
#'   - `condition`: Sample condition/name from the manifest
#'   - `sample_barcode_matched`: Logical indicating if sample barcode was found in manifest
#'
#' @examples
#' \dontrun{
#'   # Aggregate by sample barcode
#'   with_conditions <- aggregate_by_sample_barcode(
#'     mapped_reads,
#'     sample_manifest
#'   )
#' }
#'
#' @importFrom data.table merge.data.table
#' @export
aggregate_by_sample_barcode <- function(mapped_reads,
                                         sample_manifest,
                                         keep_unmatched_samples = FALSE,
                                         quiet = FALSE) {
   # Convert to data.table if needed
   if (!data.table::is.data.table(mapped_reads)) {
     mapped_reads <- data.table::as.data.table(mapped_reads)
   }
 
   manifest_dt <- data.table::as.data.table(sample_manifest)
 
   if (!quiet) {
     cli::cli_inform(c(
       "i" = "Mapping sample barcodes to conditions"
     ))
   }

   # DEBUG: Check input data
   if (!quiet) {
     n_reads <- nrow(mapped_reads)
     n_unique_sample_bc <- data.table::uniqueN(mapped_reads$sample_barcode, na.rm = TRUE)
     n_manifest <- nrow(manifest_dt)
     n_unique_manifest_bc <- data.table::uniqueN(manifest_dt$sample_barcode, na.rm = TRUE)
     
     cli::cli_inform(c(
       "=" = "aggregate_by_sample_barcode debug info:",
       "i" = "Input reads: {n_reads} rows, {n_unique_sample_bc} unique sample_barcodes",
       "i" = "Manifest: {n_manifest} rows, {n_unique_manifest_bc} unique sample_barcodes"
     ))
   }
 
   # Perform left join on sample barcode
   result <- data.table::merge.data.table(
     mapped_reads,
     manifest_dt,
     by = "sample_barcode",
     all.x = TRUE,
     sort = FALSE
   )
 
   # DEBUG: Check what happened with the merge
   if (!quiet) {
     has_condition_col <- "condition" %in% names(result)
     if (has_condition_col) {
       n_na_condition <- sum(is.na(result$condition))
       n_not_na_condition <- sum(!is.na(result$condition))
       cli::cli_inform(c(
         "i" = "Merged result has 'condition' column",
         "i" = "  NA conditions: {n_na_condition}",
         "i" = "  Non-NA conditions: {n_not_na_condition}"
       ))
     } else {
       cli::cli_inform(c(
         "x" = "ERROR: 'condition' column NOT found after merge!",
         "i" = "Columns in result: {paste(names(result), collapse=', ')}"
       ))
     }
   }
 
   # Mark which samples were matched
   result[, sample_barcode_matched := !is.na(condition)]
 
   if (!quiet) {
     n_matched_after_merge <- sum(result$sample_barcode_matched, na.rm = TRUE)
     cli::cli_inform(c(
       "i" = "Matched samples in merged data: {n_matched_after_merge}"
     ))
   }
 
    # Optionally remove unmatched samples
    if (!keep_unmatched_samples) {
      result <- result[sample_barcode_matched == TRUE]
    }
 
    # Check if result is empty
    if (nrow(result) == 0) {
      rlang::abort(c(
        "No reads matched to sample manifest.",
        "x" = "aggregate_by_sample_barcode returned empty dataset",
        "i" = "Check that sample_barcode values in mapped_reads match the sample_manifest"
      ))
    }
 
    # Summary statistics
    n_matched <- result[sample_barcode_matched == TRUE, .N]
    n_unmatched <- result[sample_barcode_matched == FALSE, .N]
 
   if (!quiet) {
     cli::cli_inform(c(
       "*" = "Sample barcode aggregation complete",
       "i" = "Reads with matched samples: {n_matched}",
       "i" = "Reads with unmatched samples: {n_unmatched}"
     ))
   }
 
   result
 }


#' Create Count Matrix
#'
#' Constructs a count matrix (samples  x  constructs/genes) from mapped and
#' aggregated reads.
#'
#' @param aggregated_reads A data.table from [aggregate_by_sample_barcode()]
#'   with construct/gene mappings and condition assignments.
#' @param level Character. Level of aggregation for the count matrix.
#'   - "construct": Rows are construct barcodes (default, most fine-grained)
#'   - "gene": Rows are genes (aggregates multiple guides per gene)
#' @param include_unmatched Logical. If `TRUE`, include reads without gene/construct
#'   assignment. Default: `FALSE`.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A list with elements:
#'   - `matrix`: A count matrix (matrix class) with dimensions (constructs/genes  x  samples).
#'     Rows are construct barcodes or genes, columns are samples/conditions.
#'   - `feature_metadata`: A data.frame with metadata for each row (barcode, construct_id, etc.)
#'   - `sample_metadata`: A data.frame with metadata for each column (condition, sample_barcode, etc.)
#'
#' @details
#' ## Count Matrix Construction
#'
#' 1. **Construct-level** (`level = "construct"`):
#'    Each row is a unique construct barcode. Multiple reads of the same
#'    barcode in the same sample are summed.
#'
#' 2. **Gene-level** (`level = "gene"`):
#'    Each row is a unique gene (identified by gene_id). Reads mapping to
#'    multiple guides (barcodes) targeting the same gene are summed across guides.
#'
#' @examples
#' \dontrun{
#'   # Create construct-level matrix
#'   construct_counts <- create_count_matrix(
#'     aggregated_reads,
#'     level = "construct"
#'   )
#'
#'   # Create gene-level matrix
#'   gene_counts <- create_count_matrix(
#'     aggregated_reads,
#'     level = "gene"
#'   )
#' }
#'
#' @importFrom data.table as.data.table dcast
#' @importFrom tidyr pivot_wider
#' @export
create_count_matrix <- function(aggregated_reads,
                                 level = c("construct", "gene"),
                                 include_unmatched = FALSE,
                                 quiet = FALSE) {
  level <- match.arg(level)

  # Convert to data.table if needed
  if (!data.table::is.data.table(aggregated_reads)) {
    aggregated_reads <- data.table::as.data.table(aggregated_reads)
  }

  # Filter data based on level and unmatched setting
  if (level == "construct") {
    if (!include_unmatched) {
      data_to_count <- aggregated_reads[!is.na(construct_id)]
    } else {
      data_to_count <- aggregated_reads
    }
    row_feature <- "construct_id"
    feature_cols <- c("barcode", "construct_id")
  } else if (level == "gene") {
    if (!include_unmatched) {
      data_to_count <- aggregated_reads[!is.na(gene_id)]
    } else {
      data_to_count <- aggregated_reads
    }
    row_feature <- "gene_id"
    feature_cols <- c("gene_symbol", "gene_id")
  }

  if (nrow(data_to_count) == 0) {
    rlang::abort(
      c("No data available for count matrix at {level} level.",
        "x" = "Try setting include_unmatched = TRUE")
    )
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Creating {level}-level count matrix"
    ))
  }

  # Count reads per feature per condition
  count_data <- data_to_count[
    !is.na(condition),
    .(read_count = .N),
    by = c(row_feature, "condition")
  ]

  # Pivot to matrix format
  count_matrix <- data.table::dcast(
    count_data,
    as.formula(paste(row_feature, "~ condition")),
    value.var = "read_count",
    fill = 0,
    fun.aggregate = sum
  )

  # Convert to matrix
  row_names <- count_matrix[[row_feature]]
  count_matrix <- as.matrix(count_matrix[, -1, with = FALSE])
  rownames(count_matrix) <- row_names

  # Create feature metadata
  if (level == "construct") {
    feature_metadata <- unique(aggregated_reads[
      !is.na(construct_id),
      ..feature_cols
    ])
    feature_metadata <- as.data.frame(feature_metadata)
    rownames(feature_metadata) <- feature_metadata$construct_id
  } else {
    feature_metadata <- unique(aggregated_reads[
      !is.na(gene_id),
      ..feature_cols
    ])
    feature_metadata <- as.data.frame(feature_metadata)
    rownames(feature_metadata) <- feature_metadata$gene_id
  }

  # Create sample metadata
  sample_metadata <- unique(aggregated_reads[
    !is.na(condition),
    .(condition, sample_barcode)
  ])
  sample_metadata <- as.data.frame(sample_metadata)
  rownames(sample_metadata) <- sample_metadata$condition

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Count matrix created: {nrow(count_matrix)} {level}s  x  {ncol(count_matrix)} samples"
    ))
  }

  list(
    matrix = count_matrix,
    feature_metadata = feature_metadata,
    sample_metadata = sample_metadata
  )
}


#' Log-Normalize Counts
#'
#' Applies log-normalization to count matrices using the PoolQ method:
#' `log2((count + 1) / (median_count_per_sample + 1))`
#'
#' This normalization accounts for differences in sequencing depth between samples
#' while avoiding log(0) by adding a pseudocount of 1.
#'
#' @param count_matrix A numeric matrix with samples as columns and features
#'   (constructs/genes) as rows.
#' @param method Character. Normalization method.
#'   - "poolq": PoolQ-style log2 normalization (default)
#'   - "log2": Simple log2(count + 1)
#'   - "cpm": Counts per million
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A numeric matrix with the same dimensions as input, containing
#'   log-normalized counts.
#'
#' @details
#' ## PoolQ Normalization Formula
#'
#' For each sample, the median read count across all features is calculated.
#' Each feature's count in that sample is then normalized as:
#'
#' ```
#' lognorm = log2((count + 1) / (median + 1))
#' ```
#'
#' This makes counts comparable across samples with different sequencing depths.
#'
#' @examples
#' \dontrun{
#'   # Create count matrix
#'   count_data <- create_count_matrix(aggregated_reads)
#'
#'   # Normalize
#'   lognorm_matrix <- lognormalize_counts(count_data$matrix)
#' }
#'
#' @export
#' @importFrom stats median
lognormalize_counts <- function(count_matrix,
                                 method = c("poolq", "log2", "cpm"),
                                quiet = FALSE) {
  method <- match.arg(method)

  if (!is.matrix(count_matrix)) {
    rlang::abort("count_matrix must be a numeric matrix")
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Applying {method} normalization"
    ))
  }

  if (method == "poolq") {
    # PoolQ method: log2((count + 1) / (median + 1))
    normalized <- count_matrix

    for (col in seq_len(ncol(count_matrix))) {
      median_count <- median(count_matrix[, col])
      normalized[, col] <- log2((count_matrix[, col] + 1) / (median_count + 1))
    }
  } else if (method == "log2") {
    # Simple log2 transformation
    normalized <- log2(count_matrix + 1)
  } else if (method == "cpm") {
    # Counts per million
    normalized <- sweep(
      count_matrix,
      2,
      colSums(count_matrix),
      FUN = "/"
    ) * 1e6
    normalized <- log2(normalized + 1)
  }

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Normalization complete"
    ))
  }

  normalized
}
