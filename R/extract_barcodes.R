#' Parse Barcode Policy String
#'
#' Parses a barcode policy string in the format "PREFIX:SEQUENCE@POSITION"
#' to extract the flanking sequence and barcode position.
#'
#' @param policy Character string specifying the barcode extraction policy.
#'   Format: "PREFIX:SEQUENCE@POSITION" where:
#'   - PREFIX identifies the flanking sequence type (currently only "PREFIX" supported)
#'   - SEQUENCE is the flanking DNA sequence (e.g., "CACCG")
#'   - POSITION is the 1-based position where the barcode starts after the sequence
#'
#' @return A list with elements:
#'   - `flanking_seq`: The flanking DNA sequence
#'   - `position_offset`: The position offset relative to the flanking sequence
#'
#' @examples
#' parse_barcode_policy("PREFIX:CACCG@6")
#' # Returns: list(flanking_seq = "CACCG", position_offset = 6)
#'
#' @keywords internal
parse_barcode_policy <- function(policy) {
  # Match pattern: PREFIX:SEQUENCE@POSITION
  pattern <- "^(PREFIX):([ACGT]+)@(\\d+)$"

  if (!stringr::str_detect(policy, pattern)) {
    rlang::abort(c(
      "Invalid barcode policy format.",
      "x" = glue::glue("Policy: {policy}"),
      "i" = "Expected format: 'PREFIX:SEQUENCE@POSITION' (e.g., 'PREFIX:CACCG@6')"
    ))
  }

  match_info <- stringr::str_match(policy, pattern)

  list(
    flanking_seq = match_info[, 3],
    position_offset = as.integer(match_info[, 4])
  )
}


#' Extract Construct and Sample Barcodes from FASTQ Files
#'
#' Extracts construct (20bp) and sample (8bp) barcodes from FASTQ reads.
#' The construct barcode position is defined by a barcode policy specifying
#' a flanking sequence and position offset.
#'
#' @param fastq_files Character vector of paths to FASTQ files (gzip compressed or not).
#' @param barcode_policy Character string specifying barcode extraction policy.
#'   Default: "PREFIX:CACCG@6" (Broad Institute standard).
#'   Format: "PREFIX:SEQUENCE@POSITION"
#' @param construct_barcode_length Integer. Length of construct barcode in bp.
#'   Default: 20 (standard for CRISPR guides).
#' @param min_quality_score Integer. Minimum mean quality score threshold.
#'   If provided, reads below this threshold are filtered. Default: NULL (no filtering).
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A [data.table::data.table] with columns:
#'   - `fastq_file`: Source FASTQ file name
#'   - `read_id`: Read identifier from FASTQ header
#'   - `construct_barcode`: Extracted 20bp construct barcode
#'   - `sample_barcode`: Extracted 8bp sample barcode (if found, else NA)
#'   - `extraction_status`: Status of barcode extraction
#'     * "success": Both barcodes extracted
#'     * "no_construct_flank": Flanking sequence not found
#'     * "construct_too_short": Read too short for barcode extraction
#'     * "no_sample_barcode": Construct extracted but sample barcode not parseable
#'
#' @details
#' ## Barcode Extraction Logic
#'
#' 1. **Construct barcode**: Searches for the flanking sequence (e.g., "CACCG") in the read.
#'    Once found, extracts the specified number of bases (default 20) starting from the
#'    position offset (e.g., position 6 after "CACCG").
#'
#' 2. **Sample barcode**: Extracts the first 8bp of the read as the sample barcode.
#'    This is based on PoolQ's dual-barcode design where the first 8bp identify the sample.
#'
#' @examples
#' \dontrun{
#'   # Extract from single FASTQ file
#'   barcodes <- extract_barcodes_from_fastq(
#'     "path/to/sample.fastq.gz"
#'   )
#'
#'   # Extract from multiple files
#'   barcodes <- extract_barcodes_from_fastq(
#'     c("file1.fastq.gz", "file2.fastq.gz")
#'   )
#' }
#'
#' @importFrom ShortRead readFastq
#' @importFrom Biostrings subseq start
#' @importFrom data.table data.table
#' @export
extract_barcodes_from_fastq <- function(fastq_files,
                                        barcode_policy = "PREFIX:CACCG@6",
                                        construct_barcode_length = 20,
                                        min_quality_score = NULL,
                                        quiet = FALSE) {
  # Validate inputs
  if (!all(file.exists(fastq_files))) {
    missing_files <- fastq_files[!file.exists(fastq_files)]
    rlang::abort(c(
      "Some FASTQ files not found.",
      "x" = paste(missing_files, collapse = ", ")
    ))
  }

  # Parse barcode policy
  policy <- parse_barcode_policy(barcode_policy)

  # Initialize results list
  all_barcodes <- list()

  # Process each FASTQ file
  for (fastq_file in fastq_files) {
    if (!quiet) {
      cli::cli_inform(c(
        "i" = "Processing: {basename(fastq_file)}"
      ))
    }

    # Read FASTQ file
    fastq_data <- ShortRead::readFastq(fastq_file)

    # Extract sequences and qualities
    sequences <- as.character(fastq_data@sread)
    quality_strings <- as.character(fastq_data@quality)
    read_ids <- as.character(fastq_data@id)

    # Convert quality strings to numeric scores
    if (!is.null(min_quality_score)) {
      quality_scores <- vapply(
        quality_strings,
        function(q) mean(as.integer(charToRaw(q)) - 33),
        numeric(1)
      )
      keep_reads <- quality_scores >= min_quality_score
    } else {
      keep_reads <- rep(TRUE, length(sequences))
    }

    sequences <- sequences[keep_reads]
    quality_strings <- quality_strings[keep_reads]
    read_ids <- read_ids[keep_reads]

    # Extract barcodes from each read
    barcodes_list <- mapply(
      function(seq, q_str, read_id) {
        extract_barcodes_from_read(
          seq,
          q_str,
          read_id,
          policy,
          construct_barcode_length
        )
      },
      sequences,
      quality_strings,
      read_ids,
      SIMPLIFY = FALSE
    )

    # Combine results
    barcodes_df <- data.table::rbindlist(barcodes_list, fill = TRUE)
    barcodes_df$fastq_file <- basename(fastq_file)

    all_barcodes[[fastq_file]] <- barcodes_df
  }

  # Combine all results
  result <- data.table::rbindlist(all_barcodes, fill = TRUE)

  if (!quiet) {
    n_success <- sum(result$extraction_status == "success", na.rm = TRUE)
    n_total <- nrow(result)
    pct <- round(100 * n_success / n_total, 2)
    cli::cli_inform(c(
      "✓" = "Extracted barcodes from {n_total} reads ({pct}% successful)"
    ))
  }

  result
}


#' Extract Barcodes from a Single FASTQ Read
#'
#' Internal function to extract construct and sample barcodes from a single read sequence.
#'
#' @param sequence Character. DNA sequence of the read.
#' @param quality_string Character. Quality string for the read.
#' @param read_id Character. Read identifier.
#' @param policy List. Parsed barcode policy (from `parse_barcode_policy()`).
#' @param construct_barcode_length Integer. Expected construct barcode length.
#'
#' @return A data.frame with one row containing extraction results.
#'
#' @keywords internal
extract_barcodes_from_read <- function(sequence,
                                        quality_string,
                                        read_id,
                                        policy,
                                        construct_barcode_length) {
  # Extract sample barcode (first 8bp)
  sample_barcode <- if (nchar(sequence) >= 8) {
    substr(sequence, 1, 8)
  } else {
    NA_character_
  }

  # Find flanking sequence
  flanking_pos <- stringr::str_locate(sequence, policy$flanking_seq)[1, 1]

  if (is.na(flanking_pos)) {
    return(data.table::data.table(
      read_id = read_id,
      construct_barcode = NA_character_,
      sample_barcode = sample_barcode,
      extraction_status = "no_construct_flank"
    ))
  }

  # Calculate construct barcode start position
  construct_start <- flanking_pos + policy$position_offset - 1

  # Check if read is long enough
  construct_end <- construct_start + construct_barcode_length - 1
  if (construct_end > nchar(sequence)) {
    return(data.table::data.table(
      read_id = read_id,
      construct_barcode = NA_character_,
      sample_barcode = sample_barcode,
      extraction_status = "construct_too_short"
    ))
  }

  # Extract construct barcode
  construct_barcode <- substr(sequence, construct_start, construct_end)

  # Check sample barcode validity (should be 8bp DNA)
  sample_status <- if (is.na(sample_barcode) || nchar(sample_barcode) != 8) {
    "no_sample_barcode"
  } else {
    "success"
  }

  data.table::data.table(
    read_id = read_id,
    construct_barcode = construct_barcode,
    sample_barcode = sample_barcode,
    extraction_status = sample_status
  )
}
