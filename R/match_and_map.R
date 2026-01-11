#' Match Extracted Barcodes to Reference Library
#'
#' Maps extracted construct barcodes to the reference library of known constructs.
#' Supports exact matching or approximate matching with a specified Hamming distance
#' tolerance.
#'
#' @param extracted_barcodes A data.table from [extract_barcodes_from_fastq()].
#' @param reference_library A tibble/data.frame from [load_reference_library()] with
#'   columns `barcode` and `construct_id`.
#' @param max_mismatches Integer. Maximum number of mismatches allowed for a barcode match.
#'   Default: 0 (exact matching only). Set >0 for approximate matching using Hamming distance.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A [data.table::data.table] with columns from the input plus:
#'   - `construct_id`: Mapped construct ID (NA if no match found)
#'   - `match_type`: Type of match ("exact", "mismatch_1bp", etc., or "no_match")
#'   - `match_quality`: Quality of the match (0-1 scale, with 1 being perfect)
#'
#' @details
#' ## Matching Strategy
#'
#' 1. **Exact matching** (default, `max_mismatches = 0`):
#'    Fast lookup using hash table. All extracted barcodes are checked against
#'    the reference library.
#'
#' 2. **Approximate matching** (`max_mismatches > 0`):
#'    Slower but catches barcodes with sequencing errors. Uses Hamming distance
#'    to find the closest matching barcode in the reference library.
#'
#' @examples
#' \dontrun{
#'   # Exact matching (fastest)
#'   matched <- match_constructs(
#'     extracted_barcodes,
#'     reference_library
#'   )
#'
#'   # Allow up to 1 mismatch
#'   matched <- match_constructs(
#'     extracted_barcodes,
#'     reference_library,
#'     max_mismatches = 1
#'   )
#' }
#'
#' @importFrom data.table as.data.table setkey merge.data.table
#' @importFrom dplyr mutate
#' @importFrom stringdist stringdist
#' @export
match_constructs <- function(extracted_barcodes,
                             reference_library,
                             max_mismatches = 0,
                             quiet = FALSE) {
  # Convert to data.table if needed
  if (!data.table::is.data.table(extracted_barcodes)) {
    extracted_barcodes <- data.table::as.data.table(extracted_barcodes)
  }

  # Filter to only successful extractions
  successfully_extracted <- extracted_barcodes[
    extraction_status %in% c("success", "no_sample_barcode")
  ]

  if (nrow(successfully_extracted) == 0) {
    rlang::abort("No successfully extracted barcodes found.")
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Matching {nrow(successfully_extracted)} barcodes to reference library"
    ))
  }

  if (max_mismatches == 0) {
    # Exact matching - fast
    result <- match_constructs_exact(
      successfully_extracted,
      reference_library,
      quiet = quiet
    )
  } else {
    # Approximate matching - slower
    result <- match_constructs_approximate(
      successfully_extracted,
      reference_library,
      max_mismatches = max_mismatches,
      quiet = quiet
    )
  }

   # Combine with failed extractions (mark as no_match)
   failed_extractions <- extracted_barcodes[
     extraction_status %nin% c("success", "no_sample_barcode")
   ]
   failed_extractions[, `:=`(
     construct_id = NA_character_,
     match_type = "no_match",
     match_quality = 0
   )]

  result <- data.table::rbindlist(
    list(result, failed_extractions),
    fill = TRUE
  )

  # Calculate summary statistics
  match_summary <- result[!is.na(construct_id), .(count = .N), by = .(match_type)]

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Matched {sum(!is.na(result$construct_id))} of {nrow(result)} reads to constructs",
      "i" = "Match summary:",
      paste0("  ", match_summary$match_type, ": ", match_summary$count, collapse = "\n")
    ))
  }

  result
}


#' Exact Barcode Matching (Internal)
#'
#' Fast exact matching using hash table lookup.
#' @param extracted_barcodes Data.table. Extracted barcode information.
#' @param reference_library Data.table. Reference construct library.
#' @param quiet Logical. Suppress messages.
#'
#' @param extracted_barcodes Data.table of extracted barcodes
#' @param reference_library Data.table of reference constructs
#' @param quiet Logical. Suppress messages.
#' @export
match_constructs_exact <- function(extracted_barcodes,
                                   reference_library,
                                   quiet = FALSE) {
  # Convert reference to data.table with key
   ref_dt <- data.table::as.data.table(reference_library)
   data.table::setkey(ref_dt, barcode)

   # Rename construct_barcode to barcode for merging
   extracted_barcodes <- extracted_barcodes[, barcode := construct_barcode]

   # Perform the merge
   result <- data.table::merge.data.table(
     extracted_barcodes,
     ref_dt,
     by = "barcode",
     all.x = TRUE,
     sort = FALSE
   )

   # Ensure sample_barcode is preserved (it should be from all.x = TRUE)
   if (!"sample_barcode" %in% names(result)) {
     result[, sample_barcode := extracted_barcodes[match(barcode, extracted_barcodes$barcode), sample_barcode]]
   }

  # Add match information
  result[, match_type := ifelse(
    is.na(construct_id),
    "no_match",
    "exact"
  )]
  result[, match_quality := ifelse(
    is.na(construct_id),
    0,
    1
  )]

  result
}


#' Approximate Barcode Matching (Internal)
#' @param extracted_barcodes Data.table. Extracted barcode information.
#' @param reference_library Data.table. Reference construct library.
#' @param max_mismatches Integer. Maximum hamming distance.
#' @param quiet Logical. Suppress messages.
#'
#' Slower approximate matching allowing for mismatches.
#'
#' @export
match_constructs_approximate <- function(extracted_barcodes,
                                          reference_library,
                                          max_mismatches = 1,
                                          quiet = FALSE) {
  # First try exact matching
  exact_matches <- match_constructs_exact(
    extracted_barcodes,
    reference_library,
    quiet = TRUE
  )

  # Find unmatched barcodes
  unmatched <- exact_matches[is.na(construct_id)]

  if (nrow(unmatched) == 0) {
    # All exact matches
    return(exact_matches)
  }

  # For unmatched, try approximate matching
  ref_barcodes <- unique(reference_library$barcode)
  approx_matches <- lapply(
    seq_len(nrow(unmatched)),
    function(i) {
      query_barcode <- unmatched[[i, "construct_barcode"]]

      # Calculate distances to all reference barcodes
      distances <- stringdist::stringdist(
        query_barcode,
        ref_barcodes,
        method = "hamming"
      )

      min_dist <- min(distances)

      if (min_dist <= max_mismatches) {
        # Found a match
        best_idx <- which.min(distances)
        matched_barcode <- ref_barcodes[best_idx]
        matched_construct <- reference_library$construct_id[
          reference_library$barcode == matched_barcode
        ][1]

        list(
          construct_id = matched_construct,
          match_type = paste0("mismatch_", min_dist, "bp"),
          match_quality = 1 - (min_dist / nchar(query_barcode))
        )
      } else {
        list(
          construct_id = NA_character_,
          match_type = "no_match",
          match_quality = 0
        )
      }
    }
  )

  # Apply approximate match results to unmatched
  for (i in seq_along(approx_matches)) {
    match_info <- approx_matches[[i]]
    unmatched[i, `:=`(
      construct_id = match_info$construct_id,
      match_type = match_info$match_type,
      match_quality = match_info$match_quality
    )]
  }

  # Combine exact and approximate
  result <- data.table::rbindlist(
    list(
      exact_matches[!is.na(construct_id)],
      unmatched
    ),
    fill = TRUE
  )

  result
}


#' Map Constructs to Genes
#'
#' Maps construct barcodes to their target genes using a chip annotation file.
#' Note that one barcode can map to multiple genes (e.g., in lax matching).
#'
#' @param matched_barcodes A data.table from [match_constructs()] with construct mappings.
#' @param chip_file A tibble/data.frame from [load_chip_file()] with columns
#'   `barcode`, `gene_symbol`, and `gene_id`.
#' @param annotation_type Character. Type of annotation to select if multiple are available.
#'   Default: "all" (keep all mappings). Other options: "strict", "lax", etc.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.
#'
#' @return A [data.table::data.table] with input columns plus:
#'   - `gene_symbol`: Target gene symbol (NA if barcode has no gene annotation)
#'   - `gene_id`: Target gene ID (NA if barcode has no gene annotation)
#'   - `barcode_to_genes`: Number of genes mapped to this barcode
#'
#' @details
#' ## Gene Mapping Strategy
#'
#' Each barcode is left-joined to the chip file. If a barcode maps to multiple genes
#' (common in lax chip files), each barcode-gene pair becomes a separate row.
#' This allows downstream analysis to keep guide-level information separate from
#' gene-level aggregation.
#'
#' @examples
#' \dontrun{
#'   # Map to genes
#'   with_genes <- map_constructs_to_genes(
#'     matched_barcodes,
#'     chip_file
#'   )
#'
#'   # Check how many genes each barcode targets
#'   with_genes %>% distinct(barcode, gene_symbol) %>% count(barcode)
#' }
#'
#' @importFrom data.table merge.data.table
#' @importFrom dplyr group_by mutate
#' @export
map_constructs_to_genes <- function(matched_barcodes,
                                     chip_file,
                                     annotation_type = "all",
                                     quiet = FALSE) {
  # Convert to data.table if needed
  if (!data.table::is.data.table(matched_barcodes)) {
    matched_barcodes <- data.table::as.data.table(matched_barcodes)
  }

  # Convert chip file to data.table
  chip_dt <- data.table::as.data.table(chip_file)

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Mapping barcodes to genes using chip file"
    ))
  }

   # Perform left join on barcode
   result <- data.table::merge.data.table(
     matched_barcodes,
     chip_dt,
     by = "barcode",
     all.x = TRUE,
     allow.cartesian = TRUE,  # Allow one-to-many joins
     sort = FALSE
   )

   # Count genes per barcode using data.table for efficiency
   # (avoid dplyr pipe which might drop unexpected columns)
   result <- result[, barcode_to_genes := dplyr::n_distinct(gene_id, na.rm = TRUE), by = barcode]

  # Summary statistics
  n_with_genes <- result[!is.na(gene_id), .N]
  n_no_genes <- result[is.na(gene_id), .N]
  n_1_to_many <- result[barcode_to_genes > 1, .N]

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Gene mapping complete",
      "i" = "Reads with gene assignment: {n_with_genes}",
      "i" = "Reads without gene assignment: {n_no_genes}",
      "i" = "Reads mapping to multiple genes: {n_1_to_many}"
    ))
  }

  result
}
