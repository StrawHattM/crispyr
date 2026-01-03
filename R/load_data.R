#' Load Reference Library from CSV
#'
#' Loads a CSV file containing construct barcodes and their associated construct IDs.
#' This is typically the reference file used for barcode-to-construct mapping during
#' PoolQ deconvolution.
#'
#' @param ref_file Path to the reference CSV file. Expected columns: 
#'   barcode sequence, construct ID, and a duplicate barcode sequence column.
#' @param quiet Logical. If `TRUE`, suppress informative messages. Default is `FALSE`.
#'
#' @return A [tibble::tibble] with columns:
#'   - `barcode`: 20bp construct barcode sequence
#'   - `construct_id`: Unique construct identifier (e.g., BRDN0001468998)
#'
#' @examples
#' \dontrun{
#'   ref_lib <- load_reference_library("path/to/CP0045_reference_20160120.csv")
#'   head(ref_lib)
#' }
#'
#' @importFrom readr read_csv
#' @importFrom dplyr select rename
#' @export
load_reference_library <- function(ref_file, quiet = FALSE) {
  if (!file.exists(ref_file)) {
    rlang::abort(c(
      "Reference file not found.",
      "x" = glue::glue("File: {ref_file}")
    ))
  }

  ref_lib <- readr::read_csv(
    ref_file,
    col_names = c("barcode", "construct_id", "barcode_dup"),
    col_types = "ccc",
    progress = !quiet
  ) %>%
    dplyr::select(barcode, construct_id)

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Loaded reference library with {nrow(ref_lib)} construct barcodes"
    ))
  }

  ref_lib
}


#' Load Gene Annotation (Chip) File
#'
#' Loads a tab-separated chip file containing the mapping between construct barcodes
#' and their target genes. Multiple annotations (strict, lax, etc.) can be loaded.
#'
#' @param chip_file Path to the chip TSV file. Expected columns:
#'   barcode sequence, gene symbol, and gene ID.
#' @param quiet Logical. If `TRUE`, suppress informative messages. Default is `FALSE`.
#'
#' @return A [tibble::tibble] with columns:
#'   - `barcode`: 20bp construct barcode sequence
#'   - `gene_symbol`: Target gene symbol (e.g., "Defb34")
#'   - `gene_id`: NCBI gene ID (e.g., "360211")
#'
#' @details
#' Note that one barcode can map to multiple genes in some chip files (lax matching).
#' Each barcode-gene pair will be a separate row in the output.
#'
#' @examples
#' \dontrun{
#'   chip <- load_chip_file("path/to/CP0045_GRCm38_NCBI_CRISPRko_strict_gene_20221201.chip")
#'   head(chip)
#' }
#'
#' @importFrom readr read_tsv
#' @export
load_chip_file <- function(chip_file, quiet = FALSE) {
  if (!file.exists(chip_file)) {
    rlang::abort(c(
      "Chip file not found.",
      "x" = glue::glue("File: {chip_file}")
    ))
  }

  chip <- readr::read_tsv(
    chip_file,
    col_types = "ccc",
    progress = !quiet
  ) %>%
    dplyr::rename(
      barcode = "Barcode Sequence",
      gene_symbol = "Gene Symbol",
      gene_id = "Gene ID"
    )

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Loaded chip file with {nrow(chip)} barcode-gene mappings"
    ))
  }

  chip
}


#' Load Sample Manifest
#'
#' Loads a CSV file containing the sample barcodes and their associated condition names.
#' This file maps the 8bp sample barcodes found in reads to experiment conditions/samples.
#'
#' @param manifest_file Path to the sample manifest CSV file. Expected columns:
#'   sample barcode (8bp) and condition/sample name.
#' @param col_names Character vector of length 2 specifying column names.
#'   Default is `c("sample_barcode", "condition")`.
#' @param quiet Logical. If `TRUE`, suppress informative messages. Default is `FALSE`.
#'
#' @return A [tibble::tibble] with columns:
#'   - `sample_barcode`: 8bp sample barcode sequence
#'   - `condition`: Sample condition or name (e.g., "Day0", "UT_r1")
#'
#' @examples
#' \dontrun{
#'   manifest <- load_sample_manifest(
#'     "path/to/GPP-6350_Martin Gonzalez Fernandez_conditions.csv"
#'   )
#'   head(manifest)
#' }
#'
#' @importFrom readr read_csv
#' @export
load_sample_manifest <- function(manifest_file,
                                  col_names = c("sample_barcode", "condition"),
                                  quiet = FALSE) {
  if (!file.exists(manifest_file)) {
    rlang::abort(c(
      "Manifest file not found.",
      "x" = glue::glue("File: {manifest_file}")
    ))
  }

  if (length(col_names) != 2) {
    rlang::abort("col_names must be a character vector of length 2")
  }

  manifest <- readr::read_csv(
    manifest_file,
    col_names = col_names,
    col_types = "cc",
    progress = !quiet
  )

  if (!quiet) {
    cli::cli_inform(c(
      "*" = "Loaded sample manifest with {nrow(manifest)} samples"
    ))
  }

  manifest
}
