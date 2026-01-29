
#' Extract data from Nine Squares plot object and assign it
#'
#' @param graph_list A Nine Squares plot object or list of Nine Squares plot objects
#' @param which Character vector specifying which data to extract. Options include:
#'  \code{"data"} (the full data frame used to generate the plot),
#'  \code{"bottom_center"}, \code{"bottom_left"},
#'  \code{"bottom_right"}, \code{"top_left"}, \code{"top_right"},
#'  \code{"top_center"}, \code{"middle_right"}, \code{"middle
#'  _left"}, \code{"center"}, \code{"neutral_slope"}.
#'  If one of the square positions is selected, the function will return
#'  the IDs of the genes in that square. If \code{"data"} is
#'  selected, the full data frame used to generate the plot will be returned.
#'
#' @returns the data frame used to generate the plot
#' @export

extractNSdata <- function(graph_list, which = "data") {

  valid_choices <- c("data", "bottom_center", "bottom_left",
                     "bottom_right", "top_left", "top_right",
                     "top_center", "middle_right", "middle_left",
                     "center", "neutral_slope")

  which <- match.arg(which, choices = valid_choices, several.ok = TRUE)

  single_graph <- inherits(graph_list, "gg")

  if (single_graph) {
    graph_list <- list(graph_list)
  }

  names <- names(graph_list)


  if ("data" %in% which) {
    if (length(which) > 1) {
      stop("If 'data' is selected, it must be the only option selected.")
    }
    result <-
      purrr::map(graph_list, function(graph) {

        graph$data

      }) %>%
      purrr::set_names(names)


  } else {
    result <-
      purrr::map(graph_list, function(graph) {
        graph$data %>%
          dplyr::filter(.data$square %in% which) %>%
          dplyr::arrange(.data$square, .data$rank) %>%
          dplyr::pull(.data$id)
      }) %>%
      purrr::set_names(names)

  }

  if (single_graph) {
    return(result[[1]])
  } else {
    return(result)
  }

}



#' Locate genes in grouped data
#'
#' Find which groups or squares contain specified genes. Accepts either a list
#' of named character vectors (e.g., from NSextractdata results) or NineSquares
#' dataframe output.
#'
#' @param genes Character vector of gene identifiers to search for.
#' @param data One of the following:
#'   \itemize{
#'     \item A named list of character vectors (each containing gene IDs)
#'     \item A dataframe with columns `id` and `square` (NineSquares format)
#'     \item A ggplot object (data will be extracted from the `data` slot)
#'   }
#' @param return_full Logical. If `TRUE` and `data` is a list, returns the full
#'   vectors containing each gene rather than just the vector names. Ignored for
#'   NineSquares data. Default: `FALSE`.
#'
#' @return For single gene input: a character vector of group/square names where
#'   the gene was found, or `character(0)` if not found. For multiple genes: a
#'   named list where each element contains the character vector of groups/squares
#'   for that gene.
#'
#' @export

where_is <- function(genes, data, return_full = FALSE) {

  # Extract data from ggplot if needed
  if (inherits(data, "gg")) {
    data <- data$data
  }

  .find_one <- function(gene) {
    # NineSquares: dataframe with square column
    if (is.data.frame(data) && "square" %in% names(data)) {
      data$square[data$id == gene]
    } else {
      # List of vectors
      hits <- vapply(data, function(x) gene %in% x, logical(1))
      found_in <- names(data)[hits]
      if (return_full) data[found_in] else found_in
    }
  }

  results <- purrr::map(purrr::set_names(genes), .find_one)

  if (length(genes) == 1L) results[[1L]] else results
}
