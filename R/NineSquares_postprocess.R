
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

NSextractdata <- function(graph_list, which = "data") {

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

  result <-
    if ("data" %in% which) {
      if (length(which) > 1) {
        stop("If 'data' is selected, it must be the only option selected.")
      }

      purrr::map(graph_list, function(graph) {

        graph$data

      }) %>%
        purrr::set_names(names)


    } else {

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


