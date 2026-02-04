#' Import RRA Results
#'
#' @param RRA_dir Path to the directory containing RRA results folders.
#' @param extra_prefix Optional prefix to add to the names of the imported dataframes.
#'
#' @returns a list of imported RRA gene summaries with names extra_prefix<comparison>
#' @export

ImportRRA <- function(RRA_dir, extra_prefix = NULL) {

  # Per generated RRA script structure, each folder in RRA_dir specifies which comparison is inside
  rra_comparisons <- list.files(RRA_dir) %>% stringr::str_replace_all("-", "_")

  dirs <- list.dirs(RRA_dir)[-1]

  purrr::map(dirs, function(dir) {

    temp_gene_route <-
      list.files(dir, pattern = ".*gene_summary.txt$", full.names = TRUE)

    temp_gene_data <-
      readr::read_delim(temp_gene_route,
                        delim= "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE,
                        show_col_types = FALSE)

    ## This step is where we clean the gene names from olfactory receptors and vomeronasal receptors;
    ##also eliminate the annoying | from column names

    temp_gene_data <-
      temp_gene_data %>%
      dplyr::filter(
        stringr::str_detect(.data$id,
                            "^Gm[:digit:]{4,5}|^Olfr[:digit:]*|^Vmn[:digit:]*",
                            negate = TRUE)
        ) %>%
      dplyr::rename_with(.cols = dplyr::everything(),
                         .fn = ~ stringr::str_replace_all(., "\\||\\-", "_"))

    return(temp_gene_data)

  }) %>%
    purrr::set_names(paste0(extra_prefix, rra_comparisons))

}


#' Build dataframe of RRA Day 0 Comparisons
#'
#' @param rra_list List of RRA dataframes as imported with ImportRRA()
#' @param pattern Regex pattern to identify RRA day 0 comparisons (or desired
#'   comparisons) in the rra_list names. Default matches names starting with "d"
#'   or "D" or containing "day0_", as `"^([dD]|[dD]ay)0_"`.
#' @param order Optional character vector or factor specifying the order of
#'   comparisons in the output dataframe.
#'
#' @returns a data frame with the LFC, minimum pvalue and FDR for each gene
#'   across all provided RRA day 0 objects
#' @export


BuildRRAdz <- function(rra_list, pattern = "^([dD]|[dD]ay)0_", order = NULL) {

  if (is.null(names(rra_list))) {
    stop("rra_list must be a named list")
  }

  if(!is.null(pattern)) {

    matching_names <- stringr::str_subset(names(rra_list), pattern = pattern)

    summ_list <- rra_list[matching_names]

    # Use the pattern to extract prefixes by removing it from names
    prefixes <- stringr::str_replace_all(names(summ_list), pattern, "")

  } else {
    summ_list <- rra_list

    # If no pattern, use full names as prefixes
    prefixes <- names(summ_list)

  }

  if (length(summ_list) == 0) {
    stop(paste0("No RRA objects found matching pattern ", pattern))
  }


  # Generate a list of built dataframes with minimum pvalue and fdr for each gene

  df_list <-
    purrr::map2(summ_list, prefixes, function(df, prefix) {

      df %>%
        as.data.frame() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(min_pval = min(dplyr::c_across(dplyr::contains("p_value"))),
                      min_fdr = min(dplyr::c_across(dplyr::contains("fdr")))) %>%
        dplyr::select("id", "num", "pos_lfc", "min_pval", "min_fdr") %>%
        dplyr::rename_with(
          .cols = 2:5,
          .fn = ~ stringr::str_replace_all(., "pos|min", prefix)
        )
    })

  # Now we merge them all together by id and num

  df_collapsed <- purrr::reduce(df_list, dplyr::full_join, by = c("id", "num")) %>% tidyr::drop_na()


  # Reorder columns if order parameter is provided
  if (!is.null(order)) {

    if(is.factor(order)) {
      order_levels <- levels(order)
    }
    else if(is.character(order)) {
      order_levels <- order
    }
    else {
      stop("order parameter must be a character vector or factor to order by.")
    }

    # Reorder columns
    df_collapsed <-
      df_collapsed %>%
      dplyr::select(
        "id",
        "num",
        !!!purrr::map(order_levels, ~ dplyr::starts_with(paste0(.x, "_")))
      )

  }

  return(df_collapsed)

}
