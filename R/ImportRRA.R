#' Import RRA Results
#'
#' @param RRA_dir Path to the directory containing RRA results folders.
#' @param extra_prefix Optional prefix to add to the names of the imported dataframes.
#'
#' @returns nothing, but assigns dataframes to the global environment with names rra_<comparison>
#' @export


ImportRRA <- function(RRA_dir, extra_prefix = NULL) {

  # Per generated RRA script structure, each folder in RRA_dir specifies which comparison is inside
  rra_comparisons <- list.files(RRA_dir) %>% stringr::str_replace_all("-", "_")

  # Now we just cycle over the folders and import each gene summary file

  for (i in 1:length(rra_comparisons)) {

    temp_gene_route <-
      list.files(list.dirs(RRA_dir)[i+1], pattern = ".*gene_summary.txt$", full.names = TRUE)

    temp_gene_data <-
      readr::read_delim(temp_gene_route, delim= "\t", escape_double = FALSE, trim_ws = TRUE)

    ## This step is where we clean the gene names from olfactory receptors and vomeronasal receptors;
    ##also eliminate the annoying | from column names

    temp_gene_data <-
      temp_gene_data %>%
      dplyr::filter(stringr::str_detect(.data$id, "^Gm[:digit:]{4,5}|^Olfr[:digit:]*|^Vmn[:digit:]*", negate = TRUE)) %>%
      dplyr::rename_with(.cols = dplyr::everything(), .fn = ~ stringr::str_replace_all(., "\\||\\-", "_"))

    ## Assign giving the user the possibility of adding a prefix to avoid overwriting previous imports.
    ## If no prefix is given, NULL is used and the name is just rra_<comparison>

    assign(x = paste0("rra_", NULL, rra_comparisons[i]),
           value = temp_gene_data,
           envir = .GlobalEnv)

  }

}




#' Build dataframe of RRA Day 0 Comparisons
#'
#' @param objects Optional character vector of RRA day 0 object names to include.
#'
#' @returns a data frame with the LFC, minimum pvalue and FDR for each gene across all provided RRA day 0 objects
#' @export


BuildRRAdz <- function(objects = NULL, order = NULL) {

  if(is.null(objects)) {
    # Get all RRA day 0 objects
    day0_objs <- stringr::str_subset(ls(envir = .GlobalEnv), pattern = "^rra_.*[dD]|[dD]ay0_.+^")
  } else {
    day0_objs <- objects
  }

  if (length(day0_objs) == 0) {
    stop("No RRA day 0 objects found matching pattern '^rra_.*[dD]|[dD]ay0_.*'")
  }

  # Obtain label from each comparison for further prefix
  prefixes <-
    stringr::str_replace_all(day0_objs, "^rra_.*([dD]|[dD]ay)0_", "")

  # Generate a list of built dataframes with minimum pvalue and fdr for each gene

  df_list<-
    lapply(seq_along(day0_objs), function(i) {

      temp_df <-
        get(day0_objs[i], envir = .GlobalEnv) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(min_pval = min(dplyr::c_across(dplyr::contains("p_value"))),
                      min_fdr = min(dplyr::c_across(dplyr::contains("fdr")))) %>%
        dplyr::select(.data$id, .data$num, .data$pos_lfc, .data$min_pval, .data$min_fdr) %>%
        dplyr::rename_with(.cols = 2:5, .fn = ~ stringr::str_replace_all(., "pos|min", prefixes[i]))

      return(temp_df)
    })

  # Now we merge them all together by id and num
  dz_df <- purrr::reduce(df_list, dplyr::full_join, by = c("id", "num"))


  # Reorder columns if order parameter is provided
  if(!is.null(order)) {

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
    dz_df <- dz_df %>%
      dplyr::select(
        .data$id,
        .data$num,
        !!!purrr::map(order_levels, ~ dplyr::starts_with(paste0(.x, "_")))
      )

  }

  return(dz_df)

}
