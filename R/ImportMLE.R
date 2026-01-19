
#' Import MLE Results from MAGeCK-MLE
#'
#' @param MLE_dir Path to the directory containing MLE results folders.
#' @param prefix Prefix pattern to identify MLE result folders.
#' @param extra_prefix Optional prefix to add to the names of the imported dataframes.
#'
#' @returns a list of imported MLE gene summaries with names based on folder names
#' @export

ImportMLE <- function(MLE_dir = ".",
                      prefix = "mageck_MLE_",
                      extra_prefix = NULL) {

  if(MLE_dir == ".") {
    message("Searching in the top layer of working directory. To specify another directory, use argument 'MLE_dir'.")
  }

  if(prefix == "mageck_MLE_") {
    message("MLE results folders will be identified with the prefix 'mageck_MLE_'. To specify another prefix, use argument 'prefix'.")
  }

  mle_dirs <- list.dirs(path = MLE_dir, recursive = FALSE) %>%
    stringr::str_subset(prefix)

  mle_names <- basename(mle_dirs) %>%
    stringr::str_replace_all(prefix, "")

  # Now we cycle over the folders and import each gene summary file

  purrr::map(mle_dirs, function(dir) {

    temp_gene_route <- list.files(dir, pattern = ".*gene_summary.txt$", full.names = TRUE)

    temp_gene_data <-
      readr::read_delim(temp_gene_route,
                        delim= "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE,
                        show_col_types = FALSE)

    ## This step is where we clean the gene names from olfactory receptors and vomeronasal receptors;
    ## also eliminate the annoying | from column names

    temp_gene_data <-
      temp_gene_data %>%
      dplyr::filter(stringr::str_detect(.data$Gene, "^Gm[:digit:]{4,5}|^Olfr[:digit:]*|^Vmn[:digit:]*", negate = TRUE)) %>%
      dplyr::rename_with(.cols = dplyr::everything(), .fn = ~ stringr::str_replace_all(., "\\||\\-", "_"))

  return(temp_gene_data)

  }) %>%
    purrr::set_names(paste0(extra_prefix, mle_names))

}
