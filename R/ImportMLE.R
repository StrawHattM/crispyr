

#' Import MLE Results from MAGeCK-MLE
#'
#' @param MLE_dir Path to the directory containing MLE results folders.
#' @param prefix Prefix pattern to identify MLE result folders.
#' @param extra_prefix Optional prefix to add to the names of the imported dataframes.
#'
#' @returns nothing, but assigns dataframes to the global environment with names mle_<comparison>
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

  MLE_summaries <-
    list.dirs(path = MLE_dir, recursive = FALSE) %>%
    stringr::str_subset(prefix) %>%
    list.files(pattern = "gene_summary", full.names = TRUE)

  # Now we cycle over the folders and import each gene summary file

  for (i in seq_along(MLE_summaries)) {

    temp_gene_data <-
      readr::read_delim(MLE_summaries[i], delim= "\t", escape_double = FALSE, trim_ws = TRUE)

    ## This step is where we clean the gene names from olfactory receptors and vomeronasal receptors;
    ## also eliminate the annoying | from column names

    temp_gene_data <-
      temp_gene_data %>%
      dplyr::filter(stringr::str_detect(.data$Gene, "^Gm[:digit:]{4,5}|^Olfr[:digit:]*|^Vmn[:digit:]*", negate = TRUE)) %>%
      dplyr::rename_with(.cols = dplyr::everything(), .fn = ~ stringr::str_replace_all(., "\\||\\-", "_"))

    ## Assign giving the user the possibility of adding a prefix to avoid overwriting previous imports.
    ## If no prefix is given, NULL is used and the name is just mle_<comparison>

    assign(x = paste0("mle_", extra_prefix, stringr::str_replace_all(basename(dirname(MLE_summaries[i])), prefix, "")),
           value = temp_gene_data,
           envir = .GlobalEnv)

  }

}
