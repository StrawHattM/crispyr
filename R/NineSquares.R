


#' Create and annotate a Nine Squares Plot
#'
#' @param data input data frame containing genes, number of sgRNAs per gene,
#'   enrichment scores and p values / p.adjusted / FDRs for all conditions you
#'   want to compare. First column should be by default the gene id, and second
#'   column should be number of sgRNAs.
#' @param control <[`tidy-select`][dplyr_tidy_select]> control enrichment score
#'   column, will be plotted on x axis. Must be numerical.
#' @param treatment <[`tidy-select`][dplyr_tidy_select]> treatment enrichment
#'   score column, will be plotted on y axis. Must be numerical.
#' @param ctrl_pval <[`tidy-select`][dplyr_tidy_select]> control pvalue /
#'   p.adjusted / FDR column. Must be numerical.
#' @param treat_pval <[`tidy-select`][dplyr_tidy_select]> treatment pvalue /
#'   p.adjusted / FDR column. Must be numerical.
#' @param min_sgrna numeric, minimum number of sgRNAs per gene to include in the
#'   analysis. Default is 3.
#' @param min_pval numeric, p-value threshold to apply to whether genes appear
#'   in the graph or not.
#' @param scale numeric, standard deviation multiplier to calculate thresholds.
#'   Default is 2.
#' @param alpha numeric, transparency of the points in the plot. Default is 0.4.
#' @param shape numeric, shape of the points in the plot. Default is 21.
#' @param top_labeled numeric, indicates how many top genes (by euclidean
#'   distance) in each group/square to label. Default is 5.
#' @param force_zero_center character, one of "none", "both", "control", or
#'   "treatment". Indicates if cutoffs should be forced to center around zero
#'   instead of the mean.
#' @param xlab character, label for the x axis. Default is "Control Enrichment
#'   Score".
#' @param ylab character, label for the y axis. Default is "Treatment Enrichment
#'   Score".
#' @param title character, graph title
#' @param xcut 2-length numerical vector indicating where cutoffs should be
#'   established for the x axis. Replaces default calculation, which is mean ±
#'   scale * standard deviation of control
#' @param ycut 2-length numerical vector indicating where cutoffs should be
#'   established for the y axis. Replaces default calculation, which is mean ±
#'   scale * standard deviation of treatment
#' @param slopecut 2-length numerical vector indicating where cutoffs should be
#'   established for the x axis. Replaces default calculation, which is mean ±
#'   scale * standard deviation of treatment - control
#' @param legend logical, indicate whether legend should be plotted. Default FALSE
#' @param filename character, path/name of the file to be written containing the data.
#' @param groups_labeled groups to have their top genes labeled. Defaults to "top_center", "bottom_center", "middle_right" and "middle_left"
#' @param goi character vector of *g*enes *o*f *i*nterest to be highlighted.
#' @param goi_auto logical, indicates whether highlighted genes of interest are coloured according to the square they are in.
#' @param goi_shape numeric, shape of the highlighted genes of interest in the plot. Default is 21.
#' @param goi_color character, color aesthetic of the highlighted genes of interest in the plot. Default is black.
#' @param goi_fill character, fill aesthetic of the highlighted genes of interest in the plot. Default is black. goi_auto needs to be set to FALSE.
#' @param goi_size size, shape of the highlighted genes of interest in the plot. Default is 4.
#' @param goi_label_type character, one of "label" or "text". Indicates whether the annotation is text or a label. Defaults to "label".
#' @param goi_label_color character, color of the label text, default is "black".
#' @param goi_label_size character, size of the label text, default is 4.
#'
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#'
#' data <- data.frame(
#' gene = paste0("Gene", 1:20000),
#'   num = sample(1:10, 20000, replace = TRUE),
#'   untreated_LFC = stats::rnorm(20000, mean = 0, sd = 3),
#'   treated_LFC = stats::rnorm(20000, mean = 0, sd = 3),
#'   untreated_pval = stats::runif(20000, min = 0, max = 1),
#'   treated_pval = stats::runif(20000, min = 0, max = 1)
#' )
#'
#' NineSquares(data,
#'   untreated_LFC,
#'   treated_LFC,
#'   untreated_pval,
#'   treated_pval,
#'   min_pval = 0.05,
#'   scale = 2,
#'   shape = 15,
#'   top_labeled = 10,
#'   xlab = "x axis label",
#'   ylab = "x axis label",
#'   title = "Example Title",
#'   groups_labeled = c("top_center", "bottom_center"),
#'   goi = c("Gene2994", "Gene405", "Gene1450", "Gene592", "Gene14261"))
#'   }

NineSquares <- function(data,
                        control,
                        treatment,
                        ctrl_pval,
                        treat_pval,
                        min_sgrna = 3,
                        min_pval,
                        scale = 2,
                        alpha = 0.4,
                        shape = 21,
                        top_labeled = 5,
                        force_zero_center = c("none", "both", "control", "treatment"),
                        xlab,
                        ylab,
                        title,
                        xcut,
                        ycut,
                        slopecut,
                        legend = FALSE,
                        filename,
                        groups_labeled = c("top_center", "bottom_center", "middle_right", "middle_left"),
                        goi,
                        goi_auto = TRUE,
                        goi_shape = 21,
                        goi_color = "black",
                        goi_fill = "black",
                        goi_size = 2.5,
                        goi_label_type = c("label", "text"),
                        goi_label_color = "black",
                        goi_label_size = 4) {

  # Build base data frame
  base_data <-
    NSbasedf(data, {{control}}, {{treatment}}, {{ctrl_pval}}, {{treat_pval}}, min_sgrna)

  # Generate cutoffs before p-val filtering, that way they're based on the whole population's distribution and not the filtered one
  ## Because it's quite a clunky arg otherwise, I'll just default it to none here
  if (missing(force_zero_center)) {
    force_zero_center <- "none"
  }

  cutoffs <-
    NSgencutoff(base_data, scale, force_zero_center)

  if (!missing(xcut)) {  ## check for user xcut input
    if (!is.numeric(cutoffs$x_cutoff) | length(cutoffs$x_cutoff) != 2) {
      stop("Argument 'xcut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }
    cutoffs$x_cutoff <- xcut
  }

  if (!missing(ycut)) { ## check for user ycut input
    if (!is.numeric(cutoffs$y_cutoff) | length(cutoffs$y_cutoff) != 2) {
      stop("Argument 'ycut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }

    cutoffs$y_cutoff <- ycut
  }

  if (!missing(slopecut)) { ## check for user slopecut input
    if (!is.numeric(cutoffs$slope_cutoff) | length(cutoffs$slope_cutoff) != 2) {
      stop("Argument 'slopecut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }

    cutoffs$slope_cutoff <- slopecut
  }


  # p-value filtering before squares so rank isn't affected
  if (!missing(min_pval)) { ## trigger is optional argument min_pval existing

    base_data <-
      NSfilterpval(base_data, min_pval = min_pval)

  }

  # Determine squares and assign in new col

  base_data <- NSsquares(base_data,
                         x_cutoff = cutoffs$x_cutoff,
                         y_cutoff = cutoffs$y_cutoff,
                         slope_cutoff = cutoffs$slope_cutoff)

  # Save data

  if(!missing(filename)){

    if(is.character(filename)) {

      if(!stringr::str_detect(filename, pattern = "\\.[:alpha:]{3}$")){
        warning("No filename extension detected - default to '.txt'")

        filename <- paste0(filename, ".txt")

      }

      if(tolower(stringr::str_extract(filename, pattern = "\\.[:alpha:]{3}$")) != ".txt"){
        warning("Your filename extension was detected to be different than .txt -
                the file generated will be a tab-separated file encoded as a .txt")

        filename <- stringr::str_replace(filename, pattern = "\\.[:alpha:]{3}$", replacement = ".txt")

      }

      filename <- filename

    } else if (filename & !missing(title)){

      filename <- paste0(getwd(), "/", title, ".txt")

    } else {stop("filename needs to be character")}

    readr::write_delim(x = base_data,
                       file = filename,
                       quote = "none",
                       delim = "\t",
                       col_names = TRUE)
  }


  # Generate base graph
  base_graph <-
    NSbasegraph(
      data = base_data,
      alpha = alpha,
      shape = shape,
      legend = legend,
      x_cutoff = cutoffs$x_cutoff,
      y_cutoff = cutoffs$y_cutoff,
      slope_cutoff = cutoffs$slope_cutoff
    )

  # Add axis labels with defaults
  if(missing(xlab)){
    xlab <- "Control Enrichment Score"
  }
  if(missing(ylab)){
    ylab <- "Treatment Enrichment Score"
  }

  graph <-
    NSaxislabels(base_graph, xlab, ylab)

  # Add title
  if(!missing(title)){
    graph <-
      NStitle(graph, title)
  }

  # Add labels for top genes in each square

  if(!is.character(groups_labeled)) {
    stop("'groups_labeled' needs to be a character vector.")
  }

  graph <-
    NSaddtoplabels(
      graph = graph ,
      data = base_data,
      groups_labeled = groups_labeled,
      top_labeled = top_labeled
    )

  # Add genes of interest highlights

  if (!missing(goi)) {

    graph <-
      NSaddgoi(
        data = base_data,
        graph = graph,
        goi_list = goi,
        goi_auto = goi_auto,
        goi_shape = goi_shape,
        goi_color = goi_color,
        goi_fill = goi_fill,
        goi_size = goi_size,
        goi_label_type = match.arg(goi_label_type),
        goi_label_color = goi_label_color,
        goi_label_size = goi_label_size
      )

  }

  return(graph)

}
