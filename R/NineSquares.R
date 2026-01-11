


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
                        top_labeled = 10,
                        force_zero_center = c("none", "both", "control", "treatment"),
                        xlab,
                        ylab,
                        title,
                        xcut,
                        ycut,
                        slopecut,
                        legend,
                        filename,
                        groups_labeled = c("top_center", "bottom_center", "middle_right", "middle_left"),
                        goi,
                        goi_auto,
                        goi_color,
                        goi_size,
                        goi_label_color,
                        goi_label_size) {

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

      if(tolower(stringr::str_extract(filename, pattern = ".[:alpha:]{3}$")) != ".txt"){
        warning("Your filename extension was detected to be different than .txt -
                the file generated will be a tab-separated file that should be encoded as a .txt")
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

  base_graph <-
    NSaxislabels(base_graph, xlab, ylab)

  # Add title
  if(!missing(title)){
    base_graph <-
      NStitle(base_graph, title)
  }

  # Add labels for top genes in each square

  if(!is.character(groups_labeled)) {
    stop("'groups_labeled' needs to be a character vector.")
  }


  # Add genes of interest highlights




}
