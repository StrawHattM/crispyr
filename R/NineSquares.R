


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
                        slopecut) {

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
    if (class(cutoffs$x_cutoff) != "numeric" | length(cutoffs$x_cutoff) != 2) {
      stop("Argument 'xcut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }
    cutoffs$x_cutoff <- xcut
  }

  if (!missing(ycut)) { ## check for user ycut input
    if (class(cutoffs$y_cutoff) != "numeric" | length(cutoffs$y_cutoff) != 2) {
      stop("Argument 'ycut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }

    cutoffs$y_cutoff <- ycut
  }

  if (!missing(slopecut)) { ## check for user slopecut input
    if (class(cutoffs$slope_cutoff) != "numeric" | length(cutoffs$slope_cutoff) != 2) {
      stop("Argument 'slopecut' must be a numerical vector of length two. Leave empty for default calculation using 'scale'")
    }

    cutoffs$slope_cutoff <- slopecut
  }


  # p-value filtering
  if (!missing(min_pval)) { ## trigger is optional argument min_pval existing
    base_data <-
      NSfilterpval(base_data, min_pval = min_pval)
  }

  # Generate base graph
  base_graph <-
    NSbasegraph(
      data = base_data,
      x_cutoff = cutoffs$x_cutoff,
      y_cutoff = cutoffs$y_cutoff,
      slope_cutoff = cutoffs$slope_cutoff,
      alpha = alpha,
      shape = shape
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
}
