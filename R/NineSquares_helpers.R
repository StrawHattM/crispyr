#' Generate cutoffs for Nine Square Plots
#'
#' @param x numerical column / vector.
#' @param scale standard deviation multiplier.
#' @param force_zero_center indicate if you want the cutoffs to be forcefully centered around zero.
#'
#' @returns A vector of length 2 with the lower and upper cutoff
#' @export
#'
#' @examples
#'
#' data <- data.frame(
#'   gene = paste0("Gene", 1:100),
#'   LFC = stats::rnorm(1000, mean = 0, sd = 2)
#' )
#'
#' cutoffs <- NScutoff(data$LFC, scale = 1)
#' cutoffs
NScutoff <- function(x, scale, force_zero_center = FALSE) {
  if (!is.numeric(x)) {
    stop("x has to be numeric")
  }
  mu <- mean(x, na.rm = TRUE)
  stdv <- stats::sd(x, na.rm = TRUE)

  if (force_zero_center == TRUE) {
    mu <- 0
  }

  return(c(mu - stdv * scale, mu + stdv * scale))
}


#' generate cutoffs using NScutoff
#'
#' @param data data frame
#' @param scale standard deviation multiplier.
#' @param force_zero_center indicate which cutoffs to be forcefully centered
#'   around zero. Options are "control", "treatment", "both" or "none"
#'   (default).
#'
#' @returns A list with x_cutoff, y_cutoff and slope_cutoff
#' @export
#'
#' @examples
#'
#' data <- data.frame(
#'   gene = paste0("Gene", 1:100),
#'   control = stats::rnorm(1000, mean = 0, sd = 2),
#'   treatment = stats::rnorm(1000, mean = 0, sd = 2)
#' )
#' cutoffs <- NSgencutoff(data, scale = 1, force_zero_center = "both")
#' cutoffs
NSgencutoff <- function(data, scale, force_zero_center = c("control", "treatment", "both", "none")) {
  if (missing(force_zero_center)) force_zero_center <- "none"

  if (force_zero_center %nin% c("control", "treatment", "both", "none")) {
    stop("force_zero_center has to be one of 'control', 'treatment', 'both' or 'none'")
  }

  if (force_zero_center == "control") {
    x_cutoff <- NScutoff(data$control, scale, force_zero_center = TRUE)
    y_cutoff <- NScutoff(data$treatment, scale)
  } else if (force_zero_center == "treatment") {
    x_cutoff <- NScutoff(data$control, scale)
    y_cutoff <- NScutoff(data$treatment, scale, force_zero_center = TRUE)
  } else if (force_zero_center == "both") {
    x_cutoff <- NScutoff(data$control, scale, force_zero_center = TRUE)
    y_cutoff <- NScutoff(data$treatment, scale, force_zero_center = TRUE)
  } else {
    # none
    x_cutoff <- NScutoff(data$control, scale)
    y_cutoff <- NScutoff(data$treatment, scale)
  }

  slope_cutoff <- NScutoff(data$treatment - data$control, scale)

  return(list(
    x_cutoff = x_cutoff,
    y_cutoff = y_cutoff,
    slope_cutoff = slope_cutoff
  ))
}


#' Generate temporary data frame using input
#'
#' @param data input data frame containing genes, number of sgRNAs per gene,
#'   enrichment scores and p values / p.adjusted / FDRs for all conditions you
#'   want to compare. First column should be by default the gene id, and second
#'   column should be number of sgRNAs
#' @param control <[`tidy-select`][dplyr_tidy_select]> control enrichment score
#'   column, will be plotted on x axis. Must be numerical.
#' @param treament <[`tidy-select`][dplyr_tidy_select]> treatment enrichment
#'   score column, will be plotted on y axis. Must be numerical.
#' @param ctrl_pval <[`tidy-select`][dplyr_tidy_select]> control pvalue /
#'   p.adjusted / FDR column. Must be numerical.
#' @param treat_pval <[`tidy-select`][dplyr_tidy_select]> treatment pvalue /
#'   p.adjusted / FDR column. Must be numerical.
#' @param min_sgrna minimum number of sgRNAs per gene to include in the
#'   analysis. Default is 3.
#'
#' @returns a data frame filtered by min_sgrna with selected columns
#' @export
#'
#' @examples
#' data <- data.frame(
#'   gene = paste0("Gene", 1:100),
#'   num = 4,
#'   untreated_LFC = stats::rnorm(100, mean = 0, sd = 2),
#'   treated_LFC = stats::rnorm(100, mean = 0, sd = 2),
#'   untreated_pval = stats::runif(100, min = 0, max = 1),
#'   treated_pval = stats::runif(100, min = 0, max = 1)
#' )
#'
#' NSbasedf(
#'   data = data,
#'   control = untreated_LFC,
#'   treament = treated_LFC,
#'   ctrl_pval = untreated_pval,
#'   treat_pval = treated_pval,
#'   min_sgrna = 3
#' )

NSbasedf <-
  function(data, control, treament, ctrl_pval, treat_pval, min_sgrna = 3) {
    if (!is.numeric(dplyr::pull(data, {{ control }})) |
      !is.numeric(dplyr::pull(data, {{ control }})) |
      !missing(ctrl_pval) && !is.numeric(dplyr::pull(data, {{ ctrl_pval }})) |
      !missing(treat_pval) && !is.numeric(dplyr::pull(data, {{ treat_pval }}))
    ) {
      stop("non-numeric control, treament, ctrl_pval and treat_pval have to be numeric columns")
    }

    tempdf <-
      data %>%
      dplyr::select(
        id = 1,
        num = 2,
        control = {{ control }},
        treatment = {{ treament }},
        ctrl_pval = {{ ctrl_pval }},
        treat_pval = {{ treat_pval }}
      ) %>%
      dplyr::filter(.data$num >= min_sgrna)

    return(tempdf)
  }


#' Filter data using minimum p value
#'
#' @param data temporary data data frame generated with NSbasedf function
#' @param min_pval minimum p value to filter genes. Default is 0.05.
#'
#' @returns a data frame filtered by min_pval
#' @export
#'
#' @examples
NSfilterpval <- function(data, min_pval = 0.05) {
  if (min_pval < 0 | min_pval > 1) {
    stop("min_pval has to be between 0 and 1")
  }

  tempdf <- data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(min_pv = min(dplyr::c_across(c("ctrl_pval", "treat_pval")), na.rm = TRUE)) %>%
    dplyr::filter(.data$min_pv < min_pval) %>%
    dplyr::select(-c("ctrl_pval", "treat_pval"))

  if (nrow(tempdf) == 0) {
    warning("No genes passed the p-value filter")
  }

  return(tempdf)
}



NSbasegraph <- function(data,
                        x_cutoff,
                        y_cutoff,
                        slope_cutoff,
                        alpha = 0.4,
                        shape = 21) {

  if (missing(x_cutoff) | missing(y_cutoff) | missing(slope_cutoff)) {
    stop("x_cutoff, y_cutoff and slope_cutoff have to be provided")
  }

  tempdf <- data %>%
    dplyr::mutate(
      diff = treatment - control,
      square = dplyr::case_when(
        dplyr::between(treatment - control, slope_cutoff[1], slope_cutoff[2]) ~ "neutral_slope",
        control < x_cutoff[1] & treatment < y_cutoff[1] ~ "bottom_left",
        control > x_cutoff[2] & treatment < y_cutoff[1] ~ "bottom_right",
        dplyr::between(control, x_cutoff[1], x_cutoff[2]) & treatment < y_cutoff[1] ~ "bottom_center",
        dplyr::between(control, x_cutoff[1], x_cutoff[2]) & treatment > y_cutoff[2] ~ "top_center",
        control < x_cutoff[1] & treatment > y_cutoff[2] ~ "top_left",
        control > x_cutoff[2] & treatment > y_cutoff[2] ~ "top_right",
        control < x_cutoff[1] & dplyr::between(treatment, y_cutoff[1], y_cutoff[2]) ~ "middle_left",
        control > x_cutoff[2] & dplyr::between(treatment, y_cutoff[1], y_cutoff[2]) ~ "middle_right",
        dplyr::between(control, x_cutoff[1], x_cutoff[2]) & dplyr::between(treatment, y_cutoff[1], y_cutoff[2]) ~ "center",
        TRUE ~ "center"
      )
    ) %>%
    dplyr::group_by(square) %>%
    dplyr::arrange(desc(abs(diff))) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()

  graph <-
    ggplot2::ggplot(tempdf, ggplot2::aes(x = control, y = treatment)) +
    ggplot2::geom_point(ggplot2::aes(color = square, fill = square), alpha = alpha, shape = shape) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = x_cutoff, linetype = "dashed", color = "grey") +
    ggplot2::geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "grey") +
    ggplot2::geom_abline(slope = 1, intercept = slope_cutoff[1], linetype = "dashed", color = "grey") +
    ggplot2::geom_abline(slope = 1, intercept = slope_cutoff[2], linetype = "dashed", color = "grey") +
    ggplot2::labs(
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::scale_color_manual(
      aesthetics = c("color", "fill"),
      guide = "none",
      values = c(
        "center" = "#D3D3D3",
        "neutral_slope" = "#D3D3D3",
        "bottom_left" = "#004488",
        "bottom_center" = "#0077BB",
        "bottom_right" = "#33BBEE",
        "middle_left" = "#CCAA00",
        "middle_right" = "#00AA88",
        "top_left" = "#EE7733",
        "top_center" = "#CC3377",
        "top_right" = "#882255"
      )
    )

  return(graph)
}


NSaxislabels <- function(graph, xlab, ylab) {

  graph <- graph +
    ggplot2::labs(
      x = xlab,
      y = ylab
    )
  return(graph)
}


NStitle <- function(graph, title) {

  graph <- graph +
    ggplot2::ggtitle(title)

  return(graph)
}




#' Add Gene of Interest highlights
#'
#' @param data temporary data data frame generated with NSbasedf function
#' @param gene_list vector of gene ids to filter the data frame
#' @param graph base graph to which add highlights
#'
#' @returns a graph with added highlights
#' @export
#'
#' @examples
NSaddgoi <- function(data, graph, gene_list) {

  tempdf <- data %>%
    dplyr::filter(.data$id %in% gene_list)

  return(tempdf)
}
