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
#' @param min_pval numerical, minimum p value to filter genes. Default is 0.05.
#'
#' @returns a data frame filtered by min_pval
#' @export

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

#' Create and assign Squares in a dataset using cutoffs
#'
#' @param data data frame containing treatment & control columns, as generated
#'   by NSbasegraph
#' @param x_cutoff numerical length-2 vector with x axis cutoff values, as from
#'   NSgencutoff
#' @param y_cutoff numerical length-2 vector with y axis cutoff values, as from
#'   NSgencutoff
#' @param slope_cutoff numerical length-2 vector with diagonal axis cutoff
#'   values, as from NSgencutoff
#'
#' @returns a data frame containing three new columns, diff, square, and rank (=
#'   highest diff per square)
#' @export

NSsquares <- function(data, x_cutoff, y_cutoff, slope_cutoff){

  if (missing(x_cutoff) | missing(y_cutoff) | missing(slope_cutoff)) {
    stop("x_cutoff, y_cutoff and slope_cutoff have to be provided")
  }

  tempdf <- data %>%
    dplyr::mutate(
      diff = .data$treatment - .data$control,
      square = dplyr::case_when(
        dplyr::between(.data$treatment - .data$control, slope_cutoff[1], slope_cutoff[2]) ~ "neutral_slope", # important to be first
        .data$control < x_cutoff[1] & .data$treatment < y_cutoff[1] ~ "bottom_left",
        .data$control > x_cutoff[2] & .data$treatment < y_cutoff[1] ~ "bottom_right",
        dplyr::between(.data$control, x_cutoff[1], x_cutoff[2]) & .data$treatment < y_cutoff[1] ~ "bottom_center",
        dplyr::between(.data$control, x_cutoff[1], x_cutoff[2]) & .data$treatment > y_cutoff[2] ~ "top_center",
        .data$control < x_cutoff[1] & .data$treatment > y_cutoff[2] ~ "top_left",
        .data$control > x_cutoff[2] & .data$treatment > y_cutoff[2] ~ "top_right",
        .data$control < x_cutoff[1] & dplyr::between(.data$treatment, y_cutoff[1], y_cutoff[2]) ~ "middle_left",
        .data$control > x_cutoff[2] & dplyr::between(.data$treatment, y_cutoff[1], y_cutoff[2]) ~ "middle_right",
        dplyr::between(.data$control, x_cutoff[1], x_cutoff[2]) & dplyr::between(.data$treatment, y_cutoff[1], y_cutoff[2]) ~ "center",
        TRUE ~ "center" # I don't think this exists but just in case
      )
    ) %>%
    dplyr::group_by(.data$square) %>%
    dplyr::mutate(euclid_dist = sqrt(.data$control^2 + .data$treatment^2), #calculate euclidean distance to 0
                  perp_dist = abs(.data$diff) / sqrt(2), # perpendicular distance to the diagonal
                  score = .data$perp_dist * (1 + .data$euclid_dist) # composite score
                  ) %>%
    dplyr::arrange(dplyr::desc(.data$score)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>% # Because it's grouped and sorted, row order = rank
    dplyr::ungroup()

  return(tempdf)
}


#' Generate base graph from preprocessed data frame for Nine Squares plot
#'
#' @param data data frame as generated from NSbasegraph and processed with
#'   NSsquares, containing square information and rank per square, as well as
#'   control and treatment enrichment value columns. Can be filtered by
#'   min-pval.
#' @param alpha numerical, transparency of geom_point, default is 0.4
#' @param shape numerical, shape of geom_point, default is 21
#' @param size numerical, size of geom_point, default is 2
#' @param legend logical, indicates whether you want square to be displayed in
#'   the legend
#' @param x_cutoff numerical length-2 vector with x axis cutoff values, as from
#'   NSgencutoff
#' @param y_cutoff numerical length-2 vector with y axis cutoff values, as from
#'   NSgencutoff
#' @param slope_cutoff numerical length-2 vector with diagonal axis cutoff
#'   values, as from NSgencutoff
#'
#' @returns a base Nine Squares graph
#' @export

NSbasegraph <- function(data,
                        x_cutoff,
                        y_cutoff,
                        slope_cutoff,
                        size = 2,
                        alpha = 0.4,
                        shape = 21,
                        legend = FALSE) {

  if(!is.logical(legend)){stop("legend needs to be TRUE or FALSE")}

  graph <-
    ggplot2::ggplot(data, ggplot2::aes(x = .data$control, y = .data$treatment)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$square, fill = .data$square),
                        alpha = alpha,
                        shape = shape,
                        size = size) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = x_cutoff, linetype = "dashed", color = "grey60") +
    ggplot2::geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "grey60") +
    ggplot2::geom_abline(slope = 1, intercept = slope_cutoff[1], linetype = "dashed", color = "grey60") +
    ggplot2::geom_abline(slope = 1, intercept = slope_cutoff[2], linetype = "dashed", color = "grey60") +
    ggplot2::labs(
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::scale_color_manual(
      aesthetics = c("color", "fill"),
      guide = ifelse(legend, "legend", "none"),
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


#' Add labels to a graph
#'
#' @param graph a graph
#' @param xlab x-axis label
#' @param ylab y-axis label
#'
#' @returns a x and y labeled graph
#' @export

NSaxislabels <- function(graph, xlab, ylab) {

  graph <- graph +
    ggplot2::labs(
      x = xlab,
      y = ylab
    )
  return(graph)
}


#' Add title to graph
#'
#' @param graph a graph
#' @param title a title
#'
#' @returns a graph with a title
#' @export

NStitle <- function(graph, title) {

  graph <- graph +
    ggplot2::ggtitle(title)

  return(graph)
}



#' Add soft labels depicting the top performing (by euclidian distance) genes
#' for each square
#'
#' @param graph graph output of basegraph
#' @param data processed dataframe, needs to be the same fed into NSbasegraph
#' @param groups_labeled square groups where top n genes will be labeled.
#'   Default is c("top_center", "bottom_center", "middle_right", "middle_left")
#' @param top_labeled number of genes from the top by euclidian distance to be labeled
#'
#' @returns a ggplot2 graph with labeled points
#' @export

NSaddtoplabels <- function(graph,
                           data,
                           groups_labeled = c("top_center", "bottom_center", "middle_right", "middle_left"),
                           top_labeled = 10) {

  nudge_distance <-
    mean(
      range(data$control, na.rm = TRUE)[2] - range(data$control, na.rm = TRUE)[1],
      range(data$treatment, na.rm = TRUE)[2] - range(data$treatment, na.rm = TRUE)[1]
    ) / 12

  label_df <-
    data %>%
    dplyr::filter(.data$square %in% groups_labeled) %>%
    dplyr::filter(.data$rank <= top_labeled) %>%
    dplyr::mutate(
      # Determine if point is above or below diagonal
      above_diagonal = .data$treatment > .data$control,

      # Nudge perpendicular to y=x based on square
      nudge_x = dplyr::case_when(
        .data$square == "top_center" ~ -nudge_distance * 0.5,  # upper-left
        .data$square == "bottom_center" ~ nudge_distance * 0.5,  # lower-right
        .data$square == "middle_right" ~ nudge_distance * 0.5,  # lower-right
        .data$square == "middle_left" ~ -nudge_distance * 0.5,  # upper-left
        .data$square == "top_left" ~ -nudge_distance * 0.5,  # upper-left
        .data$square == "top_right" & above_diagonal ~ -nudge_distance * 0.5,  # upper-left
        .data$square == "top_right" & !above_diagonal ~ nudge_distance * 0.5,  # lower-right
        .data$square == "bottom_left" & above_diagonal ~ -nudge_distance * 0.5,  # upper-left
        .data$square == "bottom_left" & !above_diagonal ~ nudge_distance * 0.5,  # lower-right
        .data$square == "bottom_right" ~ nudge_distance * 0.5,  # lower-right
        TRUE ~ ifelse(above_diagonal, -nudge_distance * 0.5, nudge_distance * 0.5)
      ),

      nudge_y = dplyr::case_when(
        .data$square == "top_center" ~ nudge_distance * 0.5,  # upper-left
        .data$square == "bottom_center" ~ -nudge_distance * 0.5,  # lower-right
        .data$square == "middle_right" ~ -nudge_distance * 0.5,  # lower-right
        .data$square == "middle_left" ~ nudge_distance * 0.5,  # upper-left
        .data$square == "top_left" ~ nudge_distance * 0.5,  # upper-left
        .data$square == "top_right" & above_diagonal ~ nudge_distance * 0.5,  # upper-left
        .data$square == "top_right" & !above_diagonal ~ -nudge_distance * 0.5,  # lower-right
        .data$square == "bottom_left" & above_diagonal ~ nudge_distance * 0.5,  # upper-left
        .data$square == "bottom_left" & !above_diagonal ~ -nudge_distance * 0.5,  # lower-right
        .data$square == "bottom_right" ~ -nudge_distance * 0.5,  # lower-right
        TRUE ~ ifelse(above_diagonal, nudge_distance * 0.5, -nudge_distance * 0.5)
      )
    )

  graph <-
    graph +
    ggrepel::geom_text_repel(
      data = label_df,
      mapping = ggplot2::aes(
        x = .data$control,
        y = .data$treatment,
        label = .data$id,
        color = .data$square
      ),
      nudge_x = label_df$nudge_x,
      nudge_y = label_df$nudge_y,
      size = 3,
      inherit.aes = FALSE,
      max.overlaps = Inf,
      min.segment.length = 0.1,
      direction = "both"
    )

  return(graph)
}

# NSaddtoplabels <- function(graph,
#                            data,
#                            groups_labeled = c("top_center", "bottom_center", "middle_right", "middle_left"),
#                            top_labeled = 10) {
#   label_df <-
#     data %>%
#     dplyr::filter(.data$square %in% groups_labeled) %>%
#     dplyr::filter(.data$rank <= top_labeled)
#
#   nudge_distance <-
#     mean(
#         range(data$control, na.rm = TRUE)[2] - range(data$control, na.rm = TRUE)[1],
#         range(data$treatment, na.rm = TRUE)[2] - range(data$treatment, na.rm = TRUE)[1]
#       ) / 12
#
#   graph <-
#     graph +
#     ggrepel::geom_text_repel(
#       data = label_df,
#       mapping = ggplot2::aes(
#         x = .data$control,
#         y = .data$treatment,
#         label = .data$id,
#         color = .data$square
#       ),
#       size = 3,
#       inherit.aes = FALSE,
#       max.overlaps = Inf,
#       min.segment.length = 0.1,
#       position = ggpp::position_nudge_center(x = nudge_distance,
#                                              y = nudge_distance,
#                                              center_x = 0,
#                                              center_y = 0,
#                                              direction = "radial",
#                                              obey_grouping = FALSE)
#     )
#
#   return(graph)
# }





#' Add Genes of Interest (goi) highlights
#'
#' @param data processed dataframe, needs to be the same fed into NSbasegraph
#' @param graph graph output of basegraph
#' @param goi_list vector list of genes of interest
#' @param goi_auto logical, indicate default processing of points with slightly
#'   darker colors than the ones in the graph. Pretttyyyyy
#' @param goi_color highlighted point color, if goi_shape is between 21 and 28 it
#'   refers to the outline.
#' @param goi_fill if goi_shape is between 21 and 28, inside part of the point
#' @param goi_shape numerical (ggplot2 aesthetics), shape of the highlighted point
#' @param goi_size numerical, size of the highlighted point, default is 2
#' @param goi_label_color color of the label text, default is black
#' @param goi_label_size size of the label text, default is 4
#' @param goi_label_type type of label, either "text" or "label"
#'
#' @returns highlighted ggplot graph
#' @export

NSaddgoi <- function(data,
                     graph,
                     goi_list,
                     goi_auto = FALSE,
                     goi_shape = 21,
                     goi_color = "black",
                     goi_fill = "black",
                     goi_size = 2.5,
                     goi_label_type = c("text", "label"),
                     goi_label_color = "black",
                     goi_label_size = 4) {

  if(missing(goi_label_type)) {
    goi_label_type <- "label"
  }

  if(missing(data)) {
    data <- graph$data
  }

  tempdf <- data %>%
    dplyr::filter(.data$id %in% goi_list)

  nudge_distance <-
    mean(
      range(data$control, na.rm = TRUE)[2] - range(data$control, na.rm = TRUE)[1],
      range(data$treatment, na.rm = TRUE)[2] - range(data$treatment, na.rm = TRUE)[1]
    ) / 12

  if(goi_auto) {

    graph <-
      graph +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(
        data = tempdf,
        mapping = ggplot2::aes(x = .data$control,
                               y = .data$treatment,
                               fill = .data$square),
        shape = 21,
        size = goi_size,
        color = goi_color,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        aesthetics = c("fill"),
        guide = "none",
        values = c(
          "center" = "#B8B8B8",
          "neutral_slope" = "#B8B8B8",
          "bottom_left" = "#003366",
          "bottom_center" = "#005599",
          "bottom_right" = "#2299CC",
          "middle_left" = "#CC9922",
          "middle_right" = "#007755",
          "top_left" = "#DD6622",
          "top_center" = "#DD2266",
          "top_right" = "#992266"
        )
      )
  } else {

    graph <-
      graph +
      ggplot2::geom_point(
        data = tempdf,
        mapping = ggplot2::aes(x = .data$control,
                               y = .data$treatment),
        shape = goi_shape,
        size = goi_size,
        color = goi_color,
        fill = goi_fill,
        inherit.aes = FALSE
      )
  }

  if(goi_label_type == "text") {

  graph <-
    graph +
    ggrepel::geom_text_repel(
      data = tempdf,
      mapping = ggplot2::aes(x = "control",
                             y = "treatment",
                             label = "id"),
      color = goi_label_color,
      size = goi_label_size,
      inherit.aes = FALSE,
      max.overlaps = Inf,
      min.segment.length = 0.4,
      show.legend = FALSE,
      position = ggpp::position_nudge_center(x = nudge_distance,
                                             y = nudge_distance,
                                             center_x = 0,
                                             center_y = 0,
                                             direction = "radial",
                                             obey_grouping = FALSE)
    )

  } else if (goi_label_type == "label") {

    graph <-
      graph +
      ggrepel::geom_label_repel(
        data = tempdf,
        mapping = ggplot2::aes(x = .data$control,
                               y = .data$treatment,
                               label = .data$id),
        color = goi_label_color,
        size = goi_label_size,
        inherit.aes = FALSE,
        max.overlaps = Inf,
        min.segment.length = 0.4,
        label.padding = 0.2,
        show.legend = FALSE,
        position = ggpp::position_nudge_center(x = nudge_distance,
                                               y = nudge_distance,
                                               center_x = 0,
                                               center_y = 0,
                                               direction = "radial",
                                               obey_grouping = FALSE)
      )

  }


  return(graph)
}
