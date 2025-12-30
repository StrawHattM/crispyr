# NScutoff - Generate Cutoff values ---------------------------------------

#' Generate cutoffs for Nine Square Plots
#'
#' @param x numerical column / vector
#' @param scale standard deviation multiplier
#' @param force_zero_center indicate if you want the cutoffs to be forcefully centered around zero
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

  if(!is.numeric(x)) {
    stop("x has to be numeric")
  }
  mu <- mean(x, na.rm = TRUE)
  stdv <- stats::sd(x, na.rm = TRUE)

  if (force_zero_center == TRUE) {
    mu <- 0
  }

  return(c(mu - stdv * scale, mu + stdv * scale))
}



# NSbasedf ----------------------------------------------------------------

NSbasedf <-
  function(data, control, treament, ctrl_pval, treat_pval, min_sgrna) {

    data %>% dplyr::select(
      id = 1,
      num = 2,
      control = {{ control }},
      treatment = {{ treament }},
      ctrl_pval = {{ ctrl_pval }},
      treat_pval = {{ treat_pval }}
    ) %>%
      dplyr::filter(.data$num >= min_sgrna)

}
