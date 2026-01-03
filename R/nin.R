#' Not In operator
#'
#' Negation of the %in% operator for checking if elements are NOT in a vector.
#'
#' @param x Character vector to check belonging
#' @param y Character vector to check x against
#'
#' @returns Boolean vector of the same length as x
#' @export
`%nin%` <- function(x, y) {
  !(x %in% y)
}
