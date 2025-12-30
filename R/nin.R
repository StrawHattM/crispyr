#' /Not in/ operator
#'
#' @param x Character vector to check belonging.
#' @param y Character vector to check x against.
#'
#' @returns Boolean vector of the same length as x indicating if each element of x is not in y.
#' @export
#'
#' @examples
#' x <- c("A", "B", "C", "D", "E", "F")
#' y <- c("B", "A", "D", "E")
#'
#' x %nin% y
#' c("FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "TRUE")
#'

`%nin%` <- function(x, y) { return(!x %in% y) }
