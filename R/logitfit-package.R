#' Logistic regression by iteratively reweighted least squares
#'
#' `logitfit` fits binomial logistic regression with the logit link. The
#' coefficient update is Fisher scoring (iteratively reweighted least
#' squares), implemented in R and in C++. It is not a wrapper around
#' [stats::glm()].
#'
#' @keywords internal
#' @import stats
#' @importFrom bench mark
#' @importFrom Rcpp evalCpp
#' @useDynLib logitfit, .registration = TRUE
"_PACKAGE"
