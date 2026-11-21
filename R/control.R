#' Control parameters for logistic regression
#'
#' @param epsilon Positive convergence tolerance. Iteration stops when the
#'   absolute change in deviance, divided by `0.1 + abs(deviance)`, falls
#'   below `epsilon`. The default matches [stats::glm.control()].
#' @param maxit Maximum number of Fisher scoring iterations. Must be a
#'   positive integer.
#' @param trace If `TRUE`, print the deviance at each iteration.
#'
#' @return A list with components `epsilon`, `maxit`, and `trace`.
#'
#' @examples
#' logit_control(epsilon = 1e-10, maxit = 50)
#'
#' @export
#' @seealso [logit_fit()]
logit_control <- function(epsilon = 1e-8, maxit = 25, trace = FALSE) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("'epsilon' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) || maxit < 1) {
    stop("'maxit' must be a positive number", call. = FALSE)
  }
  if (abs(maxit - round(maxit)) > 1e-8) {
    stop("'maxit' must be a whole number", call. = FALSE)
  }
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE", call. = FALSE)
  }
  list(epsilon = epsilon, maxit = as.integer(round(maxit)), trace = trace)
}
