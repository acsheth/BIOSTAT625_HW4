#' Compare logitfit with stats::glm
#'
#' Fits the same logistic regression with [stats::glm()], with the pure R
#' engine, and with the C++ engine. Agreement is checked with [all.equal()].
#' When `n_bench` is positive, [bench::mark()] times the three fits.
#'
#' The C++ engine is not required to be faster than `glm()`. `glm()` is
#' itself implemented in C. The timing is there to show the cost of the
#' implementation in this package.
#'
#' @param formula A two-sided model formula.
#' @param data A data frame. Required so the benchmark refits the same sample.
#' @param ... Optional arguments passed to both fitters: `weights`, `subset`,
#'   `na.action`, `offset`, `start`, and `contrasts`. `control`,
#'   `singular.ok`, `model`, `x`, and `y` are passed only to [logit_fit()].
#' @param tolerance Positive tolerance forwarded to [all.equal()].
#' @param n_bench Minimum number of timing iterations. Use `0` to skip timing.
#'
#' @return An object of class `logitfit_comparison` with components
#'   `agreement` (a data frame of `all.equal` results), `max_coef_diff`,
#'   `benchmark`, `glm`, `logitfit_r`, and `logitfit_rcpp`.
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' x <- rnorm(n)
#' yy <- rbinom(n, 1, stats::plogis(-0.2 + 0.7 * x))
#' cmp <- compare_with_glm(yy ~ x, data = data.frame(x = x, yy = yy), n_bench = 1)
#' cmp$agreement
#'
#' @export
#' @seealso [logit_fit()]
compare_with_glm <- function(formula, data, ..., tolerance = 1e-6, n_bench = 10) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula", call. = FALSE)
  }
  if (missing(data)) {
    stop("'data' is required", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("'tolerance' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(n_bench) || length(n_bench) != 1L || !is.finite(n_bench) || n_bench < 0) {
    stop("'n_bench' must be a non-negative number", call. = FALSE)
  }
  n_bench <- as.integer(round(n_bench))

  dots <- list(...)
  if (length(dots) > 0L) {
    if (is.null(names(dots)) || any(names(dots) == "")) {
      stop("arguments in ... must be named", call. = FALSE)
    }
  }
  if ("engine" %in% names(dots)) {
    warning("'engine' is ignored; compare_with_glm fits both engines", call. = FALSE)
    dots$engine <- NULL
  }
  allowed <- c(
    "weights", "subset", "na.action", "offset", "start", "contrasts",
    "control", "singular.ok", "model", "x", "y"
  )
  unknown <- setdiff(names(dots), allowed)
  if (length(unknown) > 0L) {
    stop("unknown arguments: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  glm_dots <- dots[intersect(names(dots), c("weights", "subset", "na.action", "offset", "start", "contrasts"))]

  fit_glm <- do.call(glm, c(list(formula = formula, data = data, family = binomial()), glm_dots))
  fit_r <- do.call(logit_fit, c(list(formula = formula, data = data, engine = "r"), dots))
  fit_cpp <- do.call(logit_fit, c(list(formula = formula, data = data, engine = "rcpp"), dots))

  agreement <- rbind(
    agreement_rows("r", fit_r, fit_glm, tolerance),
    agreement_rows("rcpp", fit_cpp, fit_glm, tolerance)
  )
  rownames(agreement) <- NULL

  benchmark <- NULL
  if (n_bench > 0L) {
    run_glm <- function() {
      do.call(glm, c(list(formula = formula, data = data, family = binomial()), glm_dots))
    }
    run_r <- function() {
      do.call(logit_fit, c(list(formula = formula, data = data, engine = "r"), dots))
    }
    run_cpp <- function() {
      do.call(logit_fit, c(list(formula = formula, data = data, engine = "rcpp"), dots))
    }
    benchmark <- mark(
      glm = run_glm(),
      logitfit_r = run_r(),
      logitfit_rcpp = run_cpp(),
      min_iterations = n_bench,
      check = FALSE
    )
  }

  structure(
    list(
      agreement = agreement,
      max_coef_diff = c(
        r = max_abs_diff(coef(fit_r), coef(fit_glm)),
        rcpp = max_abs_diff(coef(fit_cpp), coef(fit_glm))
      ),
      benchmark = benchmark,
      glm = fit_glm,
      logitfit_r = fit_r,
      logitfit_rcpp = fit_cpp,
      tolerance = tolerance,
      call = match.call()
    ),
    class = "logitfit_comparison"
  )
}

#' @export
print.logitfit_comparison <- function(x, ...) {
  cat("Correctness relative to stats::glm (binomial logit)\n")
  cat("all.equal tolerance:", format(x$tolerance), "\n\n")
  print(x$agreement, row.names = FALSE, digits = 4)
  cat(
    "\nMaximum absolute coefficient difference (R):   ",
    format(x$max_coef_diff[["r"]], digits = 4, scientific = TRUE), "\n",
    sep = ""
  )
  cat(
    "Maximum absolute coefficient difference (Rcpp): ",
    format(x$max_coef_diff[["rcpp"]], digits = 4, scientific = TRUE), "\n",
    sep = ""
  )
  if (!is.null(x$benchmark)) {
    cat("\nTiming (bench::mark)\n\n")
    print(x$benchmark)
  }
  invisible(x)
}

#' @noRd
agreement_rows <- function(engine, fit, fit_glm, tolerance) {
  pairs <- list(
    coefficients = list(coef(fit), coef(fit_glm)),
    fitted = list(as.numeric(fitted(fit)), as.numeric(fitted(fit_glm))),
    linear_predictors = list(as.numeric(fit$linear.predictors), as.numeric(fit_glm$linear.predictors)),
    deviance = list(deviance(fit), deviance(fit_glm)),
    null_deviance = list(fit$null.deviance, fit_glm$null.deviance),
    aic = list(AIC(fit), AIC(fit_glm))
  )
  rows <- lapply(names(pairs), function(name) {
    a <- pairs[[name]][[1L]]
    b <- pairs[[name]][[2L]]
    check <- all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)
    data.frame(
      engine = engine,
      quantity = name,
      all_equal = isTRUE(check),
      max_abs_diff = max_abs_diff(a, b),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @noRd
max_abs_diff <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b)) {
    return(NA_real_)
  }
  both_na <- is.na(a) & is.na(b)
  a[both_na] <- 0
  b[both_na] <- 0
  if (anyNA(a) || anyNA(b)) {
    return(NA_real_)
  }
  if (length(a) == 0L) {
    return(0)
  }
  max(abs(a - b))
}
