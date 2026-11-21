# Turn the model-frame response into the (proportion, weight, trials) form
# used by binomial Fisher scoring. A two-column matrix is successes and failures.

#' @param Y Response from the model frame.
#' @param prior_weights Numeric weights, or `NULL` for equal weights.
#' @return A list with `y`, `weights`, `n_vec`, and `warnings`.
#' @noRd
prepare_binomial_response <- function(Y, prior_weights) {
  warns <- character()
  if (is.logical(Y)) {
    Y <- as.numeric(Y)
  }
  if (is.matrix(Y) && ncol(Y) == 1L) {
    Y <- as.numeric(Y)
  }

  nobs <- NROW(Y)
  if (is.null(prior_weights)) {
    prior_weights <- rep(1, nobs)
  }
  prior_weights <- as.numeric(prior_weights)
  if (length(prior_weights) != nobs) {
    stop("weights and the response have different lengths", call. = FALSE)
  }
  if (anyNA(prior_weights) || any(!is.finite(prior_weights))) {
    stop("weights must be finite", call. = FALSE)
  }
  if (any(prior_weights < 0)) {
    stop("negative weights are not allowed", call. = FALSE)
  }

  if (is.factor(Y)) {
    if (nlevels(Y) != 2L) {
      stop("a factor response must have exactly two levels", call. = FALSE)
    }
    # The first level is failure, matching stats::binomial().
    Y <- as.numeric(Y != levels(Y)[1L])
  }

  if (is.matrix(Y)) {
    if (ncol(Y) != 2L) {
      stop(
        "the response must be a binary vector or a two-column matrix of successes and failures",
        call. = FALSE
      )
    }
    if (any(Y < 0, na.rm = TRUE)) {
      stop("negative counts in a binomial response", call. = FALSE)
    }
    if (any(abs(Y - round(Y)) > 0.001, na.rm = TRUE)) {
      warns <- c(warns, "non-integer counts in a binomial logitfit model")
    }
    trials <- Y[, 1L] + Y[, 2L]
    y <- ifelse(trials == 0, 0, Y[, 1L] / trials)
    weights <- prior_weights * trials
    n_vec <- trials
  } else {
    if (!is.numeric(Y)) {
      stop(
        "the response must be numeric, logical, a two-level factor, or a two-column matrix",
        call. = FALSE
      )
    }
    y <- as.numeric(Y)
    if (any(y < 0 | y > 1, na.rm = TRUE)) {
      stop("y values must be 0 <= y <= 1", call. = FALSE)
    }
    scored <- prior_weights * y
    if (any(abs(scored - round(scored)) > 0.001, na.rm = TRUE)) {
      warns <- c(warns, "non-integer #successes in a binomial logitfit model")
    }
    weights <- prior_weights
    n_vec <- rep(1, nobs)
  }

  if (anyNA(y)) {
    stop("the response contains missing values; use na.action = na.omit", call. = FALSE)
  }

  list(y = y, weights = weights, n_vec = n_vec, warnings = warns)
}

#' Factor levels among the predictors, for `predict()` on new data.
#'
#' @noRd
predictor_levels <- function(mt, mf) {
  fac <- attr(mt, "factors")
  if (is.null(fac) || length(fac) == 0L) {
    return(NULL)
  }
  vars <- intersect(rownames(fac), names(mf))
  levels_list <- lapply(mf[vars], function(col) {
    if (is.factor(col)) levels(col) else NULL
  })
  levels_list <- levels_list[!vapply(levels_list, is.null, logical(1))]
  if (length(levels_list) == 0L) NULL else levels_list
}

#' Align a user starting vector to the model-matrix columns.
#'
#' @noRd
align_start <- function(start, xnames) {
  if (!is.numeric(start)) {
    stop("'start' must be numeric", call. = FALSE)
  }
  if (is.null(names(start))) {
    if (length(start) != length(xnames)) {
      stop(
        "'start' has length ", length(start), " but the model has ",
        length(xnames), " coefficients",
        call. = FALSE
      )
    }
    beta <- as.numeric(start)
  } else {
    missing <- setdiff(xnames, names(start))
    if (length(missing) > 0L) {
      stop(
        "'start' is missing coefficients: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    beta <- as.numeric(start[xnames])
  }
  if (any(!is.finite(beta))) {
    stop("'start' must be finite", call. = FALSE)
  }
  names(beta) <- xnames
  beta
}

#' Starting linear predictor. With no `start`, this is the binomial default
#' `(weight * y + 0.5) / (weight + 1)` mapped through the logit.
#'
#' @noRd
initial_eta <- function(X, y, weights, offset, start) {
  xnames <- colnames(X)
  if (!is.null(start)) {
    beta <- align_start(start, xnames)
    eta <- as.vector(X %*% beta) + offset
    return(list(eta = eta, beta = beta, has_start = TRUE))
  }
  mustart <- (weights * y + 0.5) / (weights + 1)
  list(
    eta = log(mustart / (1 - mustart)),
    beta = setNames(rep(0, ncol(X)), xnames),
    has_start = FALSE
  )
}
