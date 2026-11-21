# Pure R Fisher scoring for the binomial logit model.
#
# At iteration t the working response and weight are
#   z_i = eta_i - offset_i + (y_i - mu_i) / (mu_i (1 - mu_i))
#   w_i = sqrt(prior_i * mu_i * (1 - mu_i))
# and the next coefficient is the weighted least-squares estimate
#   beta = argmin || w * (z - X beta) ||.
# For the canonical logit link this is Newton's method on the log-likelihood.
# The QR tolerance and the deviance convergence test follow stats::glm.fit().

#' @param X Numeric model matrix, including a column of ones when there is an intercept.
#' @param y Success proportions in \[0, 1\].
#' @param weights Prior weights, already multiplied by the number of trials.
#' @param offset Numeric offset of length `nrow(X)`.
#' @param start `NULL` or a numeric vector of length `ncol(X)`.
#' @param control List from [logit_control()].
#' @param singular.ok If `FALSE`, a rank-deficient model is an error.
#' @return A list of coefficients, fitted values, and iteration flags.
#' @noRd
logit_irls_r <- function(X, y, weights, offset, start, control, singular.ok) {
  nobs <- nrow(X)
  nvars <- ncol(X)
  xnames <- colnames(X)
  ynames <- names(y)
  msgs <- character()

  finish <- function(beta, eta, mu, dev, iter, conv, boundary) {
    names(beta) <- xnames
    names(mu) <- ynames
    names(eta) <- ynames
    list(
      coefficients = beta,
      fitted.values = mu,
      linear.predictors = eta,
      deviance = dev,
      iter = as.integer(iter),
      converged = conv,
      boundary = boundary,
      warnings = msgs,
      refit = FALSE
    )
  }

  if (nvars == 0L) {
    eta <- offset
    mu <- logit_linkinv(eta)
    dev <- sum(binomial_dev_resids(y, mu, weights))
    return(finish(
      beta = setNames(numeric(0), character(0)),
      eta = eta, mu = mu, dev = dev, iter = 0L, conv = TRUE, boundary = FALSE
    ))
  }

  tol <- min(1e-7, control$epsilon / 1000)
  conv <- FALSE
  boundary <- FALSE
  beta <- setNames(rep(NA_real_, nvars), xnames)
  coefold <- NULL
  aliased <- rep(FALSE, nvars)

  if (is.null(start)) {
    mustart <- (weights * y + 0.5) / (weights + 1)
    eta <- log(mustart / (1 - mustart))
  } else {
    beta <- align_start(start, xnames)
    coefold <- beta
    eta <- as.vector(X %*% beta) + offset
  }
  mu <- logit_linkinv(eta)
  dev <- sum(binomial_dev_resids(y, mu, weights))
  devold <- dev
  iter <- 0L

  for (iter in seq_len(control$maxit)) {
    mu_eta <- mu * (1 - mu)
    good <- (weights > 0) & (mu_eta != 0)
    if (!any(good)) {
      conv <- FALSE
      msgs <- c(msgs, paste0("no observations informative at iteration ", iter))
      break
    }

    # Working response on the coefficient scale, and square-root weights.
    z <- (eta - offset)[good] + (y - mu)[good] / mu_eta[good]
    w <- sqrt(weights[good] * mu_eta[good])
    Xw <- X[good, , drop = FALSE] * w
    qrX <- qr(Xw, tol = tol)
    if (!isTRUE(singular.ok) && qrX$rank < nvars) {
      stop("singular fit encountered", call. = FALSE)
    }
    est <- qr.coef(qrX, z * w)
    if (any(!is.finite(est[!is.na(est)]))) {
      conv <- FALSE
      msgs <- c(msgs, paste0("non-finite coefficients at iteration ", iter))
      break
    }

    aliased_new <- is.na(est)
    beta_new <- numeric(nvars)
    beta_new[!aliased_new] <- est[!aliased_new]
    eta_new <- as.vector(X %*% beta_new) + offset
    mu_new <- logit_linkinv(eta_new)
    dev_new <- sum(binomial_dev_resids(y, mu_new, weights))
    if (isTRUE(control$trace)) {
      cat("Deviance = ", dev_new, " Iterations - ", iter, "\n", sep = "")
    }

    boundary <- FALSE
    if (!is.finite(dev_new)) {
      # The full Newton step left the region where the deviance is finite.
      # Halve the step back toward the last finite coefficient, as glm.fit() does.
      if (is.null(coefold)) {
        stop(
          "no valid set of coefficients has been found: please supply starting values",
          call. = FALSE
        )
      }
      msgs <- c(msgs, "step size truncated due to divergence")
      ii <- 1L
      old <- coefold
      old[is.na(old)] <- 0
      while (!is.finite(dev_new)) {
        if (ii > control$maxit) {
          stop("inner loop 1; cannot correct step size", call. = FALSE)
        }
        ii <- ii + 1L
        beta_new <- (beta_new + old) / 2
        beta_new[aliased_new] <- 0
        eta_new <- as.vector(X %*% beta_new) + offset
        mu_new <- logit_linkinv(eta_new)
        dev_new <- sum(binomial_dev_resids(y, mu_new, weights))
      }
      boundary <- TRUE
      if (isTRUE(control$trace)) {
        cat("Step halved: new deviance = ", dev_new, "\n", sep = "")
      }
    }

    aliased <- aliased_new
    beta <- beta_new
    eta <- eta_new
    mu <- mu_new
    dev <- dev_new
    if (abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
      conv <- TRUE
      break
    }
    devold <- dev
    coefold <- beta
    coefold[aliased] <- 0
  }

  beta[aliased] <- NA_real_
  names(beta) <- xnames
  finish(beta, eta, mu, dev, iter, conv, boundary)
}

#' Whether the compiled scoring routine can be called.
#'
#' `logitfit.force_no_rcpp = TRUE` is used by the unit tests to exercise the
#' pure R fallback while the DLL is still loaded.
#'
#' @noRd
rcpp_engine_available <- function() {
  if (isTRUE(getOption("logitfit.force_no_rcpp"))) {
    return(FALSE)
  }
  fn <- get0("logit_irls_cpp", envir = asNamespace("logitfit"), inherits = FALSE)
  is.function(fn)
}

#' Choose `"r"` or `"rcpp"` from the user request and what is compiled.
#'
#' @noRd
resolve_engine <- function(engine) {
  engine <- match.arg(engine, c("auto", "rcpp", "r"))
  if (identical(engine, "r")) {
    return("r")
  }
  if (!rcpp_engine_available()) {
    if (identical(engine, "rcpp")) {
      stop(
        "The Rcpp engine is not available. Reinstall logitfit with a C++ toolchain, or set engine = \"r\".",
        call. = FALSE
      )
    }
    warning("The Rcpp engine is not available; using the pure R engine.", call. = FALSE)
    return("r")
  }
  "rcpp"
}
