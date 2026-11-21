# Binomial logit helpers shared by the fitting routine, inference, and residuals.
# The deviance and log-likelihood follow the binomial exponential family used
# by stats::binomial(), so fitted objects can be checked against stats::glm().

#' Stable inverse logit.
#'
#' Large positive and negative linear predictors are handled separately so the
#' probabilities stay inside (0, 1) and do not overflow.
#'
#' @param eta Numeric linear predictor.
#' @return Probabilities on the response scale.
#' @noRd
logit_linkinv <- function(eta) {
  nms <- names(eta)
  eta <- as.numeric(eta)
  out <- rep(NA_real_, length(eta))
  finite <- is.finite(eta)
  pos <- finite & eta >= 0
  neg <- finite & eta < 0
  out[pos] <- 1 / (1 + exp(-eta[pos]))
  out[neg] <- exp(eta[neg]) / (1 + exp(eta[neg]))
  out[is.infinite(eta) & eta > 0] <- 1
  out[is.infinite(eta) & eta < 0] <- 0
  names(out) <- nms
  out
}

#' One observation's contribution to the binomial deviance.
#'
#' Uses the convention 0 * log(0) = 0. `mu` is the mean and `wt` is the
#' binomial prior weight (trial count times any extra weight).
#'
#' @param y,mu,wt Numeric vectors of equal length.
#' @return A vector of deviance contributions, already multiplied by 2 * wt.
#' @noRd
binomial_dev_resids <- function(y, mu, wt) {
  ylogy <- function(yy, mm) {
    out <- numeric(length(yy))
    pos <- which(yy > 0)
    out[pos] <- yy[pos] * log(yy[pos] / mm[pos])
    out
  }
  2 * wt * (ylogy(y, mu) + ylogy(1 - y, 1 - mu))
}

#' Binomial log-likelihood used for AIC.
#'
#' Matches `stats::binomial()$aic` for a binary vector and for a two-column
#' matrix of successes and failures. `n_vec` is the number of trials in each
#' row before prior weights are applied. `wt` is the weight after multiplying
#' by the number of trials.
#'
#' @param y Success proportion.
#' @param mu Fitted probability.
#' @param wt Prior weights used in the scoring iteration.
#' @param n_vec Trial counts, or 1 for a binary response.
#' @return The log-likelihood as a single number.
#' @noRd
binomial_loglik <- function(y, mu, wt, n_vec) {
  m <- if (any(n_vec > 1)) n_vec else wt
  scale <- ifelse(m > 0, wt / m, 0)
  mu <- pmin(pmax(mu, 0), 1)
  ll <- dbinom(round(m * y), size = round(m), prob = mu, log = TRUE)
  sum(scale * ll)
}

#' Expected information matrix at the fitted probabilities.
#'
#' Columns with `NA` coefficients are treated as aliased and left as `NA` in
#' the covariance matrix. The active columns are scored with a weighted QR,
#' which is the same information matrix Fisher scoring uses at convergence.
#'
#' @param X Model matrix.
#' @param mu Fitted probabilities.
#' @param weights Prior weights.
#' @param beta Coefficient vector, with `NA` for aliased terms.
#' @return A p by p covariance matrix.
#' @noRd
vcov_at_fit <- function(X, mu, weights, beta) {
  p <- ncol(X)
  xnames <- colnames(X)
  V <- matrix(NA_real_, p, p, dimnames = list(xnames, xnames))
  ok <- !is.na(beta)
  if (!any(ok)) {
    return(V)
  }
  mu_eta <- mu * (1 - mu)
  good <- which(weights > 0 & mu_eta > 0 & is.finite(mu))
  if (length(good) == 0L) {
    return(V)
  }
  w <- sqrt(weights[good] * mu_eta[good])
  Xw <- X[good, ok, drop = FALSE] * w
  qrX <- qr(Xw, tol = 1e-10)
  rank <- qrX$rank
  if (rank < sum(ok)) {
    return(V)
  }
  Rmat <- qr.R(qrX)
  R1 <- Rmat[seq_len(rank), seq_len(rank), drop = FALSE]
  cov_pivot <- chol2inv(R1)
  idx <- qrX$pivot[seq_len(rank)]
  V_ok <- matrix(0, sum(ok), sum(ok))
  V_ok[idx, idx] <- cov_pivot
  V[ok, ok] <- V_ok
  V
}

#' Null mean for the intercept-only or offset-only binomial model.
#'
#' @noRd
null_mean <- function(y, weights, offset, intercept) {
  if (intercept) {
    total <- sum(weights)
    if (total == 0) {
      return(rep(NA_real_, length(y)))
    }
    return(rep(sum(weights * y) / total, length(y)))
  }
  logit_linkinv(offset)
}
