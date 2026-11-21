#' @export
print.logitfit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  coefs <- coef(x)
  if (length(coefs)) {
    cat("Coefficients:\n")
    print.default(format(coefs, digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat(
    "\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n",
    sep = " "
  )
  cat("Null Deviance:    ", format(x$null.deviance, digits = digits), "\n", sep = "")
  cat(
    "Residual Deviance:", format(x$deviance, digits = digits),
    "\tAIC:", format(x$aic, digits = digits), "\n",
    sep = " "
  )
  cat("Number of Fisher scoring iterations:", x$iter, "\n", sep = " ")
  cat("Engine:", x$engine, "\n", sep = " ")
  if (!isTRUE(x$converged)) {
    cat("Warning: algorithm did not converge\n")
  }
  invisible(x)
}

#' @export
summary.logitfit <- function(object, ...) {
  est <- coef(object)
  se <- sqrt(diag(vcov(object)))
  zval <- est / se
  pval <- 2 * pnorm(-abs(zval))
  coefmat <- cbind(est, se, zval, pval)
  dimnames(coefmat) <- list(
    names(est),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  )
  structure(
    list(
      call = object$call,
      coefficients = coefmat,
      cov.unscaled = vcov(object),
      deviance.resid = object$deviance_resid,
      deviance = object$deviance,
      null.deviance = object$null.deviance,
      df.residual = object$df.residual,
      df.null = object$df.null,
      aic = object$aic,
      iter = object$iter,
      converged = object$converged,
      boundary = object$boundary,
      engine = object$engine,
      aliased = is.na(est)
    ),
    class = "summary.logitfit"
  )
}

#' @export
print.summary.logitfit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Deviance Residuals:\n")
  resid_q <- quantile(x$deviance.resid, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  names(resid_q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resid_q, digits = digits)

  if (any(x$aliased)) {
    cat(
      "\nCoefficients: (", sum(x$aliased),
      " not defined because of singularities)\n",
      sep = ""
    )
  } else {
    cat("\nCoefficients:\n")
  }
  printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, na.print = "NA", ...)
  cat("\n(Dispersion parameter for binomial family taken to be 1)\n\n")
  cat(
    "    Null deviance: ", format(x$null.deviance, digits = digits),
    " on ", format(x$df.null), " degrees of freedom\n",
    sep = ""
  )
  cat(
    "Residual deviance: ", format(x$deviance, digits = digits),
    " on ", format(x$df.residual), " degrees of freedom\n",
    sep = ""
  )
  cat("AIC: ", format(x$aic, digits = digits), "\n\n", sep = "")
  cat("Number of Fisher scoring iterations: ", x$iter, "\n", sep = "")
  cat("Engine: ", x$engine, "\n", sep = "")
  if (!isTRUE(x$converged)) {
    cat("Warning: algorithm did not converge\n")
  }
  invisible(x)
}

#' @export
coef.logitfit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.logitfit <- function(object, ...) {
  object$vcov
}

#' @export
fitted.logitfit <- function(object, ...) {
  object$fitted.values
}

#' @export
deviance.logitfit <- function(object, ...) {
  object$deviance
}

#' @export
nobs.logitfit <- function(object, ...) {
  object$nobs
}

#' @export
logLik.logitfit <- function(object, ...) {
  val <- object$loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- object$rank
  class(val) <- "logLik"
  val
}

#' @export
formula.logitfit <- function(x, ...) {
  x$formula
}

#' @export
family.logitfit <- function(object, ...) {
  binomial()
}

#' @export
weights.logitfit <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  if (type == "prior") object$prior.weights else object$weights
}

#' @export
model.frame.logitfit <- function(formula, ...) {
  if (is.null(formula$model)) {
    stop("the model frame was not stored; refit with model = TRUE", call. = FALSE)
  }
  formula$model
}

#' @export
model.matrix.logitfit <- function(object, ...) {
  if (!is.null(object$x)) {
    return(object$x)
  }
  if (is.null(object$model)) {
    stop("refit with model = TRUE or x = TRUE to recover the model matrix", call. = FALSE)
  }
  model.matrix(object$terms, object$model, contrasts.arg = object$contrasts)
}

#' Residuals from a logistic regression
#'
#' @param object A `logitfit` object.
#' @param type Residual definition. `"deviance"` and `"pearson"` are the
#'   GLM residuals. `"working"` is the IRLS residual on the link scale.
#'   `"response"` is `y - fitted`.
#' @param ... Unused.
#'
#' @return A numeric vector with one residual per observation used in the fit.
#' @export
residuals.logitfit <- function(object,
                               type = c("deviance", "pearson", "working", "response"),
                               ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  if (is.null(y)) {
    stop("the response was not stored; refit with y = TRUE", call. = FALSE)
  }
  wt <- object$prior.weights
  res <- switch(type,
    deviance = sign(y - mu) * sqrt(pmax(0, binomial_dev_resids(y, mu, wt))),
    pearson = (y - mu) * sqrt(wt) / sqrt(pmax(mu * (1 - mu), 0)),
    working = (y - mu) / (mu * (1 - mu)),
    response = y - mu
  )
  names(res) <- names(mu)
  res
}

#' Predictions from a logistic regression
#'
#' @param object A `logitfit` object.
#' @param newdata Optional data frame of new predictor values. If omitted,
#'   fitted values for the estimation sample are returned.
#' @param type `"link"` returns log-odds. `"response"` returns probabilities.
#' @param se.fit If `TRUE`, also return standard errors. On the response
#'   scale the standard error uses the delta method.
#' @param na.action Missing-value handler used when `newdata` is supplied.
#' @param ... Unused.
#'
#' @return A numeric vector, or a list with components `fit`, `se.fit`, and
#'   `residual.scale` when `se.fit = TRUE`.
#' @export
predict.logitfit <- function(object, newdata = NULL,
                             type = c("link", "response"),
                             se.fit = FALSE, na.action = na.pass, ...) {
  type <- match.arg(type)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be TRUE or FALSE", call. = FALSE)
  }

  if (is.null(newdata)) {
    eta <- object$linear.predictors
    se <- if (se.fit) eta_se(model_matrix_from_fit(object), object) else NULL
    eta <- restore_na(object, eta)
    if (se.fit) {
      se <- restore_na(object, se)
    }
  } else {
    Terms <- delete.response(terms(object))
    mf <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
    X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
    off <- model.offset(mf)
    if (is.null(off)) {
      off <- rep(0, nrow(X))
    }
    beta <- coef(object)
    beta[is.na(beta)] <- 0
    eta <- as.vector(X %*% beta) + off
    names(eta) <- rownames(X)
    se <- if (se.fit) eta_se(X, object) else NULL
  }

  pred <- if (type == "link") eta else logit_linkinv(eta)
  if (se.fit && type == "response") {
    mu <- logit_linkinv(eta)
    se <- se * mu * (1 - mu)
  }
  if (!se.fit) {
    return(pred)
  }
  list(fit = pred, se.fit = se, residual.scale = 1)
}

#' Wald confidence intervals for logistic coefficients
#'
#' @param object A `logitfit` object.
#' @param parm Coefficients to include. A character vector of names or a
#'   numeric index. Defaults to every coefficient.
#' @param level Confidence level between 0 and 1.
#' @param ... Unused.
#'
#' @return A matrix with one row per coefficient and two columns of limits.
#' @export
confint.logitfit <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    if (any(parm < 1 | parm > length(pnames))) {
      stop("'parm' subscripts are out of range", call. = FALSE)
    }
    parm <- pnames[parm]
  }
  if (anyNA(match(parm, pnames))) {
    stop("unknown coefficient names in 'parm'", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1", call. = FALSE)
  }
  tail <- (1 - level) / 2
  probs <- c(tail, 1 - tail)
  fac <- qnorm(probs)
  ci <- cbind(cf[parm] + se[parm] * fac[1L], cf[parm] + se[parm] * fac[2L])
  colnames(ci) <- paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = 3), "%")
  rownames(ci) <- parm
  ci
}

#' @noRd
model_matrix_from_fit <- function(object) {
  if (!is.null(object$x)) {
    return(object$x)
  }
  if (!is.null(object$model)) {
    return(model.matrix(object$terms, object$model, contrasts.arg = object$contrasts))
  }
  stop(
    "standard errors need the model frame or the model matrix; refit with model = TRUE or x = TRUE",
    call. = FALSE
  )
}

#' @noRd
eta_se <- function(X, object) {
  beta <- coef(object)
  ok <- !is.na(beta) & is.finite(beta)
  if (!any(ok)) {
    return(rep(NA_real_, nrow(X)))
  }
  V <- vcov(object)[ok, ok, drop = FALSE]
  if (anyNA(V)) {
    return(rep(NA_real_, nrow(X)))
  }
  Xs <- X[, ok, drop = FALSE]
  se2 <- rowSums((Xs %*% V) * Xs)
  se2[se2 < 0] <- NA_real_
  sqrt(se2)
}

#' @noRd
restore_na <- function(object, values) {
  if (is.null(object$na.action)) {
    return(values)
  }
  napredict(object$na.action, values)
}
