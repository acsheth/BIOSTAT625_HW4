#' Fit a logistic regression by Fisher scoring
#'
#' Fits a binomial generalized linear model with the logit link. The
#' coefficient update is iteratively reweighted least squares, written out
#' in this package rather than delegated to [stats::glm()]. `engine = "r"`
#' solves each weighted least-squares step with a QR decomposition.
#' `engine = "rcpp"` builds the same normal equations in C++ and factors
#' them with a Cholesky decomposition. `engine = "auto"` uses C++ when the
#' package DLL is loaded.
#'
#' The response may be a binary 0/1 vector, a logical vector, a factor with
#' two levels (the first level is failure), or a two-column matrix of
#' successes and failures. Prior `weights` and an `offset` are allowed, as
#' in [stats::glm()].
#'
#' Wald standard errors come from the expected information matrix
#' \eqn{X'WX}{X'WX} at the fitted probabilities. The dispersion is fixed at 1.
#'
#' @param formula A two-sided model formula.
#' @param data Optional data frame (or list or environment) containing the
#'   variables in `formula`.
#' @param weights Optional numeric case weights. For a binary response these
#'   are prior weights. They must be non-negative and finite.
#' @param subset Optional expression selecting rows of `data`.
#' @param na.action Missing-value handler passed to [stats::model.frame()].
#'   The default is the global `na.action` option, usually `na.omit`.
#' @param start Optional numeric starting coefficients, of length equal to
#'   the number of columns in the model matrix. Names, if present, are
#'   matched to those columns.
#' @param offset Optional numeric offset. An `offset()` term in `formula`
#'   is also accepted.
#' @param control List produced by [logit_control()].
#' @param engine `"auto"`, `"rcpp"`, or `"r"`. See the description.
#' @param singular.ok If `FALSE`, a rank-deficient model matrix is an error.
#'   If `TRUE`, aliased coefficients are `NA`.
#' @param model If `TRUE`, store the model frame.
#' @param x If `TRUE`, store the model matrix.
#' @param y If `TRUE`, store the response proportions used in the fit.
#' @param contrasts Optional list of contrasts for factor terms.
#'
#' @return An object of class `logitfit`. Components include
#'   `coefficients`, `vcov`, `fitted.values`, `linear.predictors`,
#'   `deviance`, `null.deviance`, `aic`, `loglik`, `iter`, `converged`,
#'   `rank`, `df.residual`, `df.null`, `prior.weights`, `y`, `terms`,
#'   `formula`, `call`, and `engine`.
#'
#' @examples
#' fit <- logit_fit(died ~ age + severity + treatment + elective, data = icu_sim)
#' summary(fit)
#' predict(fit,
#'   newdata = data.frame(age = 70, severity = 3, treatment = 1, elective = 0),
#'   type = "response", se.fit = TRUE
#' )
#' fit_glm <- stats::glm(
#'   died ~ age + severity + treatment + elective,
#'   data = icu_sim, family = stats::binomial()
#' )
#' all.equal(stats::coef(fit), stats::coef(fit_glm), tolerance = 1e-6)
#'
#' @export
#' @seealso [logit_control()], [compare_with_glm()]
logit_fit <- function(formula, data, weights, subset, na.action, start = NULL,
                      offset, control = logit_control(),
                      engine = c("auto", "rcpp", "r"),
                      singular.ok = TRUE, model = TRUE, x = FALSE, y = TRUE,
                      contrasts = NULL) {
  engine <- match.arg(engine)
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula", call. = FALSE)
  }
  if (!is.logical(singular.ok) || length(singular.ok) != 1L || is.na(singular.ok)) {
    stop("'singular.ok' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'x' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(y) || length(y) != 1L || is.na(y)) {
    stop("'y' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.list(control) || is.null(control$epsilon) || is.null(control$maxit)) {
    stop("'control' must be a list from logit_control()", call. = FALSE)
  }
  if (is.null(control$trace)) {
    control$trace <- FALSE
  }

  keep_model <- model
  keep_x <- x
  keep_y <- y
  call <- match.call()

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (is.null(Y)) {
    stop("the formula needs a response", call. = FALSE)
  }
  prior <- model.weights(mf)
  prepared <- prepare_binomial_response(Y, prior)
  for (msg in prepared$warnings) {
    warning(msg, call. = FALSE)
  }
  y_vec <- prepared$y
  weights_vec <- prepared$weights
  n_vec <- prepared$n_vec
  names(y_vec) <- rownames(mf)

  offset_vec <- model.offset(mf)
  has_offset <- !is.null(offset_vec)
  if (!has_offset) {
    offset_vec <- rep(0, nrow(mf))
  }
  offset_vec <- as.numeric(offset_vec)
  if (length(offset_vec) == 1L) {
    offset_vec <- rep(offset_vec, nrow(mf))
  }
  if (length(offset_vec) != nrow(mf)) {
    stop("offset and the number of observations differ", call. = FALSE)
  }
  if (anyNA(offset_vec) || any(!is.finite(offset_vec))) {
    stop("offset must be finite", call. = FALSE)
  }
  names(offset_vec) <- rownames(mf)

  X <- model.matrix(mt, mf, contrasts.arg = contrasts)
  if (anyNA(X) || any(!is.finite(X))) {
    stop("the model matrix contains missing or non-finite values", call. = FALSE)
  }
  if (nrow(X) == 0L) {
    stop("no observations", call. = FALSE)
  }
  contrasts_used <- attr(X, "contrasts")
  xlevels <- predictor_levels(mt, mf)
  intercept <- isTRUE(attr(mt, "intercept") == 1)

  chosen <- resolve_engine(engine)
  fit <- run_engine(
    chosen = chosen,
    X = X,
    y = y_vec,
    weights = weights_vec,
    offset = offset_vec,
    start = start,
    control = control,
    singular.ok = singular.ok
  )

  beta <- fit$coefficients
  names(beta) <- colnames(X)
  # Rebuild the linear predictor from the reported coefficients so eta, mu,
  # deviance, and the information matrix all refer to the same point.
  if (ncol(X) == 0L) {
    eta <- offset_vec
  } else {
    beta_fill <- beta
    beta_fill[is.na(beta_fill)] <- 0
    eta <- as.vector(X %*% beta_fill) + offset_vec
  }
  mu <- logit_linkinv(eta)
  names(eta) <- names(y_vec)
  names(mu) <- names(y_vec)
  dev <- sum(binomial_dev_resids(y_vec, mu, weights_vec))
  # With an intercept and an offset, glm() refits the intercept-only model
  # so the null deviance includes the offset. A model with no offset uses
  # the weighted mean of y.
  if (intercept && has_offset) {
    X0 <- X[, "(Intercept)", drop = FALSE]
    null_fit <- logit_irls_r(
      X = X0, y = y_vec, weights = weights_vec, offset = offset_vec,
      start = NULL, control = control, singular.ok = TRUE
    )
    if (!isTRUE(null_fit$converged)) {
      warning("fitting to calculate the null deviance did not converge", call. = FALSE)
    }
    beta0 <- null_fit$coefficients
    beta0[is.na(beta0)] <- 0
    mu0 <- logit_linkinv(as.vector(X0 %*% beta0) + offset_vec)
  } else {
    mu0 <- null_mean(y_vec, weights_vec, offset_vec, intercept)
  }
  null_dev <- sum(binomial_dev_resids(y_vec, mu0, weights_vec))
  ll <- binomial_loglik(y_vec, mu, weights_vec, n_vec)
  rank <- sum(!is.na(beta))
  n_ok <- sum(weights_vec > 0)
  V <- vcov_at_fit(X, mu, weights_vec, beta)

  mu_eta <- mu * (1 - mu)
  working <- rep(0, length(mu))
  usable <- is.finite(mu_eta) & weights_vec > 0
  working[usable] <- weights_vec[usable] * mu_eta[usable]
  names(working) <- names(y_vec)
  names(weights_vec) <- names(y_vec)
  dev_resid <- sign(y_vec - mu) * sqrt(pmax(0, binomial_dev_resids(y_vec, mu, weights_vec)))
  names(dev_resid) <- names(y_vec)

  msgs <- fit$warnings
  if (!isTRUE(fit$converged)) {
    msgs <- c(msgs, "logit_fit: algorithm did not converge")
  }
  if (isTRUE(fit$boundary)) {
    msgs <- c(msgs, "logit_fit: algorithm stopped at boundary value")
  }
  eps <- 10 * .Machine$double.eps
  if (any(mu > 1 - eps | mu < eps, na.rm = TRUE)) {
    msgs <- c(msgs, "logit_fit: fitted probabilities numerically 0 or 1 occurred")
  }
  for (msg in unique(msgs)) {
    warning(msg, call. = FALSE)
  }

  structure(
    list(
      coefficients = beta,
      vcov = V,
      fitted.values = mu,
      linear.predictors = eta,
      residuals = y_vec - mu,
      rank = rank,
      df.residual = n_ok - rank,
      df.null = n_ok - as.integer(intercept),
      deviance = dev,
      null.deviance = null_dev,
      aic = -2 * ll + 2 * rank,
      loglik = ll,
      iter = fit$iter,
      converged = isTRUE(fit$converged),
      boundary = isTRUE(fit$boundary),
      prior.weights = weights_vec,
      weights = working,
      y = if (keep_y) y_vec else NULL,
      deviance_resid = dev_resid,
      n_trials = n_vec,
      nobs = length(y_vec),
      terms = mt,
      formula = formula(mt),
      call = call,
      contrasts = contrasts_used,
      xlevels = xlevels,
      model = if (keep_model) mf else NULL,
      x = if (keep_x) X else NULL,
      na.action = attr(mf, "na.action"),
      offset = offset_vec,
      control = control,
      engine = fit$engine,
      intercept = intercept,
      singular.ok = singular.ok
    ),
    class = "logitfit"
  )
}

#' Run the selected engine, falling back to R if C++ cannot factor X'WX.
#'
#' @noRd
run_engine <- function(chosen, X, y, weights, offset, start, control, singular.ok) {
  # Aliased columns are a property of X. The QR engine reports them as NA.
  # The Cholesky engine is used for the full-rank problems that are timed.
  if (!identical(chosen, "r") && ncol(X) > 0L && qr(X)$rank < ncol(X)) {
    chosen <- "r"
  }
  if (identical(chosen, "r") || ncol(X) == 0L) {
    out <- logit_irls_r(X, y, weights, offset, start, control, singular.ok)
    out$engine <- "r"
    return(out)
  }

  init <- initial_eta(X, y, weights, offset, start)
  beta0 <- init$beta
  out <- tryCatch(
    logit_irls_cpp(
      X = X,
      y = y,
      weights = weights,
      offset = offset,
      eta = init$eta,
      beta0 = beta0,
      has_start = init$has_start,
      maxit = control$maxit,
      epsilon = control$epsilon,
      trace = isTRUE(control$trace)
    ),
    error = function(e) list(refit = TRUE, cpp_error = conditionMessage(e))
  )
  if (isTRUE(out$refit)) {
    detail <- out$cpp_error
    if (is.null(detail)) {
      warning(
        "Rcpp engine could not factor the weighted cross-product matrix; refitting with the R engine.",
        call. = FALSE
      )
    } else {
      warning(
        "Rcpp engine failed (", detail, "); refitting with the R engine.",
        call. = FALSE
      )
    }
    out <- logit_irls_r(X, y, weights, offset, start, control, singular.ok)
    out$engine <- "r"
    return(out)
  }
  out$engine <- "rcpp"
  out
}
