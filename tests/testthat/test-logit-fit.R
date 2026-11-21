# Correctness against stats::glm, plus the edge cases the scorer can hit.

sim_binary <- function(n = 120, seed = 625) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  grp <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  off <- rnorm(n, sd = 0.15)
  eta <- -0.2 + 0.7 * x1 - 0.4 * x2 + off
  y <- rbinom(n, 1, plogis(eta))
  data.frame(y = y, x1 = x1, x2 = x2, grp = grp, off = off)
}

# Collect every warning from expr. expect_warning() only claims the first one,
# and a second warning still shows up in the test log.
capture_warnings <- function(expr) {
  msgs <- character()
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  msgs
}

expect_matches_glm <- function(formula, data, tolerance = 1e-6, extra = list()) {
  fit <- do.call(
    logit_fit,
    c(list(formula = formula, data = data, engine = "r"), extra)
  )
  glm_extra <- extra[intersect(names(extra), c("weights", "subset", "na.action", "offset", "start", "contrasts"))]
  g <- do.call(
    glm,
    c(list(formula = formula, data = data, family = binomial()), glm_extra)
  )
  expect_s3_class(fit, "logitfit")
  expect_equal(coef(fit), coef(g), tolerance = tolerance)
  expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(g)), tolerance = tolerance)
  expect_equal(fit$deviance, g$deviance, tolerance = tolerance)
  expect_equal(fit$null.deviance, g$null.deviance, tolerance = tolerance)
  expect_equal(unname(fit$aic), unname(g$aic), tolerance = tolerance)
  expect_equal(as.numeric(fit$linear.predictors), as.numeric(g$linear.predictors), tolerance = tolerance)
  list(fit = fit, glm = g)
}

test_that("the pure R engine matches glm on a main-effects model", {
  d <- sim_binary()
  expect_matches_glm(y ~ x1 + x2, d)
})

test_that("auto and rcpp engines match glm when C++ is available", {
  d <- sim_binary()
  skip_if_not(logitfit:::rcpp_engine_available(), "Rcpp engine is not compiled")
  fit_r <- logit_fit(y ~ x1 + x2, data = d, engine = "r")
  fit_cpp <- logit_fit(y ~ x1 + x2, data = d, engine = "rcpp")
  fit_auto <- logit_fit(y ~ x1 + x2, data = d, engine = "auto")
  g <- glm(y ~ x1 + x2, data = d, family = binomial())
  expect_identical(fit_cpp$engine, "rcpp")
  expect_identical(fit_auto$engine, "rcpp")
  expect_equal(coef(fit_cpp), coef(g), tolerance = 1e-6)
  expect_equal(coef(fit_auto), coef(fit_cpp), tolerance = 1e-8)
  expect_equal(fit_cpp$deviance, g$deviance, tolerance = 1e-6)
  expect_equal(vcov(fit_cpp), vcov(fit_r), tolerance = 1e-6)
})

test_that("factors, interactions, and no-intercept formulas match glm", {
  d <- sim_binary()
  expect_matches_glm(y ~ x1 * grp, d)
  expect_matches_glm(y ~ 0 + x1 + x2, d)
  expect_matches_glm(y ~ 1, d)
})

test_that("weights, offset, subset, and starting values match glm", {
  d <- sim_binary()
  d$w <- sample(1:4, nrow(d), replace = TRUE)
  expect_matches_glm(y ~ x1 + x2, d, extra = list(weights = d$w))
  expect_matches_glm(y ~ x1 + offset(off), d)
  expect_matches_glm(y ~ x1, d, extra = list(offset = d$off))
  expect_matches_glm(y ~ x1, d, extra = list(subset = d$x1 > 0))
  g <- glm(y ~ x1, data = d, family = binomial())
  out <- expect_matches_glm(y ~ x1, d, extra = list(start = coef(g)))
  expect_lt(out$fit$iter, 5L)
})

test_that("grouped binomial data and a factor response match glm", {
  d <- sim_binary(n = 80, seed = 7)
  trials <- 6
  eta <- -0.1 + 0.5 * d$x1
  d$success <- rbinom(nrow(d), trials, plogis(eta))
  d$failure <- trials - d$success
  expect_matches_glm(cbind(success, failure) ~ x1 + x2, d)

  d$outcome <- factor(ifelse(d$y == 1, "died", "lived"), levels = c("lived", "died"))
  expect_matches_glm(outcome ~ x1, d)
})

test_that("the bundled ICU sample matches glm", {
  fit <- logit_fit(died ~ age + severity + treatment + elective, data = icu_sim, engine = "r")
  g <- glm(died ~ age + severity + treatment + elective, data = icu_sim, family = binomial())
  expect_equal(coef(fit), coef(g), tolerance = 1e-6)
  expect_equal(fit$deviance, g$deviance, tolerance = 1e-6)
  expect_false(anyNA(coef(fit)))
  expect_true(fit$converged)
  expect_true(fit$iter >= 1L)
})

test_that("missing rows are omitted the same way as glm", {
  d <- sim_binary()
  d$x1[c(1, 5, 8)] <- NA
  out <- expect_matches_glm(y ~ x1 + x2, d)
  expect_lt(nobs(out$fit), nrow(d))
  expect_equal(nobs(out$fit), nobs(out$glm))
})

test_that("na.exclude puts the missing prediction back", {
  d <- sim_binary(n = 40, seed = 3)
  d$x1[4] <- NA
  fit <- logit_fit(y ~ x1, data = d, na.action = na.exclude, engine = "r")
  pred <- predict(fit, type = "response")
  expect_equal(length(pred), nrow(d))
  expect_true(is.na(pred[4]))
})

test_that("variables can live in the calling environment", {
  n <- 50
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.25 * x))
  fit <- logit_fit(y ~ x, engine = "r")
  g <- glm(y ~ x, family = binomial())
  expect_equal(coef(fit), coef(g), tolerance = 1e-6)
})

test_that("a rank-deficient design sets aliased coefficients to NA", {
  d <- sim_binary()
  d$x2 <- d$x1
  fit <- logit_fit(y ~ x1 + x2, data = d, engine = "rcpp")
  g <- glm(y ~ x1 + x2, data = d, family = binomial())
  expect_identical(fit$engine, "r")
  expect_true(anyNA(coef(fit)))
  expect_equal(coef(fit), coef(g), tolerance = 1e-5)
  expect_error(
    logit_fit(y ~ x1 + x2, data = d, singular.ok = FALSE, engine = "r"),
    "singular"
  )
})

test_that("an empty model matches glm", {
  d <- data.frame(y = c(0, 1, 0, 1, 1, 0, 1))
  fit <- logit_fit(y ~ 0, data = d, engine = "rcpp")
  g <- glm(y ~ 0, data = d, family = binomial())
  expect_identical(fit$engine, "r")
  expect_equal(length(coef(fit)), 0L)
  expect_equal(fit$deviance, g$deviance, tolerance = 1e-8)
  expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(g)), tolerance = 1e-8)
})

test_that("failed convergence and separation warn", {
  d <- sim_binary(n = 40, seed = 11)
  expect_warning(
    logit_fit(y ~ x1, data = d, control = logit_control(maxit = 1), engine = "r"),
    "did not converge"
  )
  x <- seq(-3, 3, length.out = 24)
  y <- as.numeric(x > 0)
  msgs <- capture_warnings(logit_fit(y ~ x, engine = "r"))
  expect_true(any(grepl("fitted probabilities|did not converge", msgs)))
})

test_that("zero weights and a divergent start are reported", {
  d <- sim_binary(n = 30, seed = 4)
  msgs <- capture_warnings(
    logit_fit(y ~ x1, data = d, weights = rep(0, nrow(d)), engine = "r")
  )
  expect_true(any(grepl("no observations informative|did not converge", msgs)))
  out <- logitfit:::logit_irls_r(
    X = cbind(`(Intercept)` = 1, x = d$x1),
    y = d$y,
    weights = rep(1, nrow(d)),
    offset = rep(0, nrow(d)),
    start = c(0, 1e8),
    control = logit_control(maxit = 10),
    singular.ok = TRUE
  )
  expect_true(any(grepl("step size", out$warnings)) || isTRUE(out$boundary) || !isTRUE(out$converged))
})

test_that("input validation stops bad calls", {
  d <- sim_binary(n = 20, seed = 2)
  expect_error(logit_fit("y ~ x1", data = d), "formula")
  expect_error(logit_control(epsilon = 0), "epsilon")
  expect_error(logit_control(maxit = 1.5), "whole number")
  expect_error(logit_control(maxit = 0), "maxit")
  expect_error(logit_control(trace = "yes"), "trace")
  expect_error(logit_fit(y ~ x1, data = d, control = list(epsilon = 1)), "control")
  expect_error(logit_fit(y ~ x1, data = d, singular.ok = NA), "singular")
  expect_error(logit_fit(y ~ x1, data = d, model = "yes"), "model")
  expect_error(logit_fit(y ~ x1, data = d, x = NA), "'x'")
  expect_error(logit_fit(y ~ x1, data = d, y = "no"), "'y'")
  expect_error(logit_fit(y ~ x1, data = d, engine = "fortran"), "should be one of")
  bad <- d
  bad$y[1] <- 1.5
  expect_error(logit_fit(y ~ x1, data = bad, engine = "r"), "0 <= y <= 1")
  expect_error(logit_fit(y ~ x1, data = d, weights = -d$x1^2, engine = "r"), "negative")
  expect_error(logit_fit(y ~ x1, data = d, start = 1, engine = "r"), "length")
  expect_error(logit_fit(y ~ x1, data = d, start = c(a = 0), engine = "r"), "missing")
  expect_error(
    logit_fit(cbind(y, x1, x2) ~ grp, data = d, engine = "r"),
    "two-column"
  )
  d$three <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_error(logit_fit(three ~ x1, data = d, engine = "r"), "two levels")
  expect_error(logitfit:::prepare_binomial_response(c("alive", "dead"), NULL), "numeric")
  expect_warning(
    logit_fit(y ~ x1, data = d, weights = rep(0.5, nrow(d)), engine = "r"),
    "non-integer"
  )
})

test_that("the Rcpp engine falls back when it is turned off", {
  d <- sim_binary(n = 30, seed = 8)
  old <- options(logitfit.force_no_rcpp = TRUE)
  on.exit(options(old), add = TRUE)
  expect_false(logitfit:::rcpp_engine_available())
  expect_warning(
    fit <- logit_fit(y ~ x1, data = d, engine = "auto"),
    "not available"
  )
  expect_identical(fit$engine, "r")
  expect_error(logit_fit(y ~ x1, data = d, engine = "rcpp"), "not available")
})

test_that("trace prints the deviance", {
  d <- sim_binary(n = 40, seed = 9)
  expect_output(
    logit_fit(y ~ x1, data = d, engine = "r", control = logit_control(trace = TRUE)),
    "Deviance"
  )
})

test_that("logical responses and bad binomial inputs are handled", {
  d <- sim_binary(n = 40, seed = 16)
  d$y <- as.logical(d$y)
  expect_matches_glm(y ~ x1, d)
  expect_error(logitfit:::prepare_binomial_response(cbind(c(1, -1), c(0, 1)), NULL), "negative")
  expect_error(logitfit:::prepare_binomial_response(c(0, NA_real_), NULL), "missing")
  expect_error(logitfit:::prepare_binomial_response(c(0, 1), c(1, Inf)), "finite")
  expect_error(logitfit:::prepare_binomial_response(c(0, 1), 1), "different lengths")
})

test_that("a divergent Newton step is halved in both engines", {
  set.seed(1)
  d <- data.frame(x = rnorm(40), y = rbinom(40, 1, 0.3))
  start <- c(`(Intercept)` = 0, x = 20)
  expect_warning(
    expect_output(
      fit <- logit_fit(
        y ~ x, data = d, start = start, engine = "r",
        control = logit_control(maxit = 15, trace = TRUE)
      ),
      "Step halved"
    ),
    "step size truncated"
  )
  expect_true(fit$converged)
  expect_error(
    logit_fit(y ~ x, data = d, start = start, engine = "r", control = logit_control(maxit = 1)),
    "cannot correct step size"
  )

  skip_if_not(logitfit:::rcpp_engine_available(), "Rcpp engine is not compiled")
  expect_warning(
    expect_output(
      fit_cpp <- logit_fit(
        y ~ x, data = d, start = start, engine = "rcpp",
        control = logit_control(maxit = 15, trace = TRUE)
      ),
      "Step halved"
    ),
    "step size truncated"
  )
  expect_identical(fit_cpp$engine, "rcpp")
  expect_true(fit_cpp$converged)
  expect_warning(
    expect_error(
      logit_fit(y ~ x, data = d, start = start, engine = "rcpp", control = logit_control(maxit = 1)),
      "cannot correct step size"
    ),
    "Rcpp engine failed"
  )

  X <- model.matrix(y ~ x, d)
  eta <- as.vector(X %*% start)
  expect_error(
    logitfit:::logit_irls_cpp(
      X, d$y, rep(1, nrow(d)), rep(0, nrow(d)), eta, c(0, 0), FALSE, 5L, 1e-8, FALSE
    ),
    "no valid set of coefficients"
  )
})

test_that("the C++ engine hands singular and empty problems back to R", {
  skip_if_not(logitfit:::rcpp_engine_available(), "Rcpp engine is not compiled")
  d <- data.frame(y = c(1, 0, 1, 0), x = c(1e200, -1e200, 1, -1))
  expect_warning(
    fit <- logit_fit(y ~ x, data = d, engine = "rcpp"),
    "refitting with the R engine"
  )
  expect_identical(fit$engine, "r")

  d0 <- sim_binary(n = 25, seed = 14)
  msgs <- character()
  fit0 <- withCallingHandlers(
    logit_fit(y ~ x1, data = d0, weights = rep(0, nrow(d0)), engine = "rcpp"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("refitting with the R engine", msgs)))
  expect_identical(fit0$engine, "r")

  expect_output(
    logit_fit(y ~ x1, data = d0, engine = "rcpp", control = logit_control(trace = TRUE)),
    "Deviance"
  )

  X <- cbind(`(Intercept)` = 1, x = d0$x1)
  eta <- rep(0, nrow(d0))
  bad_y <- logitfit:::logit_irls_cpp(
    X, d0$y[-1], rep(1, nrow(d0)), rep(0, nrow(d0)), eta, c(0, 0), FALSE, 2L, 1e-8, FALSE
  )
  expect_true(isTRUE(bad_y$refit))
  empty <- logitfit:::logit_irls_cpp(
    matrix(numeric(0), nrow = 0, ncol = 1),
    numeric(0), numeric(0), numeric(0), numeric(0), numeric(1),
    FALSE, 1L, 1e-8, FALSE
  )
  expect_true(isTRUE(empty$refit))
  eta_inf <- eta
  eta_inf[1] <- Inf
  eta_inf[2] <- -Inf
  finite_tail <- logitfit:::logit_irls_cpp(
    X, d0$y, rep(1, nrow(d0)), rep(0, nrow(d0)), eta_inf, c(0, 0), TRUE, 4L, 1e-8, FALSE
  )
  expect_false(isTRUE(finite_tail$refit))
})

test_that("a non-finite offset or model matrix is rejected", {
  d <- sim_binary(n = 20, seed = 13)
  d$off[1] <- Inf
  expect_error(logit_fit(y ~ x1 + offset(off), data = d, engine = "r"), "offset must be finite")
  d <- sim_binary(n = 12, seed = 2)
  d$x1[1] <- Inf
  expect_error(logit_fit(y ~ x1, data = d, engine = "r"), "model matrix")
  d <- sim_binary(n = 12, seed = 2)
  expect_error(
    logit_fit(y ~ x1, data = d, subset = rep(FALSE, nrow(d)), engine = "r"),
    "no observations"
  )
  expect_error(logitfit:::align_start("a", c("(Intercept)")), "numeric")
  expect_error(logitfit:::align_start(c(a = NA_real_), "a"), "finite")
})

test_that("zero trial counts and non-integer grouped counts are accepted", {
  d <- sim_binary(n = 16, seed = 5)
  d$success <- c(0, d$y[-1])
  d$failure <- c(0, 1 - d$y[-1])
  fit <- logit_fit(cbind(success, failure) ~ x1, data = d, engine = "r")
  expect_s3_class(fit, "logitfit")
  d$success[2] <- 1.4
  expect_warning(
    logit_fit(cbind(success, failure) ~ x1, data = d, engine = "r"),
    "non-integer"
  )
})

test_that("an intercept-only null fit can fail to converge", {
  d <- sim_binary(n = 30, seed = 6)
  msgs <- capture_warnings(
    logit_fit(
      y ~ x1 + offset(off), data = d, engine = "r",
      control = logit_control(maxit = 1)
    )
  )
  expect_true(any(grepl("null deviance did not converge", msgs)))
})

test_that("a one-column response matrix is treated as a vector", {
  d <- sim_binary(n = 40, seed = 12)
  d$ym <- cbind(d$y)
  fit <- logit_fit(ym ~ x1, data = d, engine = "r")
  g <- glm(y ~ x1, data = d, family = binomial())
  expect_equal(coef(fit), coef(g), tolerance = 1e-6)
})
