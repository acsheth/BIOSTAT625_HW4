sim_binary <- function(n = 80, seed = 21) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  city <- factor(sample(c("A", "B"), n, replace = TRUE))
  y <- rbinom(n, 1, plogis(-0.15 + 0.8 * x - 0.3 * z + 0.2 * (city == "B")))
  data.frame(y = y, x = x, z = z, city = city)
}

test_that("summary, coefficients, likelihood, and residuals agree with glm", {
  d <- sim_binary()
  fit <- logit_fit(y ~ x + z + city, data = d, engine = "r")
  g <- glm(y ~ x + z + city, data = d, family = binomial())
  sm <- summary(fit)
  expect_s3_class(sm, "summary.logitfit")
  expect_equal(coef(sm)[, "Estimate"], coef(g), tolerance = 1e-6)
  expect_equal(coef(sm)[, "Std. Error"], sqrt(diag(vcov(g))), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(AIC(fit), AIC(g), tolerance = 1e-6)
  expect_equal(deviance(fit), deviance(g), tolerance = 1e-6)
  expect_equal(nobs(fit), nobs(g))
  expect_equal(deparse1(formula(fit)), deparse1(formula(g)))
  expect_equal(family(fit)$family, "binomial")
  expect_equal(residuals(fit, type = "response"), residuals(g, type = "response"), tolerance = 1e-6)
  expect_equal(residuals(fit, type = "deviance"), residuals(g, type = "deviance"), tolerance = 1e-5)
  expect_equal(residuals(fit, type = "pearson"), residuals(g, type = "pearson"), tolerance = 1e-5)
  expect_equal(as.numeric(weights(fit, type = "prior")), as.numeric(weights(g, type = "prior")))
  expect_output(print(fit), "Coefficients")
  expect_output(print(sm), "Pr\\(>\\|z\\|\\)")
  expect_output(print(sm), "Fisher scoring")
})

test_that("predictions and Wald intervals match glm on the link and response scales", {
  d <- sim_binary()
  fit <- logit_fit(y ~ x + city, data = d, x = TRUE, engine = "r")
  g <- glm(y ~ x + city, data = d, family = binomial())
  expect_equal(predict(fit, type = "link"), predict(g, type = "link"), tolerance = 1e-6)
  expect_equal(predict(fit, type = "response"), predict(g, type = "response"), tolerance = 1e-6)
  nd <- data.frame(x = c(-1, 0, 1), city = factor(c("A", "B", "A"), levels = levels(d$city)))
  expect_equal(
    as.numeric(predict(fit, newdata = nd, type = "response")),
    as.numeric(predict(g, newdata = nd, type = "response")),
    tolerance = 1e-6
  )
  se_fit <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
  se_glm <- predict(g, newdata = nd, type = "link", se.fit = TRUE)
  expect_equal(as.numeric(se_fit$fit), as.numeric(se_glm$fit), tolerance = 1e-6)
  expect_equal(as.numeric(se_fit$se.fit), as.numeric(se_glm$se.fit), tolerance = 1e-4)
  se_resp <- predict(fit, type = "response", se.fit = TRUE)
  expect_equal(length(se_resp$se.fit), nrow(d))
  ci <- confint(fit, level = 0.9)
  ci_g <- confint.default(g, level = 0.9)
  expect_equal(ci, ci_g[rownames(ci), , drop = FALSE], tolerance = 1e-4)
  expect_equal(rownames(confint(fit, parm = 2)), names(coef(fit))[2])
  expect_equal(rownames(confint(fit, parm = "x")), "x")
  expect_error(confint(fit, parm = "missing"), "unknown")
  expect_error(confint(fit, parm = 99), "out of range")
  expect_error(confint(fit, level = 1.2), "level")
})

test_that("stored frames can be dropped and recovered", {
  d <- sim_binary(n = 40, seed = 6)
  fit <- logit_fit(y ~ x, data = d, model = FALSE, x = TRUE, y = FALSE, engine = "r")
  expect_null(fit$model)
  expect_null(fit$y)
  expect_error(model.frame(fit), "model frame")
  expect_error(residuals(fit), "response was not stored")
  expect_s3_class(summary(fit), "summary.logitfit")
  expect_equal(ncol(model.matrix(fit)), 2L)
  bare <- logit_fit(y ~ x, data = d, model = FALSE, x = FALSE, engine = "r")
  expect_error(model.matrix(bare), "model = TRUE")
  expect_error(predict(bare, se.fit = TRUE), "standard errors")
  expect_error(predict(fit, se.fit = "maybe"), "se.fit")
  expect_error(residuals(logit_fit(y ~ x, data = d, engine = "r"), type = "partial"), "should be one of")
  expect_error(weights(logit_fit(y ~ x, data = d, engine = "r"), type = "other"), "should be one of")
})

test_that("print covers an empty model, aliasing, and non-convergence", {
  empty <- logit_fit(y ~ 0, data = data.frame(y = c(0, 1, 0, 1, 1)), engine = "r")
  expect_output(print(empty), "No coefficients")
  d <- sim_binary(n = 40, seed = 6)
  d$x2 <- d$x
  aliased <- logit_fit(y ~ x + x2, data = d, engine = "r")
  expect_output(print(summary(aliased)), "singularities")
  stalled <- suppressWarnings(
    logit_fit(y ~ x, data = d, control = logit_control(maxit = 1), engine = "r")
  )
  expect_output(print(stalled), "did not converge")
  expect_output(print(summary(stalled)), "did not converge")
  expect_true(all(is.na(predict(empty, se.fit = TRUE)$se.fit)))
})

test_that("working residuals and working weights are finite on a regular fit", {
  d <- sim_binary(n = 50, seed = 15)
  fit <- logit_fit(y ~ x + z, data = d, engine = "r")
  expect_true(all(is.finite(residuals(fit, type = "working"))))
  expect_true(all(weights(fit, type = "working") > 0))
  expect_equal(as.numeric(fitted(fit)), as.numeric(predict(fit, type = "response")))
})

test_that("compare_with_glm reports agreement and can skip the timer", {
  d <- sim_binary(n = 60, seed = 18)
  cmp <- compare_with_glm(y ~ x + z, data = d, n_bench = 0)
  expect_s3_class(cmp, "logitfit_comparison")
  expect_true(all(cmp$agreement$all_equal))
  expect_null(cmp$benchmark)
  expect_lt(cmp$max_coef_diff[["r"]], 1e-6)
  expect_output(print(cmp), "Correctness")
  expect_warning(
    compare_with_glm(y ~ x, data = d, engine = "r", n_bench = 0),
    "ignored"
  )
  expect_error(compare_with_glm(y ~ x, data = d, 1, n_bench = 0), "named")
  expect_error(compare_with_glm(y ~ x, data = d, not_an_arg = 1, n_bench = 0), "unknown")
  expect_error(compare_with_glm(y ~ x, n_bench = 0), "data")
  expect_error(compare_with_glm("y ~ x", data = d, n_bench = 0), "formula")
  expect_error(compare_with_glm(y ~ x, data = d, tolerance = 0, n_bench = 0), "tolerance")
  expect_error(compare_with_glm(y ~ x, data = d, n_bench = -1), "n_bench")
})

test_that("compare_with_glm times both engines", {
  d <- sim_binary(n = 40, seed = 19)
  skip_if_not(requireNamespace("bench", quietly = TRUE))
  cmp <- compare_with_glm(y ~ x, data = d, n_bench = 1)
  expect_s3_class(cmp$benchmark, "bench_mark")
  expect_equal(nrow(cmp$benchmark), 3L)
  expect_output(print(cmp), "bench::mark")
})
