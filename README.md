# logitfit

`logitfit` fits binomial logistic regression with the logit link by Fisher scoring (iteratively reweighted least squares). The scoring update is implemented in this package, in R and in C++. It is not a wrapper around `stats::glm()`.

The model is

```text
P(Y = 1 | x) = 1 / (1 + exp(-(x'beta + offset)))
```

Each iteration forms a working response and a weight `mu * (1 - mu)`, solves a weighted least-squares problem for the next coefficient, and stops when the relative change in deviance is small. Wald standard errors come from the expected information matrix at the fitted probabilities. The dispersion is fixed at 1.

## Installation

From this repository, once it is public:

```r
remotes::install_github("acsheth/BIOSTAT625_HW4", build_vignettes = TRUE)
```

From a local copy of the package directory:

```r
install.packages("path/to/logitfit", repos = NULL, type = "source")
```

Windows needs [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching the installed R version, because the package compiles C++.

## Example

```r
library(logitfit)

fit <- logit_fit(died ~ age + severity + treatment + elective, data = icu_sim)
summary(fit)

predict(
  fit,
  newdata = data.frame(age = 70, severity = 3, treatment = 1, elective = 0),
  type = "response",
  se.fit = TRUE
)
```

`icu_sim` is a 250-row simulated intensive-care sample included in the package. The data-generating coefficients are documented in `?icu_sim`.

`coef()`, `vcov()`, `confint()`, `fitted()`, `residuals()`, `predict()`, `logLik()`, `AIC()`, and `deviance()` all have methods. `predict(..., type = "response")` returns probabilities. `type = "link"` returns log-odds.

## Check against `glm`

```r
cmp <- compare_with_glm(
  died ~ age + severity + treatment + elective,
  data = icu_sim
)
cmp$agreement
cmp$benchmark
```

`agreement` is a table of `all.equal()` results for coefficients, fitted values, linear predictors, deviance, null deviance, and AIC. `benchmark` is a `bench::mark()` timing of `glm`, the pure R engine, and the C++ engine. The C++ code does not have to be faster than `glm()`, which is itself written in C. It should usually be faster than the pure R loop.

The same comparison is written up in the vignette:

```r
vignette("logistic-regression", package = "logitfit")
```

## Engines

| `engine` | What it does |
| --- | --- |
| `"r"` | QR decomposition of the weighted model matrix, in R |
| `"rcpp"` | Cholesky factorization of `X'WX`, in C++ via Rcpp |
| `"auto"` | C++ when the DLL is loaded, otherwise R |

Rank-deficient model matrices are fit with the R engine so aliased coefficients can be returned as `NA`. Pass `singular.ok = FALSE` to make that an error.

## Documentation and tests

Help pages: `?logit_fit`, `?logit_control`, `?compare_with_glm`, `?icu_sim`.

```r
testthat::test_local()   # from a source checkout, with testthat installed
# or, after installation:
testthat::test_package("logitfit")
```

GitHub Actions runs `R CMD check` on Ubuntu, Windows, and macOS, and a coverage job that fails if R line coverage is below 90%.
