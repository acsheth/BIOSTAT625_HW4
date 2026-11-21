#' Simulated ICU mortality data
#'
#' A simulated sample of 250 intensive-care patients. In-hospital death was
#' drawn from a logistic model so the data set has a known data-generating
#' process. It is included so examples and the vignette can fit a model
#' without downloading an external file.
#'
#' The linear predictor used to simulate `died` is
#'
#' `-3.1 + 0.03 * age + 0.50 * (severity - 1) - 0.70 * treatment - 0.40 * elective`.
#'
#' @format A data frame with 250 rows and 5 variables:
#' \describe{
#'   \item{died}{In-hospital death, 1 = died and 0 = survived.}
#'   \item{age}{Age in years, rounded to one decimal place and restricted to 18-95.}
#'   \item{severity}{Integer severity score from 1 (lowest) to 4 (highest).}
#'   \item{treatment}{Treatment indicator, 1 = treated.}
#'   \item{elective}{Elective admission indicator, 1 = elective.}
#' }
#' @source Simulated inside this package. See `data-raw/icu_sim.R` for the
#'   generator. The seed is 625.
#'
#' @examples
#' head(icu_sim)
#' logit_fit(died ~ age + severity + treatment + elective, data = icu_sim)
#'
#' @docType data
#' @keywords datasets
#' @name icu_sim
NULL
