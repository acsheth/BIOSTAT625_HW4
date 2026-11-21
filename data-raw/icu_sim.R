# Generate the bundled ICU sample. Run from the package root:
#   Rscript data-raw/icu_sim.R

set.seed(625)
n <- 250
age <- round(pmin(pmax(stats::rnorm(n, mean = 64, sd = 14), 18), 95), 1)
severity <- sample(1:4, n, replace = TRUE, prob = c(0.35, 0.30, 0.22, 0.13))
treatment <- stats::rbinom(n, 1, 0.5)
elective <- stats::rbinom(n, 1, 0.4)
eta <- -3.1 + 0.03 * age + 0.50 * (severity - 1) - 0.70 * treatment - 0.40 * elective
died <- stats::rbinom(n, 1, stats::plogis(eta))

icu_sim <- data.frame(
  died = died,
  age = age,
  severity = severity,
  treatment = treatment,
  elective = elective
)

dir.create("data", showWarnings = FALSE)
save(icu_sim, file = file.path("data", "icu_sim.rda"), version = 2, compress = "gzip")
message("Wrote data/icu_sim.rda with ", nrow(icu_sim), " rows.")
