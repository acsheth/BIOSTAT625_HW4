// Fisher scoring (IRLS) for binomial logistic regression.
//
// Each iteration builds the expected information X' W X and the score
// equation X' W z, with W_ii = prior_i * mu_i * (1 - mu_i) and
// z_i = eta_i - offset_i + (y_i - mu_i) / (mu_i (1 - mu_i)).
// The linear system is solved by a Cholesky factorization after column
// scaling. If the information matrix is numerically singular, the routine
// asks the R wrapper to refit with the QR implementation, which can mark
// aliased coefficients explicitly.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;

namespace {

double logit_mu(double eta) {
  if (!std::isfinite(eta)) {
    if (eta > 0) return 1.0;
    if (eta < 0) return 0.0;
    return NA_REAL;
  }
  if (eta >= 0) {
    double e = std::exp(-eta);
    return 1.0 / (1.0 + e);
  }
  double e = std::exp(eta);
  return e / (1.0 + e);
}

double ylogy(double y, double mu) {
  if (!(y > 0.0)) return 0.0;
  return y * std::log(y / mu);
}

double dev_obs(double y, double mu, double wt) {
  return 2.0 * wt * (ylogy(y, mu) + ylogy(1.0 - y, 1.0 - mu));
}

// In-place lower Cholesky factor of a row-major p by p matrix.
// Returns false when a pivot is not safely positive.
bool cholesky_decompose(std::vector<double>& A, int p) {
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = A[i * p + j];
      for (int k = 0; k < j; ++k) {
        sum -= A[i * p + k] * A[j * p + k];
      }
      if (i == j) {
        if (!(sum > 1e-12 * (1.0 + std::fabs(A[i * p + i])))) {
          return false;
        }
        A[i * p + i] = std::sqrt(sum);
      } else {
        double pivot = A[j * p + j];
        if (!(pivot > 0.0)) return false;
        A[i * p + j] = sum / pivot;
      }
    }
    for (int j = i + 1; j < p; ++j) {
      A[i * p + j] = 0.0;
    }
  }
  return true;
}

// Solve L L' x = b. L is the lower factor stored row-major. b is overwritten.
void chol_solve(const std::vector<double>& L, int p, std::vector<double>& b) {
  for (int i = 0; i < p; ++i) {
    double sum = b[i];
    for (int k = 0; k < i; ++k) sum -= L[i * p + k] * b[k];
    b[i] = sum / L[i * p + i];
  }
  for (int i = p - 1; i >= 0; --i) {
    double sum = b[i];
    for (int k = i + 1; k < p; ++k) sum -= L[k * p + i] * b[k];
    b[i] = sum / L[i * p + i];
  }
}

}  // namespace

// [[Rcpp::export]]
List logit_irls_cpp(NumericMatrix X,
                    NumericVector y,
                    NumericVector weights,
                    NumericVector offset,
                    NumericVector eta,
                    NumericVector beta0,
                    bool has_start,
                    int maxit,
                    double epsilon,
                    bool trace) {
  int n = X.nrow();
  int p = X.ncol();
  if (n == 0 || p == 0 ||
      y.size() != n || weights.size() != n ||
      offset.size() != n || eta.size() != n ||
      beta0.size() != p) {
    return List::create(Named("refit") = true);
  }

  std::vector<double> beta(p, 0.0);
  std::vector<double> beta_old(p, 0.0);
  bool have_old = has_start;
  if (has_start) {
    for (int j = 0; j < p; ++j) beta[j] = beta0[j];
    beta_old = beta;
  }

  std::vector<double> eta_vec(n), mu(n);
  for (int i = 0; i < n; ++i) eta_vec[i] = eta[i];

  auto refresh = [&]() {
    for (int i = 0; i < n; ++i) {
      double dot = offset[i];
      for (int j = 0; j < p; ++j) dot += X(i, j) * beta[j];
      eta_vec[i] = dot;
      mu[i] = logit_mu(dot);
    }
  };
  // The first linear predictor is supplied by R (starting values or the
  // binomial mustart). Later iterations recompute it from beta.
  for (int i = 0; i < n; ++i) mu[i] = logit_mu(eta_vec[i]);

  auto deviance = [&]() {
    double dev = 0.0;
    for (int i = 0; i < n; ++i) dev += dev_obs(y[i], mu[i], weights[i]);
    return dev;
  };

  double dev = deviance();
  double devold = dev;
  bool conv = false;
  bool boundary = false;
  int iter_done = 0;
  std::vector<std::string> warnings;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_done = iter;
    std::vector<double> scale(p, 0.0);
    std::vector<double> A(p * p, 0.0);
    std::vector<double> rhs(p, 0.0);
    bool any_good = false;

    for (int i = 0; i < n; ++i) {
      double mu_eta = mu[i] * (1.0 - mu[i]);
      if (!(weights[i] > 0.0) || !(mu_eta > 0.0) || !std::isfinite(mu[i])) {
        continue;
      }
      any_good = true;
      double wi = weights[i] * mu_eta;
      double zi = (eta_vec[i] - offset[i]) + (y[i] - mu[i]) / mu_eta;
      for (int j = 0; j < p; ++j) {
        double xij = X(i, j);
        scale[j] += xij * xij * wi;
        rhs[j] += xij * wi * zi;
        for (int k = 0; k <= j; ++k) {
          A[j * p + k] += xij * X(i, k) * wi;
        }
      }
    }

    if (!any_good) {
      // The R engine owns this edge case so the warning matches logit_irls_r().
      return List::create(Named("refit") = true);
    }

    bool singular = false;
    for (int j = 0; j < p; ++j) {
      scale[j] = std::sqrt(scale[j]);
      if (!(scale[j] > 1e-12)) singular = true;
    }
    if (singular) return List::create(Named("refit") = true);

    // Scale columns so the normal equations are better conditioned, then
    // mirror the lower triangle into the upper triangle for Cholesky.
    for (int j = 0; j < p; ++j) {
      rhs[j] /= scale[j];
      for (int k = 0; k <= j; ++k) {
        A[j * p + k] /= scale[j] * scale[k];
      }
    }
    for (int j = 0; j < p; ++j) {
      for (int k = j + 1; k < p; ++k) A[j * p + k] = A[k * p + j];
    }

    if (!cholesky_decompose(A, p)) {
      return List::create(Named("refit") = true);
    }
    std::vector<double> delta = rhs;
    chol_solve(A, p, delta);
    bool finite_step = true;
    for (int j = 0; j < p; ++j) {
      delta[j] /= scale[j];
      if (!std::isfinite(delta[j])) finite_step = false;
    }
    if (!finite_step) {
      return List::create(Named("refit") = true);
    }

    beta.swap(delta);
    refresh();
    double dev_new = deviance();
    if (trace) {
      Rcout << "Deviance = " << dev_new << " Iterations - " << iter << std::endl;
    }

    boundary = false;
    if (!std::isfinite(dev_new)) {
      if (!have_old) {
        stop("no valid set of coefficients has been found: please supply starting values");
      }
      warnings.push_back("step size truncated due to divergence");
      int guard = 1;
      while (!std::isfinite(dev_new)) {
        if (guard > maxit) {
          stop("inner loop 1; cannot correct step size");
        }
        ++guard;
        for (int j = 0; j < p; ++j) beta[j] = 0.5 * (beta[j] + beta_old[j]);
        refresh();
        dev_new = deviance();
      }
      boundary = true;
      if (trace) {
        Rcout << "Step halved: new deviance = " << dev_new << std::endl;
      }
    }

    dev = dev_new;
    if (std::fabs(dev - devold) / (0.1 + std::fabs(dev)) < epsilon) {
      conv = true;
      break;
    }
    devold = dev;
    beta_old = beta;
    have_old = true;
  }

  NumericVector coef(p);
  NumericVector fitted(n);
  NumericVector linear(n);
  for (int j = 0; j < p; ++j) coef[j] = beta[j];
  for (int i = 0; i < n; ++i) {
    fitted[i] = mu[i];
    linear[i] = eta_vec[i];
  }
  coef.attr("names") = colnames(X);
  fitted.attr("names") = y.attr("names");
  linear.attr("names") = y.attr("names");

  CharacterVector warn_r(warnings.size());
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(warnings.size()); ++i) {
    warn_r[i] = warnings[i];
  }

  return List::create(
    Named("coefficients") = coef,
    Named("fitted.values") = fitted,
    Named("linear.predictors") = linear,
    Named("deviance") = dev,
    Named("iter") = iter_done,
    Named("converged") = conv,
    Named("boundary") = boundary,
    Named("warnings") = warn_r,
    Named("refit") = false
  );
}
