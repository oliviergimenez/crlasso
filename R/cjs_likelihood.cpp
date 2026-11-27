#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

using namespace Rcpp;

// Inverse-logit function
inline double inv_logit_cpp(const double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

//' CJS negative log-likelihood evaluated in C++
//'
//' Survival is modelled using X_phi and detection is intercept-only. The
//' likelihood is conditional on first capture and evaluated with a two-state
//' (alive/dead) forward recursion.
//'
//' @param par Survival regression coefficients followed by logit(p)
//' @param histories Binary encounter-history matrix, individuals in rows
//' @param X_phi Survival design matrix, one row per interval
//' @param first_capture First-capture occasions
//' @param frequencies Frequencies for encounter histories
//' @return The negative log-likelihood
// [[Rcpp::export]]
double cjs_nll_rcpp(
    const NumericVector& par,
    const IntegerMatrix& histories,
    const NumericMatrix& X_phi,
    const IntegerVector& first_capture,
    const NumericVector& frequencies) {

  const int n_histories = histories.nrow();
  const int n_occasions = histories.ncol();
  const int n_intervals = n_occasions - 1;
  const int n_phi = X_phi.ncol();

  if (X_phi.nrow() != n_intervals) {
    stop("`X_phi` must have one row per interval.");
  }
  if (par.size() != n_phi + 1) {
    stop("The parameter vector has the wrong length.");
  }
  if (first_capture.size() != n_histories || frequencies.size() != n_histories) {
    stop("Data vectors have incompatible lengths.");
  }

  // Compute interval-specific survival probabilities once per likelihood call
  NumericVector phi(n_intervals);
  for (int t = 0; t < n_intervals; ++t) {
    double eta = 0.0;
    for (int k = 0; k < n_phi; ++k) {
      eta += X_phi(t, k) * par[k];
    }
    phi[t] = inv_logit_cpp(eta);
  }
  const double p = inv_logit_cpp(par[n_phi]);

  double total_log_likelihood = 0.0;

  for (int i = 0; i < n_histories; ++i) {
    const int first = first_capture[i] - 1; // convert to zero-based indexing

    if (first < 0 || first >= n_occasions) {
      stop("Invalid first-capture occasion.");
    }
    if (first == n_occasions - 1) {
      continue;
    }

    // Conditional on first capture, the individual is alive
    double alive = 1.0;
    double dead = 0.0;
    double log_scale = 0.0;

    for (int interval = first; interval < n_intervals; ++interval) {
      // State transition
      const double alive_previous = alive;
      alive = alive_previous * phi[interval];
      dead = dead + alive_previous * (1.0 - phi[interval]);

      // Observation
      const int observation = histories(i, interval + 1);
      if (observation == 1) {
        alive *= p;
        dead = 0.0;
      } else if (observation == 0) {
        alive *= (1.0 - p);
        // A dead individual is necessarily not detected
      } else {
        stop("`histories` must contain only zeros and ones.");
      }

      // Scale the forward probabilities to avoid numerical underflow
      const double scale = alive + dead;
      if (!R_finite(scale) || scale <= 0.0) {
        return R_PosInf;
      }
      alive /= scale;
      dead /= scale;
      log_scale += std::log(scale);
    }

    total_log_likelihood += frequencies[i] * log_scale;
  }

  return -total_log_likelihood;
}
