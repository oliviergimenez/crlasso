# Functions for simulating Cormack-Jolly-Seber (CJS) capture histories

#' Simulate CJS capture histories
#'
#' Individuals are released on occasions and followed until death or
#' the final capture occasion. Survival is defined for intervals and detection
#' for occasions after first release.
#'
#' @param phi Matrix of survival probabilities. Rows correspond to individuals
#'   and columns to intervals.
#' @param p Matrix of detection probabilities with the same dimensions as
#'   `phi`.
#' @param marked Integer vector giving the number of newly marked individuals
#'   on each occasion. Its length must equal the number of capture occasions.
#' @return A binary matrix with individuals in rows and capture occasions in
#'   columns.
simulate_cjs_histories <- function(phi, p, marked) {
  phi <- as.matrix(phi)
  p <- as.matrix(p)
  marked <- as.integer(marked)

  n_occasions <- ncol(phi) + 1L
  n_individuals <- sum(marked)

  if (!length(marked) %in% c(n_occasions - 1L, n_occasions)) {
    stop("`marked` must have one entry per occasion, optionally excluding the final occasion.")
  }
  if (length(marked) == n_occasions - 1L) marked <- c(marked, 0L)
  if (!identical(dim(phi), dim(p))) {
    stop("`phi` and `p` must have identical dimensions.")
  }
  if (nrow(phi) != n_individuals) {
    stop("The number of rows in `phi` and `p` must equal sum(marked).")
  }
  if (any(phi < 0 | phi > 1) || any(p < 0 | p > 1)) {
    stop("All survival and detection probabilities must lie in [0, 1].")
  }

  histories <- matrix(0L, nrow = n_individuals, ncol = n_occasions)
  release_occasion <- rep(seq_len(n_occasions), marked)

  for (i in seq_len(n_individuals)) {
    first <- release_occasion[i]
    histories[i, first] <- 1L

    if (first == n_occasions) next

    for (occasion in seq.int(first + 1L, n_occasions)) {
      interval <- occasion - 1L
      alive <- stats::rbinom(1L, size = 1L, prob = phi[i, interval])
      if (alive == 0L) break

      histories[i, occasion] <- stats::rbinom(
        1L, size = 1L, prob = p[i, interval]
      )
    }
  }

  histories
}


#' Simulate one CJS dataset
#'
#' @param n_occasions Number of capture occasions.
#' @param marked Number newly marked per occasion.
#' @param beta_phi Regression coefficients for survival, including intercept.
#' @param beta_p Intercept-only detection on the logit scale.
#' @param covariate_correlation Pairwise correlation used when simulating more
#'   than one temporal covariate.
#' @param seed Optional random seed.
#' @return A list containing histories, design matrix, covariates and true
#'   parameters.
simulate_cjs_dataset <- function(
    n_occasions = 15L,
    marked = 15L,
    beta_phi = c(stats::qlogis(0.8), -0.6, 0),
    beta_p = stats::qlogis(0.9),
    covariate_correlation = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  n_occasions <- as.integer(n_occasions)
  n_intervals <- n_occasions - 1L

  if (length(marked) == 1L) {
    marked <- c(rep(as.integer(marked), n_occasions - 1L), 0L)
  }
  if (!length(marked) %in% c(n_occasions - 1L, n_occasions)) {
    stop("`marked` must be a scalar or have length `n_occasions - 1` or `n_occasions`.")
  }
  if (length(marked) == n_occasions - 1L) marked <- c(marked, 0L)

  n_covariates <- length(beta_phi) - 1L
  if (n_covariates < 0L) stop("`beta_phi` must include an intercept.")

  if (n_covariates == 0L) {
    covariates <- matrix(numeric(0), nrow = n_intervals, ncol = 0L)
  } else {
    sigma <- matrix(covariate_correlation, n_covariates, n_covariates)
    diag(sigma) <- 1
    chol_sigma <- chol(sigma)
    covariates <- matrix(
      stats::rnorm(n_intervals * n_covariates),
      nrow = n_intervals,
      ncol = n_covariates
    ) %*% chol_sigma
    covariates <- scale(covariates)
    colnames(covariates) <- paste0("x", seq_len(n_covariates))
  }

  X_phi <- cbind("(Intercept)" = 1, covariates)
  phi_by_interval <- stats::plogis(drop(X_phi %*% beta_phi))
  p_by_interval <- rep(stats::plogis(beta_p), n_intervals)

  n_individuals <- sum(marked)
  phi <- matrix(phi_by_interval, n_individuals, n_intervals, byrow = TRUE)
  p <- matrix(p_by_interval, n_individuals, n_intervals, byrow = TRUE)

  histories <- simulate_cjs_histories(phi = phi, p = p, marked = marked)

  list(
    histories = histories,
    X_phi = X_phi,
    covariates = covariates,
    marked = marked,
    beta_phi = beta_phi,
    beta_p = beta_p,
    phi = phi_by_interval,
    p = p_by_interval
  )
}
