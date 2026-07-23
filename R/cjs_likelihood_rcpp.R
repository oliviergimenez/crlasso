# Likelihood functions for an intercept-only detection CJS model.
# The computationally intensive forward algorithm is evaluated in C++.

#' Compile the Rcpp CJS likelihood
#'
#' Run this once at the beginning of an R session. On a cluster, compile on a
#' login node if required by local policy, or let the batch job compile once
#' before starting the simulations.
#'
#' @param cpp_file Path to `cjs_likelihood.cpp`.
#' @param rebuild Force recompilation even when a cached shared library exists.
load_cjs_rcpp <- function(cpp_file = "cjs_likelihood.cpp", rebuild = FALSE) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package `Rcpp` is required.")
  }
  if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Package `RcppArmadillo` is required.")
  }
  if (!file.exists(cpp_file)) {
    stop("C++ source file not found: ", normalizePath(cpp_file, mustWork = FALSE))
  }

  Rcpp::sourceCpp(cpp_file, rebuild = rebuild, verbose = FALSE)
  if (!exists("cjs_nll_rcpp", mode = "function", inherits = TRUE)) {
    stop("Compilation finished, but `cjs_nll_rcpp()` is unavailable.")
  }
  invisible(TRUE)
}


#' Validate and prepare CJS data
prepare_cjs_data <- function(histories, X_phi, frequencies = NULL) {
  histories <- as.matrix(histories)
  storage.mode(histories) <- "integer"
  X_phi <- as.matrix(X_phi)
  storage.mode(X_phi) <- "double"

  if (!all(histories %in% c(0L, 1L))) {
    stop("`histories` must contain only zeros and ones.")
  }
  if (ncol(histories) < 2L) stop("At least two capture occasions are needed.")
  if (nrow(X_phi) != ncol(histories) - 1L) {
    stop("`X_phi` must have one row per interval between capture occasions.")
  }
  if (is.null(frequencies)) frequencies <- rep(1, nrow(histories))
  if (length(frequencies) != nrow(histories) || any(frequencies <= 0)) {
    stop("`frequencies` must contain one positive value per history.")
  }

  first_capture <- apply(histories, 1L, function(x) {
    detected <- which(x == 1L)
    if (length(detected) == 0L) NA_integer_ else detected[1L]
  })
  if (anyNA(first_capture)) {
    stop("Every encounter history must contain at least one capture.")
  }

  if (is.null(colnames(X_phi))) {
    colnames(X_phi) <- paste0("beta_phi_", seq_len(ncol(X_phi)))
  }

  list(
    histories = histories,
    X_phi = X_phi,
    frequencies = as.numeric(frequencies),
    first_capture = as.integer(first_capture),
    n_histories = nrow(histories),
    n_occasions = ncol(histories),
    n_intervals = ncol(histories) - 1L,
    n_phi_parameters = ncol(X_phi),
    effective_sample_size = sum(frequencies)
  )
}


# Pure-R implementation retained for validation and benchmarking.
cjs_nll_R <- function(par, cjs_data) {
  n_phi <- cjs_data$n_phi_parameters
  if (length(par) != n_phi + 1L) {
    stop("The parameter vector has the wrong length.")
  }

  phi <- stats::plogis(drop(cjs_data$X_phi %*% par[seq_len(n_phi)]))
  p <- stats::plogis(par[n_phi + 1L])
  total_log_likelihood <- 0

  for (i in seq_len(cjs_data$n_histories)) {
    first <- cjs_data$first_capture[i]
    if (first == cjs_data$n_occasions) next

    alpha <- c(alive = 1, dead = 0)
    log_scale <- 0

    for (interval in seq.int(first, cjs_data$n_intervals)) {
      alive_previous <- alpha[1L]
      alpha <- c(
        alive = alive_previous * phi[interval],
        dead = alpha[2L] + alive_previous * (1 - phi[interval])
      )

      observation <- cjs_data$histories[i, interval + 1L]
      if (observation == 1L) {
        alpha <- c(alpha[1L] * p, 0)
      } else {
        alpha[1L] <- alpha[1L] * (1 - p)
      }

      scale <- sum(alpha)
      if (!is.finite(scale) || scale <= 0) return(Inf)
      alpha <- alpha / scale
      log_scale <- log_scale + log(scale)
    }

    total_log_likelihood <- total_log_likelihood +
      cjs_data$frequencies[i] * log_scale
  }

  -total_log_likelihood
}


#' Negative conditional log-likelihood for a CJS model
#'
#' This wrapper keeps the same interface as the original R function, so the
#' adaptive-Lasso code does not need to change.
cjs_nll <- function(par, cjs_data) {
  if (!exists("cjs_nll_rcpp", mode = "function", inherits = TRUE)) {
    stop(
      "The compiled function `cjs_nll_rcpp()` is unavailable. ",
      "Call `load_cjs_rcpp()` before fitting models."
    )
  }

  cjs_nll_rcpp(
    par = as.numeric(par),
    histories = cjs_data$histories,
    X_phi = cjs_data$X_phi,
    first_capture = cjs_data$first_capture,
    frequencies = cjs_data$frequencies
  )
}


#' Penalised CJS objective
cjs_penalised_nll <- function(par, cjs_data, penalty_factor) {
  if (length(penalty_factor) != length(par)) {
    stop("`penalty_factor` must have the same length as `par`.")
  }
  if (any(penalty_factor < 0)) stop("Penalty factors must be non-negative.")

  cjs_nll(par, cjs_data) + sum(penalty_factor * abs(par))
}


#' Fit an unpenalised CJS model
fit_cjs <- function(
    cjs_data,
    initial = NULL,
    method = "BFGS",
    control = list(),
    hessian = FALSE) {

  n_parameters <- cjs_data$n_phi_parameters + 1L
  if (is.null(initial)) initial <- rep(0, n_parameters)

  fit <- stats::optim(
    par = initial,
    fn = cjs_nll,
    cjs_data = cjs_data,
    method = method,
    hessian = hessian,
    control = control
  )

  parameter_names <- c(colnames(cjs_data$X_phi), "p_(Intercept)")
  names(fit$par) <- parameter_names
  fit$logLik <- -fit$value
  fit$n_parameters <- n_parameters
  fit$AIC <- 2 * fit$value + 2 * n_parameters
  fit$BIC <- 2 * fit$value + log(cjs_data$effective_sample_size) * n_parameters
  class(fit) <- c("cjs_fit", class(fit))
  fit
}
