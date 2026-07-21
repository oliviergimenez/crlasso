# Likelihood functions for an intercept-only detection CJS model.

#' Validate and prepare CJS data
#'
#' @param histories Binary capture-history matrix, individuals in rows.
#' @param X_phi Survival design matrix, one row per interval.
#' @param frequencies Optional frequencies for aggregated capture histories.
#' @return A validated list used by the fitting functions.
prepare_cjs_data <- function(histories, X_phi, frequencies = NULL) {
  histories <- as.matrix(histories)
  X_phi <- as.matrix(X_phi)

  if (!all(histories %in% c(0, 1))) {
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
    first_capture = first_capture,
    n_histories = nrow(histories),
    n_occasions = ncol(histories),
    n_intervals = ncol(histories) - 1L,
    n_phi_parameters = ncol(X_phi),
    effective_sample_size = sum(frequencies)
  )
}


#' Negative conditional log-likelihood for a CJS model
#'
#' Survival is modelled through `X_phi`; detection is constant. The likelihood
#' is conditional on first capture.
#'
#' @param par Parameter vector: survival coefficients followed by logit(p).
#' @param cjs_data Output from `prepare_cjs_data()`.
#' @return Negative log-likelihood.
cjs_nll <- function(par, cjs_data) {
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

    # Conditional on first capture, the individual is alive at `first`.
    alpha <- c(alive = 1, dead = 0)

    for (interval in seq.int(first, cjs_data$n_intervals)) {
      transition <- matrix(
        c(phi[interval], 1 - phi[interval], 0, 1),
        nrow = 2L,
        byrow = TRUE
      )
      alpha <- drop(alpha %*% transition)

      observation <- cjs_data$histories[i, interval + 1L]
      emission <- if (observation == 1L) c(p, 0) else c(1 - p, 1)
      alpha <- alpha * emission
    }

    history_probability <- sum(alpha)
    if (!is.finite(history_probability) || history_probability <= 0) {
      return(Inf)
    }
    total_log_likelihood <- total_log_likelihood +
      cjs_data$frequencies[i] * log(history_probability)
  }

  -total_log_likelihood
}


#' Penalised CJS objective
#'
#' @param par Model parameters.
#' @param cjs_data Prepared CJS data.
#' @param penalty_factor Non-negative parameter-specific penalty factors.
#' @return Negative penalised log-likelihood.
cjs_penalised_nll <- function(par, cjs_data, penalty_factor) {
  if (length(penalty_factor) != length(par)) {
    stop("`penalty_factor` must have the same length as `par`.")
  }
  if (any(penalty_factor < 0)) stop("Penalty factors must be non-negative.")

  cjs_nll(par, cjs_data) + sum(penalty_factor * abs(par))
}


#' Fit an unpenalised CJS model
fit_cjs <- function(cjs_data, initial = NULL, method = "BFGS", control = list()) {
  n_parameters <- cjs_data$n_phi_parameters + 1L
  if (is.null(initial)) initial <- rep(0, n_parameters)

  fit <- stats::optim(
    par = initial,
    fn = cjs_nll,
    cjs_data = cjs_data,
    method = method,
    hessian = TRUE,
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
