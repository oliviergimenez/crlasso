# Adaptive-Lasso fitting for the CJS likelihood.

#' Compute adaptive-Lasso weights
#'
#' @param full_fit Unpenalised fit from `fit_cjs()`.
#' @param unpenalised Integer indices of parameters that must not be penalised.
#' @param gamma Exponent in the adaptive weights, usually 1.
#' @param epsilon Lower bound preventing infinite weights.
adaptive_lasso_weights <- function(
    full_fit,
    unpenalised = c(1L, length(full_fit$par)),
    gamma = 1,
    epsilon = 1e-6) {

  weights <- 1 / pmax(abs(full_fit$par), epsilon)^gamma
  weights[unpenalised] <- 0
  weights
}


#' Fit an adaptive-Lasso regularisation path
#'
#' Information criteria are computed from the unpenalised likelihood evaluated
#' at the penalised estimates. Coefficients smaller than `zero_tolerance` are
#' treated as zero when counting effective parameters.
fit_adaptive_lasso <- function(
    cjs_data,
    lambda = NULL,
    lambda_min = 1e-3,
    lambda_max = 100,
    n_lambda = 150L,
    gamma = 1,
    unpenalised = c(1L, cjs_data$n_phi_parameters + 1L),
    zero_tolerance = 1e-3,
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-10)) {

  full_fit <- fit_cjs(cjs_data, method = method, control = control)
  weights <- adaptive_lasso_weights(
    full_fit,
    unpenalised = unpenalised,
    gamma = gamma
  )

  if (is.null(lambda)) {
    lambda <- exp(seq(log(lambda_min), log(lambda_max), length.out = n_lambda))
  }
  lambda <- sort(unique(as.numeric(lambda)))
  if (any(lambda <= 0)) stop("All values of `lambda` must be positive.")

  n_parameters <- length(full_fit$par)
  coefficients <- matrix(
    NA_real_, nrow = n_parameters, ncol = length(lambda),
    dimnames = list(names(full_fit$par), paste0("lambda_", seq_along(lambda)))
  )
  penalised_objective <- unpenalised_nll <- rep(NA_real_, length(lambda))
  convergence <- rep(NA_integer_, length(lambda))

  current <- full_fit$par
  for (j in seq_along(lambda)) {
    penalty_factor <- lambda[j] * weights
    fit <- stats::optim(
      par = current,
      fn = cjs_penalised_nll,
      cjs_data = cjs_data,
      penalty_factor = penalty_factor,
      method = method,
      control = control
    )

    current <- fit$par
    coefficients[, j] <- fit$par
    penalised_objective[j] <- fit$value
    unpenalised_nll[j] <- cjs_nll(fit$par, cjs_data)
    convergence[j] <- fit$convergence
  }

  penalised_indices <- setdiff(seq_len(n_parameters), unpenalised)
  n_selected <- colSums(abs(coefficients[penalised_indices, , drop = FALSE]) > zero_tolerance)
  n_effective <- length(unpenalised) + n_selected
  sample_size <- cjs_data$effective_sample_size

  AIC <- 2 * unpenalised_nll + 2 * n_effective
  BIC <- 2 * unpenalised_nll + log(sample_size) * n_effective

  out <- list(
    full_fit = full_fit,
    lambda = lambda,
    coefficients = coefficients,
    adaptive_weights = weights,
    penalised_objective = penalised_objective,
    unpenalised_nll = unpenalised_nll,
    n_selected = n_selected,
    n_effective = n_effective,
    AIC = AIC,
    BIC = BIC,
    convergence = convergence,
    unpenalised = unpenalised,
    zero_tolerance = zero_tolerance
  )
  class(out) <- "adaptive_lasso_cjs"
  out
}


#' Extract the model selected along a Lasso path
select_lasso_fit <- function(path_fit, criterion = c("BIC", "AIC")) {
  criterion <- match.arg(criterion)
  values <- path_fit[[criterion]]
  index <- which.min(values)

  list(
    criterion = criterion,
    index = index,
    lambda = path_fit$lambda[index],
    coefficients = path_fit$coefficients[, index],
    selected = abs(path_fit$coefficients[, index]) > path_fit$zero_tolerance,
    criterion_value = values[index],
    unpenalised_nll = path_fit$unpenalised_nll[index],
    convergence = path_fit$convergence[index]
  )
}
