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

  # penalised_indices <- setdiff(seq_len(n_parameters), unpenalised)
  # n_selected <- colSums(abs(coefficients[penalised_indices, , drop = FALSE]) > zero_tolerance)
  # n_effective <- length(unpenalised) + n_selected
  # sample_size <- cjs_data$effective_sample_size
  # 
  # AIC <- 2 * unpenalised_nll + 2 * n_effective
  # BIC <- 2 * unpenalised_nll + log(sample_size) * n_effective

  penalised_indices <- setdiff(
    seq_len(n_parameters),
    unpenalised
  )
  
  n_selected <- colSums(
    abs(coefficients[penalised_indices, , drop = FALSE]) >
      zero_tolerance
  )
  
  n_effective <- length(unpenalised) + n_selected
  sample_size <- cjs_data$effective_sample_size
  
  # Information criteria evaluated directly at penalised estimates
  AIC_penalised <- 2 * unpenalised_nll + 2 * n_effective
  BIC_penalised <- 2 * unpenalised_nll +
    log(sample_size) * n_effective
  
  
  # ----------------------------------------------------------------------
  # Post-Lasso refitting
  # ----------------------------------------------------------------------
  
  refit_nll <- rep(NA_real_, length(lambda))
  refit_AIC <- rep(NA_real_, length(lambda))
  refit_BIC <- rep(NA_real_, length(lambda))
  refit_convergence <- rep(NA_integer_, length(lambda))
  
  refit_coefficients <- vector("list", length(lambda))
  supports <- vector("list", length(lambda))
  
  for (j in seq_along(lambda)) {
    
    # Penalised survival parameters retained at this lambda
    selected_penalised <- penalised_indices[
      abs(coefficients[penalised_indices, j]) >
        zero_tolerance
    ]
    
    # Always include unpenalised parameters
    selected_indices <- sort(unique(c(
      unpenalised,
      selected_penalised
    )))
    
    supports[[j]] <- selected_indices
    
    refit <- refit_cjs_support(
      cjs_data = cjs_data,
      selected_indices = selected_indices,
      initial = coefficients[, j],
      method = method,
      control = control
    )
    
    refit_nll[j] <- refit$value
    refit_convergence[j] <- refit$convergence
    refit_coefficients[[j]] <- refit$par
    
    k <- length(refit$par)
    
    refit_AIC[j] <-
      2 * refit$value +
      2 * k
    
    refit_BIC[j] <-
      2 * refit$value +
      log(sample_size) * k
  }
  
  out <- list(
    full_fit = full_fit,
    lambda = lambda,
    coefficients = coefficients,
    adaptive_weights = weights,
    penalised_objective = penalised_objective,
    unpenalised_nll = unpenalised_nll,
    n_selected = n_selected,
    n_effective = n_effective,
    
    AIC_penalised = AIC_penalised,
    BIC_penalised = BIC_penalised,
    
    AIC = refit_AIC,
    BIC = refit_BIC,
    refit_nll = refit_nll,
    refit_coefficients = refit_coefficients,
    refit_convergence = refit_convergence,
    supports = supports,
    
    convergence = convergence,
    unpenalised = unpenalised,
    zero_tolerance = zero_tolerance
  )
  
  class(out) <- "adaptive_lasso_cjs"
  
  out
}


#' Extract the model selected along a Lasso path
select_lasso_fit <- function(
    path_fit,
    criterion = c("BIC", "AIC")) {
  
  criterion <- match.arg(criterion)
  
  values <- path_fit[[criterion]]
  index <- which.min(values)
  
  support <- path_fit$supports[[index]]
  refit_coef <- path_fit$refit_coefficients[[index]]
  
  # Full-length coefficient vector
  coefficients <- rep(
    0,
    nrow(path_fit$coefficients)
  )
  
  names(coefficients) <-
    rownames(path_fit$coefficients)
  
  coefficients[support] <- refit_coef
  
  selected <- rep(
    FALSE,
    length(coefficients)
  )
  
  names(selected) <- names(coefficients)
  selected[support] <- TRUE
  
  list(
    criterion = criterion,
    index = index,
    lambda = path_fit$lambda[index],
    
    # Estimates along the penalised path
    penalised_coefficients =
      path_fit$coefficients[, index],
    
    # Unpenalised post-Lasso estimates,
    # expanded back to the full parameter vector
    coefficients = coefficients,
    
    support = support,
    selected = selected,
    
    criterion_value = values[index],
    unpenalised_nll = path_fit$refit_nll[index],
    convergence = path_fit$refit_convergence[index]
  )
}

#' Refit an unpenalised CJS model on a selected Lasso support
#'
#' @param cjs_data Original prepared CJS data.
#' @param selected_indices Indices of selected survival parameters.
#' @param detection_index Index of the detection intercept.
#' @param initial Optional initial values for the reduced model.
#' @param method Optimisation method.
#' @param control Control arguments passed to optim().
#' @return An unpenalised CJS fit for the selected model.
refit_cjs_support <- function(
    cjs_data,
    selected_indices,
    detection_index = cjs_data$n_phi_parameters + 1L,
    initial = NULL,
    method = "BFGS",
    control = list()) {
  
  # Survival parameters retained in the reduced model
  survival_indices <- selected_indices[
    selected_indices <= cjs_data$n_phi_parameters
  ]
  
  # Always retain the survival intercept
  survival_indices <- sort(unique(c(1L, survival_indices)))
  
  # Reduced survival design matrix
  X_phi_reduced <- cjs_data$X_phi[
    ,
    survival_indices,
    drop = FALSE
  ]
  
  reduced_data <- prepare_cjs_data(
    histories = cjs_data$histories,
    X_phi = X_phi_reduced,
    frequencies = cjs_data$frequencies
  )
  
  # Initial values: survival parameters + detection intercept
  if (!is.null(initial)) {
    initial <- c(
      initial[survival_indices],
      initial[detection_index]
    )
  }
  
  fit_cjs(
    cjs_data = reduced_data,
    initial = initial,
    method = method,
    control = control
  )
}