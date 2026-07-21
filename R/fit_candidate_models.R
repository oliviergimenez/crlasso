# Exhaustive fitting of candidate CJS models using the same likelihood.

#' Fit every subset of survival covariates
#'
#' The first column of `X_phi` is assumed to be the survival intercept and is
#' retained in every model. Detection remains intercept-only.
fit_candidate_models <- function(
    histories,
    X_phi,
    frequencies = NULL,
    criterion = c("AIC", "BIC"),
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-10)) {

  criterion <- match.arg(criterion)
  X_phi <- as.matrix(X_phi)
  n_covariates <- ncol(X_phi) - 1L

  subsets <- if (n_covariates == 0L) {
    list(integer(0))
  } else {
    unlist(
      lapply(0:n_covariates, function(size) {
        combn(seq_len(n_covariates), size, simplify = FALSE)
      }),
      recursive = FALSE
    )
  }

  results <- vector("list", length(subsets))
  for (m in seq_along(subsets)) {
    selected_covariates <- subsets[[m]]
    columns <- c(1L, selected_covariates + 1L)
    model_data <- prepare_cjs_data(
      histories = histories,
      X_phi = X_phi[, columns, drop = FALSE],
      frequencies = frequencies
    )
    fit <- fit_cjs(model_data, method = method, control = control)

    results[[m]] <- data.frame(
      model = if (length(selected_covariates) == 0L) {
        "intercept"
      } else {
        paste(colnames(X_phi)[selected_covariates + 1L], collapse = " + ")
      },
      n_covariates = length(selected_covariates),
      logLik = fit$logLik,
      AIC = fit$AIC,
      BIC = fit$BIC,
      convergence = fit$convergence,
      stringsAsFactors = FALSE
    )
    results[[m]]$fit <- list(fit)
    results[[m]]$selected_columns <- list(selected_covariates + 1L)
  }

  table <- do.call(rbind, results)
  rownames(table) <- NULL
  table$delta_AIC <- table$AIC - min(table$AIC)
  table$delta_BIC <- table$BIC - min(table$BIC)
  table$selected <- FALSE
  table$selected[which.min(table[[criterion]])] <- TRUE
  table[order(table[[criterion]]), ]
}
