# White stork case study:
# adaptive-Lasso covariate selection in a CJS model

rm(list = ls())

# -------------------------------------------------------------------------
# Load functions
# -------------------------------------------------------------------------

source("R/cjs_likelihood.R")
source("R/fit_adaptive_lasso.R")
source("R/plotting.R")

# -------------------------------------------------------------------------
# Load and prepare data
# -------------------------------------------------------------------------

load("dat/whitestork.RData")

# The object `data` contains:
# - one row per unique encounter history;
# - capture occasions in all columns except the last;
# - frequencies in the last column.
#
# The object `cov` contains one row per survival interval and one column
# per climatic covariate.

encounter_histories <- as.matrix(data[, -ncol(data), drop = FALSE])
frequencies <- data[, ncol(data)]

covariates <- as.matrix(cov)

if (is.null(colnames(covariates))) {
  colnames(covariates) <- paste0("station_", seq_len(ncol(covariates)))
}

# Survival design matrix: intercept and climatic covariates
X_phi <- cbind(
  phi_intercept = 1,
  covariates
)

# Check compatibility between encounter histories and covariates
if (nrow(X_phi) != ncol(encounter_histories) - 1L) {
  stop(
    "The number of rows in `cov` must equal the number of ",
    "capture occasions minus one."
  )
}

cjs_data <- prepare_cjs_data(
  histories = encounter_histories,
  X_phi = X_phi,
  frequencies = frequencies
)

# -------------------------------------------------------------------------
# Fit the full, unpenalised model
# -------------------------------------------------------------------------

full_fit <- fit_cjs(
  cjs_data = cjs_data,
  control = list(
    maxit = 5000,
    reltol = 1e-10
  )
)

if (full_fit$convergence != 0L) {
  warning(
    "The full model did not converge. optim() convergence code: ",
    full_fit$convergence
  )
}

full_model_results <- data.frame(
  parameter = names(full_fit$par),
  estimate = unname(full_fit$par),
  stringsAsFactors = FALSE
)

print(full_model_results)

# -------------------------------------------------------------------------
# Fit the adaptive-Lasso path
# -------------------------------------------------------------------------

# The survival intercept and the detection intercept are not penalised.
unpenalised_parameters <- c(
  1L,
  ncol(X_phi) + 1L
)

# lambda_grid <- exp(
#   seq(
#     log(1e-3),
#     log(2),
#     length.out = 100
#   )
# )

lambda_grid <- exp(
  seq(
    log(1e-3),
    log(100),
    length.out = 150
  )
)

# lasso_path <- fit_adaptive_lasso(
#   cjs_data = cjs_data,
#   lambda = lambda_grid,
#   gamma = 1,
#   unpenalised = unpenalised_parameters,
#   zero_tolerance = 1e-4,
#   control = list(
#     maxit = 5000,
#     reltol = 1e-10
#   )
# )

lasso_path <- fit_adaptive_lasso(
  cjs_data = cjs_data,
  lambda = lambda_grid,
  gamma = 1,
  unpenalised = c(1L, ncol(X_phi) + 1L),
  zero_tolerance = 1e-3,
  control = list(
    maxit = 5000,
    reltol = 1e-12
  )
)

if (any(lasso_path$convergence != 0L)) {
  warning(
    sum(lasso_path$convergence != 0L),
    " fits along the regularisation path did not converge."
  )
}

# -------------------------------------------------------------------------
# Select lambda using BIC
# -------------------------------------------------------------------------

selected_fit <- select_lasso_fit(
  path_fit = lasso_path,
  criterion = "BIC"
)

cat(
  "\nSelected lambda:",
  selected_fit$lambda,
  "\nSelected BIC:",
  selected_fit$criterion_value,
  "\n"
)

# -------------------------------------------------------------------------
# Display selected coefficients
# -------------------------------------------------------------------------

coefficient_table <- data.frame(
  parameter = names(selected_fit$coefficients),
  estimate = unname(selected_fit$coefficients),
  penalised = !seq_along(selected_fit$coefficients) %in%
    unpenalised_parameters,
  selected = unname(selected_fit$selected),
  stringsAsFactors = FALSE
)

# Intercepts are always retained
coefficient_table$selected[unpenalised_parameters] <- TRUE

print(coefficient_table)

selected_covariates <- subset(
  coefficient_table,
  penalised & selected
)

cat("\nSelected climatic covariates:\n")

if (nrow(selected_covariates) == 0L) {
  cat("None\n")
} else {
  print(
    selected_covariates[, c("parameter", "estimate")],
    row.names = FALSE
  )
}

# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------

plot_lasso_criterion(
  path_fit = lasso_path,
  criterion = "BIC"
)

plot_lasso_path(
  path_fit = lasso_path,
  criterion = "BIC",
  include_unpenalised = FALSE,
  legend_position = "right"
)

ggsave(filename = "path_stork.png", dpi = 600)

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

saveRDS(
  list(
    cjs_data = cjs_data,
    full_fit = full_fit,
    lasso_path = lasso_path,
    selected_fit = selected_fit,
    coefficient_table = coefficient_table
  ),
  file = "whitestork_lasso_results.rds"
)

utils::write.csv(
  coefficient_table,
  file = "whitestork_lasso_coefficients.csv",
  row.names = FALSE
)
