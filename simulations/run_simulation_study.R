# run simulations

options(error = NULL)

library(here)
library(Rcpp)
library(RcppArmadillo)

source(here("R", "simulate_cjs.R"))
source(here("R", "cjs_likelihood.R"))
source(here("R", "fit_adaptive_lasso.R"))
source(here("R", "summarise_simulations.R"))
source(here("R", "run_simulations.R"))
load_cjs_rcpp(here("R", "cjs_likelihood.cpp"))

# Simulation design
scenarios <- make_simulation_scenarios(
  n_occasions  = c(20),
  marked        = c(20, 50),
  mean_phi      = 0.8,
  p             = c(0.5, 0.9),
  effect_size   = c(-0.3, -0.6),
  n_covariates  = c(5, 10),
  n_active      = 2, 
  correlation   = c(0, 0.5, 0.8)
)

# Run simulations
results <- run_simulation_study(
  scenarios = scenarios,
  n_replicates = 10,
  seed = 20260720,
  lambda_min = 1e-3,
  lambda_max = 1e4,
  n_lambda = 120,
  zero_tolerance = 1e-3,
  save_file = "scenario_grid_partial.rds"
)

# Check completeness
nrow(results)
table(results$scenario_id)

# Check errors
attr(results, "errors")

# Convergence diagnostics
table(results$full_fit_convergence)
table(results$path_all_converged)
table(results$null_model_at_upper)

# Overall selection performance
colMeans(
  results[, c("exact_selection", "sensitivity", "specificity")],
  na.rm = TRUE
)


library(dplyr)

factors <- c(
  "marked",
  "p",
  "effect_size",
  "n_covariates",
  "correlation"
)

performance_by_factor <- lapply(factors, function(f) {
  results |>
    group_by(across(all_of(f))) |>
    summarise(
      exact_selection = mean(exact_selection),
      sensitivity = mean(sensitivity),
      specificity = mean(specificity),
      false_positives = mean(false_positives),
      false_negatives = mean(false_negatives),
      .groups = "drop"
    )
})
names(performance_by_factor) <- factors
performance_by_factor

# Save results
write.csv(results, here("simulations", "scenario_grid_results.csv"), row.names = FALSE)
summary_results <- summarise_simulation_study(results)
write.csv(summary_results, here("simulations", "scenario_grid_summary.csv"), row.names = FALSE)
covariate_summary <- summarise_covariate_selection(results)
write.csv(covariate_summary, here("simulations", "simulation_covariate_selection.csv"), row.names = FALSE)


