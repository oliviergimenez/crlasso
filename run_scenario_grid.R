# factorial simulation design

library(Rcpp)
library(RcppArmadillo)

source("R/simulate_cjs.R")
source("R/cjs_likelihood_rcpp.R")
source("R/fit_adaptive_lasso_rcpp.R")
source("R/summarise_simulations.R")
source("R/run_simulations.R")

load_cjs_rcpp("R/cjs_likelihood.cpp")

scenarios <- make_simulation_scenarios(
  n_occasions  = c(10, 20),
  marked        = c(20, 50),
  mean_phi      = 0.8,
  p             = c(0.5, 0.9),
  effect_size   = c(-0.3, -0.6),
  n_covariates  = c(10, 20),
  n_active      = 2, # Que se passe-t-il lorsque le nombre de variables candidates augmente 
  # alors que seulement deux sont réellement informatives ?
  correlation   = c(0, 0.5, 0.8)
)

results <- run_simulation_study(
  scenarios = scenarios,
  n_replicates = 100,
  seed = 20260720,
  lambda_min = 1e-3,
  lambda_max = 100,
  n_lambda = 150,
  zero_tolerance = 1e-3,
  save_file = "scenario_grid_partial.rds"
)

saveRDS(results, "scenario_grid_results.rds")
utils::write.csv(results, "scenario_grid_results.csv", row.names = FALSE)

summary_results <- summarise_simulation_study(results)
utils::write.csv(summary_results, "scenario_grid_summary.csv", row.names = FALSE)
print(summary_results)

covariate_summary <- summarise_covariate_selection(results)
utils::write.csv(covariate_summary, "simulation_covariate_selection.csv", row.names = FALSE)

#plot_selection_performance(summary_results)
