# Adaptive Lasso for CJS models

This folder contains the simulations and case-study code for adaptive-Lasso 
covariate selection in Cormack--Jolly--Seber models.

## Files

- `R/simulate_cjs.R`: simulation of capture histories and temporal covariates.
- `R/cjs_likelihood.R`: conditional CJS likelihood and unpenalised fitting.
- `R/fit_adaptive_lasso.R`: adaptive weights, regularisation path and selection.
- `R/fit_candidate_models.R`: exhaustive candidate-model fitting by AIC/BIC.
- `R/summarise_simulations.R`: replicate-level and scenario-level summaries.
- `R/run_simulations.R`: one-replicate and complete simulation-study engine.
- `R/plotting.R`: diagnostic and coefficient-path plots.
- `run_scenario_grid.R`: example factorial scenario grid.
- `run_casestudy.R`: case study.

## Run

From the project directory:

```r
source("run_example.R")
```

The current implementation uses only base R and `stats`.

## Simulation workflow

Run the original scenario with:

```r
source("run_simulations.R")
```

The simulation output has one row per scenario, replicate and method. The four
default methods are exhaustive AIC, exhaustive BIC, adaptive-Lasso AIC and
adaptive-Lasso BIC. Raw output, summaries and marginal covariate-selection
probabilities are saved separately.

`run_simulation_study()` generates and stores a separate random seed for every
replicate, records convergence and timing, and can save an incremental RDS file
after every replicate. Failed replicates are recorded rather than silently
discarded.
