# Adaptive Lasso for capture-recapture models

This repository contains the `R` code used to implement and evaluate adaptive-Lasso covariate selection in Cormack-Jolly-Seber (CJS) capture-recapture models. It includes a simulation study assessing covariate-selection performance and a case study investigating the relationship between environmental covariates and survival in white storks (*Ciconia ciconia*).

## Repository structure

The repository is organised into three main directories:

- `R/`: core functions for data simulation, likelihood evaluation, adaptive-Lasso fitting, model selection and visualisation.
- `simulations/`: scripts to run the simulation study and analyse its results.
- `case_study/`: code and data for the white stork case study.

## Getting started

### Simulation study

To reproduce the simulation study, run:

`simulations/run_simulation_study.R`

This script defines the simulation scenarios, generates the capture-recapture datasets, fits the models and saves the results.

To summarise and visualise the simulation results, run:

`simulations/analyse_simulation_outputs.R`

### Case study

To reproduce the white stork analysis, run:

`case_study/run_casestudy.R`

This script loads the capture-recapture and environmental data, fits the adaptive Lasso, selects the model using BIC and produces the corresponding results and figures.

## Core R functions

The `R/` directory contains the functions used by both the simulation study and the case study:

- `simulate_cjs.R`: simulation of capture histories and temporal covariates.
- `cjs_likelihood.R`: CJS likelihood and unpenalised fitting.
- `fit_adaptive_lasso.R`: adaptive weights, regularisation path, unpenalised refitting and model selection.
- `run_simulations.R`: simulation engine for individual replicates and complete simulation studies.
- `summarise_simulations.R`: replicate-level and scenario-level summaries.
- `plotting.R`: visualisation of regularisation paths.
