# Functions to summarise simulation results.

#' Identify selected covariates from a coefficient vector
selected_covariates <- function(
    coefficients,
    covariate_names,
    zero_tolerance = 1e-4) {
  coefficients <- coefficients[covariate_names]
  names(coefficients)[abs(coefficients) > zero_tolerance]
}


#' Combine replicate-level simulation outputs
bind_simulation_results <- function(results) {
  if (!is.list(results) || length(results) == 0L) {
    stop("`results` must be a non-empty list.")
  }
  all_names <- unique(unlist(lapply(results, names), use.names = FALSE))
  aligned <- lapply(results, function(x) {
    missing <- setdiff(all_names, names(x))
    for (name in missing) x[[name]] <- NA
    x[all_names]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}


#' Summarise variable-selection performance by scenario and method
#'
#' @param results Output from `run_simulation_study()`.
#' @param group_columns Columns defining groups in the summary.
#' @return A data frame of selection and convergence metrics.
summarise_simulation_study <- function(
    results,
    group_columns = c("scenario_id", "method")) {
  
  results <- as.data.frame(results)
  
  if (!all(group_columns %in% names(results))) {
    stop("One or more grouping columns are absent from `results`.")
  }
  
  group_id <- interaction(
    results[group_columns],
    drop = TRUE,
    lex.order = TRUE
  )
  
  groups <- split(results, group_id)
  
  summaries <- lapply(groups, function(x) {
    
    keys <- x[1L, group_columns, drop = FALSE]
    
    converged <- is.na(x$convergence) | x$convergence == 0L
    
    data.frame(
      keys,
      n_replicates = nrow(x),
      n_converged = sum(converged),
      convergence_rate = mean(converged),
      exact_selection_probability =
        mean(x$exact_selection[converged], na.rm = TRUE),
      mean_model_size =
        mean(x$n_selected[converged], na.rm = TRUE),
      mean_true_positives =
        mean(x$true_positives[converged], na.rm = TRUE),
      mean_false_positives =
        mean(x$false_positives[converged], na.rm = TRUE),
      sensitivity =
        mean(x$sensitivity[converged], na.rm = TRUE),
      specificity =
        mean(x$specificity[converged], na.rm = TRUE),
      median_elapsed_seconds =
        stats::median(x$elapsed_seconds[converged], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  
  ## Robust ordering
  ord <- do.call(order, unname(as.list(out[group_columns])))
  out <- out[ord, , drop = FALSE]
  
  out
}

#' Summarise marginal selection probabilities for each covariate
summarise_covariate_selection <- function(
    results,
    group_columns = c("scenario_id", "method")) {

  selection_columns <- grep("^selected_x[0-9]+$", names(results), value = TRUE)
  if (!length(selection_columns)) stop("No `selected_x*` columns were found.")

  long <- do.call(rbind, lapply(selection_columns, function(column) {
    data.frame(
      results[group_columns],
      covariate = sub("^selected_", "", column),
      selected = results[[column]],
      stringsAsFactors = FALSE
    )
  }))

  group_id <- interaction(long[c(group_columns, "covariate")], drop = TRUE, lex.order = TRUE)
  groups <- split(long, group_id)
  out <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      x[1L, c(group_columns, "covariate"), drop = FALSE],
      selection_probability = mean(x$selected, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}


#' Monte Carlo standard error for a probability estimate
mcse_probability <- function(probability, n) {
  sqrt(probability * (1 - probability) / n)
}
