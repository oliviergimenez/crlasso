# Analysis of adaptive-Lasso CJS simulations
#
# This version assumes that the number of capture occasions is fixed
# at 20 in all simulation scenarios. Consequently, n_occasions is not
# included as an explanatory factor in the GLMs or in the emmeans figures.

library(tidyverse)
library(emmeans)

bic_file <- "BIC/scenario_grid_results_BIC.csv"
aic_file <- "AIC/scenario_grid_results_AIC.csv"

read_results <- function(path, expected_criterion) {
  if (!file.exists(path)) stop("File not found: ", path)
  
  out <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required <- c(
    "scenario_id", "replicate", "n_occasions", "marked", "p",
    "effect_size", "n_covariates", "correlation", "exact_selection",
    "n_selected", "true_positives", "false_positives", "false_negatives",
    "true_negatives", "sensitivity", "specificity", "convergence"
  )
  
  absent <- setdiff(required, names(out))
  
  if (length(absent)) {
    stop(path, " is missing: ", paste(absent, collapse = ", "))
  }
  
  if ("criterion" %in% names(out)) {
    observed <- unique(toupper(as.character(out$criterion)))
    observed <- observed[!is.na(observed)]
    
    if (length(observed) && !all(observed == expected_criterion)) {
      warning(
        path,
        " contains criterion value(s): ",
        paste(observed, collapse = ", "),
        "; expected ",
        expected_criterion,
        "."
      )
    }
  }
  
  out <- out %>%
    filter(is.na(convergence) | convergence == 0L) %>%
    mutate(
      individuals = factor(marked, levels = c(20, 50)),
      detectability = factor(p, levels = c(0.5, 0.9)),
      effect = factor(effect_size, levels = c(-0.3, -0.6)),
      candidates = factor(n_covariates, levels = c(5, 10)),
      correlation_factor = factor(correlation, levels = c(0, 0.5, 0.8))
    )
  
  observed_occasions <- sort(unique(out$n_occasions))
  
  if (length(observed_occasions) != 1L || observed_occasions != 20) {
    warning(
      path,
      " contains the following n_occasions value(s): ",
      paste(observed_occasions, collapse = ", "),
      ". This script assumes a single value equal to 20."
    )
  }
  
  out
}

bic <- read_results(bic_file, "BIC")
aic <- read_results(aic_file, "AIC")

# n_occasions is omitted because it is fixed at 20.
factor_variables <- c(
  "individuals",
  "detectability",
  "effect",
  "candidates",
  "correlation_factor"
)

factor_titles <- c(
  individuals = "Individuals marked per occasion",
  detectability = "Detection probability",
  effect = "Effect size",
  candidates = "Candidate covariates",
  correlation_factor = "Covariate correlation"
)

model_formula <- paste(factor_variables, collapse = " + ")

# Robust extraction across emmeans versions. The grouping-column name can vary
# after summary/confint, so it is identified as the first column.
extract_emmeans <- function(model, variable, title, metric = NA_character_) {
  emm <- emmeans::emmeans(
    model,
    specs = stats::as.formula(paste("~", variable))
  )
  
  out <- as.data.frame(
    summary(emm, type = "response", infer = c(TRUE, FALSE))
  )
  
  level_values <- as.character(out[[1L]])
  
  estimate_column <- intersect(
    c("prob", "response", "rate", "emmean"),
    names(out)
  )[1L]
  
  lower_column <- intersect(
    c("asymp.LCL", "lower.CL"),
    names(out)
  )[1L]
  
  upper_column <- intersect(
    c("asymp.UCL", "upper.CL"),
    names(out)
  )[1L]
  
  required_columns <- c(
    estimate_column,
    lower_column,
    upper_column
  )
  
  if (anyNA(required_columns)) {
    stop(
      "Could not identify the response-scale estimate or confidence limits for `",
      variable,
      "`. emmeans returned: ",
      paste(names(out), collapse = ", ")
    )
  }
  
  if (length(level_values) != nrow(out)) {
    stop("Unexpected emmeans output for `", variable, "`.")
  }
  
  data.frame(
    variable = rep(variable, nrow(out)),
    factor = rep(title, nrow(out)),
    level = level_values,
    level_order = seq_len(nrow(out)),
    metric = rep(metric, nrow(out)),
    estimate = out[[estimate_column]],
    lower = out[[lower_column]],
    upper = out[[upper_column]],
    stringsAsFactors = FALSE
  )
}

prepare_plot_axis <- function(data) {
  factor_order <- unname(factor_titles)
  
  data <- data %>%
    mutate(
      factor = factor(factor, levels = factor_order),
      level_key = paste(variable, level_order, sep = "__")
    )
  
  key_order <- data %>%
    distinct(factor, level_order, level_key) %>%
    arrange(factor, level_order) %>%
    pull(level_key)
  
  key_labels <- data %>%
    distinct(level_key, level) %>%
    { setNames(.$level, .$level_key) }
  
  data$level_key <- factor(
    data$level_key,
    levels = key_order
  )
  
  attr(data, "key_labels") <- key_labels
  data
}

# Figure 1: exact model selection
exact_model <- glm(
  as.formula(paste("exact_selection ~", model_formula)),
  family = binomial(link = "logit"),
  data = bic
)

exact_effects <- bind_rows(
  lapply(factor_variables, function(x) {
    extract_emmeans(
      exact_model,
      x,
      factor_titles[[x]]
    )
  })
) %>%
  prepare_plot_axis()

exact_labels <- attr(exact_effects, "key_labels")

figure_exact <- ggplot(
  exact_effects,
  aes(
    x = level_key,
    y = estimate,
    group = 1
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.10,
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ factor,
    scales = "free_x",
    nrow = 2
  ) +
  scale_x_discrete(labels = exact_labels) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Probability of exact selection"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(figure_exact)

ggsave(
  filename = "figure_exact_selection_adjusted.png",
  plot = figure_exact,
  width = 9,
  height = 6,
  dpi = 400
)

# Figure 2: sensitivity and specificity
sensitivity_model <- glm(
  as.formula(
    paste(
      "cbind(true_positives, false_negatives) ~",
      model_formula
    )
  ),
  family = quasibinomial(link = "logit"),
  data = bic
)

specificity_model <- glm(
  as.formula(
    paste(
      "cbind(true_negatives, false_positives) ~",
      model_formula
    )
  ),
  family = quasibinomial(link = "logit"),
  data = bic
)

performance_effects <- bind_rows(
  lapply(factor_variables, function(x) {
    bind_rows(
      extract_emmeans(
        sensitivity_model,
        x,
        factor_titles[[x]],
        "Sensitivity"
      ),
      extract_emmeans(
        specificity_model,
        x,
        factor_titles[[x]],
        "Specificity"
      )
    )
  })
) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Sensitivity", "Specificity")
    )
  ) %>%
  prepare_plot_axis()

performance_labels <- attr(
  performance_effects,
  "key_labels"
)

pd <- position_dodge(width = 0.16)

figure_performance <- ggplot(
  performance_effects,
  aes(
    x = level_key,
    y = estimate,
    group = metric,
    linetype = metric,
    shape = metric
  )
) +
  geom_line(
    linewidth = 0.55,
    position = pd
  ) +
  geom_point(
    size = 2.3,
    position = pd
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.08,
    linewidth = 0.4,
    position = pd
  ) +
  facet_wrap(
    ~ factor,
    scales = "free_x",
    nrow = 2
  ) +
  scale_x_discrete(labels = performance_labels) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Selection performance",
    linetype = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(figure_performance)

ggsave(
  filename = "figure_sensitivity_specificity_adjusted.png",
  plot = figure_performance,
  width = 9,
  height = 6,
  dpi = 400
)

# Overall AIC/BIC comparison table
criterion_comparison <- bind_rows(
  BIC = bic,
  AIC = aic,
  .id = "Criterion"
) %>%
  group_by(Criterion) %>%
  summarise(
    n_simulations = n(),
    exact_selection_probability = mean(
      exact_selection,
      na.rm = TRUE
    ),
    sensitivity = mean(
      sensitivity,
      na.rm = TRUE
    ),
    specificity = mean(
      specificity,
      na.rm = TRUE
    ),
    mean_selected_model_size = mean(
      n_selected,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Criterion = factor(
      Criterion,
      levels = c("BIC", "AIC")
    )
  ) %>%
  arrange(Criterion)

print(criterion_comparison)

# write.csv(
#   criterion_comparison,
#   "table_AIC_BIC_comparison.csv",
#   row.names = FALSE
# )
