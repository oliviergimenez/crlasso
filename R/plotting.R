library(ggplot2)

#' Plot an information criterion along the regularisation path
plot_lasso_criterion <- function(path_fit,
                                 criterion = c("BIC", "AIC")) {
  
  criterion <- match.arg(criterion)
  
  df <- data.frame(
    lambda = path_fit$lambda,
    value = path_fit[[criterion]]
  )
  
  best <- which.min(df$value)
  
  p <- ggplot(df, aes(x = lambda, y = value)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(
      xintercept = df$lambda[best],
      linetype = "dashed",
      colour = "red"
    ) +
    scale_x_log10() +
    labs(
      x = expression(lambda),
      y = criterion
    ) +
    theme_classic()
  
  print(p)
  
  invisible(best)
}


#' Plot coefficient paths
plot_lasso_path <- function(
    path_fit,
    criterion = c("BIC", "AIC"),
    include_unpenalised = FALSE,
    legend_position = "right") {
  
  criterion <- match.arg(criterion)
  
  rows <- seq_len(nrow(path_fit$coefficients))
  if (!include_unpenalised)
    rows <- setdiff(rows, path_fit$unpenalised)
  
  path <- path_fit$coefficients[rows, , drop = FALSE]
  
  df <- data.frame(
    lambda = rep(path_fit$lambda, each = nrow(path)),
    parameter = factor(
      rep(rownames(path), times = ncol(path)),
      levels = rownames(path)
    ),
    estimate = as.vector(path)
  )
  
  best <- which.min(path_fit[[criterion]])
  
  p <- ggplot(df,
              aes(x = lambda,
                  y = estimate,
                  colour = parameter)) +
    geom_line(linewidth = 0.9) +
    geom_hline(
      yintercept = 0,
      colour = "grey70"
    ) +
    geom_vline(
      xintercept = path_fit$lambda[best],
      linetype = "dashed",
      colour = "red"
    ) +
    scale_x_log10() +
    labs(
      x = expression(lambda),
      y = "Coefficient estimate",
      colour = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = legend_position
    )
  
  print(p)
  
  invisible(best)
}


#' Plot exact-selection probability across scenarios
plot_selection_performance <- function(
    summary,
    x = "scenario_id",
    metric = "exact_selection_probability",
    methods = unique(summary$method),
    legend_position = "none") {
  
  summary <- summary[
    summary$method %in% methods,
    ,
    drop = FALSE
  ]
  
  summary$method <- factor(
    summary$method,
    levels = c(
      "adaptive_lasso_AIC",
      "adaptive_lasso_BIC",
      "exhaustive_AIC",
      "exhaustive_BIC"
    ),
    labels = c(
      "Adaptive Lasso\n(AIC)",
      "Adaptive Lasso\n(BIC)",
      "Exhaustive search\n(AIC)",
      "Exhaustive search\n(BIC)"
    )
  )
  
  p <- ggplot(
    summary,
    aes(
      x = .data[[x]],
      y = .data[[metric]],
      colour = method,
      group = 1
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~method, ncol = 2) +
    labs(
      x = NULL,
      y = "Exact selection probability"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  print(p)
  
  invisible(summary)
}
