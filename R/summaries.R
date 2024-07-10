# Summarise results of simulations:
# - (relative) bias
# - (relative) empirical standard error
# - (relative) model-based standard error
# - confidence interval coverage
summarise_sim_reps <- function(sim_reps) {
  sim_reps %>%
    mutate(across(c(scenario_name, treatment_effect, compliance, dose_response,
                    sample_size, missingness_mechanism, outcome_missingness,
                    missingness_name, estimator_name),
                  fct_inorder)) %>%
    drop_na(error, std.error) %>%
    mutate(
      rel_error = if_else(true_effect < 1e-6, error, error / true_effect),
      rel_std_error = if_else(true_effect < 1e-6, std.error, std.error / true_effect)
    ) %>%
    group_by(
      scenario_name, 
      treatment_effect, compliance, dose_response,
      sample_size, missingness_mechanism, outcome_missingness,
      missingness_name, estimator_name
    ) %>%
    summarise(
      n = n(),
      bias = mean(error),
      bias_mcse = sd(error) / sqrt(n()),
      se_empirical = sd(error),
      se_empirical_mcse = se_empirical / sqrt(2 * (n() - 1)),
      se_model = sqrt(mean(std.error^2)),
      se_model_mcse = sqrt(sd(std.error^2) / n()),
      rmse = sqrt(mean(error^2)),
      rmse_mcse = sqrt(sd(error^2) / n()),
      rel_bias = mean(rel_error),
      rel_bias_mcse = sd(rel_error) / sqrt(n()),
      rel_se_empirical = sd(rel_error),
      rel_se_empirical_mcse = rel_se_empirical / sqrt(2 * (n() - 1)),
      rel_se_model = sqrt(mean(rel_std_error^2)),
      rel_se_model_mcse = sqrt(sd(rel_std_error^2) / n()),
      rel_rmse = sqrt(mean(rel_error^2)),
      rel_rmse_mcse = sqrt(sd(rel_error^2) / n()),
      ci_coverage = mean(ci_includes_truth),
      ci_coverage_mcse = sqrt(ci_coverage * (1 - ci_coverage) / n()),
    ) %>%
    ungroup() %>%
    mutate(
      missingness_name = 
        fct_relabel(missingness_name, \(x) str_replace_all(x, "_", " ")),
      treatment_effect = 
        fct_relabel(treatment_effect, \(x) str_replace_all(x, "_", " ")),
      treatment_effect = 
        fct_relabel(treatment_effect, \(x) str_replace_all(x, "compliance", "compl")),
      missingness_mechanism = 
        fct_relabel(missingness_mechanism, \(x) str_replace_all(x, "_", " ")),
    )
}

# Giant results table
full_results_table <- function(sim_reps_summary) {
  stop("FIXME")
}

mcse_summary_df <- function(sim_reps_summary) {
  sim_reps_summary %>%
    rename_with(
      \(x) if_else(str_ends(x, "_mcse"), x, glue("{x}_estimate")),
      .cols = bias:ci_coverage_mcse
    ) %>%
    pivot_longer(
      cols = bias_estimate:ci_coverage_mcse,
      names_pattern = "^([a-z_]+)_(estimate|mcse)$",
      names_to = c("parameter", ".value"),
    ) %>%
    mutate(parameter = fct_inorder(parameter)) %>%
    group_by(parameter) %>%
    summarise(
      min_mcse = min(mcse, na.rm = TRUE),
      q1_mcse = quantile(mcse, 0.25, na.rm = TRUE),
      med_mcse = median(mcse, na.rm = TRUE),
      q3_mcse = quantile(mcse, 0.75, na.rm = TRUE),
      max_mcse = max(mcse, na.rm = TRUE)
    )
}

# Plot theme - maybe this will be moved into a cpmisc package one day
theme_cp <- function(panel_border = TRUE, grid = "none") {
  t <- theme_cowplot(font_size = 11, rel_small = 9/11,
                     rel_tiny = 8/11, rel_large = 13/11) +
    theme(plot.background = element_rect(fill = "white")) +
    background_grid(major = grid)
  if (panel_border) {
    t <- t + panel_border()
  }
}

# Plots for scenarios with no missing data
make_plot_results_no_missing <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(missingness_mechanism == "none") %>%
    pivot_longer(c(rel_bias, rel_rmse, rel_se_empirical, ci_coverage),
                 names_to = "parameter",
                 values_to = "value") %>%
    mutate(
      parameter = fct_inorder(parameter),
      parameter = fct_recode(
        parameter,
        "Bias" = "rel_bias",
        "Root-mean-square error" = "rel_rmse",
        "Standard error (empirical)" = "rel_se_empirical",
        "Coverage of\n95% confidence intervals" = "ci_coverage"
      )
    )
  vline_dat <- tribble(
    ~xintercept, ~parameter,
    0, "Bias",
    0.95, "Coverage of\n95% confidence intervals"
  ) %>%
    mutate(parameter = fct_inorder(parameter)) %>%
    expand_grid(sample_size = unique(sim_reps_summary$sample_size))
  blank_dat <- tribble(
    ~value, ~parameter,
    0, "Root-mean-square error",
    1, "Coverage of\n95% confidence intervals"
  ) %>%
    mutate(parameter = fct_inorder(parameter)) %>%
    expand_grid(sample_size = unique(sim_reps_summary$sample_size),
                treatment_effect = unique(sim_reps_summary$treatment_effect),
                estimator_name = unique(sim_reps_summary$estimator_name))
  ggplot(dat, aes(y = treatment_effect, x = value,
                  colour = estimator_name, shape = estimator_name)) +
    geom_vline(data = vline_dat, aes(xintercept = xintercept),
               linetype = "dashed", colour = "grey60") +
    geom_blank(data = blank_dat) +
    geom_point(size = 2.5, position = position_dodge2(0.5, reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    facet_grid(cols = vars(parameter),
               rows = vars(sample_size),
               scales = "free_x") +
    labs(x = NULL, y = NULL, colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
}

# Save plot of results for no-missing-data scenarios
save_plot_results_no_missing <- function(p) {
  ggsave(
    here("plots", "no_missing.png"),
    plot = p,
    width = 7, height = 4, dpi = 300
  )
}

single_parameter_plot <- function(
    dat,
    parameter,
    parameter_label,
    facet_cols = NULL,
    facet_rows = NULL,
    facet_scales = "fixed",
    xintercept = NULL,
    xexpand = NULL
) {
  if (nrow(dat) < 1) {
    return(ggplot())
  }
  p <- ggplot(dat, aes(y = missingness_name, x = {{ parameter }},
                       colour = estimator_name, shape = estimator_name))
  if (!is.null(xintercept)) {
    p <- p +
      geom_vline(xintercept = xintercept,
                 linetype = "dashed", colour = "grey60")
  }
  if (!is.null(xexpand)) {
    p <- p + expand_limits(x = xexpand)
  }
  p <- p +
    geom_point(size = 2.5, position = position_dodge2(0.5, reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    facet_grid(cols = facet_cols,
               rows = facet_rows,
               scales = facet_scales) +
    labs(x = parameter_label, y = NULL,
         colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
}

# Plots for scenarios with null treatment effects (excl. no-missing-data)
make_plot_results_null <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect == "null",
           missingness_mechanism != "none") %>%
    mutate(outcome_missingness = fct_recode(
      outcome_missingness,
      "complete outcome" = "no",
      "incomplete outcome" = "yes"
    ))
  do_plot <- function(...) {
    single_parameter_plot(
      dat,
      facet_cols = vars(outcome_missingness, sample_size),
      facet_rows = vars(missingness_mechanism),
      ...
    )
  }
  list(
    bias = do_plot(rel_bias, "Bias", xintercept = 0),
    rmse = do_plot(rel_rmse, "Root-mean-square error", xexpand = 0),
    se_empirical =
      do_plot(rel_se_empirical, "Standard error (empirical)"),
    ci_coverage =
      do_plot(ci_coverage, "Coverage of 95% confidence intervals",
              xintercept = 0.95, xexpand = 1)
  )
}

# Save plots of results for null scenarios
save_plot_results_null <- function(
    plot_list, outcome_missingness_, sample_size_
) {
  c(
    ggsave(
      here("plots", glue("null_bias.png")),
      plot = plot_list$bias,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("null_rmse.png")),
      plot = plot_list$rmse,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("null_se_empirical.png")),
      plot = plot_list$se_empirical,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("null_ci_coverage.png")),
      plot = plot_list$ci_coverage,
      width = 6, height = 4, dpi = 300
    )
  )
}

# Plots for scenarios with missing data and non-null treatment effects
make_plot_results_main <- function(
    sim_reps_summary, outcome_missingness_, sample_size_
) {
  dat <- sim_reps_summary %>%
    filter(missingness_mechanism != "none",
           treatment_effect != "null",
           estimator_name != "naive",
           outcome_missingness == outcome_missingness_,
           sample_size == sample_size_)
  do_plot <- function(...) {
    single_parameter_plot(
      dat,
      facet_cols = vars(treatment_effect),
      facet_rows = vars(missingness_mechanism),
      ...
    )
  }
  list(
    bias = do_plot(rel_bias, "Bias", xintercept = 0),
    rmse = do_plot(rel_rmse, "Root-mean-square error", xexpand = 0),
    se_empirical =
      do_plot(rel_se_empirical, "Standard error (empirical)"),
    ci_coverage =
      do_plot(ci_coverage, "Coverage of 95% confidence intervals",
              xintercept = 0.95, xexpand = 1)
  )
}

# Save plots for scenarios with missing data and non-null treatment effects
save_plot_results_main <- function(
    plot_list, outcome_missingness_, sample_size_
) {
  c(
    ggsave(
      here("plots", glue("main_bias_{sample_size_}_{outcome_missingness_}.png")),
      plot = plot_list$bias,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("main_rmse_{sample_size_}_{outcome_missingness_}.png")),
      plot = plot_list$rmse,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("main_se_empirical_{sample_size_}_{outcome_missingness_}.png")),
      plot = plot_list$se_empirical,
      width = 6, height = 4, dpi = 300
    ),
    ggsave(
      here("plots", glue("main_ci_coverage_{sample_size_}_{outcome_missingness_}.png")),
      plot = plot_list$ci_coverage,
      width = 6, height = 4, dpi = 300
    )
  )
}
