# Summarise results of simulations: (relative) bias
summarise_sim_reps <- function(sim_reps) {
  sim_reps %>%
    mutate(across(c(scenario_name, treatment_effect, compliance, dose_response,
                    sample_size, missingness_mechanism, outcome_missingness,
                    missingness_name, estimator_name),
                  fct_inorder)) %>%
    drop_na(error, std.error) %>%
    mutate(rel_error = error / true_effect) %>%
    group_by(
      scenario_name, 
      treatment_effect, compliance, dose_response,
      sample_size, missingness_mechanism, outcome_missingness,
      missingness_name, estimator_name
    ) %>%
    summarise(
      bias = mean(error),
      bias_se = sd(error) / sqrt(n()),
      se_empirical = sd(error),
      se_model = sqrt(mean(std.error^2)),
      rmse = sqrt(mean(error^2)),
      rel_bias = mean(rel_error),
      rel_bias_se = sd(rel_error) / sqrt(n()),
      rel_se_empirical = sd(rel_error),
      rel_se_model = sqrt(mean(rel_error^2)),
      rel_rmse = sqrt(mean(rel_error^2)),
      ci_coverage = mean(ci_includes_truth),
      ci_coverage_se = sqrt(ci_coverage * (1 - ci_coverage) / n()),
    ) %>%
    ungroup()
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

# Plots for scenarios with no missing data (excluding null treatment effects)
make_plot_results_no_missing <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(missingness_mechanism == "none", treatment_effect != "null")
  list(
    bias = ggplot(dat, aes(y = estimator_name, x = rel_bias,
                    xmin = rel_bias - 1.96*rel_bias_se, 
                    xmax = rel_bias + 1.96*rel_bias_se,
                    colour = sample_size, shape = sample_size)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.25)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative bias", y = NULL) +
      theme_cp(grid = "x"),
    rmse = ggplot(dat, aes(y = estimator_name, x = rel_rmse,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative root-mean-square error", y = NULL) +
      theme_cp(grid = "x"),
    se_empirical = ggplot(dat, aes(y = estimator_name, x = rel_se_empirical,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative standard error (empirical)", y = NULL) +
      theme_cp(grid = "x"),
    ci_coverage = ggplot(dat, aes(y = estimator_name, x = ci_coverage,
                    xmin = ci_coverage - 1.96*ci_coverage_se, 
                    xmax = ci_coverage + 1.96*ci_coverage_se,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 1) +
      geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Coverage of 95% confidence intervals", y = NULL) +
      theme_cp(grid = "x")
  )
}

# Save plot of results for no-missing-data scenarios
save_plot_results_no_missing <- function(plot_list) {
  p <- plot_grid(
    plotlist = plot_list,
    ncol = 1
  )
  ggsave(
    here("plots", "no_missing.png"),
    plot = p,
    width = 7, height = 9, dpi = 300
  )
}

# Plots for scenarios with null treatment effects (excl. no-missing-data)
make_plot_results_null <- function(
    sim_reps_summary, outcome_missingness_, sample_size_
) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect == "null", missingness_mechanism != "none",
           outcome_missingness == outcome_missingness_,
           sample_size == sample_size_)
  list(
    bias = ggplot(dat, aes(y = missingness_name, x = bias,
                    xmin = bias - 1.96*bias_se, xmax = bias + 1.96*bias_se,
                    colour = estimator_name, shape = estimator_name)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(missingness_mechanism)) +
      labs(x = "Bias", y = NULL) +
      theme_cp(grid = "x"),
    rmse = ggplot(dat, aes(y = missingness_name, x = rmse,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Root-mean-square error", y = NULL) +
      theme_cp(grid = "x"),
    se_empirical = ggplot(dat, aes(y = missingness_name, x = se_empirical,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Standard error (empirical)", y = NULL) +
      theme_cp(grid = "x"),
    ci_coverage = ggplot(dat, aes(y = missingness_name, x = ci_coverage,
                    xmin = ci_coverage - 1.96*ci_coverage_se, 
                    xmax = ci_coverage + 1.96*ci_coverage_se,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 1) +
      geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Coverage of 95% confidence intervals", y = NULL) +
      theme_cp(grid = "x")
  )
}

# Save plots of results for null scenarios
save_plot_results_null <- function(
    plot_list, outcome_missingness_, sample_size_
) {
  p <- plot_grid(
    plotlist = plot_list,
    ncol = 1
  )
  ggsave(
    here("plots", glue("null_{sample_size_}_{outcome_missingness_}.png")),
    plot = p,
    width = 7, height = 5, dpi = 300
  )
}

# Plots for scenarios with missing data and non-null treatment effects
make_plot_results_main <- function(
    sim_reps_summary, outcome_missingness_, sample_size_
) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect != "null", missingness_mechanism != "none",
           estimator_name != "naive",
           outcome_missingness == outcome_missingness_,
           sample_size == sample_size_) %>%
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
  list(
    bias = ggplot(dat, aes(y = missingness_name, x = rel_bias,
                           colour = estimator_name, shape = estimator_name)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_grid(cols = vars(treatment_effect),
                 rows = vars(missingness_mechanism)) +
      labs(x = "Relative bias", y = NULL,
           colour = "estimator", shape = "estimator") +
      theme_cp(grid = "x") +
      theme(strip.text = element_text(size = rel(8/11)),
            legend.position = "bottom"),
    rmse = ggplot(dat, aes(y = missingness_name, x = rel_rmse,
                           colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_grid(cols = vars(treatment_effect),
                 rows = vars(missingness_mechanism)) +
      labs(x = "Relative root-mean-square error", y = NULL,
           colour = "estimator", shape = "estimator") +
      theme_cp(grid = "x") +
      theme(strip.text = element_text(size = rel(8/11)),
            legend.position = "bottom"),
    se_empirical = ggplot(dat, aes(y = missingness_name, x = rel_se_empirical,
                                   colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.25), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_grid(cols = vars(treatment_effect),
                 rows = vars(missingness_mechanism)) +
      labs(x = "Relative standard error (empirical)", y = NULL,
           colour = "estimator", shape = "estimator") +
      theme_cp(grid = "x") +
      theme(strip.text = element_text(size = rel(8/11)),
            legend.position = "bottom"),
    ci_coverage = ggplot(dat, aes(y = missingness_name, x = ci_coverage,
                                  colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 1) +
      geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
      geom_point(position = position_dodge(0.25), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_grid(cols = vars(treatment_effect),
                 rows = vars(missingness_mechanism)) +
      labs(x = "Coverage of 95% confidence intervals", y = NULL,
           colour = "estimator", shape = "estimator") +
      theme_cp(grid = "x") +
      theme(strip.text = element_text(size = rel(8/11)),
            legend.position = "bottom")
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
