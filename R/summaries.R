
summarise_sim_reps <- function(sim_reps) {
  sim_reps %>%
    mutate(across(c(scenario_name, treatment_effect, compliance, dose_response,
                    sample_size, missingness_mechanism, outcome_missingness,
                    missingness_name, estimator_name),
                  fct_inorder)) %>%
    drop_na(error, std.error) %>%
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
      rel_bias = mean(error / true_effect),
      rel_bias_se = sd(error / true_effect) / sqrt(n()),
      rel_se_empirical = sd(error / true_effect),
      rel_se_model = sqrt(mean((std.error / true_effect)^2)),
      rel_rmse = sqrt(mean((error / true_effect)^2)),
      ci_coverage = mean(ci_includes_truth),
      ci_coverage_se = sqrt(ci_coverage * (1 - ci_coverage) / n()),
    ) %>%
    ungroup()
}

theme_cp <- function(panel_border = TRUE, grid = "none") {
  t <- theme_cowplot(font_size = 11, rel_small = 9/11,
                     rel_tiny = 8/11, rel_large = 13/11) +
    theme(plot.background = element_rect(fill = "white")) +
    background_grid(major = grid)
  if (panel_border) {
    t <- t + panel_border()
  }
}

plot_results_no_missing <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(missingness_mechanism == "none", treatment_effect != "null")
  plot_grid(
    ggplot(dat, aes(y = estimator_name, x = rel_bias,
                    xmin = rel_bias - 1.96*rel_bias_se, 
                    xmax = rel_bias + 1.96*rel_bias_se,
                    colour = sample_size, shape = sample_size)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.25)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative bias", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = estimator_name, x = rel_rmse,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative root-mean-square error", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = estimator_name, x = rel_se_empirical,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Relative standard error (empirical)", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = estimator_name, x = ci_coverage,
                    xmin = ci_coverage - 1.96*ci_coverage_se, 
                    xmax = ci_coverage + 1.96*ci_coverage_se,
                    colour = sample_size, shape = sample_size)) +
      expand_limits(x = 1) +
      geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(treatment_effect)) +
      labs(x = "Coverage of 95% confidence intervals", y = NULL) +
      theme_cp(grid = "x"),
    ncol = 1
  )
  ggsave("plots/no_missing.png", width = 7, height = 9, dpi = 300)
}

plot_results_null <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect == "null", missingness_mechanism != "none",
           outcome_missingness == "no", sample_size == "small")
  plot_grid(
    ggplot(dat, aes(y = missingness_name, x = bias,
                    xmin = bias - 1.96*bias_se, xmax = bias + 1.96*bias_se,
                    colour = estimator_name, shape = estimator_name)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      scale_y_discrete(limits = rev) +
      facet_wrap(vars(missingness_mechanism)) +
      labs(x = "Bias", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = missingness_name, x = rmse,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Root-mean-square error", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = missingness_name, x = se_empirical,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 0) +
      geom_point(position = position_dodge(0.5), size = 2.5) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Standard error (empirical)", y = NULL) +
      theme_cp(grid = "x"),
    ggplot(dat, aes(y = missingness_name, x = ci_coverage,
                    xmin = ci_coverage - 1.96*ci_coverage_se, 
                    xmax = ci_coverage + 1.96*ci_coverage_se,
                    colour = estimator_name, shape = estimator_name)) +
      expand_limits(x = 1) +
      geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
      geom_pointrange(position = position_dodge(0.5)) +
      facet_wrap(vars(missingness_mechanism)) +
      scale_y_discrete(limits = rev) +
      labs(x = "Coverage of 95% confidence intervals", y = NULL) +
      theme_cp(grid = "x"),
    ncol = 1
  )
  ggsave("plots/null.png", width = 7, height = 5, dpi = 300)
}

plot_results_main <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect != "null", missingness_mechanism != "none",
           outcome_missingness == "no", sample_size == "small")
  ggplot(dat, aes(y = missingness_name, x = rel_bias,
                  xmin = rel_bias - 1.96*rel_bias_se, 
                  xmax = rel_bias + 1.96*rel_bias_se,
                  colour = estimator_name, shape = estimator_name)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_y_discrete(limits = rev) +
    facet_wrap(vars(treatment_effect, missingness_mechanism)) +
    labs(x = "Relative bias", y = NULL,
         colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
  p1 <- ggsave("plots/main_bias.png", width = 6, height = 3, dpi = 300)
  ggplot(dat, aes(y = missingness_name, x = rel_rmse,
                  colour = estimator_name, shape = estimator_name)) +
    expand_limits(x = 0) +
    geom_point(position = position_dodge(0.5), size = 2.5) +
    scale_y_discrete(limits = rev) +
    facet_wrap(vars(treatment_effect, missingness_mechanism)) +
    labs(x = "Relative root-mean-square error", y = NULL,
         colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
  p2 <- ggsave("plots/main_rmse.png", width = 6, height = 3, dpi = 300)
  ggplot(dat, aes(y = missingness_name, x = rel_se_empirical,
                  colour = estimator_name, shape = estimator_name)) +
    expand_limits(x = 0) +
    geom_point(position = position_dodge(0.25), size = 2.5) +
    scale_y_discrete(limits = rev) +
    facet_wrap(vars(treatment_effect, missingness_mechanism)) +
    labs(x = "Relative standard error (empirical)", y = NULL,
         colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
  p3 <- ggsave("plots/main_se_empirical.png", width = 6, height = 3, dpi = 300)
  ggplot(dat, aes(y = missingness_name, x = ci_coverage,
                  colour = estimator_name, shape = estimator_name)) +
    expand_limits(x = 1) +
    geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
    geom_point(position = position_dodge(0.25), size = 2.5) +
    scale_y_discrete(limits = rev) +
    facet_wrap(vars(treatment_effect, missingness_mechanism)) +
    labs(x = "Coverage of 95% confidence intervals", y = NULL,
         colour = "estimator", shape = "estimator") +
    theme_cp(grid = "x") +
    theme(strip.text = element_text(size = rel(8/11)),
          legend.position = "bottom")
  p4 <- ggsave("plots/main_ci_coverage.png", width = 6, height = 3, dpi = 300)
  c(p1, p2, p3, p4)
}
