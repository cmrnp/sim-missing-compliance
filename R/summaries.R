
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
      ci_coverage = mean(ci_includes_truth)
    ) %>%
    ungroup()
}

plot_results_no_missing <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(missingness_mechanism == "none", treatment_effect != "null")
  ggplot(dat, aes(y = estimator_name, x = rel_bias,
                  xmin = rel_bias - 1.96*rel_bias_se, 
                  xmax = rel_bias + 1.96*rel_bias_se,
                  colour = sample_size, shape = sample_size)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(position = position_dodge(0.25)) +
    facet_grid(cols = vars(compliance), rows = vars(dose_response),
               labeller = label_both)
  ggplot(dat, aes(y = estimator_name, x = rel_rmse,
                  colour = sample_size, shape = sample_size)) +
    expand_limits(x = 0) +
    geom_point(position = position_dodge(0.25)) +
    facet_grid(cols = vars(compliance), rows = vars(dose_response),
               labeller = label_both)
  ggplot(dat, aes(y = estimator_name, x = rel_se_empirical,
                  colour = sample_size, shape = sample_size)) +
    expand_limits(x = 0) +
    geom_point(position = position_dodge(0.25)) +
    facet_grid(cols = vars(compliance), rows = vars(dose_response),
               labeller = label_both)
  ggplot(dat, aes(y = estimator_name, x = ci_coverage,
                  colour = sample_size, shape = sample_size)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey60") +
    geom_point(position = position_dodge(0.25)) +
    facet_grid(cols = vars(compliance), rows = vars(dose_response),
               labeller = label_both)
}

plot_results_null <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect == "null", missingness_mechanism != "none",
           outcome_missingness == "no")
  ggplot(dat, aes(y = estimator_name, x = bias,
                  xmin = bias - 1.96*bias_se, xmax = bias + 1.96*bias_se,
                  colour = sample_size, shape = sample_size)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(position = position_dodge(0.25)) +
    facet_grid(cols = vars(missingness_name), rows = vars(missingness_mechanism))
}

plot_results_main <- function(sim_reps_summary) {
  dat <- sim_reps_summary %>%
    filter(treatment_effect != "null", missingness_mechanism != "none",
           outcome_missingness == "no",
           estimator_name != "naive")
  ggplot(dat, aes(y = estimator_name, x = rel_bias,
                  xmin = rel_bias - 1.96*rel_bias_se, 
                  xmax = rel_bias + 1.96*rel_bias_se,
                  colour = sample_size, shape = sample_size)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(position = position_dodge(0.25)) +
    facet_grid(cols = vars(dose_response, compliance), rows = vars(missingness_name, missingness_mechanism))
}
