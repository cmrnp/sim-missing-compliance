# Summarise results of simulations:
# - (relative) bias
# - (relative) empirical standard error
# - (relative) model-based standard error
# - confidence interval coverage
summarise_sim_reps <- function(sim_reps) {
  sim_reps %>%
    mutate(across(c(scenario_name, treatment_effect, compliance,
                    sample_size, missingness_mechanism, outcome_missingness,
                    missingness_name, estimator_name),
                  fct_inorder)) %>%
    drop_na(error, std.error) %>%
    mutate(
      rel_error = if_else(abs(true_effect) < 1e-6, error, error / true_effect),
      rel_std_error = if_else(abs(true_effect) < 1e-6, std.error, std.error / true_effect)
    ) %>%
    group_by(
      scenario_name, 
      treatment_effect, compliance,
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

