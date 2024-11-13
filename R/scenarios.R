
get_scenario_list <- function() {
  expand_grid(
    treatment_effect = c(
      "null",
      "power80"
    ),
    compliance = c(
      "high",
      "low"
    ),
    sample_size = c("small", "large"),
    missingness_mechanism = c("none", "mcar", "mar_weak", "mar_strong"),
    outcome_missingness = c("no"),
  ) %>%
    filter(!(missingness_mechanism == "none" & 
               outcome_missingness == "yes")) %>%
    mutate(
      null = str_detect(treatment_effect, "null"),
    ) %>%
    mutate(
      scenario_name = glue(
        "{treatment_effect}_{compliance}_{sample_size}_",
        "{missingness_mechanism}"
      ),
      .before = 1
    )
}

scenario_list <- get_scenario_list()

get_dgp_params <- function(scenarios) {
  scenarios %>%
    distinct(sample_size, null, compliance) %>%
    mutate(
      # sample size
      n = case_when(
        sample_size == "small" ~ 100,
        sample_size == "large" ~ 1000,
      ),

      p_aux = 0.5,
      p_confounder = 0.25,

      # compliance
      compliance_prop = case_when(
        compliance == "low" ~ 0.4,
        compliance == "high" ~ 0.8,
      ),
      compliance_b_confounder = log(0.6),
      compliance_b_aux = log(1.25),

      # outcome regression
      outcome_b_confounder = -0.5
    ) %>%
    group_by(compliance) %>%
    mutate(
      compliance_intercept = 
        get_compliance_intercept(
          first(compliance_prop),
          first(compliance_b_aux), first(p_aux),
          first(compliance_b_confounder), first(p_confounder)
        ),
    ) %>%
    rowwise() %>%
    mutate(
      outcome_b_response =
        case_when(
          null ~ 0,
          TRUE ~ pwr.t.test(n = n / 2, power = 0.8)$d / compliance_prop,
        ),
    ) %>%
    ungroup() %>%
    mutate(
      true_effect = outcome_b_response
    )
}

get_scenario_params <- function(scenarios, dgp_params) {
  scenarios %>%
    left_join(dgp_params) %>%
    mutate(
      # missingness parameters, exc. intercept
      miss_compliance_enabled = missingness_mechanism != "none",
      miss_compliance_b_aux = case_when(
        missingness_mechanism == "mar_weak" ~ log(0.8),
        missingness_mechanism == "mar_strong" ~ log(0.6),
        .default = 0.0
      ),
      miss_compliance_b_outcome = case_when(
        missingness_mechanism == "mar_weak" ~ log(0.9),
        missingness_mechanism == "mar_strong" ~ log(0.8),
        .default = 0.0
      ),
      
      miss_outcome_enabled = 
        missingness_mechanism != "none" & outcome_missingness == "yes",
      miss_outcome_b_aux = case_when(
        missingness_mechanism == "mar_weak" ~ 0.2,
        missingness_mechanism == "mar_strong" ~ 0.5,
        .default = 0.0
      ),
      miss_outcome_b_confounder = case_when(
        missingness_mechanism == "mar_weak" ~ 0.2,
        missingness_mechanism == "mar_strong" ~ 0.5,
        .default = 0.0
      ),
    ) %>%
    rowwise(everything()) %>%
    summarise(
      miss_compliance_intercept =
        get_mar_compliance_intercept(cur_group(),
                                     prop_missing_compliance = 0.3)
    ) %>%
    ungroup() %>%
    rowwise(everything()) %>%
    summarise(
      miss_outcome_intercept =
        get_mar_outcome_intercept(cur_group(),
                                  prop_missing_outcome = 0.3)
    ) %>%
    ungroup()
}

run_scenario <- function(scenario_params) {
  dat <- generate_scenario_data(scenario_params)

  meth <- missingness_methods
  if (scenario_params$missingness_mechanism == "none") {
    meth <- filter(meth, missingness_name == "cc")
  }
  result <- meth %>%
    rowwise(missingness_name) %>%
    reframe(missingness_fn(dat, estimators))
  bind_cols(scenario_params, result) %>%
    mutate(
      error = estimate - true_effect,
      ci_includes_truth = true_effect >= conf.low & true_effect <= conf.high,
    )
}

test_scenario_rep <- function(scenario_name, rep, global_seed) {
  rep_seed <- tar_seed_create(
    paste0("sim_reps_", scenario_name, " ", rep),
    global_seed = global_seed
  )
  name <- scenario_name
  params <- tar_read(scenario_params) %>%
    filter(scenario_name == name)
  with_seed(rep_seed, run_scenario(params))
}
