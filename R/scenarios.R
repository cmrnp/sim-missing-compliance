
get_scenario_list <- function() {
  expand_grid(
    treatment_effect = c(
      "null",
      "high_compliance_linear",
      "high_compliance_nonlinear",
      "low_compliance_linear",
      "low_compliance_nonlinear"
    ),
    sample_size = c("small", "large"),
    missingness_mechanism = c("none", "mcar", "mar_weak", "mar_strong"),
    outcome_missingness = c("no", "yes"),
  ) %>%
    filter(!(missingness_mechanism == "none" & 
               outcome_missingness == "yes")) %>%
    mutate(
      null = str_detect(treatment_effect, "null"),
      compliance = 
        if_else(str_detect(treatment_effect, "high_compliance"), 
                "high", "low"),
      dose_response = 
        if_else(str_detect(treatment_effect, "_nonlinear"),
                "nonlinear", "linear"),
    ) %>%
    mutate(
      scenario_name = glue(
        "{treatment_effect}_{sample_size}_",
        "{missingness_mechanism}_{outcome_missingness}"
      ),
      .before = 1
    )
}

scenario_list <- get_scenario_list()

get_dgp_params <- function(scenarios) {
  dgp_params <-
    scenarios %>%
    distinct(sample_size, null, compliance, dose_response) %>%
    mutate(
      # sample size
      n = case_when(
        sample_size == "small" ~ 200,
        sample_size == "large" ~ 1000,
      ),
      
      # compliance
      compliance_intercept = case_when(
        compliance == "low" ~ get_compliance_intercept(0.5, 0.5, 0.5),
        compliance == "high" ~ get_compliance_intercept(0.8, 0.5, 0.5)
      ),
      compliance_b_confounder = 0.5,
      compliance_b_aux = 0.5,
      compliance_shape = case_when(
        compliance == "low" ~ 5,
        compliance == "high" ~ 15,
      ),
      compliance_threshold = 0.8,
      
      # dose-response
      dose_response_linear = dose_response == "linear",
      dose_response_location = get_sigmoid_params(0.5, 0.2, 0.8, 0.9)["location"],
      dose_response_shape = get_sigmoid_params(0.5, 0.2, 0.8, 0.9)["shape"],
      
      # outcome regression
      outcome_b_confounder = 0.5
    ) %>%
    rowwise(n, starts_with("compliance_"), starts_with("dose_response_"),
            outcome_b_confounder) %>%
    mutate(
      outcome_b_response = 
        case_when(
          null ~ 0,
          TRUE ~ get_outcome_b_response(cur_group(), power = 0.8)
        )
    ) %>%
    rowwise(n, starts_with("compliance_"), starts_with("dose_response_"),
            starts_with("outcome_")) %>%
    mutate(
      true_effect = get_true_effect(cur_group())
    ) %>%
    ungroup()
}

get_scenario_params <- function(scenarios, dgp_params) {
  scenarios %>%
    left_join(dgp_params) %>%
    mutate(
      # missingness parameters, exc. intercept
      miss_compliance_enabled = missingness_mechanism != "none",
      miss_compliance_b_aux = case_when(
        missingness_mechanism == "mar_weak" ~ 0.1,
        missingness_mechanism == "mar_strong" ~ 0.5,
        .default = 0.0
      ),
      miss_compliance_b_outcome = case_when(
        missingness_mechanism == "mar_weak" ~ 0.1,
        missingness_mechanism == "mar_strong" ~ 0.5,
        .default = 0.0
      ),
      
      miss_outcome_enabled = 
        missingness_mechanism != "none" & outcome_missingness == "yes",
      miss_outcome_b_aux = case_when(
        missingness_mechanism == "mar_weak" ~ 0.1,
        missingness_mechanism == "mar_strong" ~ 0.5,
        .default = 0.0
      ),
      miss_outcome_b_confounder = case_when(
        missingness_mechanism == "mar_weak" ~ 0.1,
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
