

get_compliance_intercept <- function(
  mu, b_aux, p_aux, b_confounder, p_confounder, n.mc = 1e6
) {
  res <- optim(
    par = c(intercept = qlogis(mu)),
    fn = function(intercept) {
      (mean(plogis(
        intercept + b_aux*rbinom(n.mc, 1, p_aux) + b_confounder*rbinom(n.mc, 1, p_confounder)
      )) - mu)^2
    },
    method = "Brent",
    lower = -10,
    upper = 10,
    control = list(reltol = 1/sqrt(n.mc))
  )
  if (res$convergence != 0) {
    print(res)
    stop("failed to converge")
  }
  res$par
}

generate_complete_df <- function(
    n,
    p_aux,
    p_confounder,
    compliance_intercept,
    compliance_b_confounder,
    compliance_b_aux,
    outcome_b_response,
    outcome_b_confounder,
    trt_split = 0.5
) {
  n_trt <- round(trt_split * n)
  n_ctrl <- n - n_trt
  tibble(
    # trt: 0 = control, 1 = treatment
    trt = c(rep(1, n_trt), rep(0, n_ctrl)),
    # auxiliary variable: Bern(p_aux)
    aux = rbinom(n, 1, p_aux),
    # confounder variable: Bern(p_confounder)
    confounder = rbinom(n, 1, p_confounder),
    # mean of compliance variable
    compliance_prob = plogis(
      compliance_intercept +
        compliance_b_confounder * confounder +
        compliance_b_aux * aux
    ),
    # latent compliance (0 or 1): would this person comply with the intervention
    # if assigned to it
    compliance = rbinom(n, 1, compliance_prob),
    # binary compliance: always 1 for control group; otherwise latent compliance
    compliance_binary =
      as.numeric(if_else(trt == 1, compliance == 1, TRUE)),
    # dose received (0 or 1): 0 for control group, compliance for trt group
    dose_binary = if_else(trt == 1, compliance_binary, 0),
    # mean of the outcome variable:
    #  b1 * dose + b2 * confounder
    outcome_mu =
      outcome_b_response * dose_binary + outcome_b_confounder * confounder,
    # outcome: drawn from a normal distribution
    outcome = rnorm(n, outcome_mu, 1),
  )
}


add_missingness_mar <- function(
    df,
    miss_compliance_enabled,
    miss_compliance_intercept, miss_compliance_b_aux, miss_compliance_b_outcome,
    miss_outcome_enabled,
    miss_outcome_intercept, miss_outcome_b_aux, miss_outcome_b_confounder
) {
  n <- nrow(df)
  df %>%
    mutate(
      miss_compliance_lp = miss_compliance_intercept + 
        miss_compliance_b_aux*aux + miss_compliance_b_outcome*outcome,
      miss_compliance_prob = plogis(miss_compliance_lp),
      miss_compliance =
        if (miss_compliance_enabled) {
          rbinom(n, 1, miss_compliance_prob)
        } else {
          0
        },
      miss_outcome_lp = miss_outcome_intercept + 
        miss_outcome_b_aux*aux + miss_outcome_b_confounder*confounder,
      miss_outcome_prob = plogis(miss_outcome_lp),
      miss_outcome = 
        if (miss_outcome_enabled) {
          rbinom(n, 1, miss_outcome_prob)
        } else {
          0
        },
      compliance = if_else(miss_compliance == 1, NA, compliance),
      compliance_binary = 
        if_else(miss_compliance == 1 & trt == 1, NA, compliance_binary),
      dose_binary = 
        if_else(miss_compliance == 1 & trt == 1, NA, dose_binary),
      outcome = if_else(miss_outcome == 1, NA, outcome),
    )
}

generate_scenario_data <- function(scenario_params) {
  stopifnot(nrow(scenario_params) == 1)

  # generate complete dataset, i.e. no missingness yet
  complete_dat <- generate_complete_df(
    n = scenario_params$n,
    p_aux = scenario_params$p_aux,
    p_confounder = scenario_params$p_confounder,
    compliance_intercept = scenario_params$compliance_intercept,
    compliance_b_confounder = scenario_params$compliance_b_confounder,
    compliance_b_aux = scenario_params$compliance_b_aux,
    outcome_b_response = scenario_params$outcome_b_response,
    outcome_b_confounder = scenario_params$outcome_b_confounder
  )

  # add missingness according to missingness parameters
  dat <- add_missingness_mar(
    complete_dat,
    miss_compliance_enabled = scenario_params$miss_compliance_enabled,
    miss_compliance_intercept = scenario_params$miss_compliance_intercept,
    miss_compliance_b_aux = scenario_params$miss_compliance_b_aux,
    miss_compliance_b_outcome = scenario_params$miss_compliance_b_outcome,
    miss_outcome_enabled = scenario_params$miss_outcome_enabled,
    miss_outcome_intercept = scenario_params$miss_outcome_intercept,
    miss_outcome_b_aux = scenario_params$miss_outcome_b_aux,
    miss_outcome_b_confounder = scenario_params$miss_outcome_b_confounder
  )

  dat
}

get_mar_compliance_intercept <- function(
  args, prop_missing_compliance,
  n.mc = 1e6, verbose = FALSE, interval = c(-5, 5)
) {
  if (!args$miss_compliance_enabled) {
    return(0)
  }
  if (args$miss_compliance_b_aux == 0 && args$miss_compliance_b_outcome == 0) {
    return(qlogis(prop_missing_compliance))
  }

  args$miss_outcome_enabled <- FALSE
  args$miss_outcome_intercept <- 0
  args$miss_outcome_b_aux <- 0
  args$miss_outcome_b_confounder <- 0

  args$n <- n.mc

  res <- optimize(
    f = function(par) {
      call_args <- args
      call_args$miss_compliance_intercept <- par
      dat <- generate_scenario_data(call_args)
      obs_missing_compliance <- mean(is.na(dat$compliance[dat$trt == 1]))
      if (verbose) {
        cat("par =", par, 
            "; obs_missing_compliance =", obs_missing_compliance,
            "\n")
      }
      (obs_missing_compliance - prop_missing_compliance)^2
    },
    interval = interval,
    tol = 1 / sqrt(n.mc)
  )
  res$minimum
}

get_mar_outcome_intercept <- function(
  args, prop_missing_outcome,
  n.mc = 1e6, verbose = FALSE, interval = c(-5, 5)
) {
  if (!args$miss_outcome_enabled) {
    return(0)
  }

  args$n <- n.mc

  res <- optimize(
    f = function(par) {
      call_args <- args
      call_args$miss_outcome_intercept <- par
      dat <- generate_scenario_data(call_args)
      obs_missing_outcome <- mean(is.na(dat$outcome))
      if (verbose) {
        cat("par =", par, 
            "; obs_missing_outcome =", obs_missing_outcome,
            "\n")
      }
      (obs_missing_outcome - prop_missing_outcome)^2
    },
    interval = interval,
    tol = 1 / sqrt(n.mc)
  )
  cat("miss_outcome_intercept =", res$minimum, "\n")
  res$minimum
}
