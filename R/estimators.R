# Estimators for efficacy in complete data

# Naive estimator: subset the data to those who complied with the assigned
# treatment only, and calculate the difference between treatment groups
estimator_naive <- function(dat, weights = NULL) {
  dat$weights <- weights
  dat_compliers <- filter(dat, trt == 0 | dose_binary == 1)
  if (is.null(weights)) {
    m <- lm_robust(outcome ~ trt,
      data = dat_compliers,
      se_type = "classical"
    )
  } else {
    m <- lm_robust(outcome ~ trt,
      data = dat_compliers,
      weights = dat_compliers$weights, se_type = "HC2"
    )
  }
  m %>%
    tidy(conf.int = TRUE) %>%
    as_tibble() %>%
    filter(term == "trt") %>%
    select(estimate, std.error, conf.low, conf.high, df)
}

# Standardisation estimator: fit model to treatment and confounder (with
# interaction), use marginaleffects::avg_comparisons() to obtain a
# standardised result. Uses HC2 standard errors if weights are applied
estimator_standardisation <- function(dat, weights = NULL) {
  m <- lm(
    outcome ~ dose_binary + confounder + dose_binary * confounder,
    weights = weights, data = dat
  )
  if (is.null(weights)) {
    vcov_type <- TRUE
    wts <- FALSE
  } else {
    vcov_type <- "HC2"
    wts <- weights
  }
  df <- get_df(m)
  avg_comparisons(
    m, variables = "dose_binary", df = df,
    wts = wts, vcov = vcov_type
  ) %>%
    as_tibble() %>%
    select(estimate, std.error, conf.low, conf.high) %>%
    mutate(df = df)
}

# IPTW estimate using GLM for weights, HC2 robust standard errors.
estimator_iptw_glm <- function(dat, weights = NULL) {
  # Fit weights model
  weights_model <- glm(
    dose_binary ~ trt * confounder,
    family = binomial,
    data = dat
  )

  # Extract predicted probabilities from weights model and calculate weights
  # Stabilising weights based on overage probability of receiving treatment
  p_marginal <- mean(dat$dose_binary)
  probs_out <- predict(weights_model, type = "response")
  weights_out <- case_when(
    dat$dose_binary == 1 ~ p_marginal / probs_out,
    dat$dose_binary == 0 ~ (1 - p_marginal) / (1 - probs_out)
  )

  # If weights are provided (i.e. IPMW for missingness), combine (multiply)
  # these with the IPTW weights
  if (!is.null(weights)) {
    weights_out <- weights_out * weights
  }

  # Fit outcome model with HC2 standard errors
  m <- lm_robust(
    outcome ~ dose_binary,
    weights = weights_out,
    data = dat,
    se_type = "HC2"
  )

  # Extract coefficient for treatment received
  tidy(m) %>%
    filter(term == "dose_binary") %>%
    as_tibble() %>%
    select(estimate, std.error, conf.low, conf.high, df)
}

# Instrumental variables (2SLS) estimator
estimator_ivreg <- function(dat, weights = NULL) {
  # Fit instrumental variables model
  m <- ivreg(
    outcome ~ dose_binary | trt,
    weights = weights, data = dat
  )

  # Weights parameters for marginaleffects::avg_comparisons
  if (is.null(weights)) {
    vcov_type <- TRUE
    wts <- FALSE
  } else {
    vcov_type <- "HC2"
    wts <- weights
  }
  df <- get_df(m)

  # Obtain comparisons, with HC2 standard errors if using weights
  avg_comparisons(
    m, variables = "dose_binary", df = df,
    wts = wts, vcov = vcov_type
  ) %>%
    as_tibble() %>%
    select(estimate, std.error, conf.low, conf.high) %>%
    mutate(df = df)
}

# Table of estimators to apply to all generated scenarios
estimators <- tribble(
  ~estimator_name, ~estimator_fn,
  "naive", estimator_naive,
  "standardisation", estimator_standardisation,
  "iptw", estimator_iptw_glm,
  "iv", estimator_ivreg,
)
