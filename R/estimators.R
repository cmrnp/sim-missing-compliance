# Estimators for efficacy in complete data

# naive estimator: subset the data to those who complied with the assigned
# treatment only, and calculate the difference between treatment groups
estimator_naive <- function(dat, weights = NULL) {
  dat$weights <- weights
  dat_compliers <- filter(dat, trt == 0 | dose_binary == 1)
  if (is.null(weights)) {
    m <- lm_robust(outcome ~ trt, data = dat_compliers,
                   se_type = "classical")
  } else {
    m <- lm_robust(outcome ~ trt, data = dat_compliers,
                   weights = dat_compliers$weights, se_type = "HC2")
  }
  m %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "trt") %>%
    select(estimate, std.error, conf.low, conf.high, df)
}

# standardisation estimate
estimator_standardisation <- function(dat, weights = NULL) {
  if (is.null(weights)) {
    m <- lm(
      outcome ~ dose_binary + confounder + dose_binary*confounder,
      data = dat
    )
    vcov_type <- TRUE
    wts <- FALSE
  } else {
    m <- lm(
      outcome ~ dose_binary + confounder + dose_binary*confounder,
      weights = weights, data = dat
    )
    vcov_type <- "HC2"
    wts <- weights
  }
  df = get_df(m)
  avg_comparisons(
    m, variables = "dose_binary", df = df, 
    wts = wts, vcov = vcov_type
  ) %>%
    as_tibble() %>%
    select(estimate, std.error, conf.low, conf.high) %>%
    mutate(df = df)
}

# IPTW estimate
estimator_iptw_glm <- function(dat, weights = NULL) {
  weights_model <- glm(
    dose_binary ~ trt*confounder + trt*aux,
    family = binomial,
    data = dat
  )
  probs_out <- predict(weights_model, type = "response")
  weights_out <- if_else(
    dat$dose_binary == 1, 
    1 / probs_out, 
    1 / (1 - probs_out)
  )
  if (!is.null(weights)) {
    # XXX is this legit?
    weights_out <- weights_out * weights
  }
  m <- lm_robust(
    outcome ~ dose_binary,
    weights = weights_out,
    data = dat,
    se_type = "HC2"
  )
  tidy(m) %>%
    filter(term == "dose_binary") %>%
    select(estimate, std.error, conf.low, conf.high, df)
}

# IPTW estimate using GAM for weights
estimator_iptw_gam <- function(dat, weights = NULL) {
  weights_model <- gam(
    dose_binary ~ s(confounder) + s(aux),
    family = binomial,
    data = filter(dat, trt == 1)
  )
  probs_out <- as.vector(
    predict(weights_model, type = "response", newdata = dat)
  )
  weights_out <- case_when(
    dat$trt == 0 ~ 1,
    dat$dose_binary == 1 ~ 1 / probs_out,
    TRUE ~ 1 / (1 - probs_out)
  )
  if (!is.null(weights)) {
    # XXX is this legit?
    weights_out <- weights_out * weights
  }
  m <- lm_robust(
    outcome ~ dose_binary,
    weights = weights_out,
    data = dat,
    se_type = "HC2"
  )
  tidy(m) %>%
    filter(term == "dose_binary") %>%
    select(estimate, std.error, conf.low, conf.high, df)
}

estimators <- tribble(
  ~estimator_name, ~estimator_fn,
  "naive", estimator_naive,
  "standardisation", estimator_standardisation,
#  "iptw", estimator_iptw_glm,
  "iptw_glm", estimator_iptw_glm,
  "iptw_gam", estimator_iptw_gam,
)

