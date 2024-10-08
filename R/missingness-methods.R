# Methods for dealing with missing data

# Complete case analysis: drop observations with any missing data
missingness_cc <- function(dat, estimators) {
  dat_cc <- drop_na(dat, outcome, dose_binary)
  estimators %>%
    rowwise(estimator_name) %>%
    reframe(estimator_fn(dat_cc))
}

# Run an estimator function on multiply imputed datasets and pool the results
run_estimator_mi <- function(imp_all_df, estimator) {
  rows <- max(imp_all_df$.id)
  imp_results <- imp_all_df %>%
    group_by(.imp) %>%
    reframe(estimator(pick(everything())))
  pooled <- pool.scalar(imp_results$estimate, imp_results$std.error^2,
                        n = rows, k = rows - first(imp_results$df))
  tibble(
    estimate = pooled$qbar,
    std.error = sqrt(pooled$ubar),
    conf.low = estimate - qt(0.975, pooled$df)*std.error,
    conf.high = estimate + qt(0.975, pooled$df)*std.error,
    df = pooled$df
  )
}

# Multiple imputation of binary compliance measure
missingness_mi_binary <- function(dat, estimators, m = 50, iter = 10) {
  dat_for_imp <- dat %>%
    select(trt, aux, confounder, dose_binary, outcome) %>%
    mutate(dose_binary = as.factor(dose_binary))
  
  # bodge: check to see if all participants group have the
  # same dose (as binary variable) and if so, give up!
  if (all(dat_for_imp$dose_binary == 0, na.rm = TRUE)) {
    return(tibble())
  }
  
  # split data into control and treatment groups
  dat_ctrl <- filter(dat_for_imp, trt == 0)
  dat_trt <- filter(dat_for_imp, trt == 1)
  
  # Impute treatment and control separately
  imp_ctrl <- 
    suppressWarnings(mice(dat_ctrl, m = m, maxit = iter, printFlag = FALSE))
  imp_trt <- 
    suppressWarnings(mice(dat_trt, m = m, maxit = iter, printFlag = FALSE))
  # Combine imputed treatment and control groups
  imp_all <- rbind(imp_trt, imp_ctrl)
  # Convert imputed data to a single long-form data form
  imp_all_df <- complete(imp_all, action = "long", include = FALSE) %>%
    mutate(dose_binary = as.numeric(as.character(dose_binary)))

  # Run estimators on imputed data
  estimators %>%
    rowwise(estimator_name) %>%
    reframe(run_estimator_mi(imp_all_df, estimator_fn))
}

# Inverse probability weighting for missing data
missingness_ipw <- function(dat, estimators) {
  dat <- dat %>%
    mutate(complete_case = !is.na(dose_binary) & !is.na(outcome))
  weights_model <- glm(
    complete_case ~ trt*aux + trt*outcome,
    family = binomial,
    data = dat
  )

  dat_cc <- drop_na(dat, outcome, dose_binary)
  probs_out <- predict(weights_model, type = "response", newdata = dat_cc)
  weights_out <- 1 / probs_out

  estimators %>%
    rowwise(estimator_name) %>%
    reframe(estimator_fn(dat_cc, weights = weights_out))
}

missingness_methods <- tribble(
  ~missingness_name, ~missingness_fn,
  "cc", missingness_cc,
  "mi", missingness_mi_binary,
  "ipw", missingness_ipw,
)
