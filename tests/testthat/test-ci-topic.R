# Confidence-interval smoke tests:
# - DR-Learner linear CI via StatsModelsLinearRegression
# - bootstrap_ci opt-in shared machinery + hybrid policy
# - CEVAE Pyro Predictive (ATE + ITE CI)
# - BCF MCMC diagnostics (n_chains > 1 + Rhat)
# - ci_method label propagation (summary.cf_model + metadata)
# - res$ate_ci / res$tau_ci unified location + validator.

.simulate_data <- function(n = 150, seed = 0L) {
  set.seed(seed)
  X1 <- stats::rnorm(n)
  X2 <- stats::rnorm(n)
  T  <- stats::rbinom(n, 1, 0.5)
  Y  <- 1.5 * T + 0.5 * X1 + stats::rnorm(n)
  data.frame(X1 = X1, X2 = X2, T = T, Y = Y)
}

# ----- ate_ci / tau_ci / ci_method on R-native methods -----

test_that("grf populates res$ate_ci, res$tau_ci, and ci_method = 'asymptotic Wald'", {
  skip_if_not_installed("grf")
  df <- .simulate_data(n = 200, seed = 1L)
  m <- cf_fit(Y ~ T | X1 + X2, data = df, method = "grf", random_state = 0L)
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 200L)
  expect_true(all(c("lower", "upper") %in% colnames(res$tau_ci)))
  expect_identical(res$ci_method, "asymptotic Wald")
  expect_match(format(summary(m)), "asymptotic Wald", all = FALSE)
})

# ----- hybrid policy: error paths before any fit work -----

test_that("bootstrap_ci = TRUE errors for NN methods (tarnet/dragonnet/ganite)", {
  df <- .simulate_data(n = 80, seed = 5L)
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "tarnet", bootstrap_ci = TRUE),
    class = "cfomics_bootstrap_unsupported_nn"
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "dragonnet", bootstrap_ci = TRUE),
    class = "cfomics_bootstrap_unsupported_nn"
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "ganite", bootstrap_ci = TRUE),
    class = "cfomics_bootstrap_unsupported_nn"
  )
})

# ----- bootstrap_ci accept paths (metalearners) -----

test_that("slearner with bootstrap_ci = TRUE produces ate_ci + tau_ci", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 100, seed = 6L)
  m <- cf_fit(
    Y ~ T | X1 + X2, data = df, method = "slearner",
    bootstrap_ci = TRUE, n_bootstrap = 10L, random_state = 1L
  )
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 100L)
  expect_identical(res$ci_method, "nonparametric bootstrap")
})

# ----- DR-Learner linear via StatsModelsLinearRegression -----

test_that("drlearner (model_final = 'linear') returns asymptotic Wald ATE + ITE CIs", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 200, seed = 7L)
  m <- cf_fit(
    Y ~ T | X1 + X2, data = df, method = "drlearner",
    model_final = "linear", random_state = 1L
  )
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 200L)
  expect_identical(res$ci_method, "asymptotic Wald")
})

# ----- CEVAE Pyro Predictive -----

test_that("cevae populates res$ate_ci / res$tau_ci via Pyro Predictive", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 80, seed = 8L)
  m <- cf_fit(
    Y ~ T | X1 + X2, data = df, method = "cevae",
    random_state = 1L, num_epochs = 3L, num_samples = 5L
  )
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 80L)
  expect_identical(res$ci_method, "Bayesian posterior")
})

# ----- BCF MCMC diagnostic -----

test_that("bcf with n_chains > 1 produces a Gelman-Rubin Rhat in metadata", {
  skip_if_not_installed("bcf")
  df <- .simulate_data(n = 100, seed = 9L)
  m <- suppressWarnings(cf_fit(
    Y ~ T | X1 + X2, data = df, method = "bcf",
    n_burn = 50L, n_iter = 100L, n_chains = 2L,
    random_state = 1L
  ))
  expect_equal(m$res$metadata$n_chains, 2L)
  expect_true(is.numeric(m$res$metadata$rhat))
  expect_false(is.na(m$res$metadata$rhat))
  expect_identical(m$res$ci_method, "Bayesian posterior")
})

test_that("bcf does not carry res$credible_intervals", {
  skip_if_not_installed("bcf")
  df <- .simulate_data(n = 80, seed = 10L)
  m <- suppressWarnings(cf_fit(
    Y ~ T | X1 + X2, data = df, method = "bcf",
    n_burn = 50L, n_iter = 100L, random_state = 1L
  ))
  expect_null(m$res$credible_intervals)
  expect_s3_class(m$res$tau_ci, "data.frame")
})

# ----- validator: bad ate_ci should abort -----

test_that("validate_cfomics_result rejects malformed ate_ci / tau_ci", {
  skip_if_not_installed("grf")
  df <- .simulate_data(n = 100, seed = 11L)
  m <- cf_fit(Y ~ T | X1 + X2, data = df, method = "grf", random_state = 0L)

  m_bad <- m
  m_bad$res$ate_ci <- c(2, 1)  # reversed bounds
  expect_error(validate_cfomics_result(m_bad), class = "cfomics_invalid_result")

  m_bad2 <- m
  m_bad2$res$tau_ci <- data.frame(lo = 1:10)  # wrong columns
  expect_error(validate_cfomics_result(m_bad2), class = "cfomics_invalid_result")
})
