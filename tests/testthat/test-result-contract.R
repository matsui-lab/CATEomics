# Test file for result contract validation
# Task 2.1: Verify cfomics_result contains estimand field

test_that("cfomics_result contains estimand field", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf")

  expect_true("estimand" %in% names(fit$meta))
  # All methods are CATE-capable; grf reports the CATE estimand (ATE is a
  # derived summary).
  expect_equal(fit$meta$estimand$type, "CATE")
})

test_that("create_cf_meta includes estimand with default value", {
  meta <- create_cf_meta(
    formula = Y ~ T | X1,
    n = 100,
    p = 1,
    outcome_name = "Y",
    treatment_name = "T",
    covariate_names = "X1",
    method = "grf"
  )

  expect_true("estimand" %in% names(meta))
  expect_equal(meta$estimand$type, "ATE")
})

test_that("create_cf_meta accepts custom estimand", {
  meta <- create_cf_meta(
    formula = Y ~ T | X1,
    n = 100,
    p = 1,
    outcome_name = "Y",
    treatment_name = "T",
    covariate_names = "X1",
    method = "grf",
    estimand = "ATT"
  )

  expect_true("estimand" %in% names(meta))
  expect_equal(meta$estimand$type, "ATT")
})

test_that("validate_cfomics_result passes with valid structure", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf")

  # Should not throw error

  expect_true(validate_cfomics_result(fit))
})
