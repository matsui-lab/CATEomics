# Tests for Python method structure and error handling
# These tests focus on validation and structure without requiring GPU or full Python execution

# ============================================================================

# ============================================================================
# DRLearner Error Handling Tests
# ============================================================================

test_that("cf_fit with drlearner returns a cf_model when Python is available", {
  skip_if_no_cf_python()

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  result <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")
  expect_s3_class(result, "cf_model")
  expect_equal(result$method, "drlearner")
})

test_that("predict_cf_drlearner validates type argument", {
  skip_if_no_cf_python()

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  # Test all valid types
  expect_type(predict(fit, type = "ate"), "double")
  expect_type(predict(fit, type = "tau"), "double")
  expect_type(predict(fit, type = "mu0"), "double")
  expect_type(predict(fit, type = "mu1"), "double")
  expect_type(predict(fit, type = "summary"), "list")

  # Invalid type is rejected (mirrors the grf sibling below).
  expect_error(
    predict(fit, type = "invalid_type"),
    "should be one of"
  )
})

# ============================================================================
# GANITE Error Handling Tests
# ============================================================================

test_that("cf_fit with ganite returns a cf_model when Python is available", {
  skip_if_no_cf_python()

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  result <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite",
                   iterations = 10L, verbose = FALSE)
  expect_s3_class(result, "cf_model")
  expect_equal(result$method, "ganite")
})

# ============================================================================
# CEVAE Error Handling Tests
# ============================================================================

test_that("cf_fit with cevae returns a cf_model when Python is available", {
  skip_if_no_cf_python()

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  result <- cf_fit(Y ~ T | X1 + X2, data = data, method = "cevae",
                   num_epochs = 10L, verbose = FALSE)
  expect_s3_class(result, "cf_model")
  expect_equal(result$method, "cevae")
})

# ============================================================================
# Cross-method Tests
# ============================================================================

test_that("cf_fit method auto-selection prefers R methods when Python unavailable", {
  # When Python is not available, auto should select an R method
  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Auto method should work regardless of Python availability
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "auto")

  expect_s3_class(fit, "cf_model")

  # "prefer R" contract: auto selects an R-native method (grf/bcf) here.
  expect_true(fit$method %in% c("grf", "bcf"))
})

test_that("all methods return consistent result structure", {
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Test with grf (always available when installed)
  fit_grf <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

  # Check structure (four-pillar)
  expect_true("method" %in% names(fit_grf))
  expect_true("res" %in% names(fit_grf))
  expect_true("model" %in% names(fit_grf))
  expect_true("meta" %in% names(fit_grf))

  # Check res structure
  res <- fit_grf$res
  expect_true("ate" %in% names(res))
  expect_true("tau_hat" %in% names(res))
  expect_true("mu0_hat" %in% names(res))
  expect_true("mu1_hat" %in% names(res))
  expect_true("summary" %in% names(res))

  # Check meta structure
  meta <- fit_grf$meta
  expect_true("n" %in% names(meta))
  expect_true("p" %in% names(meta))
})

test_that("Python method result structure matches R method structure", {
  skip_if_no_cf_python()

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit_r <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf", partition = "none")
  fit_py <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  # R-native partition = "none" fits do not populate `oof_ite` (no K-fold
  # cross-fitting), whereas the Python oof path does. Compare the core contract
  # fields shared across the result object, not strict `names()` equality.
  shared <- c("method", "res", "model", "meta", "data")
  expect_true(all(shared %in% names(fit_r)))
  expect_true(all(shared %in% names(fit_py)))

  # Both should have same result names
  expect_true(all(c("ate", "tau_hat", "mu0_hat", "mu1_hat") %in% names(fit_r$res)))
  expect_true(all(c("ate", "tau_hat", "mu0_hat", "mu1_hat") %in% names(fit_py$res)))

  # Both should have same meta names (at minimum)
  common_meta <- c("n", "p")
  expect_true(all(common_meta %in% names(fit_r$meta)))
  expect_true(all(common_meta %in% names(fit_py$meta)))
})

# ============================================================================
# Prediction Type Tests for All Methods
# ============================================================================

test_that("predict types work for grf method", {
  skip_if_not_installed("grf")
  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf", partition = "none")

  # Test all prediction types
  ate <- predict(fit, type = "ate")
  expect_length(ate, 1)
  expect_type(ate, "double")

  ite <- predict(fit, type = "tau")
  expect_length(ite, n)
  expect_type(ite, "double")

  y0 <- predict(fit, type = "mu0")
  expect_length(y0, n)

  y1 <- predict(fit, type = "mu1")
  expect_length(y1, n)

  summary <- predict(fit, type = "summary")
  expect_type(summary, "list")
  expect_true("ate" %in% names(summary))
})

test_that("predict validates type argument", {
  skip_if_not_installed("grf")
  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none")

  expect_error(
    predict(fit, type = "invalid_type"),
    "should be one of"
  )
})
