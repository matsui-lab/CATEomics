# =============================================================================
# test-sup-5-trim.R
# =============================================================================
#
# Unified PS trim API: `trim = c(L, U)` is accepted across methods. PS-using
# methods (drlearner / xlearner / bcf) actually clip; PS-free methods silently
# ignore it (dispatch does not forward it, so no unused-argument error).
# Python integration is gated by cf_has_python.

test_that("cf_fit validates `trim` (length-2 with 0 < L < U < 1)", {
  skip_if_not_installed("grf")
  set.seed(42)
  df <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50), X2 = rnorm(50)
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "grf",
            trim = c(0.1, 0.1)),
    class = "cfomics_invalid_trim"
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "grf",
            trim = c(-0.1, 0.9)),
    class = "cfomics_invalid_trim"
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "grf",
            trim = 0.5),
    class = "cfomics_invalid_trim"
  )
})

test_that("PS-free methods silently accept default `trim` (no error)", {
  skip_if_not_installed("grf")
  set.seed(7)
  df <- data.frame(
    Y  = rnorm(60),
    T  = rbinom(60, 1, 0.5),
    X1 = rnorm(60), X2 = rnorm(60)
  )
  # grf (PS-free CATE); should run unchanged.
  fit <- suppressMessages(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "grf",
            trim = c(0.05, 0.95), partition = "none")
  )
  expect_s3_class(fit, "cf_model")
})

test_that("bcf trim replaces the legacy hardcoded c(0.01, 0.99)", {
  skip_if_not_installed("bcf")
  skip_if_not_installed("glmnet")
  set.seed(11)
  df <- data.frame(
    Y  = rnorm(80),
    T  = rbinom(80, 1, 0.5),
    X1 = rnorm(80), X2 = rnorm(80)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "bcf",
            trim = c(0.05, 0.95), n_burn = 50L, n_iter = 100L,
            random_state = 1L)
  ))
  expect_s3_class(fit, "cf_model")
})

test_that("ClippedPropensityModel clips predict_proba to [L, U]", {
  skip_if_no_cf_python()
  np  <- reticulate::import("numpy", delay_load = FALSE)
  skl <- reticulate::import("sklearn.linear_model", delay_load = FALSE)
  set.seed(99)
  n <- 60L
  X <- matrix(rnorm(n * 2), n, 2)
  y <- rbinom(n, 1, plogis(X[, 1]))
  clipped_mod <- reticulate::import_from_path(
    module = "clipped_propensity",
    path   = system.file("python", package = "cfomics")
  )
  base <- skl$LogisticRegression()
  wrapped <- clipped_mod$ClippedPropensityModel(base, lower = 0.1, upper = 0.9)
  wrapped$fit(X, y)
  proba <- reticulate::py_to_r(wrapped$predict_proba(X))
  expect_true(all(proba >= 0.1 - 1e-12))
  expect_true(all(proba <= 0.9 + 1e-12))
  # Rows still sum to ~1.
  expect_true(all(abs(rowSums(proba) - 1) < 1e-9))
})

test_that("drlearner applies ClippedPropensityModel wrap (smoke)", {
  skip_if_no_cf_python()
  set.seed(13)
  df <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50), X2 = rnorm(50)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "drlearner",
            trim = c(0.05, 0.95), random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
})

test_that("xlearner accepts trim arg without error (smoke)", {
  skip_if_no_cf_python()
  set.seed(21)
  df <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50), X2 = rnorm(50)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "xlearner",
            trim = c(0.05, 0.95), random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
})
