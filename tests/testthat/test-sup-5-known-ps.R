# =============================================================================
# test-sup-5-known-ps.R
# =============================================================================
#
# Known PS: cf_fit(..., propensity_score = ...) handles the RCT context where
# the PS is known. Supports both scalar (50/50 RCT) and vector (stratified
# RCT). PS-using methods use it; PS-free methods silently ignore it (dispatch
# does not forward it).

test_that("cf_fit validates `propensity_score` range and length", {
  skip_if_not_installed("bcf")
  set.seed(11)
  df <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40), X2 = rnorm(40)
  )
  # Out of range.
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf",
            propensity_score = 1.5, partition = "none"),
    class = "cateomics_invalid_propensity_score"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf",
            propensity_score = 0, partition = "none"),
    class = "cateomics_invalid_propensity_score"
  )
  # Length mismatch.
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf",
            propensity_score = rep(0.5, 5), partition = "none"),
    class = "cateomics_invalid_propensity_score"
  )
})

test_that("cf_fit raises conflict when both propensity_score and model_propensity set", {
  skip_if_not_installed("bcf")
  set.seed(13)
  df <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40), X2 = rnorm(40)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf",
            propensity_score = 0.5, model_propensity = "rf",
            partition = "none"),
    class = "cateomics_propensity_args_conflict"
  )
})

test_that("bcf with vector known PS (stratified RCT)", {
  skip_if_not_installed("bcf")
  set.seed(202)
  n <- 100L
  strata <- rbinom(n, 1, 0.5)
  ps_strata <- ifelse(strata == 1L, 0.7, 0.3)
  Tr <- rbinom(n, 1, ps_strata)
  Y  <- 0.4 * Tr + 0.2 * strata + rnorm(n)
  df <- data.frame(Y = Y, T = Tr, X1 = rnorm(n), X2 = rnorm(n))
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", propensity_score = ps_strata,
            partition = "none", random_state = 42L, params = list(nburn = 50L,
                nsim = 100L))
  ))
  expect_s3_class(fit, "cf_model")
})

test_that("bcf accepts known PS via dispatch (pihat forward)", {
  skip_if_not_installed("bcf")
  set.seed(303)
  n <- 80L
  Tr <- rbinom(n, 1, 0.5)
  Y  <- 0.3 * Tr + rnorm(n)
  df <- data.frame(Y = Y, T = Tr, X1 = rnorm(n), X2 = rnorm(n))
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", propensity_score = 0.5,
            partition = "none", random_state = 1L, params = list(nburn = 50L,
                nsim = 100L))
  ))
  expect_s3_class(fit, "cf_model")
})

test_that("drlearner with known PS (smoke) — scalar only via cross-fit constraint", {
  skip_if_no_cf_python()
  set.seed(505)
  n <- 60L
  Tr <- rbinom(n, 1, 0.5)
  Y  <- 0.5 * Tr + rnorm(n)
  df <- data.frame(Y = Y, T = Tr, X1 = rnorm(n), X2 = rnorm(n))
  # scalar OK
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "drlearner",
            propensity_score = 0.5, partition = "none",
            random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
  # vector unsupported for cross-fit methods
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "drlearner",
            propensity_score = rep(0.5, n), partition = "none"),
    class = "cateomics_propensity_vector_unsupported"
  )
})

test_that("xlearner rejects vector PS (cross-fit constraint)", {
  set.seed(606)
  n <- 60L
  Tr <- rbinom(n, 1, 0.5)
  Y  <- 0.4 * Tr + rnorm(n)
  df <- data.frame(Y = Y, T = Tr, X1 = rnorm(n), X2 = rnorm(n))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "xlearner",
            propensity_score = rep(0.5, n), partition = "none"),
    class = "cateomics_propensity_vector_unsupported"
  )
})

test_that("xlearner with scalar known PS via FixedPropensityModel (smoke)", {
  skip_if_no_cf_python()
  set.seed(707)
  n <- 60L
  Tr <- rbinom(n, 1, 0.5)
  Y  <- 0.5 * Tr + rnorm(n)
  df <- data.frame(Y = Y, T = Tr, X1 = rnorm(n), X2 = rnorm(n))
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "xlearner",
            propensity_score = 0.5, partition = "none",
            random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
})
