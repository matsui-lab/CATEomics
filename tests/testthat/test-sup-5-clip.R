# =============================================================================
# test-sup-5-clip.R
# =============================================================================
#
# Unified PS clip API: `clip = c(L, U)` is accepted across methods. PS-using
# methods (drlearner / xlearner / bcf) actually clip; PS-free methods silently
# ignore it (dispatch does not forward it, so no unused-argument error).
# The Python methods are skipped when their environment is unavailable.

test_that("cf_fit validates `clip` (length-2 with 0 < L < U < 1)", {
  skip_if_not_installed("grf")
  set.seed(42)
  df <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50), X2 = rnorm(50)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf",
            clip = c(0.1, 0.1)),
    class = "cateomics_invalid_clip"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf",
            clip = c(-0.1, 0.9)),
    class = "cateomics_invalid_clip"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf",
            clip = 0.5),
    class = "cateomics_invalid_clip"
  )
})

test_that("PS-free methods silently accept default `clip` (no error)", {
  skip_if_not_installed("grf")
  set.seed(7)
  df <- data.frame(
    Y  = rnorm(60),
    T  = rbinom(60, 1, 0.5),
    X1 = rnorm(60), X2 = rnorm(60)
  )
  # grf (PS-free CATE); should run unchanged.
  fit <- suppressMessages(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf",
            clip = c(0.05, 0.95), partition = "none")
  )
  expect_s3_class(fit, "cf_model")
})

test_that("bcf clip replaces the legacy hardcoded c(0.01, 0.99)", {
  skip_if_not_installed("bcf")
  skip_if_not_installed("glmnet")
  set.seed(11)
  df <- data.frame(
    Y  = rnorm(80),
    T  = rbinom(80, 1, 0.5),
    X1 = rnorm(80), X2 = rnorm(80)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", clip = c(0.05,
            0.95), random_state = 1L, params = list(nburn = 50L,
            nsim = 100L))
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
    path   = system.file("python", package = "CATEomics")
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
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "drlearner",
            clip = c(0.05, 0.95), random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
})

test_that("xlearner accepts clip arg without error (smoke)", {
  skip_if_no_cf_python()
  set.seed(21)
  df <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50), X2 = rnorm(50)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "xlearner",
            clip = c(0.05, 0.95), random_state = 42L)
  ))
  expect_s3_class(fit, "cf_model")
})


# bcf actually uses the clip and the known propensity it is given ---------
#
# The older tests here end at expect_s3_class(fit, "cf_model"), which a bcf
# that discarded both arguments would also satisfy. These observe the
# propensity the model ended up with.

.bcf_confounded <- function(n = 80L, seed = 11L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  T_vec <- rbinom(n, 1L, plogis(2 * x1))   # extreme propensities, so a clip bites
  data.frame(Y = x2 + T_vec * (1 + x1) + rnorm(n), T = T_vec, x1 = x1, x2 = x2)
}

.bcf_short <- list(bcf = list(nburn = 10L, nsim = 10L, n_chains = 1L))

test_that("bcf's stored propensity respects the clip it was given", {
  skip_if_not_installed("bcf")
  df <- .bcf_confounded()
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "bcf",
    params = .bcf_short, clip = c(0.2, 0.8), partition = "none",
    random_state = 1L))
  ps <- fit$res$ps
  expect_true(all(ps >= 0.2 & ps <= 0.8))
  # The clip has to bind, or the assertion above would hold for any bounds.
  expect_gt(sum(ps %in% c(0.2, 0.8)), 0L)

  wider <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "bcf",
    params = .bcf_short, clip = c(0.02, 0.98), partition = "none",
    random_state = 1L))
  expect_lt(min(wider$res$ps), 0.2)
})

test_that("bcf uses the known propensity it is handed, scalar or per sample", {
  skip_if_not_installed("bcf")
  df <- .bcf_confounded()
  scalar <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "bcf",
    params = .bcf_short, propensity_score = 0.5, partition = "none",
    random_state = 1L))
  expect_true(all(scalar$res$ps == 0.5))
  expect_equal(scalar$res$ps_type, "known")

  strata <- ifelse(df$x1 > 0, 0.7, 0.3)
  vec <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "bcf",
    params = .bcf_short, propensity_score = strata, partition = "none",
    random_state = 1L))
  expect_equal(unname(vec$res$ps), unname(strata))
  expect_equal(vec$res$ps_type, "known")
})
