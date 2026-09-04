# =============================================================================
# test-sup-5-6-integration.R
# =============================================================================
#
# End-to-end chain of trim / known PS + metadata 4-bucket /
# compute_diagnostics / omics-integrity S3 methods, smoke-tested through a
# single cf_fit chain.

test_that("cf_fit end-to-end populates meta buckets + trim + propensity_score", {
  skip_if_not_installed("bcf")
  set.seed(31)
  n <- 80L
  df <- data.frame(
    Y  = rnorm(n) + 0.3 * rbinom(n, 1, 0.5),
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "bcf",
            partition = "none",
            trim = c(0.05, 0.95),
            propensity_score = 0.5,
            n_burn = 50L, n_iter = 100L,
            random_state = 42L)
  ))
  # Top-level identifiers preserved.
  expect_identical(fit$meta$method, "bcf")
  # The meta buckets present are estimand / provenance / estimator.
  expect_true(all(c("estimand", "provenance", "estimator") %in%
                    names(fit$meta)))
  expect_identical(fit$meta$estimand$type, "CATE")
})

test_that("on-demand cf_diagnostics runs after a bcf fit", {
  skip_if_not_installed("bcf")
  set.seed(32)
  n <- 60L
  df <- data.frame(
    Y  = rnorm(n) + 0.4 * rbinom(n, 1, 0.5),
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "bcf",
            partition = "none", n_burn = 50L, n_iter = 100L,
            random_state = 42L)
  ))
  dg <- cf_diagnostics(fit)
  expect_s3_class(dg, "cfomics_diagnostics")
  expect_true(isTRUE(dg$evalue$available))
})

test_that("Propensity-args conflict error fires", {
  skip_if_not_installed("bcf")
  set.seed(34)
  df <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40), X2 = rnorm(40)
  )
  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "bcf",
            propensity_score = 0.5,
            model_propensity = "rf",
            partition = "none"),
    class = "cfomics_propensity_args_conflict"
  )
})
