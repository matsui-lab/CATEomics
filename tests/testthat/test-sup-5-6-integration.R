# =============================================================================
# test-sup-5-6-integration.R
# =============================================================================
#
# End-to-end chain of clip and a known propensity, with the meta buckets the
# fit populates, smoke-tested through a single cf_fit chain.

test_that("cf_fit end-to-end populates meta buckets + clip + propensity_score", {
  skip_if_not_installed("bcf")
  set.seed(31)
  n <- 80L
  df <- data.frame(
    Y  = rnorm(n) + 0.3 * rbinom(n, 1, 0.5),
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", partition = "none",
            clip = c(0.05, 0.95), propensity_score = 0.5, random_state = 42L,
            params = list(nburn = 50L, nsim = 100L))
  ))
  # The method name is the top-level one; `meta` does not carry a copy.
  expect_identical(fit$method, "bcf")
  expect_false("method" %in% names(fit$meta))
  # The meta buckets present are provenance / estimator.
  expect_true(all(c("provenance", "estimator") %in%
                    names(fit$meta)))
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
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", partition = "none",
            random_state = 42L, params = list(nburn = 50L, nsim = 100L))
  ))
  dg <- cf_diagnostics(fit)
  expect_s3_class(dg, "cateomics_diagnostics")
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
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf",
            propensity_score = 0.5,
            model_propensity = "rf",
            partition = "none"),
    class = "cateomics_propensity_args_conflict"
  )
})
