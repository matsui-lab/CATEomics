# Test file for result contract validation
# The meta contract: which fields a result carries, and what the validator
# refuses.

test_that("the result carries the contract's meta fields", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf")

  expect_true(all(c("n", "p", "covariate_names", "provenance", "estimator") %in%
                    names(fit$meta)))
  # The variables the caller named are recorded once, under named_args.
  expect_equal(fit$meta$named_args$Y, "Y")
  expect_equal(fit$meta$named_args$T, "T")
  expect_equal(fit$meta$named_args$X_baseline, "X1")
})

test_that("create_cf_meta assembles the buckets", {
  meta <- create_cf_meta(
    n = 100,
    p = 1,
    covariate_names = "X1",
    method = "grf"
  )

  expect_true(all(c("provenance", "estimator") %in% names(meta)))
  expect_identical(meta$provenance$backend, "grf (R)")
  # Every method here estimates a conditional effect, so there is nothing for
  # an estimand field to vary over: a constant would record nothing.
  expect_null(meta$estimand)
})

test_that("a result with no per-unit estimate is refused", {
  # Every method estimates a conditional effect, so a result whose tau_hat is
  # entirely missing is not one this contract can describe. This is the check
  # that used to sit behind the estimand field.
  skip_if_not_installed("grf")
  set.seed(1)
  d <- data.frame(Y = stats::rnorm(40), T = stats::rbinom(40, 1, 0.5),
                  x1 = stats::rnorm(40), x2 = stats::rnorm(40))
  fit <- suppressMessages(cf_fit(data = d, Y = "Y", T = "T",
                                 X_baseline = c("x1", "x2"), method = "grf",
                                 preprocess = FALSE))
  expect_error(validate_cateomics_result(fit), NA)
  bad <- fit
  bad$res$tau_hat <- rep(NA_real_, length(bad$res$tau_hat))
  expect_error(validate_cateomics_result(bad),
               class = "cateomics_invalid_result")
})


test_that("validate_cateomics_result passes with valid structure", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf")

  # Should not throw error

  expect_true(validate_cateomics_result(fit))
})

# The interval that does not contain its estimate -------------------------

test_that("a fit says so when its interval excludes the point estimate", {
  warn_of <- function(res) {
    w <- NULL
    withCallingHandlers(
      CATEomics:::.cf_warn_ci_excludes_estimate(res),
      cateomics_bootstrap_ci_excludes_estimate = function(cnd) {
        w <<- cnd; invokeRestart("muffleWarning")
      })
    w
  }
  boot <- function(ate, lo, hi)
    list(res = list(ate = ate, ate_ci = c(lo, hi),
                    ci_method = "nonparametric bootstrap"))

  # The estimate above its own interval, and below it.
  w <- warn_of(boot(0.56, 0.28, 0.52))
  expect_s3_class(w, "cateomics_bootstrap_ci_excludes_estimate")
  expect_match(rlang::cnd_message(w), "systematically below", fixed = TRUE)
  expect_match(rlang::cnd_message(w), "0.56", fixed = TRUE)
  w_lo <- warn_of(boot(0.10, 0.28, 0.52))
  expect_match(rlang::cnd_message(w_lo), "systematically above", fixed = TRUE)

  # Silent when the estimate is inside, and when there is no interval.
  expect_null(warn_of(boot(0.40, 0.28, 0.52)))
  expect_null(warn_of(list(res = list(ate = 0.4, ate_ci = NULL))))
  expect_null(warn_of(list(res = list(ate = 0.4, ate_ci = c(NA_real_, NA_real_)))))

  # Both shapes a fitter can return reach the same check.
  expect_s3_class(warn_of(list(ate = 0.56, ate_ci = c(0.28, 0.52),
                               ci_method = "nonparametric bootstrap")),
                  "cateomics_bootstrap_ci_excludes_estimate")
})
