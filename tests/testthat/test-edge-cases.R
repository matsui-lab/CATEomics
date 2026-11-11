# test-edge-cases.R - Edge case tests and input validation for cf_fit()

test_that("cf_fit handles single covariate", {
  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  expect_no_error(cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"))
})

test_that("cf_fit handles small sample sizes (minimum 10)", {
  set.seed(123)
  # 20 observations should work
  data <- data.frame(Y = rnorm(20), T = c(rep(0, 10), rep(1, 10)), X1 = rnorm(20))
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
})

test_that("cf_fit rejects sample sizes below minimum", {
  set.seed(123)
  # 9 observations should fail (minimum is 10)
  data <- data.frame(Y = rnorm(9), T = c(rep(0, 5), rep(1, 4)), X1 = rnorm(9))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "Insufficient sample size"
  )
})

test_that("cf_fit validates binary treatment", {
  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rnorm(50), X1 = rnorm(50))  # T not binary
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "binary"
  )
})

test_that("cf_fit validates numeric outcome", {
  set.seed(123)
  # Create a factor outcome - won't be converted to NA like character
  data <- data.frame(
    Y = factor(sample(c("low", "high"), 50, replace = TRUE)),
    T = rbinom(50, 1, 0.5),
    X1 = rnorm(50),
    stringsAsFactors = FALSE
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "numeric"
  )
})

test_that("cf_fit requires both treatment groups", {
  set.seed(123)
  # All treated
  data <- data.frame(Y = rnorm(50), T = rep(1, 50), X1 = rnorm(50))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "treatment and control"
  )

  # All control
  data2 <- data.frame(Y = rnorm(50), T = rep(0, 50), X1 = rnorm(50))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data2, method = "grf", partition = "none"),
    "treatment and control"
  )
})

test_that("cf_fit handles multiple covariates", {
  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n),
    X4 = rnorm(n),
    X5 = rnorm(n)
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3", "X4", "X5"), data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$meta$p, 5)
})

test_that("cf_fit auto-expands factor covariates", {
  # factor covariates are auto one-hot encoded via model.matrix(~ . - 1). The
  # train levels are recorded in meta$factor_levels and re-applied at newdata
  # predict time.
  set.seed(123)
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = factor(sample(c("A", "B", "C"), 50, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
  expect_setequal(fit$meta$factor_levels$X1, c("A", "B", "C"))
})

test_that("cf_fit rejects NA values", {
  set.seed(123)
  # NA in outcome
  data <- data.frame(Y = c(NA, rnorm(49)), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "NA values"
  )

  # NA in treatment
  data2 <- data.frame(Y = rnorm(50), T = c(NA, rbinom(49, 1, 0.5)), X1 = rnorm(50))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data2, method = "grf", partition = "none"),
    "NA values"
  )

  # NA in covariate
  data3 <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = c(NA, rnorm(49)))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data3, method = "grf", partition = "none"),
    "NA values"
  )
})

test_that("cf_fit handles logical treatment (converts to integer)", {
  set.seed(123)
  data <- data.frame(
    Y = rnorm(50),
    T = sample(c(TRUE, FALSE), 50, replace = TRUE),
    X1 = rnorm(50)
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
})

test_that("cf_fit handles highly imbalanced treatment", {
  set.seed(123)
  # 95% treated, 5% control (but still valid)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = c(rep(1, 95), rep(0, 5)),
    X1 = rnorm(n)
  )
  # The design runs, but the evaluation nuisance has to be one the 5 control
  # rows can support: the default library has 7 members, and non-negative least
  # squares does not determine 7 weights from 5 rows.
  fit <- suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf",
                                 partition = "none",
                                 eval_library = c("mean", "ridge")))
  expect_s3_class(fit, "cf_model")
  expect_equal(sum(fit$meta$estimator$eval_sl_weights$mu0), 1, tolerance = 1e-8)

  # With the default library the same data is refused, and the message names
  # the size that would work. Before this was checked, the weights were
  # reported as if they meant something: measured on 3 rows against 7 members,
  # the combination fitted those rows worse than their own mean.
  expect_error(
    suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf",
                            partition = "none")),
    class = "cateomics_sl_library_too_large")
  # A single learner needs no weights and is unaffected.
  expect_s3_class(
    suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf",
                            partition = "none", eval_nuisance = "rf")),
    "cf_model")
})

test_that("cf_fit returns correct result structure", {
  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "grf", partition = "none")

  # Check class

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")

  # Check structure

  expect_true("method" %in% names(fit))
  expect_true("res" %in% names(fit))
  expect_true("model" %in% names(fit))
  expect_true("meta" %in% names(fit))

  # Check meta fields
  expect_equal(fit$meta$n, n)
  expect_equal(fit$meta$p, 2)
  expect_equal(fit$meta$named_args$Y, "Y")
  expect_equal(fit$meta$named_args$T, "T")
  expect_equal(fit$meta$covariate_names, c("X1", "X2"))

  # Check results (four-pillar)
  expect_true("ate" %in% names(fit$res))
  expect_true("tau_hat" %in% names(fit$res))
  expect_length(fit$res$tau_hat, n)
})

test_that("cf_fit handles empty data gracefully", {
  data <- data.frame(Y = numeric(0), T = integer(0), X1 = numeric(0))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "no observations|at least 10"
  )
})

test_that("cf_fit validates variable name constraints", {
  set.seed(123)
  # `make.names()`-valid identifiers are accepted.
  # `X.1` (period) is now valid; only names that fail `make.names()`
  # (leading digit, hyphen, whitespace, ...) are rejected.
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    `X.1` = rnorm(50),  # period is a valid name char
    check.names = FALSE
  )
  expect_no_error(
    suppressMessages(cf_fit(Y = "Y", T = "T", X_baseline = "X.1", data = data, method = "grf", partition = "none"))
  )
})
