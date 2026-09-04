test_that("cf_fit_drlearner basic execution", {
  skip_if_no_cf_python()

  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "drlearner")

  # The model pillar holds a serialised model path, not a live object.
  expect_true(file.exists(fit$model$model$model_path))

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  # True ATE is 1.5 in this DGP (mirrors the grf sibling in test-traditional.R).
  expect_true(ate > 0.5 && ate < 2.5)

  ite <- predict(fit, type = "tau")
  expect_equal(length(ite), n)

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("tau_mean" %in% names(summary_stats))
})

test_that("drlearner prediction on new data", {
  skip_if_no_cf_python()

  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  newdata <- data.frame(X1 = rnorm(10), X2 = rnorm(10))

  ite_new <- predict_cf_drlearner(fit, newdata = newdata, type = "tau")
  expect_equal(length(ite_new), 10)
})

test_that("drlearner y0/y1 are direct nuisance estimates (Option B)", {
  skip_if_no_cf_python()

  set.seed(456)
  n <- 200
  X1 <- rnorm(n); X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner",
                random_state = 0L)

  ite <- predict(fit, type = "tau")
  y0  <- predict(fit, type = "mu0")
  y1  <- predict(fit, type = "mu1")

  # By design: y0/y1 are μ_0, μ_1 from models_regression, while ite is the
  # orthogonalized CATE. They are NOT required to satisfy y1-y0 == ite.
  expect_length(y0, n)
  expect_length(y1, n)
  expect_true(all(is.finite(y0)))
  expect_true(all(is.finite(y1)))

  md <- fit$res$metadata
  expect_equal(md$backend, "econml.dr.DRLearner")
  expect_true(nzchar(md$note))
})
