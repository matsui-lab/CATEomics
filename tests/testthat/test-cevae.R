.skip_if_no_pyro <- function() {
  # CEVAE runs inside the CATEomics basilisk environment.
  skip_if_no_cf_python()
}

test_that("cf_fit with method='cevae' runs on continuous Y and hits contract", {
  .skip_if_no_pyro()

  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.3 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "cevae",
          params = list(num_epochs = 5L, latent_dim = 4L, hidden_dim = 16L,
              batch_size = 50L))

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "cevae")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_length(ate, 1)

  ite <- predict(fit, type = "tau")
  expect_equal(length(ite), n)

  y0 <- predict(fit, type = "mu0")
  y1 <- predict(fit, type = "mu1")
  expect_equal(length(y0), n)
  expect_equal(length(y1), n)
  expect_equal(ite, y1 - y0, tolerance = 1e-6)

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("tau_mean" %in% names(summary_stats))
  expect_true("tau_iqr" %in% names(summary_stats))
})

test_that("predict_cf_cevae supports newdata for all types", {
  .skip_if_no_pyro()

  set.seed(321)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.3 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "cevae",
          params = list(num_epochs = 5L, latent_dim = 4L, hidden_dim = 16L,
              batch_size = 50L))

  n_new <- 15
  newdata <- data.frame(X1 = rnorm(n_new), X2 = rnorm(n_new))

  y0_new <- predict(fit, newdata = newdata, type = "mu0")
  y1_new <- predict(fit, newdata = newdata, type = "mu1")
  ite_new <- predict(fit, newdata = newdata, type = "tau")
  ate_new <- predict(fit, newdata = newdata, type = "ate")

  expect_length(y0_new, n_new)
  expect_length(y1_new, n_new)
  expect_length(ite_new, n_new)
  expect_equal(ite_new, y1_new - y0_new, tolerance = 1e-6)
  expect_equal(ate_new, mean(ite_new), tolerance = 1e-8)

  # type = "summary" is an in-sample concept; newdata does not support it
  # (parity with grf / metalearners newdata predict).

  expect_error(
    predict(fit, newdata = matrix(rnorm(10), ncol = 1), type = "tau"),
    "columns"
  )
})
