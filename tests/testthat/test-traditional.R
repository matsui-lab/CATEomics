test_that("cf_fit_grf basic execution", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "grf")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 0.5 && ate < 2.5)

  ite <- predict(fit, type = "tau")
  expect_equal(length(ite), n)

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
})

test_that("cf_fit_grf y0/y1 use Y.hat - W.hat*tau form (B3)", {
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 200
  X1 <- rnorm(n); X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf",
                random_state = 0L)

  tau_forest <- fit$model$model
  Y_hat <- as.numeric(tau_forest$Y.hat)
  W_hat <- as.numeric(tau_forest$W.hat)
  ite   <- predict(fit, type = "tau")

  expected_y0 <- Y_hat - W_hat * ite
  expected_y1 <- Y_hat + (1 - W_hat) * ite

  expect_equal(predict(fit, type = "mu0"), expected_y0, tolerance = 1e-10)
  expect_equal(predict(fit, type = "mu1"), expected_y1, tolerance = 1e-10)

  # Sanity: the new form differs from the legacy observed-value form on
  # treated units (old: y0 = Y - ite; new: y0 = Y.hat - W.hat*ite).
  treated <- which(data$T == 1)
  legacy_y0 <- Y[treated] - ite[treated]
  new_y0    <- predict(fit, type = "mu0")[treated]
  expect_false(isTRUE(all.equal(new_y0, legacy_y0, tolerance = 1e-6)))
})

test_that("cf_fit_grf restores .Random.seed after fitting", {
  skip_if_not_installed("grf")

  set.seed(777)
  pre_seed <- .Random.seed

  n <- 100
  X1 <- rnorm(n); X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- T + X1 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  set.seed(777)  # reset before the call
  pre_seed <- .Random.seed
  cf_fit(Y ~ T | X1 + X2, data = data, method = "grf",
         partition = "none",                   # at n=100 the K=5 default oof
         random_state = 1L, num.trees = 200L)  # gives insufficient grf
                                                # honesty; use none. This test
                                                # (.Random.seed restore) is
                                                # independent of partition.

  # After cf_fit_grf returns, .GlobalEnv RNG state must match pre-call state.
  expect_identical(.Random.seed, pre_seed)
})
