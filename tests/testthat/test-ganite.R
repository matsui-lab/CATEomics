test_that("cf_fit with method='ganite' basic execution", {
  skip_if_no_cf_python()


  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- as.integer(X1 > 0)
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)


  # Fit with GANITE method (reduced iterations for a fast smoke test)
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite",
                iterations = 20L, batch_size = 32L, verbose = FALSE)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "ganite")

  # Predict ATE
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_length(ate, 1)

  # Predict ITE
  ite <- predict(fit, type = "tau")
  expect_type(ite, "double")
  expect_equal(length(ite), n)

  # Predict Y0 and Y1
  y0 <- predict(fit, type = "mu0")
  y1 <- predict(fit, type = "mu1")
  expect_equal(length(y0), n)
  expect_equal(length(y1), n)

  # Verify ITE = Y1 - Y0
  expect_equal(ite, y1 - y0, tolerance = 1e-6)
})

test_that("cf_fit_ganite returns correct structure with normalization metadata", {
  skip_if_no_cf_python()

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(456)
  n <- 50
  X <- matrix(rnorm(n * 3), ncol = 3)
  T <- rbinom(n, 1, 0.5)
  Y <- X[, 1] + 0.5 * X[, 2] + T + rnorm(n)

  fit <- cf_fit_ganite(
    X = X, T = T, Y = Y,
    covariate_names = c("X1", "X2", "X3"),
    random_state = 42L,
    iterations = 20L,
    batch_size = 16L,
    verbose = FALSE
  )

  expect_s3_class(fit, "cf_model_ganite")
  # The raw method returns a two-pillar (model + res) list; the
  # model pillar holds the serialised path + backend, and cf_fit() adds the
  # central summary (summary/samples are not on the raw method output).
  expect_equal(fit$model$backend, "ganite.Ganite")
  expect_true(file.exists(fit$model$model_path))
  expect_true(!is.null(fit$res$ate))
  expect_true(!is.null(fit$res$tau_hat))
  expect_true(!is.null(fit$res$mu0_hat))
  expect_true(!is.null(fit$res$mu1_hat))

  # Y normalization must persist for newdata scoring (X is standardised
  # centrally and re-applied via meta$x_zscore). The raw method records the Y
  # scale as res$y_zscore (promoted to meta$y_zscore by cf_fit()).
  ystats <- fit$res$y_zscore
  expect_false(is.null(ystats))
  expect_true(is.numeric(ystats$Y_mean))
  expect_true(is.numeric(ystats$Y_sd))
  expect_true(ystats$Y_sd > 0)

  # Hyperparameters (alpha/beta defaults) recorded.
  hp <- fit$res$metadata$hyperparameters
  expect_equal(hp$alpha, 2.0)
  expect_equal(hp$beta, 5.0)
})

test_that("predict_cf_ganite returns correct types (in-sample)", {
  skip_if_no_cf_python()

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(789)
  n <- 30
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- as.integer(X1 > 0)
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)


  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite",
                iterations = 20L, batch_size = 16L, verbose = FALSE)

  # Test summary type
  summary_result <- predict(fit, type = "summary")
  expect_type(summary_result, "list")
  expect_true("ate" %in% names(summary_result))
  expect_true("tau_mean" %in% names(summary_result))
  expect_true("tau_sd" %in% names(summary_result))
  expect_true("tau_iqr" %in% names(summary_result))

  # There is no "samples" predict token; mu0 / mu1 / tau are
  # separate types now.
  expect_equal(length(predict(fit, type = "mu0")), nrow(data))
  expect_equal(length(predict(fit, type = "mu1")), nrow(data))
  expect_equal(length(predict(fit, type = "tau")), nrow(data))
})

test_that("predict_cf_ganite supports newdata via saved inference network", {
  skip_if_no_cf_python()

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  set.seed(321)
  n <- 60
  X <- matrix(rnorm(n * 2), ncol = 2)
  T <- rbinom(n, 1, 0.5)
  Y <- X[, 1] + 0.5 * X[, 2] + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X[, 1], X2 = X[, 2])

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite",
                iterations = 20L, batch_size = 16L, verbose = FALSE,
                random_state = 7L)

  # Score fresh data via the saved Inference network.
  n_new <- 15
  X_new <- data.frame(X1 = rnorm(n_new), X2 = rnorm(n_new))

  y0_new <- predict(fit, newdata = X_new, type = "mu0")
  y1_new <- predict(fit, newdata = X_new, type = "mu1")
  ite_new <- predict(fit, newdata = X_new, type = "tau")

  expect_type(y0_new, "double")
  expect_type(y1_new, "double")
  expect_type(ite_new, "double")
  expect_length(y0_new, n_new)
  expect_length(y1_new, n_new)
  expect_length(ite_new, n_new)

  # Additivity contract within numerical tolerance.
  expect_equal(ite_new, y1_new - y0_new, tolerance = 1e-6)

  # ate on newdata is a scalar mean(ite_new).
  ate_new <- predict(fit, newdata = X_new, type = "ate")
  expect_type(ate_new, "double")
  expect_length(ate_new, 1)
  expect_equal(ate_new, mean(ite_new), tolerance = 1e-8)

  # Dimensionality mismatch must error.
  expect_error(
    predict(fit, newdata = matrix(rnorm(10), ncol = 1), type = "tau"),
    "columns"
  )
})
