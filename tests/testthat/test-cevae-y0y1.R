test_that("CEVAEExtended.predict_y returns numeric y0, y1 of correct length", {
  skip_if_no_cf_python()

  set.seed(42)
  n <- 60
  X <- matrix(rnorm(n * 2), ncol = 2)
  T <- rbinom(n, 1, 0.5)
  Y <- X[, 1] + T + rnorm(n)

  fit <- cf_fit_cevae(
    X = X, T = T, Y = Y,
    covariate_names = c("X1", "X2"),
    random_state = 7L,
    num_epochs = 3L,
    latent_dim = 4L,
    hidden_dim = 16L,
    batch_size = 30L
  )

  expect_s3_class(fit, "cf_model_cevae")
  expect_false(is.null(fit$model))
  expect_length(fit$res$mu0_hat, n)
  expect_length(fit$res$mu1_hat, n)
  expect_type(fit$res$mu0_hat, "double")
  expect_type(fit$res$mu1_hat, "double")
  expect_equal(fit$res$tau_hat, fit$res$mu1_hat - fit$res$mu0_hat, tolerance = 1e-6)
  expect_equal(fit$res$ate, mean(fit$res$tau_hat), tolerance = 1e-8)

  expect_false(is.null(fit$res$metadata$hyperparameters))
  expect_equal(fit$res$metadata$hyperparameters$outcome_dist, "normal")
})
