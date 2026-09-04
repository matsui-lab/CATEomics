test_that("cf_fit_tarnet smoke (catenets path)", {
  skip_if_no_cf_python()

  set.seed(123)
  n <- 80
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  T <- as.integer(rbinom(n, 1, 0.5))
  Y <- as.numeric(X[, 1] + 0.5 * T + rnorm(n))
  df <- data.frame(Y = Y, T = T, X)
  names(df)[3:7] <- paste0("X", 1:5)

  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4 + X5,
                data = df, method = "tarnet",
                random_state = 42L, n_iter = 20L)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_false(is.na(ate))

  ite <- predict(fit, type = "tau")
  expect_length(ite, n)
  expect_false(any(is.na(ite)))

  mu0 <- fit$res$mu0_hat
  mu1 <- fit$res$mu1_hat
  expect_length(mu0, n)
  expect_length(mu1, n)
  expect_false(any(is.na(mu0)))
  expect_false(any(is.na(mu1)))

  # Class A4: mu1 - mu0 == ite (plug-in CATE)
  expect_equal(mu1 - mu0, ite, tolerance = 1e-5)
})

test_that("cf_fit_tarnet reproducibility (RNG bridge)", {
  skip_if_no_cf_python()

  set.seed(0)
  n <- 60
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  T <- as.integer(rbinom(n, 1, 0.5))
  Y <- as.numeric(X[, 1] + 0.5 * T + rnorm(n))
  df <- data.frame(Y = Y, T = T, X)
  names(df)[3:6] <- paste0("X", 1:4)

  fit1 <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = df, method = "tarnet",
                 random_state = 7L, n_iter = 10L)
  fit2 <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = df, method = "tarnet",
                 random_state = 7L, n_iter = 10L)

  # Wrapper saves/restores global RNGs; same seed -> bit-exact same CATE.
  expect_equal(fit1$res$tau_hat, fit2$res$tau_hat, tolerance = 1e-6)
})
