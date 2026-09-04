test_that("cf_fit_bcf basic execution", {
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()

  set.seed(123)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[, 1])
  T <- rbinom(n, 1, ps)
  tau_true <- 2.0 + 0.5 * X[, 1]
  Y <- tau_true * T + 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  tree_dir <- tempfile("cfomics_bcf_test_")
  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 100, n_iter = 200, n_chains = 1L,
                save_tree_directory = tree_dir)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "bcf")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  # True ATE is approximately mean(2.0 + 0.5 * X[,1]) ~ 2.0
  expect_true(ate > 1.0 && ate < 3.5)

  ite <- predict(fit, type = "tau")
  expect_equal(length(ite), n)

  # Metadata should record save_tree_directory.
  expect_equal(fit$res$metadata$save_tree_directory, tree_dir)
})

# bcf posterior draws (res$samples) are dropped after CI computation, so there
# is no "returns posterior samples" test.

test_that("cf_fit_bcf provides credible intervals", {
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()

  set.seed(789)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[, 1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 50, n_iter = 100, n_chains = 1L,
                save_tree_directory = tempfile("cfomics_bcf_test_"))

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
  expect_true(summary_stats$ate_ci_lower < summary_stats$ate)
  expect_true(summary_stats$ate_ci_upper > summary_stats$ate)
})

test_that("cf_fit_bcf supports newdata prediction via save_tree_directory", {
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()
  skip_if_not_installed("glmnet")

  set.seed(2024)
  n <- 120
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, plogis(0.3 * X[, 1]))
  Y <- 1.5 * T + 0.5 * X[, 1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  tree_dir <- tempfile("cfomics_bcf_newdata_")
  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 50, n_iter = 100, n_chains = 1L,
                save_tree_directory = tree_dir)

  n_new <- 15
  X_new <- matrix(rnorm(n_new * p), n_new, p)
  colnames(X_new) <- cov_names
  pihat_new <- pmax(pmin(plogis(0.3 * X_new[, 1]), 0.99), 0.01)

  ite_new <- predict_cf_bcf(fit, newdata = X_new, pihat_new = pihat_new,
                            type = "tau")
  expect_length(ite_new, n_new)
  expect_true(all(is.finite(ite_new)))

  ate_new <- predict_cf_bcf(fit, newdata = X_new, pihat_new = pihat_new,
                            type = "ate")
  expect_length(ate_new, 1L)
  expect_true(is.finite(ate_new))

  y0_new <- predict_cf_bcf(fit, newdata = X_new, pihat_new = pihat_new,
                           type = "mu0")
  y1_new <- predict_cf_bcf(fit, newdata = X_new, pihat_new = pihat_new,
                           type = "mu1")
  expect_length(y0_new, n_new)
  expect_length(y1_new, n_new)
  # Plug-in identity by construction: y1 - y0 == ite on new data.
  expect_equal(y1_new - y0_new, ite_new, tolerance = 1e-10)
})

test_that("cf_fit_bcf newdata auto-computes pihat_new", {
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()

  set.seed(11)
  n <- 80
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- T + X[, 1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 40, n_iter = 80, n_chains = 1L,
                save_tree_directory = tempfile("cfomics_bcf_err_"))

  # BCF newdata predict auto-computes pihat_new from the persisted PS model.
  expect_true(!is.null(fit$model$ps_model))

  X_new <- matrix(rnorm(10 * p), 10, p)
  colnames(X_new) <- cov_names

  out <- predict_cf_bcf(fit, newdata = X_new, pihat_new = NULL, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat", "mu1_hat"))
  expect_length(out$tau_hat, 10L)
  expect_length(out$mu0_hat, 10L)
})
