## Unified newdata predict API smoke tests.
##
## Covers:
##   * type = "all" flat list contract (CATE-capable methods)
##   * type = "tau_ci" / "ate_ci" CI surface for grf
##   * CI on unsupported methods raises cfomics_ite_ci_unsupported

make_demo <- function(n = 80, p = 4, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.3 * X[, 1] + 0.4 * T + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(
    df  = data.frame(Y = Y, T = T, X),
    X_new = matrix(rnorm(20 * p), 20, p,
                   dimnames = list(NULL, paste0("X", seq_len(p))))
  )
}

test_that("predict.cf_model(type = 'all') returns flat list for grf", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = d$df, method = "grf",
                num.trees = 200L, tune.parameters = "none")
  out <- predict(fit, newdata = d$X_new, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat", "mu1_hat"))
  expect_length(out$tau_hat, 20L)
  expect_length(out$mu0_hat, 20L)
  expect_length(out$mu1_hat, 20L)
  expect_true(is.numeric(out$ate))
  # GRF reconstructs y0/y1 from forest.Y/W; plug-in identity holds.
  expect_equal(out$mu1_hat - out$mu0_hat, out$tau_hat, tolerance = 1e-8)
  # Per-point CI on newdata via IJ variance (Wald) is valid...
  expect_true(is.data.frame(out$tau_ci))
  expect_named(out$tau_ci, c("lower", "upper"))
  # ...but there is no honest ATE CI on newdata (grf has no outcomes there),
  # so ate_ci is NA rather than a sqrt(mean(var_i)) approximation.
  expect_true(all(is.na(out$ate_ci)))
})

test_that("predict.cf_model raises cfomics_ite_ci_unsupported for non-linear DR final", {
  skip_if_no_cf_python()
  d <- make_demo()
  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = d$df, method = "drlearner",
                model_final = "rf")
  expect_error(
    predict(fit, newdata = d$X_new, type = "tau_ci"),
    class = "cfomics_drlearner_ci_unsupported_nonlinear_final"
  )
})

test_that("predict.cf_model returns a data.frame ite_ci for grf newdata", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = d$df, method = "grf",
                num.trees = 200L, tune.parameters = "none")
  ite_ci <- predict(fit, newdata = d$X_new, type = "tau_ci")
  expect_s3_class(ite_ci, "data.frame")
  expect_named(ite_ci, c("lower", "upper"))
  expect_true(all(ite_ci$lower <= ite_ci$upper))
})

test_that("predict.cf_model dispatches stored fields for newdata = NULL", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4, data = d$df, method = "grf",
                num.trees = 200L, tune.parameters = "none")
  out <- predict(fit, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat", "mu1_hat"))
  expect_true(all(!is.na(out$tau_hat)))
  expect_true(all(!is.na(out$mu0_hat)))
})
