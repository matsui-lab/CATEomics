# Propensity-score exposure through predict() and the res pillar.
#
# The x-learner / dr-learner / dragonnet scores are read out of backend
# attributes (econml `propensity_models` / `models_propensity`, the catenets
# propensity head), so these tests pin the contract that extraction keeps:
# right length, strictly inside (0, 1), and an honest `ps_type`.

make_ps_data <- function(n = 80L, p = 4L, seed = 7L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))
  # Confounded assignment, so an estimated propensity is not constant.
  lin <- 0.9 * X[, 1] - 0.6 * X[, 2]
  T <- stats::rbinom(n, 1L, stats::plogis(lin))
  Y <- 1.2 * T + X[, 1] + stats::rnorm(n)
  list(df = data.frame(Y = Y, T = T, X), X = X,
       Xb = colnames(X),
       X_new = matrix(stats::rnorm(10L * p), 10L, p,
                      dimnames = list(NULL, colnames(X))))
}

expect_valid_ps <- function(ps, n) {
  expect_true(is.numeric(ps))
  expect_length(ps, n)
  expect_false(anyNA(ps))
  expect_true(all(ps > 0 & ps < 1))
}


test_that("predict() always carries ps and ps_type, NA for PS-free methods", {
  skip_if_not_installed("grf")
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "grf",
    partition = "none", params = list(num.trees = 100L,
                                      tune.parameters = "none")))
  out <- predict(fit, type = "all")
  expect_true(all(c("ps", "ps_type") %in% names(out)))
  expect_valid_ps(out$ps, nrow(d$df))
  expect_identical(out$ps_type, "oob")
})


test_that("partition = 'oof' assembles the propensity from the held-out folds", {
  skip_if_not_installed("grf")
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "grf",
    partition = "oof", n_folds = 3L, random_state = 1L,
    params = list(num.trees = 100L, tune.parameters = "none")))
  expect_identical(fit$res$ps_type, "oof")
  expect_valid_ps(fit$res$ps, nrow(d$df))
  # The held-out propensity is a different vector from the full-model one.
  full <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "grf",
    partition = "none", random_state = 1L,
    params = list(num.trees = 100L, tune.parameters = "none")))
  expect_false(isTRUE(all.equal(fit$res$ps, full$res$ps)))
})


test_that("a known propensity is reported as known, not as an estimate", {
  skip_if_not_installed("grf")
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "grf",
    partition = "none", propensity_score = 0.5,
    params = list(num.trees = 100L, tune.parameters = "none")))
  expect_identical(fit$res$ps_type, "known")
  expect_true(all(fit$res$ps == 0.5))
})


test_that("the clip bounds the propensity the model reports", {
  skip_if_not_installed("grf")
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "grf",
    partition = "none", clip = c(0.2, 0.8),
    params = list(num.trees = 100L, tune.parameters = "none")))
  expect_true(all(fit$res$ps >= 0.2 & fit$res$ps <= 0.8))
  ps_new <- predict(fit, newdata = d$X_new, type = "ps")
  expect_true(all(ps_new >= 0.2 & ps_new <= 0.8))
})


test_that("x-learner surfaces its econml propensity model", {
  skip_if_no_cf_python()
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "xlearner",
    partition = "none", random_state = 3L))
  expect_identical(fit$res$ps_type, "in-sample")
  expect_valid_ps(fit$res$ps, nrow(d$df))
  # Estimated, not the 0.5 placeholder the field used to hold.
  expect_true(stats::sd(fit$res$ps) > 0)

  out <- predict(fit, newdata = d$X_new, type = "all")
  expect_valid_ps(out$ps, 10L)
  expect_identical(out$ps_type, "out-of-sample")
})


test_that("dr-learner averages its cross-fit propensity models", {
  skip_if_no_cf_python()
  d <- make_ps_data()
  fit <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = "drlearner",
    partition = "none", random_state = 3L))
  expect_identical(fit$res$ps_type, "cross-fit (approx)")
  expect_valid_ps(fit$res$ps, nrow(d$df))
  expect_true(stats::sd(fit$res$ps) > 0)
})


test_that("the neural methods report no propensity", {
  skip_if_no_cf_python()
  d <- make_ps_data()
  # catenets gives DragonNet (and TARNet) a propensity head, but it is an
  # auxiliary training term: TARNet sets prop_loss_multiplier = 0 so the head
  # never learns, and DragonNet's head is neither regularised nor calibrated
  # (on a randomised design it puts a third of held-out rows outside
  # [0.025, 0.975]). Neither is published as a propensity estimate.
  for (m in c("dragonnet", "tarnet")) {
    fit <- suppressMessages(cf_fit(
      data = d$df, Y = "Y", T = "T", X_baseline = d$Xb, method = m,
      partition = "none", random_state = 3L, params = list(n_iter = 200L)))
    expect_true(all(is.na(predict(fit, type = "ps"))))
    expect_true(is.na(predict(fit, type = "all")$ps_type))
  }
})
