# Central covariate z-score: binary/one-hot skip, continuous standardisation,
# leakage-free stats, and the grf/bcf scale-invariance the centralisation relies
# on.

test_that(".cf_compute_x_zscore standardises continuous cols and skips binary", {
  set.seed(1)
  X <- cbind(
    g1  = rnorm(50, mean = 3, sd = 2),
    g2  = rnorm(50, mean = -1, sd = 0.5),
    bin = rbinom(50, 1, 0.5),          # binary 0/1
    onehot = rep(c(0, 1), each = 25),  # one-hot style
    const = rep(2, 50)                 # constant column
  )
  st <- cfomics:::.cf_compute_x_zscore(X)

  # continuous columns get their true mean / sd
  expect_equal(unname(st$X_mean["g1"]), mean(X[, "g1"]))
  expect_true(st$X_sd["g1"] > 1)
  # binary / one-hot columns kept native (mean 0, sd 1 => identity)
  expect_equal(unname(st$X_mean["bin"]), 0)
  expect_equal(unname(st$X_sd["bin"]), 1)
  expect_equal(unname(st$X_mean["onehot"]), 0)
  expect_equal(unname(st$X_sd["onehot"]), 1)
  # constant column: sd falls back to 1 (no divide-by-zero)
  expect_equal(unname(st$X_sd["const"]), 1)
})

test_that(".cf_apply_x_zscore yields unit-scale continuous cols, native binary", {
  set.seed(2)
  X <- cbind(
    g1  = rnorm(80, 5, 3),
    bin = rbinom(80, 1, 0.4)
  )
  st <- cfomics:::.cf_compute_x_zscore(X)
  Z  <- cfomics:::.cf_apply_x_zscore(X, st)
  expect_equal(mean(Z[, "g1"]), 0, tolerance = 1e-8)
  expect_equal(stats::sd(Z[, "g1"]), 1, tolerance = 1e-8)
  # binary column unchanged (skip)
  expect_equal(Z[, "bin"], X[, "bin"])
})

test_that("central z-score is applied once in cf_fit and stored in meta", {
  set.seed(3)
  n <- 60; p <- 6
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("g", seq_len(p))
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X)
  f  <- stats::as.formula(paste("Y ~ T |", paste(paste0("g", seq_len(p)), collapse = " + ")))
  fit <- cf_fit(f, data = df, method = "grf", partition = "none")
  expect_false(is.null(fit$meta$x_zscore))
  expect_equal(length(fit$meta$x_zscore$X_mean), p)
  # X_preprocessed columns should already be standardised (unit sd)
  Xp <- fit$data$X_preprocessed
  if (!is.null(Xp)) {
    expect_equal(mean(apply(Xp, 2, stats::sd)), 1, tolerance = 0.1)
  }
})

test_that("grf and bcf are invariant to the central z-score", {
  skip_if_not_installed("grf")
  set.seed(1); n <- 120; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- rnorm(n) + X[, 1] * W
  Xz <- scale(X)

  set.seed(2); f1 <- grf::causal_forest(X,  Y, W, num.trees = 500)
  set.seed(2); f2 <- grf::causal_forest(Xz, Y, W, num.trees = 500)
  expect_equal(predict(f1)$predictions, predict(f2)$predictions, tolerance = 1e-8)

  skip_if_not_installed("bcf")
  ps <- rep(0.5, n)
  set.seed(4)
  b1 <- suppressWarnings(bcf::bcf(y = Y, z = W, x_control = X, x_moderate = X,
                                  pihat = ps, nburn = 100, nsim = 100, verbose = FALSE))
  set.seed(4)
  b2 <- suppressWarnings(bcf::bcf(y = Y, z = W, x_control = Xz, x_moderate = Xz,
                                  pihat = ps, nburn = 100, nsim = 100, verbose = FALSE))
  expect_equal(colMeans(b1$tau), colMeans(b2$tau), tolerance = 1e-6)
})
