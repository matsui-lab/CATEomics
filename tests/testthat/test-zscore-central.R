# Central covariate z-score: every column standardised (binary indicators
# included), leakage-free stats, and the grf/bcf scale-invariance the
# centralisation relies on.

test_that(".cf_compute_x_zscore standardises every column, binary included", {
  set.seed(1)
  X <- cbind(
    g1  = rnorm(50, mean = 3, sd = 2),
    g2  = rnorm(50, mean = -1, sd = 0.5),
    bin = rbinom(50, 1, 0.5),          # binary 0/1
    onehot = rep(c(0, 1), each = 25),  # one-hot style
    const = rep(2, 50)                 # constant column
  )
  st <- CATEomics:::.cf_compute_x_zscore(X)

  # continuous columns get their true mean / sd
  expect_equal(unname(st$X_mean["g1"]), mean(X[, "g1"]))
  expect_true(st$X_sd["g1"] > 1)
  # Binary indicators are standardised too: on their native scale a penalised
  # learner shrinks them far harder than the z-scored continuous columns.
  expect_equal(unname(st$X_mean["bin"]), mean(X[, "bin"]))
  expect_equal(unname(st$X_sd["bin"]), stats::sd(X[, "bin"]))
  expect_equal(unname(st$X_mean["onehot"]), mean(X[, "onehot"]))
  expect_equal(unname(st$X_sd["onehot"]), stats::sd(X[, "onehot"]))
  # constant column: sd falls back to 1 (no divide-by-zero)
  expect_equal(unname(st$X_sd["const"]), 1)
})

test_that(".cf_apply_x_zscore yields unit-scale columns, binary included", {
  set.seed(2)
  X <- cbind(
    g1  = rnorm(80, 5, 3),
    bin = rbinom(80, 1, 0.4)
  )
  st <- CATEomics:::.cf_compute_x_zscore(X)
  Z  <- CATEomics:::.cf_apply_x_zscore(X, st)
  expect_equal(mean(Z[, "g1"]), 0, tolerance = 1e-8)
  expect_equal(stats::sd(Z[, "g1"]), 1, tolerance = 1e-8)
  expect_equal(mean(Z[, "bin"]), 0, tolerance = 1e-8)
  expect_equal(stats::sd(Z[, "bin"]), 1, tolerance = 1e-8)
})

test_that("central z-score is applied once in cf_fit and stored in meta", {
  set.seed(3)
  n <- 60; p <- 6
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("g", seq_len(p))
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X)
  fit <- cf_fit(data = df, Y = "Y", T = "T", X_baseline = paste0("g", seq_len(p)),
                method = "grf", partition = "none")
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
  # bcf's own defaults write the MCMC tree draws and its progress log to the
  # working directory, which during a test run is tests/testthat (one run left
  # 28 MB there). CATEomics redirects both into a tempdir on every call, but this
  # reference calls bcf directly, so it has to do the same.
  bdir <- file.path(tempdir(), "zscore_bcf_trees")
  dir.create(bdir, showWarnings = FALSE, recursive = TRUE)
  blog <- tempfile("bcf_log_", tmpdir = bdir, fileext = ".txt")
  b1 <- suppressWarnings(bcf::bcf(y = Y, z = W, x_control = X, x_moderate = X,
                                  pihat = ps, nburn = 100, nsim = 100, verbose = FALSE,
                                  save_tree_directory = bdir, log_file = blog))
  set.seed(4)
  b2 <- suppressWarnings(bcf::bcf(y = Y, z = W, x_control = Xz, x_moderate = Xz,
                                  pihat = ps, nburn = 100, nsim = 100, verbose = FALSE,
                                  save_tree_directory = bdir, log_file = blog))
  expect_equal(colMeans(b1$tau), colMeans(b2$tau), tolerance = 1e-6)
})
