# Confidence-interval smoke tests:
# - DR-Learner linear CI via StatsModelsLinearRegression
# - bootstrap_ci opt-in shared machinery + hybrid policy
# - CEVAE Pyro Predictive (ATE + ITE CI)
# - BCF MCMC diagnostics (n_chains > 1 + Rhat)
# - ci_method label propagation (summary.cf_model + metadata)
# - res$ate_ci / res$tau_ci unified location + validator.

.simulate_data <- function(n = 150, seed = 0L) {
  set.seed(seed)
  X1 <- stats::rnorm(n)
  X2 <- stats::rnorm(n)
  T  <- stats::rbinom(n, 1, 0.5)
  Y  <- 1.5 * T + 0.5 * X1 + stats::rnorm(n)
  data.frame(X1 = X1, X2 = X2, T = T, Y = Y)
}

# ----- ate_ci / tau_ci / ci_method on R-native methods -----

test_that("grf populates res$ate_ci, res$tau_ci, and ci_method = 'asymptotic Wald'", {
  skip_if_not_installed("grf")
  df <- .simulate_data(n = 200, seed = 1L)
  m <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf", random_state = 0L)
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 200L)
  expect_true(all(c("lower", "upper") %in% colnames(res$tau_ci)))
  expect_identical(res$ci_method, "asymptotic Wald")
  expect_match(format(summary(m)), "asymptotic Wald", all = FALSE)
})

# ----- hybrid policy: error paths before any fit work -----

test_that("bootstrap_ci = TRUE errors for NN methods (tarnet/dragonnet/ganite)", {
  df <- .simulate_data(n = 80, seed = 5L)
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "tarnet", bootstrap_ci = TRUE),
    class = "cateomics_bootstrap_unsupported_nn"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "dragonnet", bootstrap_ci = TRUE),
    class = "cateomics_bootstrap_unsupported_nn"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "ganite", bootstrap_ci = TRUE),
    class = "cateomics_bootstrap_unsupported_nn"
  )
})

test_that("bootstrap_ci = TRUE errors for the backends that cannot resample", {
  # grf / bcf / cevae have no `bootstrap_ci` formal, so the request never
  # reached them: it is refused at the entry instead of being dropped.
  # The check reads the fitter's formals, so these three are not a hard-coded
  # list and the test also guards that resolution.
  df <- .simulate_data(n = 80, seed = 5L)
  for (meth in c("grf", "bcf", "cevae")) {
    expect_error(
      cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = meth, bootstrap_ci = TRUE),
      class = "cateomics_bootstrap_unsupported",
      info = meth
    )
  }
  # The metalearner path does accept it, so it must not be caught here.
  expect_false("bootstrap_ci" %in% names(formals(CATEomics:::cf_fit_grf)))
  expect_true("bootstrap_ci" %in% names(formals(CATEomics:::cf_fit_tlearner)))
})

# ----- bootstrap_ci accept paths (metalearners) -----

test_that("slearner with bootstrap_ci = TRUE produces ate_ci + tau_ci", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 100, seed = 6L)
  m <- cf_fit(
    Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "slearner",
    bootstrap_ci = TRUE, n_bootstrap = 10L, random_state = 1L
  )
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 100L)
  expect_identical(res$ci_method, "nonparametric bootstrap")
})

# ----- DR-Learner linear via StatsModelsLinearRegression -----

test_that("drlearner (model_final = 'linear') returns asymptotic Wald ATE + ITE CIs", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 200, seed = 7L)
  m <- cf_fit(
    Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "drlearner",
    model_final = "linear", random_state = 1L
  )
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 200L)
  expect_identical(res$ci_method, "asymptotic Wald")
})

# ----- CEVAE Pyro Predictive -----

test_that("cevae populates res$ate_ci / res$tau_ci via Pyro Predictive", {
  skip_if_no_cf_python()
  df <- .simulate_data(n = 80, seed = 8L)
  m <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "cevae",
          random_state = 1L, params = list(num_epochs = 3L,
              num_samples = 5L))
  res <- m$res
  expect_true(is.numeric(res$ate_ci) && length(res$ate_ci) == 2L)
  expect_lte(res$ate_ci[1], res$ate_ci[2])
  expect_s3_class(res$tau_ci, "data.frame")
  expect_equal(nrow(res$tau_ci), 80L)
  expect_identical(res$ci_method, "Bayesian posterior")
})

# ----- BCF MCMC diagnostic -----

test_that("bcf with n_chains > 1 produces a Gelman-Rubin Rhat in metadata", {
  skip_if_not_installed("bcf")
  df <- .simulate_data(n = 100, seed = 9L)
  m <- suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", random_state = 1L,
          params = list(nburn = 50L, nsim = 100L, n_chains = 2L)))
  expect_equal(m$res$metadata$n_chains, 2L)
  expect_true(is.numeric(m$res$metadata$rhat))
  expect_false(is.na(m$res$metadata$rhat))
  expect_identical(m$res$ci_method, "Bayesian posterior")
})

test_that("bcf does not carry res$credible_intervals", {
  skip_if_not_installed("bcf")
  df <- .simulate_data(n = 80, seed = 10L)
  m <- suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "bcf", random_state = 1L,
          params = list(nburn = 50L, nsim = 100L)))
  expect_null(m$res$credible_intervals)
  expect_s3_class(m$res$tau_ci, "data.frame")
})

# ----- validator: bad ate_ci should abort -----

test_that("validate_cateomics_result rejects malformed ate_ci / tau_ci", {
  skip_if_not_installed("grf")
  df <- .simulate_data(n = 100, seed = 11L)
  m <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf", random_state = 0L)

  m_bad <- m
  m_bad$res$ate_ci <- c(2, 1)  # reversed bounds
  expect_error(validate_cateomics_result(m_bad), class = "cateomics_invalid_result")

  m_bad2 <- m
  m_bad2$res$tau_ci <- data.frame(lo = 1:10)  # wrong columns
  expect_error(validate_cateomics_result(m_bad2), class = "cateomics_invalid_result")
})


# The half-sample bootstrap is calibrated -------------------------------

test_that("the half-sample SE estimates the full-sample SE, not sqrt(2) times it", {
  # Half-sample resampling draws floor(n/2) units without replacement. The
  # smaller sample inflates the spread by about sqrt(2), and the
  # finite-population correction cancels it; that cancellation is why the
  # plain sd of the replicates is used as the SE. Checked on a statistic
  # whose SE is known in closed form.
  set.seed(101)
  n <- 400L
  x <- rnorm(n, sd = 2)
  se_hs <- CATEomics:::.cf_halfsample_se(n, function(idx) mean(x[idx]),
                                       n_bootstrap = 4000L, random_state = 1L)
  se_analytic <- stats::sd(x) / sqrt(n)
  expect_equal(se_hs / se_analytic, 1, tolerance = 0.05)
  # Without the correction this ratio would sit near sqrt(2).
  expect_lt(se_hs / se_analytic, 1.2)
})

test_that("the half-sample SE is reproducible and needs a random_state", {
  set.seed(7)
  x <- rnorm(200L)
  f <- function(idx) mean(x[idx])
  a <- CATEomics:::.cf_halfsample_se(200L, f, n_bootstrap = 200L, random_state = 3L)
  b <- CATEomics:::.cf_halfsample_se(200L, f, n_bootstrap = 200L, random_state = 3L)
  cc <- CATEomics:::.cf_halfsample_se(200L, f, n_bootstrap = 200L, random_state = 4L)
  expect_identical(a, b)
  expect_false(identical(a, cc))
  expect_gt(a, 0)
})

test_that("a method-side bootstrap interval is reproducible and brackets the point", {
  skip_if_no_cf_python()
  set.seed(4)
  n <- 120L
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n), T = rbinom(n, 1L, 0.5))
  df$Y <- df$X1 + df$T * (1 + df$X2) + rnorm(n, sd = 0.5)

  fit <- function(seed) suppressMessages(cf_fit(
    Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "tlearner", bootstrap_ci = TRUE,
    n_bootstrap = 20L, random_state = seed))
  a <- fit(11L); b <- fit(11L); cc <- fit(12L)

  expect_equal(a$res$ate_ci, b$res$ate_ci)
  expect_equal(a$res$tau_ci, b$res$tau_ci)
  expect_equal(a$res$ate, b$res$ate)
  # A different seed has to give a different resampling. Seeding the stream
  # once per bootstrap rather than once per replicate is what makes this
  # hold: `seed + b` had neighbouring seeds sharing all but one resample.
  expect_false(isTRUE(all.equal(a$res$ate_ci, cc$res$ate_ci)))
  expect_lte(a$res$ate_ci[1], a$res$ate)
  expect_gte(a$res$ate_ci[2], a$res$ate)
})


test_that("every opt-in bootstrap path is reproducible and brackets its point", {
  skip_if_no_cf_python()
  # s / t / x-learner and a non-linear drlearner are the four paths whose
  # interval comes from resampling; grf, bcf, cevae and a linear drlearner
  # have native intervals. The older tests here checked the shape of the
  # returned object, which holds whatever numbers the resampling produced.
  set.seed(9)
  n <- 80L
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n), T = rbinom(n, 1L, 0.5))
  df$Y <- df$X1 + df$T * (1 + df$X2) + rnorm(n, sd = 0.5)

  for (meth in c("slearner", "xlearner", "drlearner")) {
    extra <- if (meth == "drlearner") list(model_final = "ridge") else list()
    run <- function(rs) suppressMessages(do.call(cf_fit, c(list(
      data = df, Y = "Y", T = "T", X_baseline = c("X1", "X2"),
      method = meth, bootstrap_ci = TRUE,
      n_bootstrap = 20L, random_state = rs), extra)))
    a <- run(21L); b <- run(21L)
    expect_equal(a$res$ate_ci, b$res$ate_ci,
                 info = paste(meth, "is not reproducible"))
    expect_equal(a$res$tau_ci, b$res$tau_ci, info = meth)
    expect_identical(a$res$ci_method, "nonparametric bootstrap",
                     info = meth)
    expect_lte(a$res$ate_ci[1], a$res$ate)
    expect_gte(a$res$ate_ci[2], a$res$ate)
    expect_equal(nrow(a$res$tau_ci), n, info = meth)
    # Reproducibility alone does not separate one seeded stream from a
    # fresh seed per replicate: both repeat exactly. The defect that was
    # here made neighbouring seeds share all but one resample, so require
    # an adjacent seed to give a different interval.
    cc <- run(22L)
    expect_false(isTRUE(all.equal(a$res$ate_ci, cc$res$ate_ci)),
                 info = paste(meth, "gives seed 21 and 22 the same interval"))
  }
})
