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
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "ganite",
          params = list(iterations = 20L, batch_size = 32L))

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

  # `covariate_names` is not one of this fitter's arguments: nothing in it
  # reads the names.
  fit <- cf_fit_ganite(
    X = X, T = T, Y = Y,
    random_state = 42L,
    iterations = 20L,
    batch_size = 16L
  )

  expect_s3_class(fit, "cf_model_ganite")
  # The raw method returns a two-pillar (model + res) list; the
  # model pillar holds the serialised path + backend, and cf_fit() adds the
  # central summary (summary/samples are not on the raw method output).
  expect_equal(fit$model$backend, "CATEomics.ganite_impl.GaniteModel")
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


  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "ganite",
          params = list(iterations = 20L, batch_size = 16L))

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

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "ganite",
          random_state = 7L, params = list(iterations = 20L,
              batch_size = 16L))

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


# Structure of the implementation ------------------------------------------
#
# These do not check accuracy. They check that the four networks are wired the
# way the paper states: the effect is the difference of the two inference
# heads, both outcome heads move away from their initialisation, every
# hyperparameter reaches the backend, and the trace covers both blocks.

.gn_data <- function(n = 120L, p = 8L, seed = 4L, binary = FALSE) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  T <- stats::rbinom(n, 1L, 0.5)
  tau <- 1 + 2 * X[, 1]
  Y <- if (binary) stats::rbinom(n, 1L, stats::plogis(X[, 2] + T * tau))
       else X[, 2] + T * tau + stats::rnorm(n)
  # `prog` is the prognostic part of the outcome, which both heads carry and
  # the effect does not, so it separates "the head learned something" from
  # "the effect was recovered".
  list(df = data.frame(Y = Y, T = T, X), cv = colnames(X), tau = tau,
       prog = X[, 2])
}

.gn_fit <- function(d, ..., iterations = 150L) {
  cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv, method = "ganite",
         preprocess = FALSE, random_state = 0L, partition = "none",
         params = c(list(iterations = iterations), list(...)))
}

test_that("the effect stays the difference of the two ITE heads", {
  skip_if_no_cf_python()
  d <- .gn_data()
  f <- .gn_fit(d, iterations = 40L)
  # The R side subtracts the two heads, so this identity is arithmetic rather
  # than evidence about the network; it pins that the assembly keeps working
  # that way instead of starting to carry a separately estimated effect, as
  # the residualising methods do.
  expect_equal(f$res$mu1_hat - f$res$mu0_hat, f$res$tau_hat, tolerance = 0)
  # What the backend itself reports has to agree with the assembled result:
  # `predict(newdata = )` reloads the saved ITE block and recomputes both
  # potential outcomes in Python, so this crosses the serialisation boundary.
  X_tr <- as.matrix(d$df[, d$cv, drop = FALSE])
  expect_equal(predict(f, newdata = X_tr, type = "tau"),
               predict(f, newdata = X_tr, type = "mu1") -
                 predict(f, newdata = X_tr, type = "mu0"), tolerance = 0)
})

test_that("the adversarial term responds to its weight", {
  skip_if_no_cf_python()
  d <- .gn_data()
  # Pinned: the trace keeps one row per 100 steps plus the step the block
  # ended on, so the last row is the end of a 150-step block. Shorter runs
  # would make the last row the first step, measured before the generator has
  # been updated at all, which no choice of alpha can move.
  loss <- function(a) {
    f <- .gn_fit(d, alpha = a, iterations = 150L)
    tr <- f$meta$estimator$training_trace
    cf <- tr[tr$block == "counterfactual", ]
    expect_identical(cf$iteration[nrow(cf)], 149L)
    cf$d_g[nrow(cf)]
  }
  # What is pinned: the recorded loss responds to alpha at this setting
  # (0.177 at alpha = 0 against 0.005 at alpha = 10 on this draw). The trace
  # records a training loss, which is not the discriminator's quality -- it
  # memorises the training rows -- so this is a regression pin on one draw,
  # not a property of the method.
  expect_gt(loss(0), loss(10))
})

test_that("a zero supervised weight loses the effect, a positive one keeps it", {
  skip_if_no_cf_python()
  d <- .gn_data()
  cor_at <- function(a)
    stats::cor(.gn_fit(d, alpha = a, iterations = 150L)$res$tau_hat, d$tau)
  # With no supervised term the generator only has to fool the discriminator,
  # which says nothing about the outcome scale. A regression pin on this draw
  # at this length (0.233 against 0.768), not a property of the method.
  expect_lt(cor_at(0), 0.5)
  expect_gt(cor_at(2), 0.5)
})

test_that("both outcome heads move away from their initialisation", {
  skip_if_no_cf_python()
  d <- .gn_data()
  # A head that receives no gradient is what this catches. Asking whether the
  # head varies across samples does not: an untrained head is a
  # Xavier-initialised network of X, so it varies anyway. What a frozen head
  # cannot do is move, and the seed makes the initialisation identical
  # between these two fits.
  untrained <- .gn_fit(d, iterations = 1L)$res
  trained   <- .gn_fit(d, iterations = 150L)$res
  # Measured: 2.07 for mu0 and 5.43 for mu1. A frozen head gives exactly 0.
  expect_gt(max(abs(trained$mu0_hat - untrained$mu0_hat)), 0.1)
  expect_gt(max(abs(trained$mu1_hat - untrained$mu1_hat)), 0.1)
  # And the control head has to learn the prognostic part, not just move:
  # measured cor 0.787 after 150 steps against -0.167 untrained.
  expect_gt(stats::cor(trained$mu0_hat, d$prog), 0.5)
})

test_that("the trace covers both blocks and ends where training ended", {
  skip_if_no_cf_python()
  d <- .gn_data()
  tr <- .gn_fit(d, iterations = 120L)$meta$estimator$training_trace
  expect_true(is.data.frame(tr))
  expect_true(all(c("block", "iteration", "d_g", "g", "g_factual",
                    "d_i", "i", "l_s_i") %in% names(tr)))
  expect_setequal(unique(tr$block), c("counterfactual", "ite"))
  # One row per 100 steps leaves the end of a 120-step block unrecorded unless
  # the last step is always kept, which is what a reader of a loss curve wants
  # and what the sibling helper for the variational fit already did.
  expect_setequal(unique(tr$iteration), c(0L, 100L, 119L))
  # The effect discriminator exists, so its loss is a number rather than the
  # absent column it would be without that network.
  ite <- tr[tr$block == "ite", ]
  expect_true(all(is.finite(ite$d_i)))
  expect_true(all(is.finite(ite$i)))
  expect_true(all(is.na(ite$d_g)))
})

test_that("a binary outcome keeps the estimates on the probability scale", {
  skip_if_no_cf_python()
  d <- .gn_data(binary = TRUE)
  f <- .gn_fit(d)
  # The heads are read through a logistic function, so the potential outcomes
  # are in [0, 1] and the effect is a risk difference in [-1, 1]. The bounds
  # are inclusive: a saturated head is exactly 1 in single precision, which is
  # the end of the range rather than outside it.
  expect_true(all(f$res$mu0_hat >= 0 & f$res$mu0_hat <= 1))
  expect_true(all(f$res$mu1_hat >= 0 & f$res$mu1_hat <= 1))
  # A regression on 0/1 would leave the range, so the effect would too.
  expect_true(all(abs(f$res$tau_hat) <= 1))
})

test_that("the same seed reproduces the fit and a different one moves it", {
  skip_if_no_cf_python()
  d <- .gn_data()
  # The ambient RNG is moved between the calls: CATEomics restores
  # `.Random.seed`, so two back-to-back calls would agree even if the seed were
  # not honoured.
  set.seed(111); a <- .gn_fit(d, iterations = 40L)$res$tau_hat
  set.seed(222); b <- .gn_fit(d, iterations = 40L)$res$tau_hat
  expect_equal(a, b, tolerance = 0)
  c2 <- cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
               method = "ganite", preprocess = FALSE, random_state = 1L,
               partition = "none", params = list(iterations = 40L))$res$tau_hat
  expect_gt(max(abs(a - c2)), 1e-6)
})

test_that("every hyperparameter reaches the backend", {
  skip_if_no_cf_python()
  d <- .gn_data()
  base <- .gn_fit(d, iterations = 60L)$res$tau_hat
  # Being on the allow-list is not the same as arriving, and `d_steps` -- the
  # ratio of discriminator steps to its counterpart -- is the one nothing else
  # here would notice the loss of.
  alt <- list(alpha = 10, d_steps = 5L, beta = 0, h_dim = 32L,
              num_layers = 3L, batch_size = 32L, learning_rate = 0.01,
              weight_decay = 0.1)
  for (nm in names(alt)) {
    args <- c(list(d = d, iterations = 60L), stats::setNames(alt[nm], nm))
    moved <- max(abs(do.call(.gn_fit, args)$res$tau_hat - base))
    expect_gt(moved, 1e-6)
  }
})

test_that("a hyperparameter below one is refused rather than floored", {
  skip_if_no_cf_python()
  d <- .gn_data()
  # The backend floors these to keep its loops well formed, so a zero would
  # run: an empty mini-batch leaves every network at its initialisation and
  # one step per block is not a fit. Neither says so, hence the refusal here.
  for (nm in c("h_dim", "num_layers", "iterations", "batch_size", "d_steps")) {
    pl <- stats::setNames(list(0L), nm)
    if (nm != "iterations") pl$iterations <- 20L
    expect_error(
      cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
             method = "ganite", preprocess = FALSE, random_state = 0L,
             partition = "none", params = pl),
      class = "cateomics_invalid_ganite_param")
  }
})

test_that("the hyperparameters are the implementation's own", {
  skip_if_no_cf_python()
  d <- .gn_data()
  # ganite is implemented in this package, so its names are checked here
  # rather than left to an upstream constructor.
  expect_error(.gn_fit(d, depth = 3L), class = "cateomics_unknown_param")
  expect_error(.gn_fit(d, verbose = TRUE), class = "cateomics_unknown_param")
  expect_error(.gn_fit(d, seed = 1L), class = "cateomics_rng_param")
  expect_setequal(
    .cf_params_allowed("ganite"),
    c("h_dim", "num_layers", "alpha", "beta", "iterations", "batch_size",
      "learning_rate", "weight_decay", "d_steps"))
})
