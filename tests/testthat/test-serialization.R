# Test serialization (save/load) of CATEomics results.
# R-native methods store the backend object directly. Python methods
# store a serialised-model path in the `model` pillar (no live Python object),
# so the whole cf_model is saveRDS-able and predict() reloads from that path.

test_that("R-native grf results can be saved and restored", {
  skip_if_not_installed("grf")

  set.seed(789)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T + 0.3 * data$X1

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "grf")

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  expect_s3_class(fit2, "cf_model")
  expect_equal(fit2$method, "grf")

  # ATE should be preserved
  ate1 <- predict(fit, type = "ate")
  ate2 <- predict(fit2, type = "ate")
  expect_equal(ate1, ate2, tolerance = 1e-10)

  # ITE should be preserved
  ite1 <- predict(fit, type = "tau")
  ite2 <- predict(fit2, type = "tau")
  expect_equal(ite1, ite2, tolerance = 1e-10)

  # Verify summary works
  expect_no_error(summary(fit2))
})

# Portable output is plain saveRDS() / write.csv() / sessionInfo() (see the
# R-native saveRDS round-trip tests above).

test_that("meta information is preserved after serialization", {
  skip_if_not_installed("grf")
  set.seed(222)
  n <- 40
  data <- data.frame(
    outcome = rnorm(n),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10)
  )
  data$outcome <- data$outcome + 1.5 * data$treatment

  fit <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = data, method = "grf",
                partition = "none", random_state = 42)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  # Meta should be fully preserved. The variables the caller named are in
  # `named_args`, the seed in `provenance$rng_seed` and the split in
  # `estimator`: read them where they live, so that a field moving again makes
  # this fail rather than pass on two NULLs.
  expect_equal(fit2$meta$n, fit$meta$n)
  expect_equal(fit2$meta$p, fit$meta$p)
  expect_identical(fit2$meta$named_args, fit$meta$named_args)
  expect_identical(fit2$meta$named_args$Y, "outcome")
  expect_identical(fit2$meta$named_args$T, "treatment")
  expect_equal(fit2$meta$covariate_names, fit$meta$covariate_names)
  expect_equal(fit2$meta$provenance$rng_seed, fit$meta$provenance$rng_seed)
  expect_identical(fit2$meta$estimator$partition, fit$meta$estimator$partition)
  # The whole pillar, so a field that is added or dropped is caught too.
  expect_identical(names(fit2$meta), names(fit$meta))
})

test_that("Python-backed models survive saveRDS and predict from the restored path", {
  skip_if_no_cf_python()
  set.seed(202)
  n <- 80L
  df <- data.frame(Y = rnorm(n) + 0.5 * rbinom(n, 1, 0.5),
                   T = rbinom(n, 1, 0.5),
                   X1 = rnorm(n), X2 = rnorm(n))
  nd <- data.frame(X1 = rnorm(10), X2 = rnorm(10))

  # one joblib (econml) and one torch.save (catenets) backend.
  for (m in c("drlearner", "tarnet")) {
    extra <- if (m == "tarnet") list(params = list(n_iter = 50L)) else list()
    fit <- suppressWarnings(suppressMessages(do.call(cf_fit, c(
      list(data = df, Y = "Y", T = "T", X_baseline = c("X1", "X2"),
           method = m, partition = "none",
           random_state = 1L), extra))))

    # The model pillar is a path (no live Python object), so saveRDS works.
    expect_type(fit$model$model$model_path, "character")
    tau_before <- predict(fit, newdata = nd, type = "tau")

    tmp <- tempfile(fileext = ".rds"); on.exit(unlink(tmp), add = TRUE)
    saveRDS(fit, tmp)
    rm(fit)
    restored <- readRDS(tmp)

    expect_s3_class(restored, "cf_model")
    expect_true(file.exists(restored$model$model$model_path))  # cache persists
    # predict on the restored object reloads the serialised model from its path.
    tau_after <- predict(restored, newdata = nd, type = "tau")
    expect_equal(length(tau_after), nrow(nd))
    expect_equal(tau_after, tau_before, tolerance = 1e-8)
  }
})

# ---------------------------------------------------------------------------
# What a returned object carries beyond its own content.
#
# In R an environment travels by reference, and a formula and a closure each
# hold one. A synthesised formula built with the calling frame, or a closure
# written inline, therefore carries everything that frame can reach -- which
# was the input container itself. These are the assertions that fail if that
# comes back: the point of the change is structural, so nothing about the
# numbers would notice.
# ---------------------------------------------------------------------------

test_that("a fit does not reach the objects of the frame it was called from", {
  skip_if_not_installed("grf")
  fit_in_frame <- function(df, with_big) {
    big <- if (with_big) matrix(0, 400L, 400L) else NULL   # 1.2 MB when present
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df,
           method = "grf", partition = "none", random_state = 1L,
           params = list(num.trees = 100L, tune.parameters = "none"))
  }
  set.seed(4)
  n  <- 60L
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                   X1 = rnorm(n), X2 = rnorm(n))
  lean <- fit_in_frame(df, with_big = FALSE)
  fat  <- fit_in_frame(df, with_big = TRUE)
  # The returned object is the same size whatever else the calling frame held,
  # which is what fails if anything in it captures that frame. Measured with
  # serialize, not object.size: the latter does not follow environments.
  expect_equal(length(serialize(fat, NULL)), length(serialize(lean, NULL)))
  # No environment anywhere in the fit can see the 1.2 MB matrix.
  reaches_big <- function(x, depth = 0L) {
    if (depth > 12L) return(FALSE)
    if (is.environment(x)) return("big" %in% ls(x))
    if (is.list(x)) return(any(vapply(x, reaches_big, logical(1), depth + 1L)))
    FALSE
  }
  expect_false(reaches_big(fat))
  # The canonical formula the named arguments are turned into is built in
  # `baseenv()`, which is what keeps the frame out. Assert it directly: the
  # field that used to carry it is gone, so nothing else would notice if the
  # environment came back.
  built <- CATEomics:::.build_formula_from_named(
    Y = "Y", T = "T", X_baseline = c("X1", "X2"))
  expect_identical(environment(built$formula), baseenv())
})

test_that("a treatment rule carries its cut-point, not the model", {
  skip_if_not_installed("grf")
  set.seed(5)
  n  <- 80L
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                   X1 = rnorm(n), X2 = rnorm(n))
  df$Y <- df$Y + df$T * df$X1
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df,
                method = "grf", partition = "oof", n_folds = 3L,
                random_state = 1L, eval_nuisance = "ridge",
                params = list(num.trees = 100L, tune.parameters = "none"))

  # What identifies each rule is a plain value, not a closure, so nothing in
  # the returned object can reach the fit it came from.
  rl <- cf_rule(fit, type = "threshold", threshold = 0, 
                compute_policy_value = FALSE, random_state = 1L)
  expect_identical(rl$threshold, 0)
  expect_identical(rl$cutoff, NA_real_)
  expect_identical(rl$decisions,
                   as.integer(fit$res$tau_hat_oof > rl$threshold))

  rl <- cf_rule(fit, type = "topk", topk_frac = 1 / 3,
                compute_policy_value = FALSE, random_state = 1L)
  expect_true(is.numeric(rl$cutoff) && length(rl$cutoff) == 1L)
  expect_true(is.finite(rl$cutoff))
  expect_identical(rl$threshold, NA_real_)
  # the cut-point reproduces the decisions it was read off
  expect_identical(rl$decisions,
                   as.integer(is.finite(fit$res$tau_hat_oof) &
                                fit$res$tau_hat_oof >= rl$cutoff))

  rl <- cf_rule(fit, type = "tree", compute_policy_value = FALSE,
                random_state = 1L)
  expect_identical(rl$cutoff, NA_real_)
  expect_s3_class(rl$tree, "policy_tree")

  # No component of a rule carries an environment any more, so the size of the
  # rule cannot follow the size of the fit it was built from.
  expect_false(any(vapply(rl, is.function, logical(1))))
  big <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df,
                method = "grf", partition = "oof", n_folds = 3L,
                random_state = 1L, eval_nuisance = "ridge",
                params = list(num.trees = 2000L, tune.parameters = "none"))
  expect_gt(length(serialize(big, NULL)), 5 * length(serialize(fit, NULL)))
  r_small <- cf_rule(fit, type = "topk", compute_policy_value = FALSE,
                     random_state = 1L)
  r_big   <- cf_rule(big, type = "topk", compute_policy_value = FALSE,
                     random_state = 1L)
  expect_lt(abs(length(serialize(r_big, NULL)) -
                  length(serialize(r_small, NULL))), 2000L)
})
