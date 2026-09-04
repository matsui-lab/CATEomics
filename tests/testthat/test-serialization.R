# Test serialization (save/load) of cfomics results.
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

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

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

  fit <- cf_fit(outcome ~ treatment | age, data = data, method = "grf",
                partition = "none", random_state = 42)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  # Meta should be fully preserved
  expect_equal(fit2$meta$n, fit$meta$n)
  expect_equal(fit2$meta$p, fit$meta$p)
  expect_equal(fit2$meta$outcome_name, fit$meta$outcome_name)
  expect_equal(fit2$meta$treatment_name, fit$meta$treatment_name)
  expect_equal(fit2$meta$covariate_names, fit$meta$covariate_names)
  # random_state lives in `provenance$rng_seed`; `estimand` is a list with
  # field `type`.
  expect_equal(fit2$meta$provenance$rng_seed, fit$meta$provenance$rng_seed)
  expect_equal(fit2$meta$estimand$type, fit$meta$estimand$type)
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
    extra <- if (m == "tarnet") list(n_iter = 50L) else list()
    fit <- suppressWarnings(suppressMessages(do.call(cf_fit, c(
      list(Y ~ T | X1 + X2, data = df, method = m, partition = "none",
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
