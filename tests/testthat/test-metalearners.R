test_that("S/T/X-learner smoke: fit, predict, contract compliance", {
  skip_if_no_cf_python()

  df <- sim_simple_linear(n = 200L, p = 5L, seed = 123L, tau = 1.5)
  cov_names <- paste0("X", 1:5)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  for (meth in c("slearner", "tlearner", "xlearner")) {
    fit <- suppressMessages(
      cf_fit(fml, data = df, method = meth, random_state = 42L)
    )

    expect_s3_class(fit, "cf_model")
    expect_equal(fit$method, meth)

    ate <- predict(fit, type = "ate")
    ite <- predict(fit, type = "tau")
    y0  <- predict(fit, type = "mu0")
    y1  <- predict(fit, type = "mu1")
    s   <- predict(fit, type = "summary")

    expect_length(ate, 1L)
    expect_true(is.finite(ate))
    expect_length(ite, nrow(df))
    expect_length(y0, nrow(df))
    expect_length(y1, nrow(df))
    expect_true(all(is.finite(ite)))
    expect_true(all(is.finite(y0)))
    expect_true(all(is.finite(y1)))

    expect_true(all(c("ate", "tau_mean", "tau_sd", "tau_iqr")
                    %in% names(s)))

    # Contract: 4-class character.
    # - S/T-learner (Class A): ite == mu1_hat - mu0_hat exactly
    # - X-learner (Class B'): mu1 - mu0 != ite by design (Stage 3 blend),
    #   roxygen note added; smoke only checks non-NA + finite.
    if (meth %in% c("slearner", "tlearner")) {
      expect_equal(ite, y1 - y0, tolerance = 1e-6)
    } else {
      # Class B' (xlearner): mu1 - mu0 is the Stage-1 plug-in; ite is the
      # Stage-3 propensity-weighted blend. They may differ in finite samples.
      expect_true(all(is.finite(ite)))
      expect_true(all(is.finite(y1 - y0)))
    }
  }
})

test_that("S/T/X-learner accept newdata in predict", {
  skip_if_no_cf_python()

  df <- sim_simple_linear(n = 150L, p = 4L, seed = 7L, tau = 1.0)
  cov_names <- paste0("X", 1:4)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- suppressMessages(
    cf_fit(fml, data = df, method = "tlearner", random_state = 1L)
  )

  newX <- df[1:10, cov_names, drop = FALSE]
  expect_length(predict(fit, newdata = newX, type = "tau"), 10L)
  expect_length(predict(fit, newdata = newX, type = "mu0"), 10L)
  expect_length(predict(fit, newdata = newX, type = "mu1"), 10L)
  expect_length(predict(fit, newdata = newX, type = "ate"), 1L)
})

test_that("X-learner accepts model_propensity argument", {
  skip_if_no_cf_python()

  df <- sim_simple_linear(n = 150L, p = 4L, seed = 11L)
  cov_names <- paste0("X", 1:4)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- suppressMessages(
    cf_fit(fml, data = df, method = "xlearner",
           model_outcome = "linear", model_propensity = "logistic",
           random_state = 3L)
  )

  md <- fit$res$metadata
  expect_equal(md$model_outcome, "linear")
  expect_equal(md$model_propensity, "logistic")
  expect_equal(md$backend, "econml.metalearners.XLearner")
})
