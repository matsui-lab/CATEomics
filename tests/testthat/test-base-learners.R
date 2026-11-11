# Tests for the unified base learner infrastructure (.fit_base_learner /
# .predict_base_learner). Covers:
#   * the learners this file exercises, for both outcome types
#   * reproducibility via random_state
#   * input validation
#   * predict output shape / type

make_beta <- function(p) {
  # Deterministic, non-zero-dominated coefficient vector of any length p.
  base <- c(1, -0.5, 0.3, 0.2, -0.4, 0.1)
  rep_len(base, p)
}

make_reg_data <- function(n = 80, p = 4, seed = 1L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  y <- as.numeric(x %*% make_beta(p) + stats::rnorm(n, sd = 0.5))
  list(x = x, y = y)
}

make_cls_data <- function(n = 80, p = 4, seed = 1L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  lp <- as.numeric(x %*% make_beta(p))
  pr <- 1 / (1 + exp(-lp))
  y <- as.integer(stats::runif(n) < pr)
  list(x = x, y = y)
}

learners <- c("linear", "lasso", "ridge", "elastic_net", "rf")

test_that("all learners fit and predict for regression", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")

  d <- make_reg_data()
  for (lrn in learners) {
    obj <- .fit_base_learner(d$x, d$y, type = "regression",
                             learner = lrn, random_state = 42L)
    expect_s3_class(obj, "cf_base_learner")
    expect_identical(obj$learner, lrn)
    expect_identical(obj$type, "regression")
    expect_identical(obj$feature_names, colnames(d$x))

    pred <- .predict_base_learner(obj, d$x)
    expect_type(pred, "double")
    expect_length(pred, nrow(d$x))
    expect_true(all(is.finite(pred)))
  }
})

test_that("all learners fit and predict for classification", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")

  d <- make_cls_data()
  # "logistic" is an alias for "linear"; cover it here.
  for (lrn in c("logistic", "lasso", "ridge", "elastic_net", "rf")) {
    obj <- .fit_base_learner(d$x, d$y, type = "classification",
                             learner = lrn, random_state = 42L)
    expect_s3_class(obj, "cf_base_learner")
    expect_identical(obj$type, "classification")

    pred <- .predict_base_learner(obj, d$x)
    expect_type(pred, "double")
    expect_length(pred, nrow(d$x))
    expect_true(all(is.finite(pred)))
    expect_true(all(pred >= 0 & pred <= 1))
  }
})

test_that("random_state yields reproducible fits", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")

  d <- make_reg_data(n = 120, p = 5, seed = 7L)

  for (lrn in c("lasso", "ridge", "elastic_net", "rf")) {
    f1 <- .fit_base_learner(d$x, d$y, type = "regression",
                            learner = lrn, random_state = 123L)
    f2 <- .fit_base_learner(d$x, d$y, type = "regression",
                            learner = lrn, random_state = 123L)
    p1 <- .predict_base_learner(f1, d$x)
    p2 <- .predict_base_learner(f2, d$x)
    expect_equal(p1, p2, tolerance = 1e-10,
                 info = sprintf("learner = %s", lrn))
  }
})

test_that("data.frame input is accepted", {
  skip_if_not_installed("glmnet")

  d <- make_reg_data()
  df <- as.data.frame(d$x)
  obj <- .fit_base_learner(df, d$y, type = "regression", learner = "lasso",
                           random_state = 1L)
  pred <- .predict_base_learner(obj, df)
  expect_length(pred, nrow(df))
})

test_that("classification accepts two-level factor y", {
  skip_if_not_installed("glmnet")

  d <- make_cls_data()
  yf <- factor(d$y, levels = c(0, 1))
  obj <- .fit_base_learner(d$x, yf, type = "classification",
                           learner = "lasso", random_state = 1L)
  pred <- .predict_base_learner(obj, d$x)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("invalid learner spec errors", {
  d <- make_reg_data()
  # xgboost is part of the valid vocabulary — use a genuinely bogus token here
  # instead.
  expect_error(
    .fit_base_learner(d$x, d$y, type = "regression", learner = "bogus_learner"),
    class = "cateomics_base_learner_error"
  )
  expect_error(
    .fit_base_learner(d$x, d$y, type = "regression",
                      learner = c("linear", "lasso")),
    class = "cateomics_base_learner_error"
  )
})

test_that("xgboost learner fits and predicts (smoke; skips if package absent)", {
  skip_if_not_installed("xgboost")
  d <- make_reg_data()
  obj <- .fit_base_learner(d$x, d$y, type = "regression",
                           learner = "xgboost", random_state = 1L,
                           nrounds = 10L)
  expect_s3_class(obj, "cf_base_learner")
  pred <- .predict_base_learner(obj, d$x)
  expect_length(pred, length(d$y))
  expect_true(all(is.finite(pred)))
})

test_that("non-numeric x is rejected", {
  x_bad <- data.frame(a = letters[1:5], b = 1:5)
  expect_error(
    .fit_base_learner(x_bad, 1:5, type = "regression", learner = "linear"),
    class = "cateomics_base_learner_error"
  )
})

test_that("classification rejects non-0/1 numeric y", {
  d <- make_cls_data()
  y_bad <- d$y + 2L  # values in {2, 3}
  expect_error(
    .fit_base_learner(d$x, y_bad, type = "classification", learner = "lasso"),
    class = "cateomics_base_learner_error"
  )
})

test_that("predict errors on non-cf_base_learner object", {
  expect_error(
    .predict_base_learner(list(), matrix(1:4, 2, 2)),
    class = "cateomics_base_learner_error"
  )
})

test_that("predict aligns columns by name", {
  skip_if_not_installed("glmnet")

  d <- make_reg_data()
  obj <- .fit_base_learner(d$x, d$y, type = "regression", learner = "lasso",
                           random_state = 1L)

  # Permute columns: prediction should be invariant to column order.
  perm <- d$x[, rev(colnames(d$x)), drop = FALSE]
  p_orig <- .predict_base_learner(obj, d$x)
  p_perm <- .predict_base_learner(obj, perm)
  expect_equal(p_orig, p_perm, tolerance = 1e-10)

  # Missing column -> error.
  expect_error(
    .predict_base_learner(obj, d$x[, -1, drop = FALSE]),
    class = "cateomics_base_learner_error"
  )
})
