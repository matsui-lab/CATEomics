## Regression coverage of post-fit functions against models fit through the
## named-args facade.
##
## The facade is a pure entry-point convertor: named-args are resolved into the
## canonical formula and the rest of the pipeline (Steps 1-7, base learners,
## dispatch, prediction, meta-recording) is untouched. These tests therefore
## verify that downstream API functions (predict / cf_evaluate / cf_interpret /
## cf_rule) operate correctly on a model fit with named-args (i.e., on a result
## whose `meta$named_args` is populated and whose `meta$propensity_X` may be
## set).

library(testthat)

.make_df <- function(n = 120L, p = 6L, seed = 42L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  T <- stats::rbinom(n, 1, 0.5)
  Y <- as.numeric(T) * 1.5 + rowSums(X) + 0.1 * stats::rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  colnames(df) <- c("Y", "T", paste0("X", seq_len(p)))
  df
}

.fit_named <- function(method = "grf", ...) {
  df <- .make_df()
  args <- list(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
    data = df, method = method,
    random_state = 42L, partition = "oof", n_folds = 3L
  )
  if (identical(method, "grf")) args$params <- list(tune.parameters = "none")
  args <- c(args, list(...))
  suppressWarnings(suppressMessages(do.call(cf_fit, args)))
}

test_that("predict newdata works on a named-args fitted model", {
  fit <- .fit_named("grf")
  expect_false(is.null(fit$meta$named_args))
  nd <- .make_df()[1:10, paste0("X", 1:6)]
  pred <- predict(fit, newdata = nd, type = "all")
  expect_true(is.list(pred))
  expect_true(is.finite(pred$ate))
  expect_length(pred$tau_hat, 10L)
})

test_that("cf_evaluate runs on a named-args fitted model", {
  fit <- .fit_named("grf")
  ev <- suppressWarnings(cf_evaluate(fit, metrics = "rate",
                                       target_quantile = 0.2))
  expect_true(is.list(ev))
  # Either flat metrics or a metrics sub-list — both shapes are accepted.
  m <- ev$metrics %||% ev
  expect_true(length(m) > 0L)
})

test_that("cf_interpret runs on a named-args fitted model", {
  fit <- .fit_named("grf")
  itp <- suppressWarnings(cf_interpret(fit, surrogate = "ridge"))
  expect_s3_class(itp, "cateomics_interpret")
  expect_true(!is.null(itp$vi))
})

test_that("cf_rule runs on a named-args fitted model", {
  fit <- .fit_named("grf")
  rule <- suppressWarnings(cf_rule(fit, type = "threshold"))
  expect_true(inherits(rule, "cateomics_rule") || is.list(rule))
})

test_that("propensity_X is preserved on $meta after K-fold cross-fit", {
  # `.fit_named` already defaults to partition='oof' + n_folds=3L; just
  # pass propensity_X as the only extra arg (do.call would otherwise raise
  # "matched by multiple actual arguments" for the duplicates).
  fit <- suppressWarnings(.fit_named("grf", propensity_X = c("X1", "X3")))
  expect_equal(fit$meta$propensity_X, c("X1", "X3"))
})
