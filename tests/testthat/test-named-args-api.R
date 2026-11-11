## Tests for the named-argument entry of cf_fit / cf_compare.
##
## Verifies that the covariate slots are read as documented and that the
## validation surface (missing input, X slot type checks) raises the expected
## cateomics_* error classes.

library(testthat)

# ---- Small fixture ------------------------------------------------------
.make_demo <- function(n = 100L, p = 4L, seed = 42L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  T <- stats::rbinom(n, 1, 0.5)
  Y <- as.numeric(T) * 1.5 + rowSums(X) + 0.1 * stats::rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  colnames(df) <- c("Y", "T", paste0("X", seq_len(p)))
  df
}

test_that("facade: the entry records the names it was given", {
  df <- .make_demo()
  fit <- cf_fit(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:4),
    data = df, method = "grf",
    random_state = 42L, partition = "none"
  )
  expect_s3_class(fit, "cf_model")
  expect_false(is.null(fit$meta$named_args))
  expect_equal(fit$meta$named_args$Y, "Y")
  expect_equal(fit$meta$named_args$T, "T")
  expect_equal(fit$meta$named_args$X_baseline, paste0("X", 1:4))
})

test_that("facade: neither Y nor T supplied errors", {
  df <- .make_demo()
  expect_error(
    cf_fit(data = df, method = "grf"),
    class = "cateomics_invalid_input"
  )
})

test_that("facade: invalid Y / T argument types are rejected", {
  df <- .make_demo()
  # Y as numeric vector instead of character scalar
  expect_error(
    cf_fit(Y = 1L, T = "T", X_baseline = paste0("X", 1:4),
           data = df, method = "grf"),
    class = "cateomics_invalid_Y"
  )
  # T as character vector of length 2
  expect_error(
    cf_fit(Y = "Y", T = c("T1", "T2"), X_baseline = paste0("X", 1:4),
           data = df, method = "grf"),
    class = "cateomics_invalid_T"
  )
})

test_that("facade: at least one X slot must be supplied", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T", data = df, method = "grf"),
    class = "cateomics_no_X_supplied"
  )
})

test_that("facade: X slot type checks", {
  df <- .make_demo()
  # X_baseline as numeric instead of character
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = 1:3, data = df, method = "grf"),
    class = "cateomics_invalid_X_slot"
  )
  # `omics` as an unnamed list: every element must say which layer it is.
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = paste0("X", 1:4),
           omics = list(c("ATAC1")),  # unnamed
           data = df, method = "grf"),
    class = "cateomics_invalid_omics"
  )
})

test_that("facade: Y / T cannot appear among X slots", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("Y", "X1", "X2"),
           data = df, method = "grf"),
    class = "cateomics_Y_in_X"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("T", "X1"),
           data = df, method = "grf"),
    class = "cateomics_T_in_X"
  )
})

test_that("facade: a baseline column repeated in the omics block is rejected", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("X1", "X2"),
           omics      = c("X2", "X3"),  # X2 is already a baseline covariate
           preprocess = FALSE,
           data = df, method = "grf"),
    class = "cateomics_column_conflict"
  )
})

test_that("facade: propensity_X is accepted + stored on meta", {
  df <- .make_demo()
  fit <- cf_fit(Y = "Y", T = "T",
                X_baseline = paste0("X", 1:4),
                propensity_X = c("X1", "X2"),
                data = df, method = "grf",
                random_state = 42L, partition = "none")
  expect_equal(fit$meta$propensity_X, c("X1", "X2"))
})

test_that("facade: cf_compare accepts named-args", {
  df <- .make_demo()
  cmp <- suppressWarnings(cf_compare(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:4),
    data = df, methods = c("grf"),
    random_state = 42L, n_folds = 3L
  ))
  expect_s3_class(cmp, "cateomics_compare")
  expect_true("grf" %in% cmp$scoreboard$method)
})

