## Tests for the named-args facade (dual API entry).
##
## Verifies that cf_fit / cf_compare accept both the legacy
## formula API and the new named-args API, that the two routes produce
## identical results for the same inputs, and that the validation surface
## (mode conflict, missing input, X slot type checks) raises the expected
## cfomics_* error classes.

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

test_that("facade: legacy formula API and named-args API are equivalent", {
  df <- .make_demo()
  fmla <- stats::as.formula("Y ~ T | X1 + X2 + X3 + X4")

  fit_legacy <- suppressWarnings(
    cf_fit(fmla, data = df, method = "grf",
           random_state = 42L, partition = "none")
  )
  fit_named <- cf_fit(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:4),
    data = df, method = "grf",
    random_state = 42L, partition = "none"
  )

  expect_s3_class(fit_legacy, "cf_model")
  expect_s3_class(fit_named, "cf_model")
  # Same RNG seed + same inputs → identical ATE
  expect_identical(fit_legacy$res$ate, fit_named$res$ate)
  # Legacy path leaves named_args = NULL; new path populates it
  expect_null(fit_legacy$meta$named_args)
  expect_false(is.null(fit_named$meta$named_args))
  expect_equal(fit_named$meta$named_args$Y, "Y")
  expect_equal(fit_named$meta$named_args$T, "T")
  expect_equal(fit_named$meta$named_args$X_baseline, paste0("X", 1:4))
})

test_that("facade: mode conflict (both formula and named args) errors", {
  df <- .make_demo()
  fmla <- stats::as.formula("Y ~ T | X1 + X2")
  expect_error(
    cf_fit(formula = fmla, Y = "Y", T = "T",
           X_baseline = paste0("X", 1:4),
           data = df, method = "grf"),
    class = "cfomics_input_mode_conflict"
  )
})

test_that("facade: neither formula nor Y/T supplied errors", {
  df <- .make_demo()
  expect_error(
    cf_fit(data = df, method = "grf"),
    class = "cfomics_invalid_input"
  )
})

test_that("facade: invalid Y / T argument types are rejected", {
  df <- .make_demo()
  # Y as numeric vector instead of character scalar
  expect_error(
    cf_fit(Y = 1L, T = "T", X_baseline = paste0("X", 1:4),
           data = df, method = "grf"),
    class = "cfomics_invalid_Y"
  )
  # T as character vector of length 2
  expect_error(
    cf_fit(Y = "Y", T = c("T1", "T2"), X_baseline = paste0("X", 1:4),
           data = df, method = "grf"),
    class = "cfomics_invalid_T"
  )
})

test_that("facade: at least one X slot must be supplied", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T", data = df, method = "grf"),
    class = "cfomics_no_X_supplied"
  )
})

test_that("facade: X slot type checks", {
  df <- .make_demo()
  # X_baseline as numeric instead of character
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = 1:3, data = df, method = "grf"),
    class = "cfomics_invalid_X_slot"
  )
  # X_assays as un-named list — for data.frame input, `.resolve_X_columns`
  # rejects X_assays first (data.frame has no assay concept), so the
  # actual error class is `cfomics_X_assays_unsupported_df`. The
  # `cfomics_invalid_X_assays` class fires only when X_assays unnamed is
  # supplied via SE / MAE input (where the slot is supported in principle
  # but the list shape is invalid). Both errors are user-actionable. The class
  # assertion matches the actual order of validation.
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = paste0("X", 1:4),
           X_assays = list(c("ATAC1")),  # unnamed
           data = df, method = "grf"),
    class = "cfomics_X_assays_unsupported_df"
  )
})

test_that("facade: Y / T cannot appear among X slots", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("Y", "X1", "X2"),
           data = df, method = "grf"),
    class = "cfomics_Y_in_X"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("T", "X1"),
           data = df, method = "grf"),
    class = "cfomics_T_in_X"
  )
})

test_that("facade: duplicate X column names across slots are rejected", {
  df <- .make_demo()
  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = c("X1", "X2"),
           X_rnaseq   = c("X2", "X3"),  # X2 duplicated
           data = df, method = "grf"),
    class = "cfomics_duplicate_X"
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
  expect_s3_class(cmp, "cfomics_compare")
  expect_true("grf" %in% cmp$scoreboard$method)
})

test_that("facade: legacy formula API still emits deprecation warning once", {
  df <- .make_demo()
  fmla <- stats::as.formula("Y ~ T | X1 + X2 + X3 + X4")
  # `rlang::warn(.frequency = "once")` tracks the first emission across the
  # whole R session by the .frequency_id, so an earlier test may have already
  # consumed the warning. Use `rlang::reset_warning_verbosity()` to clear the
  # tracker for the cfomics_formula_deprecation id before the assertion.
  if (utils::packageVersion("rlang") >= "1.0.0") {
    tryCatch(rlang::reset_warning_verbosity("cfomics_formula_deprecation"),
             error = function(e) NULL)
  }
  expect_warning(
    cf_fit(fmla, data = df, method = "grf",
           random_state = 42L, partition = "none"),
    class = "cfomics_formula_deprecated"
  )
})
