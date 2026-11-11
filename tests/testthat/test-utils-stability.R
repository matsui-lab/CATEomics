# Smoke tests for the oof-matrix cleaner / correlation-knee helpers: verify return-structure contracts on
# small synthetic inputs only; no PEHE / accuracy checks.

set.seed(123)
.make_test_oof_tau_matrix <- function(n = 50L, K = 4L, n_failed = 0L) {
  mat <- matrix(rnorm(n * K), nrow = n, ncol = K)
  colnames(mat) <- paste0("m", seq_len(K))
  if (n_failed > 0L) {
    # First `n_failed` columns become all-NA (mimics failed methods in
    # cf_compare).
    for (j in seq_len(n_failed)) mat[, j] <- NA_real_
  }
  mat
}


# .cateomics_clean_oof_tau_matrix --------------------------------------------

test_that(".cateomics_clean_oof_tau_matrix drops failed methods + incomplete rows", {
  mat <- .make_test_oof_tau_matrix(n = 30L, K = 5L, n_failed = 2L)
  # Introduce an NA in 1 row of a surviving column to test row drop.
  mat[5L, 3L] <- NA_real_

  out <- .cateomics_clean_oof_tau_matrix(mat)

  expect_equal(out$status, "ok")
  expect_equal(length(out$kept_methods), 3L)
  expect_equal(length(out$dropped_methods), 2L)
  expect_setequal(out$dropped_methods, c("m1", "m2"))
  expect_equal(out$n_rows_dropped, 1L)
  expect_equal(out$n_rows_kept, 29L)
  expect_equal(ncol(out$mat), 3L)
  expect_true(all(stats::complete.cases(out$mat)))
})

test_that(".cateomics_clean_oof_tau_matrix returns no_input for NULL / wrong type", {
  out1 <- .cateomics_clean_oof_tau_matrix(NULL)
  expect_equal(out1$status, "no_input")
  expect_null(out1$mat)

  out2 <- .cateomics_clean_oof_tau_matrix(c(1, 2, 3))
  expect_equal(out2$status, "no_input")
})

test_that(".cateomics_clean_oof_tau_matrix returns no_columns when all methods fail", {
  mat <- .make_test_oof_tau_matrix(n = 20L, K = 3L, n_failed = 3L)
  out <- .cateomics_clean_oof_tau_matrix(mat)
  expect_equal(out$status, "no_columns")
  expect_null(out$mat)
})


# .cateomics_tau_correlation_knee --------------------------------------------

test_that(".cateomics_tau_correlation_knee returns expected list shape", {
  mat <- .make_test_oof_tau_matrix(n = 80L, K = 5L)
  out <- .cateomics_tau_correlation_knee(mat, cor_method = "spearman")

  expect_true(is.matrix(out$cor_matrix))
  expect_equal(dim(out$cor_matrix), c(5L, 5L))
  expect_true(all(abs(diag(out$cor_matrix) - 1) < 1e-12))
  expect_equal(length(out$ordering), 5L)
  expect_setequal(out$ordering, seq_len(5L))
  expect_equal(length(out$sorted_mean_cor), 5L)
  # sorted_mean_cor must be non-increasing (descending order).
  expect_true(all(diff(out$sorted_mean_cor) <= 0 + 1e-12))
  expect_true(out$knee_k >= 2L && out$knee_k <= 5L)
  expect_equal(out$cor_method, "spearman")
})

test_that(".cateomics_tau_correlation_knee handles K < 2 gracefully", {
  mat <- matrix(rnorm(20), nrow = 20L, ncol = 1L)
  colnames(mat) <- "m1"
  out <- .cateomics_tau_correlation_knee(mat)
  expect_null(out$cor_matrix)
  expect_null(out$ordering)
  expect_true(is.na(out$knee_k))
})

test_that(".cateomics_tau_correlation_knee supports all 3 cor methods", {
  mat <- .make_test_oof_tau_matrix(n = 40L, K = 3L)
  for (cm in c("kendall", "spearman", "pearson")) {
    out <- .cateomics_tau_correlation_knee(mat, cor_method = cm)
    expect_equal(out$cor_method, cm)
    expect_equal(dim(out$cor_matrix), c(3L, 3L))
  }
})

test_that(".cateomics_locate_knee respects min_knee bound", {
  # Sharp drop at position 2.
  x <- c(0.9, 0.85, 0.1, 0.05, 0.02)
  expect_equal(.cateomics_locate_knee(x, min_knee = 2L), 2L)

  # Largest drop at position 4 (penultimate to last).
  y <- c(0.9, 0.88, 0.85, 0.8, 0.2)
  expect_equal(.cateomics_locate_knee(y, min_knee = 2L), 4L)

  # When length <= min_knee, return length.
  z <- c(0.9, 0.8)
  expect_equal(.cateomics_locate_knee(z, min_knee = 2L), 2L)
})


# cor_method selects the correlation, it is not just echoed back ----------

test_that("the requested correlation method is the one that is computed", {
  set.seed(4)
  a <- rnorm(100)
  # An exponential column is monotone in `a` but not linear, so Pearson and
  # Spearman must disagree. Computing Pearson whatever was asked for would
  # make every stored matrix identical.
  mat <- cbind(m1 = a, m2 = exp(a), m3 = a + rnorm(100, sd = 0.3))

  out <- .cateomics_tau_correlation_knee(mat, cor_method = "pearson")
  out_s <- .cateomics_tau_correlation_knee(mat, cor_method = "spearman")
  expect_equal(out$cor_matrix[1, 2], stats::cor(a, exp(a)), tolerance = 1e-10)
  expect_equal(out_s$cor_matrix[1, 2], 1, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(out$cor_matrix[1, 2], out_s$cor_matrix[1, 2])))
})
