# Smoke tests for cf_model_stability + cf_plot_stability.
# Goals (smoke only):
#   - 3 correlation matrices (Pearson / Spearman / Kendall)
#   - per-method mean tau correlation vector
#   - cf_plot_stability renders for each correlation method

set.seed(2026)
.make_mock_compare_stability <- function(n = 60L, K = 5L,
                                          include_outlier = FALSE) {
  # Build correlated tau columns to produce a non-degenerate stability output.
  base <- rnorm(n, sd = 1)
  tau_mat <- sapply(seq_len(K), function(k) base + rnorm(n, sd = 0.4))
  if (include_outlier) {
    # Replace last column with uncorrelated noise -> outlier method.
    tau_mat[, K] <- rnorm(n, sd = 1)
  }
  colnames(tau_mat) <- paste0("m", seq_len(K))

  Y <- rnorm(n)
  T_vec <- rbinom(n, 1L, 0.5)
  mu0 <- rnorm(n, sd = 0.3); mu1 <- mu0 + rnorm(n, sd = 0.2)
  pi  <- rep(0.5, n)
  models <- list()
  models[["m1"]] <- list(data = list(Y = Y, T = T_vec,
                                     X = matrix(0, n, 2)))

  structure(
    list(
      scoreboard = data.frame(method = colnames(tau_mat)),
      models     = models,
      metadata   = list(
        n_folds        = 5L,
        nuisance_oof   = list(mu0_oof = mu0, mu1_oof = mu1, pi_oof = pi),
        oof_tau_matrix = tau_mat
      )
    ),
    class = c("cateomics_compare", "list")
  )
}


# Core output contract ------------------------------------------------------

test_that("cf_model_stability returns the expected list shape", {
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  out <- suppressMessages(cf_model_stability(cmp))

  expect_s3_class(out, "cateomics_model_stability")
  # One thing is read (the K columns of out-of-fold tau-hat) and the
  # correlations between them are reported; which samples each method would
  # call high-benefit is cf_benefit_groups' question.
  expect_setequal(names(out),
                  c("cor_matrices", "per_method_tau_correlation", "meta"))
  expect_null(out[["benefit_group_jaccard"]])
  expect_null(out[["outlier_methods"]])
  expect_null(out[["summary"]])

  # 3 cor matrices requested by default.
  expect_setequal(names(out$cor_matrices),
                  c("pearson", "spearman", "kendall"))
  for (cm in names(out$cor_matrices)) {
    expect_equal(dim(out$cor_matrices[[cm]]), c(4L, 4L))
    expect_true(all(abs(diag(out$cor_matrices[[cm]]) - 1) < 1e-12))
  }

  expect_equal(length(out$per_method_tau_correlation), 4L)
  # The counts sit with the other settings; the aggregates are computed where
  # they are shown.
  expect_identical(length(out$meta$methods_used), 4L)
  expect_identical(out$meta$n_rows, 50L)
  expect_null(out$meta$cor_methods)   # = names(cor_matrices)
  expect_null(out$meta$benefit_top_q)
  expect_null(out$meta$outlier_mad)
  # The arguments behind the removed quantities are gone too.
  expect_error(suppressMessages(cf_model_stability(cmp, benefit_top_q = 0.25)),
               class = "cateomics_unknown_argument")
  expect_error(suppressMessages(cf_model_stability(cmp, outlier_mad = 2)),
               class = "cateomics_unknown_argument")
})


test_that("the least agreeing method has the lowest mean tau correlation", {
  cmp <- .make_mock_compare_stability(n = 80L, K = 5L,
                                       include_outlier = TRUE)
  out <- suppressMessages(cf_model_stability(cmp))
  # m5 is the noise column. The vector says so; the function does not decide
  # for the reader where the cut between "agreeing" and "not" lies.
  v <- out$per_method_tau_correlation
  expect_identical(names(v)[which.min(v)], "m5")
})


# Cleaning + input validation ----------------------------------------------

test_that("cf_model_stability handles failed methods via OOF-matrix cleaning", {
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  cmp$metadata$oof_tau_matrix[, "m1"] <- NA_real_
  out <- suppressMessages(cf_model_stability(cmp))
  expect_true("m1" %in% out$meta$methods_dropped)
  expect_equal(length(out$meta$methods_used), 3L)
})

test_that("cf_model_stability refuses a comparison with one usable method", {
  # Everything this function reports is an agreement between methods, so one
  # surviving column leaves nothing to estimate: the correlation and Jaccard
  # matrices would be 1 by definition and every agreement statistic NA.
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  for (m in c("m2", "m3", "m4")) cmp$metadata$oof_tau_matrix[, m] <- NA_real_
  expect_error(suppressMessages(cf_model_stability(cmp)),
               class = "cateomics_stability_single_method")
  # Two usable columns is the smallest comparison it accepts.
  cmp2 <- .make_mock_compare_stability(n = 50L, K = 4L)
  for (m in c("m3", "m4")) cmp2$metadata$oof_tau_matrix[, m] <- NA_real_
  out <- suppressMessages(cf_model_stability(cmp2))
  expect_equal(length(out$meta$methods_used), 2L)
  expect_true(all(is.finite(out$per_method_tau_correlation)))
  # Kendall on a single column is where this used to fail: R's own
  # symmetrisation sizes a matrix from the diagonal, and a self-correlation of
  # 1 - 2.2e-16 collapses it to 0 x 0.
  expect_true(all(is.finite(out$cor_matrices$kendall)))
})

test_that("cf_model_stability errors on non-cateomics_compare", {
  expect_error(
    cf_model_stability(list()),
    class = "cateomics_invalid_compare_result"
  )
})


# cf_plot_stability rendering ---------------------------------------------

test_that("cf_plot_stability renders the correlation heatmap for each cor_method", {
  skip_if_not_installed("ggplot2")
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  out <- suppressMessages(cf_model_stability(cmp))

  # The stability plot is the cor heatmap, selectable via cor_method.
  expect_s3_class(cf_plot_stability(out), "ggplot")
  for (cm in c("pearson", "spearman", "kendall")) {
    expect_s3_class(cf_plot_stability(out, cor_method = cm), "ggplot")
  }
})

test_that("cf_plot_stability errors on non-stability input", {
  expect_error(
    cf_plot_stability(list(), type = "cor"),
    class = "cateomics_invalid_input"
  )
})


# print / summary smoke ----------------------------------------------------

test_that("print and summary methods do not error", {
  cmp <- .make_mock_compare_stability(n = 40L, K = 3L)
  out <- suppressMessages(cf_model_stability(cmp))
  expect_output(print(out), "<cateomics_model_stability>")
  # The print reports what was measured; it carries no advice on how to use it.
  expect_output(print(out), "primary_cor")
  # What is reported is a correlation of the estimates, so that is what it is
  # called.
  expect_output(print(out), "tau correlation")
  expect_false(any(grepl("agreement",
                         capture.output(print(out)), fixed = TRUE)))
  expect_false(any(grepl("outlier",
                         capture.output(print(out)), fixed = TRUE)))
  s <- summary(out)
  expect_s3_class(s, "summary.cateomics_model_stability")
  expect_output(print(s), "cf_model_stability summary")
})


test_that("with two methods print gives the one number, and names it", {
  mock <- .make_mock_compare_stability(n = 60L, K = 2L)
  out <- suppressMessages(cf_model_stability(mock))
  txt <- paste(capture.output(print(out)), collapse = "\n")
  # One pair, so one number -- and on its own line it has to say which
  # correlation it is.
  expect_match(txt, "pairwise tau correlation (spearman):", fixed = TRUE)
  expect_false(grepl("per-method mean tau correlation", txt, fixed = TRUE))
  expect_false(grepl("pairwise range", txt, fixed = TRUE))
})
