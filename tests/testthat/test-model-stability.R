# Smoke tests for cf_model_stability + cf_plot_stability.
# Goals (smoke only):
#   - 3 correlation matrices (Pearson / Spearman / Kendall)
#   - Jaccard high-benefit matrix
#   - per-method agreement vector
#   - outlier detection (descriptive)
#   - cf_plot_stability 3 plot types render without error

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
    class = c("cfomics_compare", "list")
  )
}


# Core output contract ------------------------------------------------------

test_that("cf_model_stability returns the expected list shape", {
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  out <- suppressMessages(cf_model_stability(cmp))

  expect_s3_class(out, "cfomics_model_stability")
  expect_true(all(c("cor_matrices", "benefit_group_jaccard",
                    "per_method_tau_correlation", "outlier_methods",
                    "summary", "meta") %in% names(out)))

  # 3 cor matrices requested by default.
  expect_setequal(names(out$cor_matrices),
                  c("pearson", "spearman", "kendall"))
  for (cm in names(out$cor_matrices)) {
    expect_equal(dim(out$cor_matrices[[cm]]), c(4L, 4L))
    expect_true(all(abs(diag(out$cor_matrices[[cm]]) - 1) < 1e-12))
  }

  expect_equal(dim(out$benefit_group_jaccard), c(4L, 4L))
  expect_true(all(abs(diag(out$benefit_group_jaccard) - 1) < 1e-12))

  expect_equal(length(out$per_method_tau_correlation), 4L)
  expect_true(all(out$summary$n_methods == 4L))
})


test_that("cf_model_stability respects benefit_top_q bounds", {
  cmp <- .make_mock_compare_stability(n = 50L, K = 3L)
  expect_error(
    suppressMessages(cf_model_stability(cmp, benefit_top_q = 0)),
    class = "cfomics_stability_invalid_top_q"
  )
  expect_error(
    suppressMessages(cf_model_stability(cmp, benefit_top_q = 0.8)),
    class = "cfomics_stability_invalid_top_q"
  )
})


# Outlier detection (descriptive) ------------------------------------------

test_that("cf_model_stability flags an obviously-uncorrelated method", {
  cmp <- .make_mock_compare_stability(n = 80L, K = 5L,
                                       include_outlier = TRUE)
  out <- suppressMessages(cf_model_stability(cmp, outlier_mad = 0.5))
  # m5 is the noise column; should be flagged as outlier.
  expect_true("m5" %in% out$outlier_methods)
})


# Cleaning + input validation ----------------------------------------------

test_that("cf_model_stability handles failed methods via OOF-matrix cleaning", {
  cmp <- .make_mock_compare_stability(n = 50L, K = 4L)
  cmp$metadata$oof_tau_matrix[, "m1"] <- NA_real_
  out <- suppressMessages(cf_model_stability(cmp))
  expect_true("m1" %in% out$meta$methods_dropped)
  expect_equal(length(out$meta$methods_used), 3L)
})

test_that("cf_model_stability errors on non-cfomics_compare", {
  expect_error(
    cf_model_stability(list()),
    class = "cfomics_invalid_compare_result"
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
    class = "cfomics_plot_stability_invalid_input"
  )
})


# print / summary smoke ----------------------------------------------------

test_that("print and summary methods do not error", {
  cmp <- .make_mock_compare_stability(n = 40L, K = 3L)
  out <- suppressMessages(cf_model_stability(cmp))
  expect_output(print(out), "<cfomics_model_stability>")
  expect_output(print(out), "descriptive diagnostic")
  s <- summary(out)
  expect_s3_class(s, "summary.cfomics_model_stability")
  expect_output(print(s), "cf_model_stability summary")
})
