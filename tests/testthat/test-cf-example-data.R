test_that("cf_example_data returns a valid SummarizedExperiment with ground truth", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  suppressPackageStartupMessages(library(SummarizedExperiment))

  se <- cf_example_data(n = 60L, p = 200L, n_prog = 20L, n_pred = 10L, seed = 1L)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(200L, 60L))

  # Counts: integer, non-negative (valid for lognorm / edgeR).
  cts <- SummarizedExperiment::assay(se, "counts")
  expect_true(is.integer(cts))
  expect_true(all(cts >= 0))

  # Real gene symbols as row names, no duplicates.
  expect_type(rownames(se), "character")
  expect_false(anyDuplicated(rownames(se)) > 0)

  # colData carries treatment / outcome / clinical / truth.
  cd <- SummarizedExperiment::colData(se)
  expect_true(all(c("treatment", "outcome", "age", "sex", "cf_tau_true") %in%
                    colnames(cd)))
  expect_true(all(cd$treatment %in% c(0L, 1L)))
  expect_true(is.factor(cd$sex))
  expect_setequal(levels(cd$sex), c("F", "M"))
  expect_length(cd$cf_tau_true, 60L)

  # rowData truth flags: exactly n_pred predictive, n_prog prognostic, disjoint.
  rd <- SummarizedExperiment::rowData(se)
  expect_true(all(c("is_prognostic", "is_predictive") %in% colnames(rd)))
  expect_equal(sum(rd$is_predictive), 10L)
  expect_equal(sum(rd$is_prognostic), 20L)
  expect_equal(sum(rd$is_predictive & rd$is_prognostic), 0L)

  # metadata: embedded pathway + truth.
  md <- S4Vectors::metadata(se)
  expect_true(all(c("cf_pathway", "cf_predictive_genes", "cf_prognostic_genes",
                    "cf_true_ate", "cf_snr") %in% names(md)))
  expect_equal(md$cf_pathway, "HALLMARK_APOPTOSIS")
  expect_length(md$cf_predictive_genes, 10L)
  # predictive genes are exactly the is_predictive rows.
  expect_setequal(md$cf_predictive_genes, rownames(se)[rd$is_predictive])
})

test_that("cf_example_data hits the requested SNR and is reproducible", {
  skip_if_not_installed("SummarizedExperiment")
  suppressPackageStartupMessages(library(SummarizedExperiment))

  se <- cf_example_data(n = 200L, p = 300L, snr = 2.0, seed = 7L)
  md <- S4Vectors::metadata(se)
  # cf_snr is the realised Var(tau)/Var(epsilon); the noise sd targets `snr`
  # against the realised Var(tau), so the realised ratio lands near the target
  # up to finite-sample variation in Var(epsilon).
  expect_equal(md$cf_snr, 2.0, tolerance = 0.2)
  # true ATE is the mean individual effect.
  expect_equal(md$cf_true_ate,
               mean(SummarizedExperiment::colData(se)$cf_tau_true))
  # A non-trivial minority benefit negatively (rule/benefit-group analyses work).
  expect_gt(mean(SummarizedExperiment::colData(se)$cf_tau_true < 0), 0.02)

  # Reproducible for a fixed seed; different for another.
  se2 <- cf_example_data(n = 200L, p = 300L, snr = 2.0, seed = 7L)
  expect_identical(SummarizedExperiment::assay(se, "counts"),
                   SummarizedExperiment::assay(se2, "counts"))
  se3 <- cf_example_data(n = 200L, p = 300L, snr = 2.0, seed = 8L)
  expect_false(identical(SummarizedExperiment::assay(se, "counts"),
                         SummarizedExperiment::assay(se3, "counts")))

  # The caller's RNG stream is restored (no side effect on global seed).
  set.seed(123L); a <- runif(1)
  set.seed(123L); invisible(cf_example_data(n = 40L, p = 60L)); b <- runif(1)
  expect_equal(a, b)
})

test_that("cf_example_data validates arguments", {
  expect_error(cf_example_data(n = 2L), class = "cfomics_arg_error")
  expect_error(cf_example_data(p = 5L), class = "cfomics_arg_error")
  expect_error(cf_example_data(p = 25L, n_prog = 20L, n_pred = 10L),
               class = "cfomics_arg_error")
  expect_error(cf_example_data(snr = -1), class = "cfomics_arg_error")
})

test_that("cf_example_data feeds cf_fit (R-native) end to end", {
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 80L, p = 150L, seed = 3L)
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = se, Y = "outcome", T = "treatment",
           X_baseline = c("age", "sex"), method = "grf",
           partition = "none", random_state = 1L)))
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
  expect_length(fit$res$tau_hat, 80L)
})
