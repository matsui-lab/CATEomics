# =============================================================================
# MAE per-assay normalisation contract.
#
# Only the RNA-seq assay is internally log-normalised; every other assay passes
# through as-is (the user is expected to pre-normalise it). The RNA-seq assay is
# identified by experiment name (rnaseq / rna / rna_seq) or the named-arg
# `X_rnaseq` slot. When no RNA-seq assay is present, no assay is normalised.
# `feature_method = "pathway"` is not defined for the multi-assay MAE matrix.
# =============================================================================

.mae_assay <- function(name, kind, seed, n = 60L, g = 8L) {
  set.seed(seed)
  m <- if (kind == "counts") matrix(rpois(g * n, 50), g, n) else
    matrix(rnorm(g * n, 5, 1), g, n)
  rownames(m) <- paste0(name, "_g", seq_len(g))
  colnames(m) <- paste0("s", seq_len(n))
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
}

.mae_coldata <- function(n = 60L, seed = 9L) {
  set.seed(seed)
  cd <- S4Vectors::DataFrame(
    treatment = rbinom(n, 1, 0.5),
    age       = rnorm(n, 50, 5),
    row.names = paste0("s", seq_len(n))
  )
  cd$outcome <- cd$treatment * 1.2 + 0.05 * cd$age + rnorm(n)
  cd
}

test_that("MAE with an RNA-seq experiment flags only that assay for normalisation", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  cd <- .mae_coldata()
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = .mae_assay("rnaseq", "counts", 1),
         methyl = .mae_assay("methyl", "cont", 2)),
    colData = cd)

  ex <- .extract_mae_assay(mae, outcome ~ treatment | age)
  expect_true(ex$primary_is_rnaseq)
  expect_equal(ex$primary_exp, "rnaseq")

  fit <- cf_fit(outcome ~ treatment | age, data = mae, method = "grf",
                partition = "none", random_state = 1L)
  expect_s3_class(fit, "cf_model")
})

test_that("MAE without an RNA-seq experiment passes all assays through", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  cd <- .mae_coldata()
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(methyl  = .mae_assay("methyl", "cont", 2),
         protein = .mae_assay("protein", "cont", 3)),
    colData = cd)

  ex <- .extract_mae_assay(mae, outcome ~ treatment | age)
  expect_false(ex$primary_is_rnaseq)
  expect_equal(ex$primary_exp, "methyl")  # experiments[1] carrier

  fit <- cf_fit(outcome ~ treatment | age, data = mae, method = "grf",
                partition = "none", random_state = 1L)
  expect_s3_class(fit, "cf_model")
})

test_that("same count matrix is normalised as RNA-seq but not as a generic assay", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("edgeR")
  cd <- .mae_coldata()

  # Identical count assay, once designated RNA-seq (normalised) and once a
  # generic layer (passthrough). The stored feature matrices must differ.
  as_rna  <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = .mae_assay("rnaseq", "counts", 1)), colData = cd)
  as_gen  <- MultiAssayExperiment::MultiAssayExperiment(
    list(methyl = .mae_assay("methyl", "counts", 1)), colData = cd)

  f_rna <- cf_fit(outcome ~ treatment | age, data = as_rna, method = "grf",
                  partition = "none", random_state = 1L)
  f_gen <- cf_fit(outcome ~ treatment | age, data = as_gen, method = "grf",
                  partition = "none", random_state = 1L)

  x_rna <- f_rna$data$X[, grep("_g", colnames(f_rna$data$X)), drop = FALSE]
  x_gen <- f_gen$data$X[, grep("_g", colnames(f_gen$data$X)), drop = FALSE]
  expect_false(isTRUE(all.equal(unname(as.numeric(x_rna)),
                                unname(as.numeric(x_gen)))))
})

test_that("feature_method='pathway' is rejected for MAE input", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  cd <- .mae_coldata()
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = .mae_assay("rnaseq", "counts", 1)), colData = cd)
  expect_error(
    cf_fit(outcome ~ treatment | age, data = mae, method = "grf",
           feature_method = "pathway"),
    class = "cfomics_mae_pathway_unsupported")
})
