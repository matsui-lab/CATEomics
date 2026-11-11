# =============================================================================
# MAE per-assay normalisation contract.
#
# Only the RNA-seq assay is internally log-normalised; every other assay passes
# through as-is (the user is expected to pre-normalise it). The RNA-seq assay is
# identified by experiment name (rnaseq / rna / rna_seq / rna-seq) or by
# `omics_type`. When no RNA-seq assay is present, no assay is normalised.
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

  fit <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = mae, method = "grf",
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

  fit <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = mae, method = "grf",
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

  f_rna <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = as_rna, method = "grf",
                  partition = "none", random_state = 1L)
  f_gen <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = as_gen, method = "grf",
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
    cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = mae, method = "grf",
           feature_method = "pathway"),
    class = "cateomics_mae_pathway_unsupported")
})


# Layer concatenation ------------------------------------------------------

.mae_two_layer <- function(n = 60L) {
  MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = .mae_assay("rnaseq", "counts", 1, n = n),
         methyl = .mae_assay("methyl", "cont", 2, n = n)),
    colData = .mae_coldata(n = n))
}

test_that("every requested layer enters X, prefixed by its experiment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  mae <- .mae_two_layer()

  fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                X_baseline = "age", omics = "all", method = "grf",
                partition = "none", random_state = 1L)
  cn <- fit$meta$covariate_names
  expect_setequal(cn, c("age",
                        paste0("rnaseq_rnaseq_g", 1:8),
                        paste0("methyl_methyl_g", 1:8)))

  # Naming an experiment as a whole must reach the same columns.
  fit2 <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                 X_baseline = "age",
                 omics = list(rnaseq = "all", methyl = "all"),
                 method = "grf", partition = "none", random_state = 1L)
  expect_setequal(fit2$meta$covariate_names, cn)
})

test_that("named features from two layers are both kept, and prefixed", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  mae <- .mae_two_layer()

  fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                X_baseline = "age",
                omics = list(rnaseq = c("rnaseq_g1", "rnaseq_g2"),
                             methyl = "methyl_g3"),
                method = "grf", partition = "none", random_state = 1L)
  expect_setequal(fit$meta$covariate_names,
                  c("age", "rnaseq_rnaseq_g1", "rnaseq_rnaseq_g2",
                    "methyl_methyl_g3"))
})

test_that("feature selection is applied within each layer", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  mae <- .mae_two_layer()

  for (fm in c("variance", "pca")) {
    fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                  X_baseline = "age", omics = "all", method = "grf",
                  feature_method = fm, n_features = 3, partition = "none",
                  random_state = 1L)
    cn <- fit$meta$covariate_names
    # 3 features from each of the two layers, plus the clinical column.
    expect_length(cn, 7L)
    expect_equal(sum(startsWith(cn, "rnaseq_")), 3L)
    expect_equal(sum(startsWith(cn, "methyl_")), 3L)
  }
})

test_that("MAE newdata rebuilds the training design matrix layer by layer", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  mae <- .mae_two_layer()

  for (fm in c("none", "variance", "pca")) {
    fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                  X_baseline = "age", omics = "all", method = "grf",
                  feature_method = fm, n_features = 3, partition = "none",
                  random_state = 1L)
    nd <- .apply_newdata_preprocessing(fit,
            .extract_X_from_bioc_newdata(fit, mae))
    expect_identical(colnames(as.matrix(nd)),
                     colnames(fit$data$X_preprocessed))
    expect_equal(as.matrix(nd), fit$data$X_preprocessed,
                 ignore_attr = TRUE, tolerance = 1e-12)
  }
})

test_that("a missing training feature in newdata is named", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  mae <- .mae_two_layer()
  fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                X_baseline = "age", omics = "all", method = "grf",
                partition = "none", random_state = 1L)
  short <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = .mae_assay("rnaseq", "counts", 1, g = 4L),
         methyl = .mae_assay("methyl", "cont", 2)),
    colData = .mae_coldata())
  expect_error(.extract_X_from_bioc_newdata(fit, short),
               class = "cateomics_newdata_mae_layer_missing")
})
