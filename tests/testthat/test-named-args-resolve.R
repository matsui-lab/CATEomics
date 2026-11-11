## Tests for `omics` resolution per container.
##
## `.cf_resolve_omics()` decides which layer(s) and which features make up the
## omics block, and which of them count as RNA-seq for normalisation.
## Column-existence validation for the formula covariates remains the job of
## `parse_cf_formula` and is not duplicated here.

library(testthat)

.make_df <- function(n = 120L, p = 4L, seed = 42L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  T <- stats::rbinom(n, 1, 0.5)
  Y <- as.numeric(T) * 1.5 + rowSums(X) + 0.1 * stats::rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  colnames(df) <- c("Y", "T", paste0("X", seq_len(p)))
  df
}

.make_se <- function(n = 40L, ng = 6L, assays = c("counts")) {
  cd <- S4Vectors::DataFrame(
    Y   = stats::rnorm(n),
    T   = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n)
  )
  mats <- lapply(assays, function(a) {
    matrix(stats::rnorm(ng * n), nrow = ng, ncol = n,
           dimnames = list(paste0("G", seq_len(ng)), paste0("S", seq_len(n))))
  })
  names(mats) <- assays
  SummarizedExperiment::SummarizedExperiment(assays = mats, colData = cd)
}

test_that("data.frame defaults to no omics block", {
  df <- .make_df()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:4),
                data = df, method = "grf", random_state = 42L,
                partition = "none", params = list(tune.parameters = "none"))
  expect_s3_class(fit, "cf_model")
  expect_length(fit$meta$estimator$preprocess$layers, 0L)
  expect_equal(fit$meta$p, 4L)
})

test_that("data.frame omics columns become the assay block", {
  df <- .make_df()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"),
                omics = c("X3", "X4"), preprocess = FALSE,
                feature_method = "none",
                data = df, method = "grf", random_state = 42L,
                partition = "none", params = list(tune.parameters = "none"))
  expect_s3_class(fit, "cf_model")
  expect_setequal(fit$meta$named_args$X_baseline, c("X1", "X2"))
  expect_setequal(fit$meta$estimator$preprocess$layers$omics$features,
                  c("X3", "X4"))
  expect_equal(fit$meta$p, 4L)
})

test_that("a column that is not in the data is rejected", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", omics = c("nope"),
           data = df, method = "grf"),
    class = "cateomics_omics_columns_missing"
  )
})

test_that("SE: assays cannot be combined", {
  skip_if_not_installed("SummarizedExperiment")
  se <- .make_se(assays = c("counts", "logcounts"))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "age",
           omics = list(counts = "all", logcounts = "all"),
           data = se, method = "grf"),
    class = "cateomics_omics_multiple_assays"
  )
})

test_that("SE: several assays must be named explicitly", {
  skip_if_not_installed("SummarizedExperiment")
  se <- .make_se(assays = c("counts", "logcounts"))
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "age", omics = "all",
           data = se, method = "grf"),
    class = "cateomics_omics_assay_ambiguous"
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "age", omics = "counts",
                preprocess = FALSE, data = se, method = "grf",
                random_state = 42L, partition = "none",
                params = list(tune.parameters = "none"))
  # The layer is named after the assay it came from.
  expect_named(fit$meta$estimator$preprocess$layers, "counts")
})

test_that("SE: an unknown layer name is rejected", {
  skip_if_not_installed("SummarizedExperiment")
  se <- .make_se()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "age",
           omics = list(atac = "all"), data = se, method = "grf"),
    class = "cateomics_omics_unknown_layer"
  )
})

test_that("omics_type overrides the layer-name rule", {
  skip_if_not_installed("SummarizedExperiment")
  se <- .make_se()
  spec <- .cf_resolve_omics("all", omics_type = "methyl", data = se)
  expect_false(unname(spec$types[["counts"]]))
  spec2 <- .cf_resolve_omics("all", omics_type = NULL, data = se)
  expect_true(unname(spec2$types[["counts"]]))
})

test_that("MAE: layers map to experiments and only RNA-seq is normalised", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  n <- 40L
  cd_outer <- S4Vectors::DataFrame(
    Y   = stats::rnorm(n),
    T   = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n),
    row.names = paste0("S", seq_len(n))
  )
  make_exp <- function(prefix, ng) {
    m <- matrix(stats::rnorm(ng * n), nrow = ng, ncol = n,
                dimnames = list(paste0(prefix, seq_len(ng)),
                                paste0("S", seq_len(n))))
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  }
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rnaseq = make_exp("R", 5L), mirna = make_exp("M", 4L)),
    colData     = cd_outer
  )
  spec <- .cf_resolve_omics("all", data = mae)
  expect_setequal(names(spec$layers), c("rnaseq", "mirna"))
  expect_true(unname(spec$types[["rnaseq"]]))
  expect_false(unname(spec$types[["mirna"]]))
})

test_that("omics = 'none' leaves a MAE with nothing to use", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  n <- 20L
  cd <- S4Vectors::DataFrame(Y = stats::rnorm(n), T = stats::rbinom(n, 1, 0.5),
                             age = stats::rnorm(n),
                             row.names = paste0("S", seq_len(n)))
  m <- matrix(stats::rnorm(3 * n), nrow = 3,
              dimnames = list(paste0("R", 1:3), paste0("S", seq_len(n))))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rnaseq = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = m))),
    colData = cd
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "age", omics = "none",
           data = mae, method = "grf"),
    class = "cateomics_omics_none_mae"
  )
})

test_that("omics must be a character vector or a named list", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", omics = 1:3,
           data = df, method = "grf"),
    class = "cateomics_invalid_omics"
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", omics = list("X2"),
           data = df, method = "grf"),
    class = "cateomics_invalid_omics"
  )
})
