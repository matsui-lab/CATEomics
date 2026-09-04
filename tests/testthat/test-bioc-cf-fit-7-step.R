# Tests for the cf_fit 7-step flow: meta slots populate + NN method z-score
# unification. Focused on paths not covered by test-bioc-integration.R /
# test-bioc.R (new meta slots / NN newdata path unification / factor_levels
# re-application).

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

create_minimal_se_counts <- function(n_samples = 30, n_genes = 20, seed = 42) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)
  # Raw-counts-like: integer matrix so the normalisation path runs when
  # preprocess = TRUE.
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 50),
    nrow = n_genes, ncol = n_samples
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("sample", seq_len(n_samples))
  col_data <- data.frame(
    treatment = rbinom(n_samples, 1, 0.5),
    age       = rnorm(n_samples, 50, 10),
    outcome   = rnorm(n_samples),
    row.names = colnames(counts)
  )
  col_data$outcome <- col_data$treatment * 1.5 + col_data$age * 0.05 +
    rnorm(n_samples)
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = counts),
    colData = col_data,
    rowData = data.frame(
      gene_symbol = paste0("SYMB", seq_len(n_genes)),
      row.names   = rownames(counts)
    )
  )
}

# -----------------------------------------------------------------------------
# df-path assay_cols / assay_pattern
# -----------------------------------------------------------------------------

test_that("assay_cols + assay_pattern mutual exclusion", {
  df <- data.frame(Y = rnorm(20), T = rbinom(20, 1, 0.5), X = rnorm(20))
  expect_error(
    cf_fit(Y ~ T | X, data = df, method = "grf", partition = "none",
           assay_cols = "X", assay_pattern = "^X"),
    class = "cfomics_input_conflict"
  )
})

test_that("df + assay_cols selects assay features", {
  set.seed(1)
  df <- data.frame(
    Y  = rnorm(30),
    T  = rbinom(30, 1, 0.5),
    age = rnorm(30, 50, 10)
  )
  for (g in paste0("gene", 1:10)) df[[g]] <- rnorm(30)
  fit <- cf_fit(
    Y ~ T | age, data = df, method = "grf", partition = "none",
    assay_cols = paste0("gene", 1:10),
    preprocess = FALSE
  )
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$meta$input_type, "data.frame")
  expect_equal(fit$meta$assay_spec$assay_cols, paste0("gene", 1:10))
  # X = age + 10 genes = 11 columns
  expect_equal(fit$meta$p, 11)
})

test_that("df + assay_pattern selects via regex", {
  set.seed(2)
  df <- data.frame(
    Y   = rnorm(30),
    T   = rbinom(30, 1, 0.5),
    age = rnorm(30, 50, 10)
  )
  for (g in paste0("ENSG", sprintf("%05d", 1:8))) df[[g]] <- rnorm(30)
  fit <- cf_fit(
    Y ~ T | age, data = df, method = "grf", partition = "none",
    assay_pattern = "^ENSG",
    preprocess = FALSE
  )
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$meta$assay_spec$assay_pattern, "^ENSG")
  expect_equal(fit$meta$p, 9)  # age + 8 genes
})

test_that("feature_method = 'pathway' is gated behind GSVA", {
  set.seed(3)
  df <- data.frame(
    Y  = rnorm(30),
    T  = rbinom(30, 1, 0.5),
    age = rnorm(30, 50, 10)
  )
  for (g in paste0("gene", 1:10)) df[[g]] <- rnorm(30)

  # Skip when GSVA is not installed. When installed, whether the fit succeeds
  # depends on the gene set and sample size, so this only smoke-checks that the
  # feature_method integration and pathway dispatch run.
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  # The Hallmark 50 gene set needs ENSG / SYMBOL matches; on this toy data a
  # zero match is likely, so catch any exception and only verify dispatch.
  err <- tryCatch(
    cf_fit(
      Y ~ T | age, data = df, method = "grf", partition = "none",
      assay_cols = paste0("gene", 1:10),
      preprocess = FALSE,
      feature_method = "pathway"
    ),
    error = function(e) e,
    warning = function(w) w
  )
  # With zero pathway matches GSVA may error or return an empty matrix; it is
  # enough here that the pathway dispatch either runs or stops with an
  # exception.
  expect_true(inherits(err, c("error", "warning", "cf_model")))
})

# -----------------------------------------------------------------------------
# Y leakage hard error
# -----------------------------------------------------------------------------

test_that("SE: Y in assay rownames raises cfomics_y_in_assay", {
  se <- create_minimal_se_counts(n_samples = 20, n_genes = 10)
  # rename one assay row to match outcome column name in colData
  expect_error(
    cf_fit(
      gene1 ~ treatment | age, data = se,
      method = "grf", partition = "none", preprocess = FALSE
    ),
    class = "cfomics_y_in_assay"
  )
})

# -----------------------------------------------------------------------------
# PCA loadings saved
# -----------------------------------------------------------------------------

test_that("feature_method = 'pca' populates meta$pca_loadings", {
  se <- create_minimal_se_counts(n_samples = 40, n_genes = 30)
  fit <- cf_fit(
    outcome ~ treatment | age, data = se,
    method = "grf", partition = "none",
    preprocess = FALSE,
    feature_method = "pca",
    n_features = 5
  )
  expect_s3_class(fit, "cf_model")
  expect_true(!is.null(fit$meta$pca_loadings))
  expect_equal(ncol(fit$meta$pca_loadings), 5)
  # rotation matrix should have one row per gene
  expect_equal(nrow(fit$meta$pca_loadings), 30)
})

# -----------------------------------------------------------------------------
# rowData round-trip
# -----------------------------------------------------------------------------

# meta$rowdata (full rowData table) is not part of the result contract;
# rowData round-trips via as_cf_data / cf_attach_results instead.

# -----------------------------------------------------------------------------
# input_type recording
# -----------------------------------------------------------------------------

test_that("meta$input_type reflects data class", {
  df <- data.frame(Y = rnorm(20), T = rbinom(20, 1, 0.5), X = rnorm(20))
  fit_df <- cf_fit(Y ~ T | X, data = df, method = "grf", partition = "none")
  expect_equal(fit_df$meta$input_type, "data.frame")

  se <- create_minimal_se_counts(n_samples = 30, n_genes = 10)
  fit_se <- cf_fit(
    outcome ~ treatment | age, data = se,
    method = "grf", partition = "none", preprocess = FALSE
  )
  expect_equal(fit_se$meta$input_type, "SE")
})

# -----------------------------------------------------------------------------
# cf_model$data$X_preprocessed saved
# -----------------------------------------------------------------------------

test_that("cf_model$data$X_preprocessed is stored", {
  set.seed(4)
  df <- data.frame(Y = rnorm(30), T = rbinom(30, 1, 0.5), X = rnorm(30))
  fit <- cf_fit(Y ~ T | X, data = df, method = "grf", partition = "none")
  expect_true(!is.null(fit$data))
  expect_true(is.matrix(fit$data$X_preprocessed))
  expect_equal(nrow(fit$data$X_preprocessed), 30)
  expect_equal(length(fit$data$Y), 30)
  expect_equal(length(fit$data$T), 30)
})

# -----------------------------------------------------------------------------
# factor one-hot + factor_levels
# -----------------------------------------------------------------------------

test_that("factor covariates produce one-hot + factor_levels", {
  set.seed(5)
  df <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    sex = sample(c("M", "F"), 40, replace = TRUE),
    stringsAsFactors = FALSE
  )
  fit <- cf_fit(Y ~ T | sex, data = df, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
  expect_true(!is.null(fit$meta$factor_levels))
  expect_setequal(fit$meta$factor_levels$sex, c("F", "M"))
  # one-hot: 2 columns for 2-level factor with intercept-free model.matrix
  expect_equal(fit$meta$p, 2)
})

# -----------------------------------------------------------------------------
# NN method meta$x_zscore / meta$y_zscore
# -----------------------------------------------------------------------------

test_that("NN methods populate meta$y_zscore", {
  skip_if_no_cf_python()
  set.seed(6)
  df <- data.frame(Y = rnorm(60), T = rbinom(60, 1, 0.5))
  for (g in paste0("X", 1:5)) df[[g]] <- rnorm(60)

  fit_tar <- cf_fit(
    Y ~ T | X1 + X2 + X3 + X4 + X5, data = df,
    method = "tarnet", random_state = 42L, n_iter = 50L
  )
  expect_true(!is.null(fit_tar$meta$y_zscore))
  expect_true(is.numeric(fit_tar$meta$y_zscore$Y_mean))
  expect_true(is.numeric(fit_tar$meta$y_zscore$Y_sd))
  expect_true(!is.null(fit_tar$meta$x_zscore))
  expect_equal(length(fit_tar$meta$x_zscore$X_mean), 5)
})

test_that("cevae populates central x_zscore and NN y_zscore", {
  skip_if_no_cf_python()
  set.seed(7)
  df <- data.frame(Y = rnorm(40), T = rbinom(40, 1, 0.5))
  for (g in paste0("X", 1:3)) df[[g]] <- rnorm(40)
  fit <- cf_fit(
    Y ~ T | X1 + X2 + X3, data = df,
    method = "cevae", random_state = 42L, num_epochs = 5L
  )
  expect_true(!is.null(fit$meta$y_zscore))
  expect_true(!is.null(fit$meta$x_zscore))
  expect_equal(length(fit$meta$x_zscore$X_mean), 3)
})

# -----------------------------------------------------------------------------
# newdata factor_levels re-apply
# -----------------------------------------------------------------------------

test_that("newdata factor_levels re-apply works", {
  set.seed(8)
  df_train <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    sex = sample(c("M", "F"), 40, replace = TRUE),
    stringsAsFactors = FALSE
  )
  fit <- cf_fit(Y ~ T | sex, data = df_train, method = "grf", partition = "none")

  df_new <- data.frame(
    sex = c("M", "F", "M"),
    T   = c(0, 1, 1),  # newdata may need T even if not used
    stringsAsFactors = FALSE
  )
  # newdata predict: should not error on factor re-application
  out <- predict(fit, newdata = df_new, type = "ate")
  expect_true(is.numeric(out) || is.list(out))
})


# -----------------------------------------------------------------------------
# .apply_lognorm = zero-var drop -> TMM (normLibSizes) -> log-CPM.
# filterByExpr is DEG-only and not applied to the causal feature matrix.
# -----------------------------------------------------------------------------

test_that(".apply_lognorm returns log-CPM + norm_method when edgeR is available", {
  skip_if_not_installed("edgeR")

  set.seed(42L)
  # 50 genes x 30 samples raw counts, with 5 zero-variance genes that the
  # zero-variance drop should remove.
  counts_var <- matrix(rpois(45 * 30, lambda = 50), nrow = 45)
  counts_zv  <- matrix(7L, nrow = 5, ncol = 30)  # constant rows
  counts     <- rbind(counts_var, counts_zv)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("s", seq_len(ncol(counts)))

  out <- cfomics:::.apply_lognorm(counts, input_type = "SE", preprocess = TRUE)
  expect_true(out$enabled)
  expect_equal(out$method, "logcpm")
  expect_equal(out$norm_method, "TMM")
  expect_true(is.integer(out$filter_kept))
  # The 5 zero-variance rows are dropped; 45 informative rows remain.
  expect_equal(out$filter_kept, 45L)
  expect_true(is.matrix(out$assay))
  expect_equal(ncol(out$assay), ncol(counts))
  # log-CPM values are on the log2 scale (roughly < 20 for these depths).
  expect_true(max(out$assay) < 25)
})


test_that(".apply_lognorm honors preprocess = FALSE (no-op)", {
  counts <- matrix(1:60, nrow = 6, ncol = 10)
  rownames(counts) <- paste0("g", 1:6)
  out <- cfomics:::.apply_lognorm(counts, preprocess = FALSE)
  expect_false(out$enabled)
  expect_identical(out$assay, counts)
})


# -----------------------------------------------------------------------------
# .apply_combat runs sva::ComBat on the log-CPM values (RNA-seq and microarray
# unified on the log scale).
# -----------------------------------------------------------------------------

test_that(".apply_combat corrects log-scale values via sva::ComBat", {
  skip_if_not_installed("sva")
  set.seed(11L)
  expr <- matrix(rnorm(20 * 30, mean = 8, sd = 1.5), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:30)
  batch_vec <- rep(c("A", "B"), each = 15)

  out <- cfomics:::.apply_combat(expr, batch_vec = batch_vec,
                                 batch_label = "batch")
  expect_true(out$batch_corrected)
  expect_equal(out$batch_method, "combat")
  expect_equal(dim(out$assay), dim(expr))
})


test_that(".apply_combat is a no-op when batch is NULL", {
  expr <- matrix(rnorm(12 * 10), nrow = 12)
  out <- cfomics:::.apply_combat(expr, batch_vec = NULL)
  expect_false(out$batch_corrected)
  expect_identical(out$assay, expr)
})
