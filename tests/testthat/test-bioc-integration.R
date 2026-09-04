# Tests for Bioconductor integration (SummarizedExperiment, MultiAssayExperiment)

# Helper function to create a mock SummarizedExperiment
create_test_se <- function(n_samples = 50, n_genes = 100) {
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)

  # Gene expression matrix (genes x samples)
  counts <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 3),
    nrow = n_genes, ncol = n_samples
  )
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)

  # Sample metadata
  col_data <- data.frame(
    treatment = rbinom(n_samples, 1, 0.5),
    age = rnorm(n_samples, 50, 10),
    outcome = rnorm(n_samples),
    row.names = colnames(counts)
  )
  # Make outcome correlated with treatment and age
  col_data$outcome <- col_data$treatment * 1.5 + col_data$age * 0.05 + rnorm(n_samples)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
}

# Helper function to create a mock MultiAssayExperiment
create_test_mae <- function(n_samples = 30, n_genes = 50, n_proteins = 40) {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)

  # Shared sample names
  sample_names <- paste0("sample", 1:n_samples)

  # RNAseq experiment
  rnaseq <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 3),
    nrow = n_genes, ncol = n_samples
  )
  rownames(rnaseq) <- paste0("gene", 1:n_genes)
  colnames(rnaseq) <- sample_names

  # Proteomics experiment
  proteomics <- matrix(
    rnorm(n_proteins * n_samples, mean = 5, sd = 2),
    nrow = n_proteins, ncol = n_samples
  )
  rownames(proteomics) <- paste0("protein", 1:n_proteins)
  colnames(proteomics) <- sample_names

  # Sample metadata
  col_data <- data.frame(
    treatment = rbinom(n_samples, 1, 0.5),
    age = rnorm(n_samples, 50, 10),
    outcome = rnorm(n_samples),
    row.names = sample_names
  )
  col_data$outcome <- col_data$treatment * 1.5 + col_data$age * 0.05 + rnorm(n_samples)

  # Create experiment list
  exp_list <- list(
    rnaseq = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = rnaseq)
    ),
    proteomics = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = proteomics)
    )
  )

  # Sample map
  sample_map <- data.frame(
    assay = rep(c("rnaseq", "proteomics"), each = n_samples),
    primary = rep(sample_names, 2),
    colname = rep(sample_names, 2)
  )

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = exp_list,
    colData = col_data,
    sampleMap = sample_map
  )
}

# --- SummarizedExperiment tests ---

test_that("is_se correctly identifies SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 10, n_genes = 20)
  df <- data.frame(a = 1:10)


  expect_true(is_se(se))
  expect_false(is_se(df))
  expect_false(is_se(NULL))
  expect_false(is_se(list()))
})

test_that("as_cf_data.SummarizedExperiment extracts colData variables", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(se, formula = outcome ~ treatment | age)

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  # Pure converter: colData covariate + EVERY assay feature (no selection).
  # age + gene1..gene100 = 101 covariates.
  expect_true("age" %in% result$covariate_names)
  expect_true(all(paste0("gene", 1:100) %in% result$covariate_names))
  expect_equal(ncol(result$X), 101L)
  expect_equal(length(result$Y), 50)
  expect_equal(length(result$T), 50)
  expect_true(all(result$T %in% c(0, 1)))
})

test_that("as_cf_data.SummarizedExperiment extracts gene as outcome", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(se, formula = gene1 ~ treatment | age)

  expect_equal(result$outcome_name, "gene1")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(length(result$Y), 50)
})

test_that("as_cf_data.SummarizedExperiment fails for non-existent variable", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  expect_error(
    as_cf_data(se, formula = nonexistent ~ treatment | age),
    "not found"
  )
})

test_that("as_cf_data.SummarizedExperiment requires formula or outcome/treatment", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  expect_error(
    as_cf_data(se),
    "formula|outcome|treatment"
  )
})

test_that("as_cf_data.SummarizedExperiment with manual specification", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age"
  )

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  # Pure converter: age + every assay feature.
  expect_true("age" %in% result$covariate_names)
  expect_equal(ncol(result$X), 101L)
})

test_that("as_cf_data.SummarizedExperiment with gene as outcome (manual)", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "gene1",
    treatment = "treatment",
    covariates = "age"
  )

  expect_equal(result$outcome_name, "gene1")
})

# as_cf_data is a pure converter (no feature_select / PCA / n_features
# arguments). Feature selection lives in cf_fit / cf_compare and is covered by
# test-coverage-dgp-methods.R and the cf_fit SE feature-selection test below.

test_that("cf_fit works with SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("grf")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  # preprocess = FALSE: the SE path applies normalisation by default when
  # preprocess = TRUE. create_test_se returns continuous rnorm values (not raw
  # counts), so preprocess = FALSE is set explicitly (clinical-only X = age).
  fit <- cf_fit(
    outcome ~ treatment | age,
    data = se,
    method = "grf", partition = "none",
    preprocess = FALSE
  )

  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "grf")
  expect_equal(fit$meta$n, 50)
})

test_that("cf_fit works with SummarizedExperiment and feature selection", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("grf")

  # n=150 to satisfy GRF honest splitting (B9 default); smaller n triggers
  # grf's "honesty fraction too close to 1 or 0" error from tune_forest.
  se <- create_test_se(n_samples = 150, n_genes = 100)

  # preprocess = FALSE: the test data are not raw counts, so normalisation is
  # skipped. feature_method = "variance" selects the top 10 features, giving
  # X = age + 10 features = 11 columns.
  fit <- cf_fit(
    outcome ~ treatment | age,
    data = se,
    method = "grf",
    preprocess = FALSE,
    feature_method = "variance",
    n_features = 10
  )

  expect_s3_class(fit, "cf_model")
  # Should include age + 10 variance-selected features
  expect_equal(fit$meta$p, 11)
})

# --- MultiAssayExperiment tests ---

test_that("is_mae correctly identifies MultiAssayExperiment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 10, n_genes = 20, n_proteins = 15)
  df <- data.frame(a = 1:10)

  expect_true(is_mae(mae))
  expect_false(is_mae(df))
  expect_false(is_mae(NULL))
})

test_that("as_cf_data.MultiAssayExperiment extracts colData variables", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  # When all variables are in colData, parse like data.frame
  result <- as_cf_data(mae, formula = outcome ~ treatment | age)

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(result$covariate_names, "age")
})

test_that("as_cf_data.MultiAssayExperiment with multi-omics features", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  # Pure converter: every feature of every assay is included (prefixed), with
  # no per-assay selection.
  result <- as_cf_data(
    mae,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age"
  )

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_true("age" %in% result$covariate_names)
  # Features from both rnaseq and proteomics (all of them).
  expect_true(any(grepl("^rnaseq_", result$covariate_names)))
  expect_true(any(grepl("^proteomics_", result$covariate_names)))
})

test_that("as_cf_data.MultiAssayExperiment requires outcome and treatment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  expect_error(
    as_cf_data(mae, outcome = "outcome"),  # missing treatment
    "treatment"
  )
})

test_that("cf_fit works with MultiAssayExperiment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("grf")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  # preprocess = FALSE + feature_method = "none": the MAE path extracts the
  # primary experiment and applies normalisation by default. This smoke test
  # only checks the cf_model structure, with clinical-only covariates
  # (normalisation skipped, feature_method = none, all primary genes).
  fit <- cf_fit(
    outcome ~ treatment | age,
    data = mae,
    method = "grf", partition = "none",
    preprocess = FALSE,
    feature_method = "none"
  )

  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "grf")
})

# `select_features_internal` is not part of the as_cf_data feature
# selection path. Feature selection (variance / PCA) is exercised through
# cf_fit / cf_compare (`.feature_select_step`) in test-coverage-dgp-methods.R and the
# cf_fit SE feature-selection test above.


test_that("as_cf_data.data.frame accepts named args equivalent to the formula", {
  set.seed(1)
  df <- data.frame(Y = rnorm(20), T = rbinom(20, 1, 0.5),
                   X1 = rnorm(20), X2 = rnorm(20))
  by_formula <- as_cf_data(df, formula = Y ~ T | X1 + X2)
  by_named   <- as_cf_data(df, outcome = "Y", treatment = "T",
                           covariates = c("X1", "X2"))
  expect_identical(by_named$outcome_name, by_formula$outcome_name)
  expect_identical(by_named$treatment_name, by_formula$treatment_name)
  expect_identical(by_named$covariate_names, by_formula$covariate_names)
  expect_equal(by_named$X, by_formula$X)
  # covariates default to all non-outcome/treatment columns.
  by_default <- as_cf_data(df, outcome = "Y", treatment = "T")
  expect_setequal(by_default$covariate_names, c("X1", "X2"))
  # outcome / treatment are required when no formula is given.
  expect_error(as_cf_data(df, covariates = "X1"),
               class = "cfomics_invalid_input")
})


test_that("parse_cf_formula rejects unknown columns and outcome-in-covariates", {
  set.seed(1)
  df <- data.frame(Y = rnorm(15), T = rbinom(15, 1, 0.5), X1 = rnorm(15))
  # Unknown column (typo) in any of the three roles -> cfomics_unknown_column,
  # on both entry paths.
  expect_error(as_cf_data(df, outcome = "Yx", treatment = "T", covariates = "X1"),
               class = "cfomics_unknown_column")
  expect_error(as_cf_data(df, outcome = "Y", treatment = "Tx", covariates = "X1"),
               class = "cfomics_unknown_column")
  expect_error(as_cf_data(df, formula = Y ~ T | Xnope),
               class = "cfomics_unknown_column")
  # Outcome listed among covariates -> cfomics_outcome_in_covariates
  # (symmetric with the pre-existing treatment-in-covariates check).
  expect_error(as_cf_data(df, outcome = "Y", treatment = "T",
                          covariates = c("Y", "X1")),
               class = "cfomics_outcome_in_covariates")
  expect_error(as_cf_data(df, formula = Y ~ T | Y + X1),
               class = "cfomics_outcome_in_covariates")
})
