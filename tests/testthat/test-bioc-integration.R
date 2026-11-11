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







test_that("cf_fit works with SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("grf")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  # preprocess = FALSE: the SE path applies normalisation by default when
  # preprocess = TRUE. create_test_se returns continuous rnorm values (not raw
  # counts), so preprocess = FALSE is set explicitly (clinical-only X = age).
  fit <- cf_fit(
    Y = "outcome", T = "treatment", X_baseline = "age",
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
    Y = "outcome", T = "treatment", X_baseline = "age",
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
    Y = "outcome", T = "treatment", X_baseline = "age",
    data = mae,
    method = "grf", partition = "none",
    preprocess = FALSE,
    feature_method = "none"
  )

  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "grf")
})

# selection path. Feature selection (variance / PCA) is exercised through
# cf_fit / cf_compare (`.feature_select_step`) in test-coverage-dgp-methods.R and the
# cf_fit SE feature-selection test above.




test_that("parse_cf_formula rejects unknown columns and outcome-in-covariates", {
  set.seed(1)
  df <- data.frame(Y = rnorm(15), T = rbinom(15, 1, 0.5), X1 = rnorm(15))
  # Unknown column (typo) in any of the three roles -> cateomics_unknown_column.
  expect_error(parse_cf_formula(Y ~ Tx | X1, df),
               class = "cateomics_unknown_column")
  expect_error(parse_cf_formula(Yx ~ T | X1, df),
               class = "cateomics_unknown_column")
  expect_error(parse_cf_formula(Y ~ T | Xnope, df),
               class = "cateomics_unknown_column")
  # Outcome listed among covariates -> cateomics_outcome_in_covariates
  # (symmetric with the pre-existing treatment-in-covariates check).
  expect_error(parse_cf_formula(Y ~ T | Y + X1, df),
               class = "cateomics_outcome_in_covariates")
  # Both reach the user through the named-args entry as well.
  expect_error(cf_fit(data = df, Y = "Y", T = "T", X_baseline = "Xnope",
                      omics = "none", method = "grf"),
               class = "cateomics_unknown_column")
})


# Writing variable importance into rowData -------------------------------

test_that("rowData VI columns follow the overwrite rule, as a pair", {
  skip_if_not_installed("SummarizedExperiment")
  vi_df <- data.frame(feature_name = c("g1", "g2", "g3"),
                      importance   = c(0.9, 0.5, 0.1),
                      rank         = 1:3,
                      stringsAsFactors = FALSE)
  rd <- S4Vectors::DataFrame(gene_type = c("a", "b", "c"),
                             row.names = c("g1", "g2", "g3"))

  out <- CATEomics:::.cf_rowdata_with_vi(rd, vi_df)
  expect_true(all(c("cf_vi_importance", "cf_vi_rank") %in% colnames(out)))
  expect_equal(out$cf_vi_importance, vi_df$importance)

  # A single stale column blocks both writes, so the pair cannot disagree.
  rd2 <- out
  rd2$cf_vi_importance <- NULL
  expect_warning(kept <- CATEomics:::.cf_rowdata_with_vi(rd2, vi_df),
                 class = "cateomics_attach_skipped")
  expect_false("cf_vi_importance" %in% colnames(kept))

  forced <- CATEomics:::.cf_rowdata_with_vi(out, vi_df, overwrite = TRUE)
  expect_equal(forced$cf_vi_rank, vi_df$rank)

  # A prefix matches the MAE covariate naming.
  rd3 <- S4Vectors::DataFrame(row.names = c("g1", "g2", "g3"))
  vi_pref <- vi_df
  vi_pref$feature_name <- paste0("rnaseq_", vi_df$feature_name)
  out3 <- CATEomics:::.cf_rowdata_with_vi(rd3, vi_pref, prefix = "rnaseq")
  expect_equal(out3$cf_vi_importance, vi_df$importance)
})

test_that("a covariate name cannot be produced twice", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  n <- 40L
  set.seed(3)
  cd <- S4Vectors::DataFrame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                             age = rnorm(n), row.names = paste0("s", seq_len(n)))
  # experiment "a" holding "b_x" and experiment "a_b" holding "x" would both
  # become "a_b_x".
  a  <- matrix(rnorm(2L * n), 2L,
               dimnames = list(c("b_x", "z"), paste0("s", seq_len(n))))
  ab <- matrix(rnorm(2L * n), 2L,
               dimnames = list(c("x", "w"), paste0("s", seq_len(n))))
  mae <- MultiAssayExperiment::MultiAssayExperiment(list(a = a, a_b = ab),
                                                    colData = cd)
  expect_error(
    cf_fit(data = mae, Y = "Y", T = "T", X_baseline = "age", omics = "all",
           method = "grf", partition = "none", random_state = 1L),
    class = "cateomics_duplicate_covariate_names")
})

test_that("an experiment without feature names is rejected at fit time", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  n <- 40L
  set.seed(4)
  cd <- S4Vectors::DataFrame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                             age = rnorm(n), row.names = paste0("s", seq_len(n)))
  named <- matrix(rnorm(3L * n), 3L,
                  dimnames = list(paste0("r", 1:3), paste0("s", seq_len(n))))
  unnamed <- matrix(rnorm(3L * n), 3L,
                    dimnames = list(NULL, paste0("s", seq_len(n))))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = named, blank = unnamed), colData = cd)
  expect_error(
    cf_fit(data = mae, Y = "Y", T = "T", X_baseline = "age", omics = "all",
           method = "grf", partition = "none", random_state = 1L),
    class = "cateomics_mae_layer_no_rownames")
})


# Results land on the sample they belong to ------------------------------

test_that("a MultiAssayExperiment fit keeps the colData sample order", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(5)
  n <- 30L; p <- 6L
  samp <- paste0("s", seq_len(n))
  counts <- matrix(rpois(p * n, 60), p, n,
                   dimnames = list(paste0("g", seq_len(p)), samp))
  cd <- S4Vectors::DataFrame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                             marker = seq_len(n) * 1.0, row.names = samp)
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(expr = counts), colData = cd)
  fit <- suppressMessages(cf_fit(
    data = mae, Y = "Y", T = "T", X_baseline = "marker", omics = "all",
    preprocess = FALSE, method = "grf", partition = "none", random_state = 1L))
  # The sampleMap join goes through table(), whose names come back sorted:
  # s1, s10, s11, ... Anything written back row by row would then land on the
  # wrong sample.
  expect_identical(rownames(fit$data$X_preprocessed), samp)
})

test_that("cf_attach_results matches the estimate to its own sample", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(5)
  n <- 30L; p <- 6L
  samp <- paste0("s", seq_len(n))
  counts <- matrix(rpois(p * n, 60), p, n,
                   dimnames = list(paste0("g", seq_len(p)), samp))
  marker <- seq_len(n) * 1.0
  T_vec <- rbinom(n, 1L, 0.5)
  cd <- S4Vectors::DataFrame(Y = marker + T_vec * marker + rnorm(n, sd = 0.01),
                             T = T_vec, marker = marker, row.names = samp)
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(expr = counts), colData = cd)
  fit <- suppressMessages(cf_fit(
    data = mae, Y = "Y", T = "T", X_baseline = "marker", omics = "all",
    preprocess = FALSE, method = "grf", partition = "none", random_state = 1L))

  out <- suppressMessages(cf_attach_results(mae, fit, what = "tau"))
  att <- as.data.frame(MultiAssayExperiment::colData(out))
  by_name <- stats::setNames(fit$res$tau_hat,
                             rownames(fit$data$X_preprocessed))
  expect_equal(unname(att$cf_tau), unname(by_name[rownames(att)]))

  # The companion fix makes the two orders agree, so the line above is the
  # identity permutation and would also pass if the write went back to
  # being positional. Permute the model's samples and require the values to
  # follow their names.
  ord <- rev(seq_len(n))
  shuffled <- fit
  shuffled$res$tau_hat <- fit$res$tau_hat[ord]
  rownames(shuffled$data$X_preprocessed) <-
    rownames(fit$data$X_preprocessed)[ord]
  out2 <- suppressMessages(cf_attach_results(mae, shuffled, what = "tau"))
  att2 <- as.data.frame(MultiAssayExperiment::colData(out2))
  expect_equal(att2$cf_tau, att$cf_tau)
  expect_false(isTRUE(all.equal(att2$cf_tau, shuffled$res$tau_hat)))

  # A permutation of the same length used to pass unnoticed; naming the
  # samples differently is now an error rather than a silent mismatch.
  wrong <- fit
  rownames(wrong$data$X_preprocessed) <- paste0("other", seq_len(n))
  expect_error(cf_attach_results(mae, wrong, what = "tau"),
               class = "cateomics_attach_sample_mismatch")
})

test_that("the attached tau is the out-of-fold one when the fit has it", {
  set.seed(11)
  n <- 50L; p <- 5L
  samp <- paste0("s", seq_len(n))
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  T_vec <- rbinom(n, 1L, 0.5)
  Y <- X[, 1] + T_vec * (1 + X[, 2]) + rnorm(n)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = matrix(rpois(p * n, 60), p, n,
                   dimnames = list(paste0("g", seq_len(p)), samp))),
    colData = S4Vectors::DataFrame(Y = Y, T = T_vec, row.names = samp))
  fit <- suppressMessages(cf_fit(
    data = data.frame(Y = Y, T = T_vec, X), Y = "Y", T = "T",
    X_baseline = colnames(X), method = "grf", partition = "oof",
    n_folds = 3L, random_state = 1L, preprocess = FALSE))

  # The two differ by more than rounding, so the assertion distinguishes them.
  expect_gt(max(abs(fit$res$tau_hat_oof - fit$res$tau_hat)), 0.05)
  v <- SummarizedExperiment::colData(
    suppressMessages(cf_attach_results(se, fit, what = "tau")))$cf_tau
  expect_equal(unname(v), unname(fit$res$tau_hat_oof))

  # A fit with no out-of-fold estimate falls back to the in-sample one.
  fit0 <- suppressMessages(cf_fit(
    data = data.frame(Y = Y, T = T_vec, X), Y = "Y", T = "T",
    X_baseline = colnames(X), method = "grf", partition = "none",
    random_state = 1L, preprocess = FALSE))
  expect_null(fit0$res$tau_hat_oof)
  v0 <- SummarizedExperiment::colData(
    suppressMessages(cf_attach_results(se, fit0, what = "tau")))$cf_tau
  expect_equal(unname(v0), unname(fit0$res$tau_hat))
})

test_that("a repeated sample name is rejected, not resolved silently", {
  skip_if_not_installed("MultiAssayExperiment")
  set.seed(6)
  n <- 12L
  ids <- c("p1", "p1", paste0("p", 2:11))
  cnt <- matrix(rpois(8L * n, 50), 8L, n,
                dimnames = list(paste0("g", 1:8), NULL))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cnt),
    colData = S4Vectors::DataFrame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                                   age = rnorm(n)))
  colnames(se) <- ids
  rownames(SummarizedExperiment::colData(se)) <- ids
  # Indexing a container by a repeated name silently returns the first match,
  # so a reorder-by-name would use one row and drop the other; and a
  # positional reorder would repeat the sample and inflate n. Neither is
  # something the caller can see, so refuse the input.
  expect_error(
    cf_fit(data = se, Y = "Y", T = "T", X_baseline = "age", omics = "counts",
           method = "grf", partition = "none", random_state = 1L,
           n_features = 4, preprocess = FALSE),
    class = "cateomics_duplicate_samples")

  # The MAE accepts a colData whose rownames repeat while the sampleMap names
  # the primary once, which is how the repeat reaches CATEomics at all.
  ex <- matrix(rpois(6L * 3L, 50), 6L, 3L,
               dimnames = list(paste0("g", 1:6), c("a1", "a2", "a3")))
  cd <- S4Vectors::DataFrame(Y = c(1.1, 2.2, 3.3, 4.4), T = c(0L, 1L, 0L, 1L),
                             age = 1:4, row.names = NULL)
  rownames(cd) <- c("q1", "q1", "q2", "q3")
  smap <- S4Vectors::DataFrame(assay = factor(rep("e", 3L)),
                               primary = c("q1", "q2", "q3"),
                               colname = c("a1", "a2", "a3"))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(e = ex), colData = cd, sampleMap = smap)
  expect_error(
    cf_fit(data = mae, Y = "Y", T = "T", X_baseline = "age", omics = "all",
           method = "grf", partition = "none", random_state = 1L,
           preprocess = FALSE, n_features = 3),
    class = "cateomics_duplicate_samples")
})

# --------------------------------------------------------------------------
# A colData column whose name is not a syntactic R identifier. `as.data.frame()`
# on a DataFrame rewrites such a name (`HLA-B` -> `HLA.B`) and ignores
# `check.names`, so the formula -- built from the names the caller gave --
# stopped finding its own column. Clinical columns are named this way in
# practice ("BMI (kg/m2)"), and so are gene symbols.
# --------------------------------------------------------------------------

.cf_nonsyntactic_fixture <- function(nm, n = 40L) {
  set.seed(1)
  d <- data.frame(Y = stats::rnorm(n), T = stats::rbinom(n, 1, 0.5),
                  age = stats::rnorm(n), z = stats::rnorm(n),
                  check.names = FALSE)
  colnames(d)[4] <- nm
  rownames(d) <- paste0("s", seq_len(n))
  assay <- matrix(stats::rpois(8L * n, 20), 8L, n,
                  dimnames = list(paste0("g", 1:8), paste0("s", seq_len(n))))
  list(df = d, assay = assay,
       cd = S4Vectors::DataFrame(d, check.names = FALSE))
}

test_that("a non-syntactic colData name survives the SE and MAE routes", {
  skip_if_not_installed("MultiAssayExperiment")
  for (nm in c("HLA-B", "bad name", "BMI (kg/m2)", "if")) {
    fx <- .cf_nonsyntactic_fixture(nm)
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = fx$assay), colData = fx$cd)
    fit <- suppressWarnings(suppressMessages(
      cf_fit(data = se, Y = "Y", T = "T", X_baseline = c("age", nm),
             omics = "none", method = "grf", n_features = 2)))
    expect_s3_class(fit, "cf_model")
    # The name reaches the design matrix unchanged, which is what every table
    # built from it will show.
    expect_true(nm %in% colnames(fit$data$X_preprocessed))

    mae <- MultiAssayExperiment::MultiAssayExperiment(
      experiments = MultiAssayExperiment::ExperimentList(list(rna = fx$assay)),
      colData = fx$cd)
    fit_mae <- suppressWarnings(suppressMessages(
      cf_fit(data = mae, Y = "Y", T = "T", X_baseline = c("age", nm),
             omics = "rna", method = "grf", n_features = 4)))
    expect_s3_class(fit_mae, "cf_model")
    expect_true(nm %in% colnames(fit_mae$data$X_preprocessed))
  }
})

test_that("predict(newdata = SE) finds a non-syntactic column too", {
  fx <- .cf_nonsyntactic_fixture("HLA-B")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fx$assay), colData = fx$cd)
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = se, Y = "Y", T = "T", X_baseline = c("age", "HLA-B"),
           omics = "none", method = "grf", n_features = 2)))
  tau <- predict(fit, newdata = se, type = "tau")
  expect_length(tau, ncol(se))
  expect_true(all(is.finite(tau)))
})

test_that("the data.frame route was never affected", {
  fx <- .cf_nonsyntactic_fixture("HLA-B")
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = fx$df, Y = "Y", T = "T", X_baseline = c("age", "HLA-B"),
           method = "grf", n_features = 2)))
  expect_true("HLA-B" %in% colnames(fit$data$X_preprocessed))
})
