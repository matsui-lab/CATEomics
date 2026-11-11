# =============================================================================
# cf_fit partition / eval_nuisance / oof_ite regression tests.
#
# Covers both partition = 'oof' (default) and 'none' (opt-out) paths, plus
# downstream wiring, determinism, and keep_kfold_models. grf is the primary
# method (R-native, CPU-light, reproducible).
#
# Test data: n = 200 (grf default honesty.fraction 0.5 keeps fold n ~ 40),
# p = 4 numeric covariates, moderate treatment effect.
#
# MCMC methods such as bcf are skipped here (CPU-heavy).
# =============================================================================

make_partition_test_data <- function(n = 200L, p = 4L, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  T <- rbinom(n, 1, plogis(0.3 * X[, 1]))
  Y <- 0.5 * X[, 1] + 0.3 * X[, 2] + 0.8 * T + 0.4 * T * X[, 1] + rnorm(n)
  df <- as.data.frame(X)
  df$T <- T
  df$Y <- Y
  df
}

# ---------------------------------------------------------------------------
# 1. partition='oof' (default): oof_ite + nuisance populated
# ---------------------------------------------------------------------------

test_that("partition='oof' grf: oof_ite + nuisance + meta populated", {
  df <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                method = "grf", partition = "oof", random_state = 1L)
  expect_equal(length(fit$res$tau_hat_oof), 200L)
  expect_false(anyNA(fit$res$tau_hat_oof))
  expect_equal(fit$meta$estimator$partition, "oof")
  expect_true(is.list(fit$res$nuisance))
  expect_true(all(c("mu0_oof", "mu1_oof", "pi_oof") %in% names(fit$res$nuisance)))
  expect_equal(length(fit$res$nuisance$mu0_oof), 200L)
})

# ---------------------------------------------------------------------------
# 2. partition='none' (opt-out): oof_ite NULL, nuisance populated
# ---------------------------------------------------------------------------

test_that("partition='none' grf: oof_ite NULL but nuisance populated (spec L685)", {
  df <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                method = "grf", partition = "none", random_state = 1L)
  expect_null(fit$res$tau_hat_oof)
  expect_equal(fit$meta$estimator$partition, "none")
  expect_true(is.list(fit$res$nuisance))
  expect_equal(length(fit$res$nuisance$mu0_oof), 200L)
})

# ---------------------------------------------------------------------------
# 3b. res$fold_id: the split that produced the out-of-fold quantities travels
#     with them, and both entry points record it the same way.
# ---------------------------------------------------------------------------

test_that("partition='oof' records res$fold_id alongside tau_hat_oof", {
  df  <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"),
                data = df, method = "grf", partition = "oof",
                n_folds = 4L, random_state = 1L)
  fid <- fit$res$fold_id
  expect_true(is.integer(fid))
  expect_length(fid, 200L)
  expect_length(fid, length(fit$res$tau_hat_oof))
  expect_false(anyNA(fid))
  # every fold is used and every row belongs to exactly one
  expect_identical(sort(unique(fid)), 1:4)
  expect_equal(sum(table(fid)), 200L)
  # the split is stratified on T, so each fold holds both arms
  expect_true(all(vapply(split(df$T, fid),
                         function(t) length(unique(t)) == 2L, logical(1))))
})

test_that("partition='none' leaves res$fold_id NULL", {
  df  <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"),
                data = df, method = "grf", partition = "none",
                random_state = 1L)
  expect_null(fit$res$fold_id)
})

test_that("cf_compare models and the ensemble carry the shared res$fold_id", {
  df  <- make_partition_test_data()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"),
                    data = df, methods = "grf", n_folds = 4L,
                    eval_nuisance = "ridge", random_state = 1L,
                    BPPARAM = BiocParallel::SerialParam())
  expect_false(is.null(cmp$metadata$fold_id))
  for (lab in names(cmp$models)) {
    expect_identical(cmp$models[[lab]]$res$fold_id, cmp$metadata$fold_id)
  }
  ens <- suppressWarnings(cf_ensemble(cmp, method = "r_stacking"))
  expect_identical(ens$res$fold_id, cmp$metadata$fold_id)
  # the field indexes the out-of-fold estimates it sits next to
  expect_length(ens$res$fold_id, length(ens$res$tau_hat_oof))
})

# ---------------------------------------------------------------------------
# 4. partition='oof' + downstream (cf_evaluate / cf_interpret / cf_rule)
# ---------------------------------------------------------------------------

test_that("partition='oof' grf -> cf_evaluate / cf_interpret / cf_rule chain", {
  df <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                method = "grf", partition = "oof", random_state = 1L)

  ev <- cf_evaluate(fit, metrics = c("rate", "blp", "gates"),
                    random_state = 1L)
  expect_s3_class(ev, "cateomics_evaluate")
  expect_true(is.numeric(ev$rate$toc))

  # cf_interpret is surrogate-SHAP (no `type`/`blp_x`); output is per-surrogate
  # under $vi$<surrogate>$features.
  itp <- suppressWarnings(cf_interpret(fit, surrogate = "ridge",
                                       random_state = 1L))
  expect_s3_class(itp, "cateomics_interpret")
  expect_true(!is.null(itp$vi$ridge$features))

  rl <- cf_rule(fit, type = "threshold", threshold = 0)
  expect_s3_class(rl, "cateomics_rule")
  expect_equal(length(rl$decisions), 200L)
})

# ---------------------------------------------------------------------------
# 5. Determinism with random_state
# ---------------------------------------------------------------------------

test_that("partition='oof' deterministic with same random_state", {
  df <- make_partition_test_data()
  fit1 <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                 method = "grf", partition = "oof", random_state = 7L)
  fit2 <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                 method = "grf", partition = "oof", random_state = 7L)
  expect_equal(fit1$res$tau_hat_oof, fit2$res$tau_hat_oof, tolerance = 1e-10)
  expect_equal(fit1$res$nuisance$mu0_oof, fit2$res$nuisance$mu0_oof,
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 6. keep_kfold_models = TRUE populates fold_cache
# ---------------------------------------------------------------------------

test_that("keep_kfold_models = TRUE populates fold_cache, FALSE leaves NULL", {
  df <- make_partition_test_data()
  fit_keep <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                     method = "grf", partition = "oof",
                     keep_kfold_models = TRUE, n_folds = 5L, random_state = 1L)
  expect_equal(length(fit_keep$fold_cache), 5L)
  expect_true(all(vapply(fit_keep$fold_cache, inherits, logical(1),
                          what = "cf_model")))

  fit_drop <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                     method = "grf", partition = "oof", random_state = 1L)
  expect_null(fit_drop$fold_cache)
})

# ---------------------------------------------------------------------------
# 7. n_folds + eval_nuisance pass-through
# ---------------------------------------------------------------------------

test_that("n_folds = 3 reduces fold_cache length", {
  df <- make_partition_test_data()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                method = "grf", partition = "oof",
                keep_kfold_models = TRUE, n_folds = 3L, random_state = 1L)
  expect_equal(length(fit$fold_cache), 3L)
})

test_that("eval_nuisance = 'lasso' produces nuisance distinct from 'rf'", {
  df <- make_partition_test_data()
  fit_rf <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                    method = "grf", partition = "oof",
                    eval_nuisance = "rf", random_state = 1L)
  fit_la <- cf_fit(Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3", "x4"), data = df,
                    method = "grf", partition = "oof",
                    eval_nuisance = "lasso", random_state = 1L)
  expect_false(isTRUE(all.equal(fit_rf$res$nuisance$mu0_oof,
                                  fit_la$res$nuisance$mu0_oof,
                                  tolerance = 1e-6)))
})

# ---------------------------------------------------------------------------
# 8. SE / MAE partition smoke: SE / MAE inputs also work through the shared
#    preprocessing and the flatten inside `.kfold_oof_predict()`.
# ---------------------------------------------------------------------------

make_se_partition <- function(n_samples = 200L, n_genes = 10L, seed = 42L) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50),
                    nrow = n_genes, ncol = n_samples)
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

test_that("SE input + partition='oof' grf: oof_ite + nuisance populated (Q-impl-1)", {
  se <- make_se_partition()
  fit <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = se,
                method = "grf", partition = "oof",
                preprocess = FALSE, random_state = 1L)
  expect_equal(length(fit$res$tau_hat_oof), 200L)
  expect_false(anyNA(fit$res$tau_hat_oof))
  expect_equal(fit$meta$estimator$partition, "oof")
  expect_true(is.list(fit$res$nuisance))
  expect_equal(length(fit$res$nuisance$mu0_oof), 200L)
  expect_equal(fit$meta$input_type, "SE")
})

test_that("MAE input + partition='oof' grf: oof_ite + nuisance populated (Q-impl-1)", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  se1 <- make_se_partition()
  se2 <- make_se_partition(seed = 7L)
  # Use single experiment for MAE primary to keep flatten simple.
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = se1),
    colData     = SummarizedExperiment::colData(se1)
  )
  fit <- cf_fit(Y = "outcome", T = "treatment", X_baseline = "age", data = mae,
                method = "grf", partition = "oof",
                preprocess = FALSE, random_state = 1L)
  expect_equal(length(fit$res$tau_hat_oof), 200L)
  expect_false(anyNA(fit$res$tau_hat_oof))
  expect_equal(fit$meta$estimator$partition, "oof")
  expect_true(is.list(fit$res$nuisance))
  expect_equal(fit$meta$input_type, "MAE")
})
