## Unified newdata predict API smoke tests.
##
## Covers:
##   * type = "all" flat list contract (CATE-capable methods)
##   * type = "tau_ci" / "ate_ci" CI surface for grf
##   * CI on unsupported methods raises cateomics_ite_ci_unsupported

make_demo <- function(n = 80, p = 4, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.3 * X[, 1] + 0.4 * T + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(
    df  = data.frame(Y = Y, T = T, X),
    X_new = matrix(rnorm(20 * p), 20, p,
                   dimnames = list(NULL, paste0("X", seq_len(p))))
  )
}

test_that("predict.cf_model(type = 'all') returns flat list for grf", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3", "X4"), data = d$df, method = "grf",
          params = list(num.trees = 200L, tune.parameters = "none"))
  out <- predict(fit, newdata = d$X_new, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat", "mu1_hat",
                      "ps", "ps_type"))
  # The propensity forest scores rows it never saw, and the fit-time clip
  # bounds it exactly as in-sample.
  expect_length(out$ps, 20L)
  expect_identical(out$ps_type, "out-of-sample")
  expect_true(all(out$ps >= 0.025 & out$ps <= 0.975))
  expect_identical(predict(fit, newdata = d$X_new, type = "ps"), out$ps)
  expect_length(out$tau_hat, 20L)
  expect_length(out$mu0_hat, 20L)
  expect_length(out$mu1_hat, 20L)
  expect_true(is.numeric(out$ate))
  # GRF reconstructs y0/y1 from forest.Y/W; plug-in identity holds.
  expect_equal(out$mu1_hat - out$mu0_hat, out$tau_hat, tolerance = 1e-8)
  # Per-point CI on newdata via IJ variance (Wald) is valid...
  expect_true(is.data.frame(out$tau_ci))
  expect_named(out$tau_ci, c("lower", "upper"))
  # ...but there is no honest ATE CI on newdata (grf has no outcomes there),
  # so ate_ci is NA rather than a sqrt(mean(var_i)) approximation.
  expect_true(all(is.na(out$ate_ci)))
})

test_that("predict.cf_model raises cateomics_ite_ci_unsupported for non-linear DR final", {
  skip_if_no_cf_python()
  d <- make_demo()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3", "X4"), data = d$df, method = "drlearner",
                model_final = "rf")
  expect_error(
    predict(fit, newdata = d$X_new, type = "tau_ci"),
    class = "cateomics_drlearner_ci_unsupported_nonlinear_final"
  )
})

test_that("predict.cf_model returns a data.frame ite_ci for grf newdata", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3", "X4"), data = d$df, method = "grf",
          params = list(num.trees = 200L, tune.parameters = "none"))
  ite_ci <- predict(fit, newdata = d$X_new, type = "tau_ci")
  expect_s3_class(ite_ci, "data.frame")
  expect_named(ite_ci, c("lower", "upper"))
  expect_true(all(ite_ci$lower <= ite_ci$upper))
})

test_that("predict.cf_model dispatches stored fields for newdata = NULL", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3", "X4"), data = d$df, method = "grf",
          params = list(num.trees = 200L, tune.parameters = "none"))
  out <- predict(fit, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat", "mu1_hat",
                      "ps", "ps_type"))
  # partition = "oof" is the cf_fit default, so the stored PS is the one
  # assembled from the held-out fold predictions.
  expect_identical(out$ps_type, "oof")
  expect_length(out$ps, length(out$tau_hat))
  expect_true(all(!is.na(out$tau_hat)))
  expect_true(all(!is.na(out$mu0_hat)))
})

test_that("newdata works when the whole omics layer was taken (omics = 'all')", {
  # Regression: the newdata path used to rebuild the omics block only from
  # features the user had listed by name, so a fit that took the whole assay
  # left it with nothing to rebuild and predict() failed on a new container.
  skip_if_not_installed("grf")
  skip_if_not_installed("SummarizedExperiment")
  set.seed(21)
  n <- 40L
  counts <- matrix(round(abs(stats::rnorm(12L * n)) * 100), nrow = 12L,
                   dimnames = list(paste0("g", 1:12), paste0("s", seq_len(n))))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = counts),
    colData = S4Vectors::DataFrame(Y = stats::rnorm(n),
                                   T = stats::rbinom(n, 1, 0.5),
                                   age = stats::rnorm(n))
  )
  for (fm in c("none", "variance", "pca")) {
    fit <- cf_fit(data = se, Y = "Y", T = "T", X_baseline = "age",
                  omics = "all", method = "grf", partition = "none",
                  random_state = 1L, feature_method = fm, n_features = 5L,
                  params = list(tune.parameters = "none"))
    tau_new <- predict(fit, newdata = se[, 1:8], type = "tau")
    expect_length(tau_new, 8L)
    expect_false(anyNA(tau_new))
  }
  # The rows that entered preprocessing are recorded for that rebuild.
  fit <- cf_fit(data = se, Y = "Y", T = "T", X_baseline = "age",
                omics = "all", method = "grf", partition = "none",
                random_state = 1L, params = list(tune.parameters = "none"))
  lay <- fit$meta$estimator$preprocess$layers
  expect_named(lay, "expr")
  expect_setequal(lay$expr$features, rownames(counts))
})


# A factor baseline covariate survives the newdata rebuild ----------------

test_that("newdata prediction works when a baseline covariate is a factor", {
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 40, p = 60, n_prog = 6, n_pred = 3, seed = 1)
  expect_s3_class(SummarizedExperiment::colData(se)$sex, "factor")

  for (fm in c("none", "variance", "pca")) {
    fit <- suppressMessages(cf_fit(
      data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
      omics = "all", feature_method = fm, n_features = 15, method = "grf",
      partition = "none", random_state = 1L))
    # The training covariates are the expanded dummies; newdata holds the raw
    # factor, so the rebuild has to expand it rather than look for sexF / sexM.
    expect_true(all(c("sexF", "sexM") %in% fit$meta$covariate_names))
    rebuilt <- suppressMessages(CATEomics:::.apply_newdata_preprocessing(
      fit, CATEomics:::.extract_X_from_bioc_newdata(fit, se)))
    expect_identical(colnames(as.matrix(rebuilt)),
                     colnames(fit$data$X_preprocessed))
    expect_equal(as.matrix(rebuilt), fit$data$X_preprocessed,
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_length(predict(fit, newdata = se)$tau_hat, 40L)
  }
})

test_that("a clinical column missing from newdata is still named", {
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 40, p = 60, n_prog = 6, n_pred = 3, seed = 1)
  fit <- suppressMessages(cf_fit(
    data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
    omics = "all", feature_method = "none", method = "grf",
    partition = "none", random_state = 1L))
  se2 <- se
  SummarizedExperiment::colData(se2)$sex <- NULL
  expect_error(predict(fit, newdata = se2),
               class = "cateomics_newdata_missing_columns")
})


# The missing-value guard covers what reaches X, and no more --------------

test_that("an omics feature with a missing value is dropped, whatever the layer", {
  skip_if_not_installed("SummarizedExperiment")
  # The outcome must not depend on whether the layer happens to be normalised.
  se0 <- cf_example_data(n = 40, p = 60, n_prog = 6, n_pred = 3, seed = 1)
  a0 <- SummarizedExperiment::assay(se0, 1)
  gene <- names(sort(apply(a0, 1, stats::var), decreasing = TRUE))[1]
  a0[gene, 3] <- NA
  SummarizedExperiment::assay(se0, 1) <- a0

  for (pre in c(TRUE, FALSE)) {
    expect_warning(
      fit <- suppressMessages(cf_fit(
        data = se0, Y = "outcome", T = "treatment", X_baseline = "age",
        omics = "all", preprocess = pre, feature_method = "none",
        method = "grf", partition = "none", random_state = 1L)),
      class = "cateomics_assay_na_features_dropped")
    expect_false(gene %in% fit$meta$covariate_names)
    expect_false(anyNA(fit$data$X_preprocessed))
  }

  # Selection still returns the number of features that was asked for.
  expect_warning(
    fit <- suppressMessages(cf_fit(
      data = se0, Y = "outcome", T = "treatment", X_baseline = "age",
      omics = "all", feature_method = "variance", n_features = 5,
      method = "grf", partition = "none", random_state = 1L)),
    class = "cateomics_assay_na_features_dropped")
  expect_length(fit$meta$estimator$preprocess$selected_features, 5L)

  # A missing value in a named clinical covariate is an error, not a drop.
  se1 <- cf_example_data(n = 40, p = 60, n_prog = 6, n_pred = 3, seed = 1)
  SummarizedExperiment::colData(se1)$age[3] <- NA
  expect_error(
    cf_fit(data = se1, Y = "outcome", T = "treatment", X_baseline = "age",
           omics = "all", method = "grf", partition = "none",
           random_state = 1L),
    class = "cateomics_na_in_variable")
})

test_that("newdata with a missing value in a used feature is an error", {
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 40, p = 60, n_prog = 6, n_pred = 3, seed = 1)
  fit <- suppressMessages(cf_fit(
    data = se, Y = "outcome", T = "treatment", X_baseline = "age",
    omics = "all", feature_method = "none", method = "grf",
    partition = "none", random_state = 1L))
  se2 <- se
  a <- SummarizedExperiment::assay(se2, 1)
  a[fit$meta$estimator$preprocess$selected_features[1], 2] <- NA
  SummarizedExperiment::assay(se2, 1) <- a
  # Dropping a feature here would not match the training columns, and the
  # sample cannot be scored without it.
  expect_error(predict(fit, newdata = se2), class = "cateomics_assay_has_na")
})


# Pathway newdata with the default gene sets ------------------------------

test_that("pathway newdata works when the default gene sets were used", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("edgeR")
  se <- cf_example_data(n = 40, p = 200, n_prog = 6, n_pred = 3, seed = 1)
  fit <- suppressMessages(cf_fit(
    data = se, Y = "outcome", T = "treatment", X_baseline = "age",
    omics = "all", feature_method = "pathway", method = "grf",
    partition = "none", random_state = 1L))
  prep <- fit$meta$estimator$preprocess
  # The default collection is resolved inside the fit, so it has to be the
  # resolved sets that are stored, not the NULL the caller passed.
  expect_equal(prep$gene_sets_source, "hallmark_msigdbr")
  expect_gt(length(prep$gene_sets), 0L)
  expect_length(suppressMessages(predict(fit, newdata = se))$tau_hat, 40L)
})


# The omics block of newdata goes through the training chain ---------------

test_that("newdata missing an omics feature the fit used is an error", {
  set.seed(41)
  n <- 60L; p <- 40L
  genes <- paste0("G", seq_len(p))
  cnt <- matrix(rpois(p * n, 60), n, p,
                dimnames = list(paste0("s", seq_len(n)), genes))
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1L, 0.5), age = rnorm(n),
                   cnt, check.names = FALSE)
  fit <- suppressMessages(suppressWarnings(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = "age", omics = genes,
    preprocess = TRUE, feature_method = "variance", n_features = 8,
    method = "grf", partition = "none", random_state = 1L)))
  sel <- fit$meta$estimator$preprocess$selected_features

  # Every feature, including the ones selection left behind, sets the library
  # size, so a newdata frame that carries only the selected ones is not on
  # the training scale even though it looks complete.
  short <- df[, c("age", union(sel, genes[1:25]))]
  expect_error(predict(fit, newdata = short),
               class = "cateomics_newdata_omics_missing")

  full <- df[, c("age", genes)]
  rebuilt <- suppressMessages(CATEomics:::.apply_newdata_preprocessing(fit, full))
  expect_equal(as.matrix(rebuilt), fit$data$X_preprocessed,
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("new samples are normalised jointly with the training block", {
  skip_if_not_installed("edgeR")
  set.seed(77)
  n <- 60L; p <- 40L
  genes <- paste0("G", seq_len(p))
  depth <- runif(n, 0.6, 1.8)
  cnt <- t(vapply(seq_len(n), function(i) rpois(p, 60 * depth[i]),
                  numeric(p)))
  dimnames(cnt) <- list(paste0("s", seq_len(n)), genes)
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1L, 0.5), age = rnorm(n),
                   cnt, check.names = FALSE)
  tr <- df[1:52, ]; te <- df[53:60, ]

  fit <- suppressMessages(suppressWarnings(cf_fit(
    data = tr, Y = "Y", T = "T", X_baseline = "age", omics = genes,
    preprocess = TRUE, feature_method = "variance", n_features = 8,
    method = "grf", partition = "none", random_state = 1L)))
  sel <- fit$meta$estimator$preprocess$selected_features

  got <- suppressMessages(CATEomics:::.apply_newdata_preprocessing(
    fit, te[, c("age", genes)]))
  got <- as.matrix(got)[, sel, drop = FALSE]

  # What those eight samples would have carried had they been in the fit:
  # TMM and the library size are taken over the samples in the matrix, so
  # this only matches if the new samples were normalised with the training
  # block rather than on their own.
  lc_all <- CATEomics:::.apply_lognorm(t(cnt), input_type = "SE",
                                     preprocess = TRUE)$assay
  expected <- t(lc_all[sel, 53:60, drop = FALSE])
  xs  <- fit$meta$x_zscore
  idx <- match(sel, names(xs$X_mean))
  expected <- sweep(sweep(expected, 2, xs$X_mean[idx], "-"),
                    2, xs$X_sd[idx], "/")
  expect_equal(unname(got), unname(expected), tolerance = 1e-10)

  # Normalising the eight on their own does not land there.
  lc_alone <- CATEomics:::.cf_newdata_logcpm(t(cnt[53:60, , drop = FALSE]))
  expect_gt(max(abs(lc_alone[sel, ] - lc_all[sel, 53:60])), 1e-3)
})


test_that("an ensemble scores newdata as the weighted member combination", {
  skip_if_not_installed("grf")
  d <- make_demo()
  # Two grf candidates differing only in forest size: no Python, and it
  # exercises the `configs` entry at the same time.
  cmp <- suppressWarnings(suppressMessages(cf_compare(
    data = d$df, Y = "Y", T = "T", X_baseline = paste0("X", 1:4),
    n_folds = 3L, random_state = 0L, verbose = FALSE,
    eval_library = c("mean", "ridge"),
    configs = list(
      grf_small = list(method = "grf",
                       params = list(num.trees = 200L, num.threads = 1L,
                                     tune.parameters = "none")),
      grf_large = list(method = "grf",
                       params = list(num.trees = 400L, num.threads = 1L,
                                     tune.parameters = "none"))))))
  ens <- suppressMessages(cf_ensemble(cmp, method = "mean"))

  out <- predict(ens, newdata = d$X_new, type = "all")
  expect_named(out, c("ate", "ate_ci", "tau_hat", "tau_ci", "mu0_hat",
                      "mu1_hat", "ps", "ps_type"))
  expect_length(out$tau_hat, 20L)
  expect_true(all(is.finite(out$tau_hat)))

  # The combination is the members' own newdata effects under the stored
  # weights, so it reproduces exactly rather than approximately.
  m_tau <- vapply(names(ens$model$members), function(k)
    as.numeric(predict(ens$model$members[[k]], newdata = d$X_new,
                       type = "tau")), numeric(20L))
  manual <- as.numeric(m_tau %*% ens$weights[colnames(m_tau)])
  expect_equal(out$tau_hat, manual, tolerance = 0)
  expect_equal(predict(ens, newdata = d$X_new, type = "tau"), manual,
               tolerance = 0)
  expect_equal(out$ate, mean(manual), tolerance = 0)

  # The potential outcomes come from the shared evaluation nuisance in-sample,
  # which has no value for rows it was not fitted on.
  expect_true(all(is.na(out$mu0_hat)))
  expect_true(all(is.na(out$mu1_hat)))
  expect_length(out$mu0_hat, 20L)
  expect_true(all(is.na(out$ps)))
  expect_true(is.na(out$ps_type))
  # In-sample they are still produced, so the asymmetry is the nuisance, not
  # the ensemble forgetting how to form them.
  expect_true(all(is.finite(predict(ens, type = "all")$mu0_hat)))
})


test_that("newdata columns are picked by name whether or not a layer exists", {
  skip_if_not_installed("grf")
  # Covariates only, so the preprocessing record holds no omics layer. The
  # selecting feature methods used to skip the column pass on the strength of
  # `feature_method` alone, which is only safe when a layer is there to rebuild
  # X from; with no layer a stray column reached the backend and the user saw
  # grf's own message about a non-numeric feature matrix.
  set.seed(1)
  n <- 60L
  X <- matrix(rnorm(n * 8L), n, 8L, dimnames = list(NULL, paste0("g", 1:8)))
  df <- data.frame(Y = X[, 1] + rnorm(n), T = rbinom(n, 1L, 0.5),
                   age = rnorm(n, 60, 10), X)
  cv <- c("age", paste0("g", 1:8))
  for (fm in c("variance", "none")) {
    fit <- suppressMessages(cf_fit(
      data = df, Y = "Y", T = "T", X_baseline = cv, method = "grf",
      feature_method = fm, partition = "oof", n_folds = 3L, random_state = 1L,
      params = list(num.trees = 150L, tune.parameters = "none")))
    nd <- df
    nd$note <- sample(letters, n, TRUE)          # a column the fit never saw
    expect_equal(predict(fit, newdata = nd, type = "tau"),
                 predict(fit, newdata = df, type = "tau"), info = fm)
    expect_error(predict(fit, newdata = df[, setdiff(names(df), "g5")]),
                 class = "cateomics_newdata_missing_columns", info = fm)
  }
})

test_that("a non-syntactic colData name survives a MAE newdata", {
  skip_if_not_installed("MultiAssayExperiment")
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cd  <- SummarizedExperiment::colData(se)
  mk  <- function() {
    s <- SummarizedExperiment::SummarizedExperiment(
      assays = list(expr = SummarizedExperiment::assay(se)[1:8, ]),
      colData = S4Vectors::DataFrame(
        outcome   = cd$outcome,
        treatment = cd$treatment,
        `BMI (kg/m2)` = as.numeric(cd$age),
        check.names = FALSE,
        row.names = colnames(se)))
    MultiAssayExperiment::MultiAssayExperiment(
      experiments = list(expr = s),
      colData = SummarizedExperiment::colData(s))
  }
  mae <- mk()
  fit <- cf_fit(data = mae, Y = "outcome", T = "treatment",
                X_baseline = "BMI (kg/m2)", method = "grf",
                random_state = 0L, propensity_score = 0.5,
                eval_propensity_score = 0.5)
  # The column is real; the default `as.data.frame()` would have renamed it to
  # BMI..kg.m2. and the fit's own column would look missing.
  pr <- suppressWarnings(predict(fit, newdata = mk()))
  expect_equal(length(pr$tau_hat), ncol(mae[["expr"]]))
})

test_that("a clinical-only fit scores a SummarizedExperiment newdata", {
  se <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cd <- SummarizedExperiment::colData(se)
  df <- data.frame(Y = cd$outcome, treatment = cd$treatment, age = cd$age)
  fit <- cf_fit(data = df, Y = "Y", T = "treatment", X_baseline = "age",
                method = "grf", random_state = 0L, propensity_score = 0.5,
                eval_propensity_score = 0.5)
  # colData carries the outcome, the treatment and a factor the fit never saw.
  # Re-featuring is what would normally drop them, and there is no layer to
  # re-feature, so the columns are selected here instead.
  on_se <- suppressWarnings(predict(fit, newdata = se))
  on_df <- suppressWarnings(predict(fit, newdata = df))
  expect_equal(on_se$tau_hat, on_df$tau_hat)
})

test_that("a data.frame fit with an omics layer refuses a SE newdata by name", {
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cd  <- SummarizedExperiment::colData(se)
  cnt <- t(SummarizedExperiment::assay(se)[1:6, ])
  colnames(cnt) <- paste0("g", seq_len(6))
  df  <- cbind(data.frame(Y = cd$outcome, treatment = cd$treatment,
                          age = cd$age), cnt)
  fit <- cf_fit(data = df, Y = "Y", T = "treatment", X_baseline = "age",
                omics = paste0("g", seq_len(6)), omics_type = "proteomics",
                method = "grf", random_state = 0L, propensity_score = 0.5,
                eval_propensity_score = 0.5)
  # The layer is called "omics"; the newdata has no assay of that name, so the
  # promise to continue cannot be kept and is refused instead of dying raw.
  expect_error(suppressWarnings(predict(fit, newdata = se)),
               class = "cateomics_newdata_layer_not_in_assays")
})

test_that("a per-sample known propensity says it cannot reach newdata", {
  se <- cf_example_data(n = 60, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cd <- SummarizedExperiment::colData(se)
  df <- data.frame(Y = cd$outcome, treatment = cd$treatment, age = cd$age,
                   sex = cd$sex)
  nd <- df[1:10, ]
  vec <- cf_fit(data = df, Y = "Y", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                random_state = 0L,
                propensity_score = rep(c(0.4, 0.6), length.out = 60))
  w <- tryCatch(predict(vec, newdata = nd), warning = function(w) w)
  expect_s3_class(w, "cateomics_newdata_propensity_unavailable")
  pr <- suppressWarnings(predict(vec, newdata = nd))
  expect_true(all(is.na(pr$mu0_hat)))
  expect_false(any(is.na(pr$tau_hat)))
  # A scalar known score can be replayed, so it says nothing and keeps the mu.
  sc <- cf_fit(data = df, Y = "Y", T = "treatment",
               X_baseline = c("age", "sex"), method = "grf",
               random_state = 0L, propensity_score = 0.5,
               eval_propensity_score = 0.5)
  pr2 <- suppressWarnings(predict(sc, newdata = nd))
  expect_false(any(is.na(pr2$mu0_hat)))
})
