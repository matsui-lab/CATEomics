# Coverage test — every method's fit/predict contract, and the major
# analysis functions (compare / select / ensemble / evaluate / rule
# / interpret / aggregate_vi / benefit_groups / model_stability / diagnostics /
# attach_results / as_shapviz) plus the plot family,
# exercised across the three input containers (data.frame, SummarizedExperiment,
# MultiAssayExperiment) using the recommended named-argument API. Also covers
# propensity-score regimes (observational / known-RCT / trimmed) and binary
# outcomes. Smoke level: runs + correct class/shape/contract, not estimation
# accuracy (precision is the separate benchmark track).

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

.cov_make_df <- function(n = 60L, p = 8L, seed = 7L, binary = FALSE) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(0.4 * X[, 1]))
  if (binary) {
    Y <- rbinom(n, 1, plogis(0.3 * X[, 1] + 0.6 * T))
  } else {
    Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.4)
  }
  data.frame(Y = Y, T = T, X)
}

.cov_covnames <- function(p = 8L) paste0("X", seq_len(p))

.cov_coldata <- function(n, seed = 11L, binary = FALSE) {
  set.seed(seed)
  cd <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age       = rnorm(n, 50, 10),
    sex       = rbinom(n, 1, 0.5),          # clinical covariate
    bmi       = rnorm(n, 25, 4),            # clinical covariate
    outcome   = rnorm(n),
    row.names = paste0("s", seq_len(n))
  )
  lin <- 0.5 * cd$treatment + 0.03 * cd$age + 0.2 * cd$sex + 0.02 * cd$bmi
  cd$outcome <- if (binary) rbinom(n, 1, plogis(lin - mean(lin))) else
    lin + rnorm(n, sd = 0.4)
  cd
}

# The RNA-seq assay holds non-negative raw counts so the default
# feature_method = "variance" path (zero-var drop -> TMM -> log-CPM) runs
# cleanly. Negative "counts" would break edgeR::normLibSizes once the omics
# assay is included by default.
.cov_make_se <- function(n = 60L, n_genes = 30L, seed = 12L, binary = FALSE) {
  set.seed(seed)
  m <- matrix(rnbinom(n_genes * n, size = 10, mu = 200),
              nrow = n_genes, ncol = n)
  rownames(m) <- paste0("gene", seq_len(n_genes))
  colnames(m) <- paste0("s", seq_len(n))
  cd <- .cov_coldata(n, binary = binary)
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = m), colData = S4Vectors::DataFrame(cd))
}

# MAE with an RNA-seq experiment plus two other omics layers (methylation,
# protein). Only the RNA-seq assay is internally log-normalised
# (so it holds raw counts); the other layers are user-pre-normalised and pass
# through as-is (continuous, e.g. M-values / log-intensities).
.cov_make_mae <- function(n = 60L, seed = 13L) {
  set.seed(seed)
  sn <- paste0("s", seq_len(n))
  mk_counts <- function(prefix, k, mu) {
    a <- matrix(rnbinom(k * n, size = 10, mu = mu), nrow = k, ncol = n)
    rownames(a) <- paste0(prefix, seq_len(k)); colnames(a) <- sn
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = a))
  }
  mk_norm <- function(prefix, k, mean, sd) {
    a <- matrix(rnorm(k * n, mean = mean, sd = sd), nrow = k, ncol = n)
    rownames(a) <- paste0(prefix, seq_len(k)); colnames(a) <- sn
    SummarizedExperiment::SummarizedExperiment(assays = list(values = a))
  }
  exps <- list(
    rnaseq  = mk_counts("gene", 20L, 200),   # raw counts -> internal log-CPM
    methyl  = mk_norm("cpg",  15L, 0, 1),    # pre-normalised passthrough
    protein = mk_norm("prot", 10L, 3, 1)     # pre-normalised passthrough
  )
  cd <- .cov_coldata(n)
  sm <- do.call(rbind, lapply(names(exps), function(e)
    data.frame(assay = e, primary = sn, colname = sn)))
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(exps),
    colData = S4Vectors::DataFrame(cd), sampleMap = sm)
}

# --------------------------------------------------------------------------
# Block 1: per-method fit/predict contract (named-args, data.frame)
# --------------------------------------------------------------------------

test_that("R-native methods honour the four-pillar contract", {
  skip_if_not_installed("grf")
  df <- .cov_make_df()
  Xb <- .cov_covnames()

  for (m in c("grf", "bcf")) {
    if (m == "bcf" && !requireNamespace("bcf", quietly = TRUE)) next
    fit <- suppressWarnings(suppressMessages(
      cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
             method = m, partition = "none", random_state = 1L)))
    expect_s3_class(fit, "cf_model")
    expect_s3_class(fit, "cf_model")
    expect_true(all(c("tau_hat", "ate", "mu0_hat", "mu1_hat") %in% names(fit$res)))
    expect_equal(length(predict(fit, type = "tau")), nrow(df))
    expect_length(predict(fit, type = "ate"), 1L)
    expect_false(anyNA(predict(fit, type = "tau")))       # CATE method
    s <- predict(fit, type = "summary")
    expect_true(all(c("ate", "tau_mean", "tau_iqr") %in% names(s)))
  }

})

test_that("Python methods honour the four-pillar contract (basilisk)", {
  skip_if_no_cf_python()
  df <- .cov_make_df()
  Xb <- .cov_covnames()

  for (m in c("slearner", "tlearner", "xlearner", "drlearner",
              "cevae", "tarnet", "dragonnet")) {
    nn <- m %in% c("cevae", "tarnet", "dragonnet")
    extra <- if (m == "cevae") list(params = list(num_epochs = 10L)) else
             if (nn) list(params = list(n_iter = 50L)) else list()
    fit <- suppressWarnings(suppressMessages(do.call(cf_fit, c(
      list(data = df, Y = "Y", T = "T", X_baseline = Xb,
           method = m, partition = "none", random_state = 1L), extra))))
    expect_s3_class(fit, "cf_model")
    expect_true(all(c("tau_hat", "ate", "mu0_hat", "mu1_hat") %in% names(fit$res)))
    expect_equal(length(predict(fit, type = "tau")), nrow(df))
    expect_true(file.exists(fit$model$model$model_path))  # serialised model round-trips
  }
})

# --------------------------------------------------------------------------
# Block 2: major analysis functions across df / SE / MAE (named-args)
# --------------------------------------------------------------------------

# Runs the full analysis workflow with two R-native CATE methods (grf + bcf)
# so it exercises every major function deterministically without the Python
# backends. `data` + named-arg (Y/T/X_baseline) change per container.
.cov_run_workflow <- function(data, Y, T, X_baseline) {
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = data, Y = Y, T = T, X_baseline = X_baseline,
               methods = c("grf", "bcf"), n_folds = 3L,
               verbose = FALSE, random_state = 1L)))
  expect_s3_class(cmp, "cateomics_compare")

  sel <- cf_select(cmp, metric = "r_loss")
  expect_s3_class(sel, "cateomics_select")
  model <- sel$model

  expect_s3_class(suppressWarnings(cf_ensemble(cmp, method = "mean")),
                  "cateomics_ensemble")

  ev <- suppressWarnings(suppressMessages(cf_evaluate(model)))
  expect_s3_class(ev, "cateomics_evaluate")

  rule <- suppressWarnings(cf_rule(model, type = "topk", topk_frac = 0.3))
  expect_s3_class(rule, "cateomics_rule")

  rg <- suppressWarnings(cf_rule(model, type = "topk", grid = TRUE))
  expect_s3_class(rg, "cateomics_rule_grid")

  itp <- suppressWarnings(suppressMessages(cf_interpret(model)))
  expect_s3_class(itp, "cateomics_interpret")

  # Agreement between methods needs two of them. bcf's BART degenerates on the
  # multi-omics container and cf_compare records that and continues, so this
  # helper is reached with one surviving method on some containers and two on
  # others. Both outcomes are asserted rather than assuming either. The cause
  # is the row count the forest actually sees: the comparison fits each fold on
  # n - ceiling(n / n_folds) rows, which here is 40 against 48 covariates,
  # while the single-assay container's 33 covariates fit in the same 40 rows.
  surviving <- sum(vapply(seq_len(ncol(cmp$metadata$oof_tau_matrix)),
                          function(j) !all(is.na(cmp$metadata$oof_tau_matrix[, j])),
                          logical(1)))
  stab <- NULL
  if (surviving >= 2L) {
    stab <- suppressWarnings(cf_model_stability(cmp))
    expect_s3_class(stab, "cateomics_model_stability")
  } else {
    expect_error(suppressWarnings(cf_model_stability(cmp)),
                 class = "cateomics_stability_single_method")
  }

  grp <- suppressWarnings(cf_benefit_groups(model, method = "median"))
  expect_s3_class(grp, "cateomics_groups")

  diag <- suppressWarnings(suppressMessages(cf_diagnostics(model)))
  expect_s3_class(diag, "cateomics_diagnostics")

  if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    expect_s3_class(suppressWarnings(cf_aggregate_vi(cmp)),
                    "cateomics_aggregate_vi")
  }
  if (requireNamespace("shapviz", quietly = TRUE)) {
    expect_s3_class(suppressWarnings(as_shapviz(itp)), "shapviz")
  }

  invisible(list(cmp = cmp, model = model, ev = ev, rule = rule, rg = rg,
                 itp = itp, stab = stab, grp = grp, diag = diag))
}

test_that("major functions run on a data.frame (named-args)", {
  skip_if_not_installed("grf"); skip_if_not_installed("bcf")
  df <- .cov_make_df()
  wf <- .cov_run_workflow(df, "Y", "T", .cov_covnames())

  df2 <- cf_attach_results(df, wf$model, what = "all", groups = wf$grp)
  expect_true("cf_tau" %in% colnames(df2))
  # The average effect is one number for the fit, not a column of n copies.
  expect_false("cf_ate" %in% colnames(df2))
})

test_that("major functions run on a SummarizedExperiment (clinical + omics)", {
  skip_if_not_installed("grf"); skip_if_not_installed("bcf")
  skip_if_not_installed("SummarizedExperiment")
  se <- .cov_make_se()
  # X_baseline = clinical colData vars; the assay features are added as omics.
  wf <- .cov_run_workflow(se, "outcome", "treatment", c("age", "sex", "bmi"))

  se2 <- cf_attach_results(se, wf$model, what = "all", groups = wf$grp)
  expect_true(is_se(se2))
  expect_true("cf_tau" %in% colnames(SummarizedExperiment::colData(se2)))
})

test_that("major functions run on a multi-omics MultiAssayExperiment", {
  skip_if_not_installed("grf"); skip_if_not_installed("bcf")
  skip_if_not_installed("MultiAssayExperiment")
  mae <- .cov_make_mae()   # rnaseq + methyl + protein
  wf <- .cov_run_workflow(mae, "outcome", "treatment", c("age", "sex", "bmi"))
  expect_s3_class(wf$cmp, "cateomics_compare")
})

# --------------------------------------------------------------------------
# Block 2b: predict — in-sample and newdata
# --------------------------------------------------------------------------

test_that("predict returns in-sample and newdata estimates (R-native)", {
  skip_if_not_installed("grf")
  df  <- .cov_make_df()
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = .cov_covnames(),
           method = "grf", partition = "none", random_state = 1L)))
  nd <- .cov_make_df(n = 15L, seed = 99L)[, .cov_covnames(), drop = FALSE]

  expect_equal(length(predict(fit, type = "tau")), nrow(df))
  tau_new  <- predict(fit, newdata = nd, type = "tau")
  expect_equal(length(tau_new), nrow(nd))
  expect_equal(length(predict(fit, newdata = nd, type = "mu0")), nrow(nd))
  expect_equal(length(predict(fit, newdata = nd, type = "mu1")), nrow(nd))
  expect_length(predict(fit, newdata = nd, type = "ate"), 1L)
})

test_that("predict on newdata reloads the serialised Python model (basilisk)", {
  skip_if_no_cf_python()
  df <- .cov_make_df()
  nd <- .cov_make_df(n = 12L, seed = 98L)[, .cov_covnames(), drop = FALSE]

  # one econml (joblib) and one catenets (torch.save) backend cover both
  # serialisation paths through the newdata reload.
  for (m in c("xlearner", "tarnet")) {
    extra <- if (m == "tarnet") list(params = list(n_iter = 50L)) else list()
    fit <- suppressWarnings(suppressMessages(do.call(cf_fit, c(
      list(data = df, Y = "Y", T = "T", X_baseline = .cov_covnames(),
           method = m, partition = "none", random_state = 1L), extra))))
    tau_new <- predict(fit, newdata = nd, type = "tau")
    expect_equal(length(tau_new), nrow(nd))
    expect_false(anyNA(tau_new))
  }
})

# --------------------------------------------------------------------------
# Block 3: plot + print/summary smoke (no crash)
# --------------------------------------------------------------------------

test_that("plot family and print/summary run without error", {
  skip_if_not_installed("grf"); skip_if_not_installed("bcf")
  skip_if_not_installed("ggplot2")
  df <- .cov_make_df()
  wf <- .cov_run_workflow(df, "Y", "T", .cov_covnames())
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)

  ggsave1 <- function(p) expect_no_error(suppressMessages(suppressWarnings(
    ggplot2::ggsave(f, p, width = 4, height = 3))))
  ggsave1(cf_plot_compare(wf$cmp))
  ggsave1(cf_plot_evaluate(wf$ev, type = "gates"))
  ggsave1(cf_plot_evaluate(wf$ev, type = "rate"))
  ggsave1(cf_plot_evaluate(wf$ev, type = "blp"))
  ggsave1(cf_plot_policy_curve(wf$rg))
  ggsave1(cf_plot_stability(wf$stab))
  if (requireNamespace("shapviz", quietly = TRUE)) {
    ggsave1(cf_plot_interpret(wf$itp, type = "bar"))
    # Feature names live in the shapviz object (as_shapviz(itp)$S cols), not
    # itp$features. This block only runs when shapviz is present.
    feat <- colnames(as_shapviz(wf$itp)$S)[1]
    ggsave1(cf_plot_interpret(wf$itp, type = "dependence", feature = feat))
  }
  # `cf_plot_rule()` renders through DiagrammeR, so that is what gates it.
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    rtree <- suppressWarnings(cf_rule(wf$model, type = "tree", tree_depth = 2L))
    expect_no_error(suppressMessages(cf_plot_rule(rtree)))
  }

  ggsave1(cf_plot_cate(wf$model))                        # CATE distribution
  ggsave1(cf_plot_diagnostics(wf$diag))                  # cf_plot_diagnostics

  for (obj in list(wf$cmp, wf$model, wf$ev, wf$rule, wf$rg, wf$itp,
                   wf$stab, wf$grp, wf$diag)) {
    expect_output(print(obj))
  }
  expect_output(print(summary(wf$cmp)))
  expect_output(print(summary(wf$model)))
})

# --------------------------------------------------------------------------
# Block 4: propensity-score regimes (observational / known RCT / clip)
# --------------------------------------------------------------------------

test_that("propensity-score regimes: observational, fixed scalar/vector, clip", {
  skip_if_not_installed("grf")
  df <- .cov_make_df(); Xb <- .cov_covnames()

  # (a) observational: PS estimated from the data (default).
  expect_s3_class(suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
           method = "grf", partition = "none", random_state = 1L))), "cf_model")

  # (b) RCT / known PS as a scalar (1:1 randomised trial).
  if (requireNamespace("bcf", quietly = TRUE)) {
    f_rct <- suppressWarnings(suppressMessages(
      cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
              method = "bcf", partition = "none", propensity_score = 0.5,
              random_state = 1L, params = list(nburn = 50L, nsim = 100L))))
    expect_true(is.finite(f_rct$res$ate))

    # (c) known PS supplied as a length-n vector.
    f_vec <- suppressWarnings(suppressMessages(
      cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
              method = "bcf", partition = "none", propensity_score = plogis(0.4 *
                  df$X1), random_state = 1L, params = list(nburn = 50L,
                  nsim = 100L))))
    expect_true(is.finite(f_vec$res$ate))
  }

  # (d) PS trimming (clip extreme weights).
  expect_s3_class(suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb, method = "grf",
           partition = "none", clip = c(0.05, 0.95), random_state = 1L))),
    "cf_model")

  # (e) known PS threads through a PS-using Python method (basilisk).
  if (isTRUE(.cf_basilisk_ok())) {
    expect_s3_class(suppressWarnings(suppressMessages(
      cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb, method = "drlearner",
             partition = "none", propensity_score = 0.5, random_state = 1L))),
      "cf_model")
  }
})

# --------------------------------------------------------------------------
# Block 5: binary outcome -> risk-difference contract, all containers
# --------------------------------------------------------------------------

test_that("binary outcome is auto-detected and yields a risk-difference contract", {
  skip_if_not_installed("grf")
  df <- .cov_make_df(binary = TRUE); Xb <- .cov_covnames()

  # CATE method (grf): RD per individual within [-1, 1].
  f_grf <- suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
           method = "grf", partition = "none", random_state = 1L)))
  tau <- predict(f_grf, type = "tau")
  expect_equal(length(tau), nrow(df))
  expect_true(all(tau >= -1 & tau <= 1))

  # Binary Y through a SummarizedExperiment (clinical + omics).
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    se_b <- .cov_make_se(binary = TRUE)
    f_se <- suppressWarnings(suppressMessages(
      cf_fit(data = se_b, Y = "outcome", T = "treatment",
             X_baseline = c("age", "sex", "bmi"), method = "grf",
             partition = "none", random_state = 1L)))
    expect_s3_class(f_se, "cf_model")
  }

  # Binary Y through a Python CATE backend (basilisk).
  if (isTRUE(.cf_basilisk_ok())) {
    f_x <- suppressWarnings(suppressMessages(
      cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
             method = "xlearner", partition = "none", random_state = 1L)))
    expect_equal(length(predict(f_x, type = "tau")), nrow(df))
  }
})

# --------------------------------------------------------------------------
# Block 6: partition modes + assay feature selection
# --------------------------------------------------------------------------

test_that("partition oof vs none populate the OOF fields as specified", {
  skip_if_not_installed("grf")
  df <- .cov_make_df(); Xb <- .cov_covnames()

  f_oof <- suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
           method = "grf", partition = "oof", random_state = 1L)))
  expect_false(is.null(f_oof$res$tau_hat_oof))
  expect_true(is.list(f_oof$res$nuisance))

  f_none <- suppressWarnings(suppressMessages(
    cf_fit(data = df, Y = "Y", T = "T", X_baseline = Xb,
           method = "grf", partition = "none", random_state = 1L)))
  expect_null(f_none$res$tau_hat_oof)
})

