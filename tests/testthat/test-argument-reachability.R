# Four documented entry-point arguments that no other test passed. Each one is
# checked by its effect, not by the absence of an error: an argument that is
# accepted and then dropped would pass a "does not error" test.
#
# Found while closing the omics-only covariate route: a scan of every cf_fit /
# cf_compare call in the shipped tests, examples and vignettes showed these
# four reaching no call at all, which is the same coverage shape that let the
# omics route ship broken.

ar_se <- function(n = 60L, g = 6L, seed = 5L, batch = FALSE) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)
  cnt <- matrix(stats::rpois(n * g, 60), nrow = g,
                dimnames = list(paste0("gene", seq_len(g)),
                                paste0("s", seq_len(n))))
  arm <- rep(0:1, each = n / 2)
  cd <- data.frame(resp = stats::rnorm(n) + 0.5 * arm,
                   arm = arm, age = stats::rnorm(n),
                   row.names = colnames(cnt))
  if (batch) {
    cd$site <- factor(rep(c("a", "b"), times = n / 2))
    # Shift one batch so a correction has something to remove.
    cnt[, cd$site == "b"] <- cnt[, cd$site == "b"] + 40L
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cnt), colData = S4Vectors::DataFrame(cd))
}

ar_df <- function(n = 80L, p = 4L, seed = 6L, binary_y = FALSE) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("X", seq_len(p))))
  T <- rep(0:1, each = n / 2)
  lin <- 0.8 * X[, 1] + T * (1 + X[, 2])
  Y <- if (binary_y) stats::rbinom(n, 1, stats::plogis(lin)) else
                     lin + stats::rnorm(n, sd = 0.5)
  data.frame(Y = Y, T = T, X)
}

test_that("BPPARAM reaches BiocParallel and the parallel run matches serial", {
  se <- ar_se()
  args <- list(Y = "resp", T = "arm", X_baseline = "age", omics = "all",
               data = se, methods = c("grf", "tlearner"), preprocess = FALSE,
               random_state = 0L, eval_nuisance = "ridge")
  skip_if_no_cf_python()
  ser <- do.call(cf_compare, c(args, list(BPPARAM = BiocParallel::SerialParam())))
  # Two workers only: R CMD check caps forks at two.
  par <- do.call(cf_compare, c(args,
                 list(BPPARAM = BiocParallel::MulticoreParam(2))))
  # No config may be lost under forking (a port clash once dropped the Python
  # backends silently).
  expect_identical(nrow(par$scoreboard), nrow(ser$scoreboard))
  expect_identical(sort(par$scoreboard$method), sort(ser$scoreboard$method))
  expect_equal(par$oof_tau, ser$oof_tau, tolerance = 1e-10)
  expect_equal(par$scoreboard[order(par$scoreboard$method), "r_loss"],
               ser$scoreboard[order(ser$scoreboard$method), "r_loss"],
               tolerance = 1e-10)
  # The value is handed to BiocParallel rather than inspected by CATEomics, so
  # something that is not a BiocParallelParam has to fail.
  expect_error(do.call(cf_compare, c(args, list(BPPARAM = "serial"))))
})

test_that("propensity_params reaches the propensity learner of bcf", {
  skip_if_not_installed("bcf")
  skip_if_not_installed("glmnet")
  skip_if_skip_bcf()
  df <- ar_df()
  base <- list(Y = "Y", T = "T", X_baseline = paste0("X", 1:4), data = df,
               method = "bcf", preprocess = FALSE, random_state = 0L,
               eval_nuisance = "ridge", model_propensity = "ridge",
               params = list(bcf = list(nburn = 20L, nsim = 20L,
                                        n_threads = 1L)))
  # Two penalty grids far apart: the fitted propensity has to move.
  loose <- do.call(cf_fit, c(base, list(
    propensity_params = list(lambda = c(1e-3, 5e-4, 1e-4)))))
  tight <- do.call(cf_fit, c(base, list(
    propensity_params = list(lambda = c(1e3, 5e2, 1e2)))))
  expect_true(all(is.finite(loose$res$ps)))
  expect_gt(max(abs(loose$res$ps - tight$res$ps)), 1e-6)
  # A name the learner does not take is refused rather than ignored.
  expect_error(do.call(cf_fit, c(base,
                 list(propensity_params = list(not_a_glmnet_arg = 1)))))
  # And the argument is refused outright for methods with no pluggable model.
  expect_error(cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = df,
                      method = "grf", propensity_params = list(nfolds = 3L)),
               class = "cateomics_propensity_params_unused")
})

test_that("batch reaches the correction step and is recorded", {
  skip_if_not_installed("sva")
  se <- ar_se(batch = TRUE)
  with_b <- cf_fit(Y = "resp", T = "arm", X_baseline = "age", omics = "all",
                   data = se, batch = "site", method = "grf",
                   preprocess = TRUE, random_state = 0L,
                   eval_nuisance = "ridge")
  without <- cf_fit(Y = "resp", T = "arm", X_baseline = "age", omics = "all",
                    data = se, method = "grf", preprocess = TRUE,
                    random_state = 0L, eval_nuisance = "ridge")
  rec <- with_b$meta$estimator$preprocess
  # Whether the correction ran belongs to the layer it ran on; the column and
  # the method are settings of the call, so they stay beside the layers.
  expect_true(isTRUE(rec$layers[[1L]]$batch_corrected))
  expect_identical(rec$batch_method, "combat")
  expect_identical(rec$batch_column, "site")
  expect_false(isTRUE(
    without$meta$estimator$preprocess$layers[[1L]]$batch_corrected))
  # The correction has to change the design matrix, not just the record.
  genes <- paste0("gene", 1:6)
  expect_gt(max(abs(with_b$data$X_preprocessed[, genes] -
                    without$data$X_preprocessed[, genes])), 1e-8)
  # A column that is not there is named in the error.
  expect_error(cf_fit(Y = "resp", T = "arm", X_baseline = "age",
                      omics = "all", data = se, batch = "no_such_column",
                      method = "grf", preprocess = TRUE, random_state = 0L))
})

test_that("outcome_type overrides the automatic detection", {
  df <- ar_df(binary_y = TRUE)
  base <- list(Y = "Y", T = "T", X_baseline = paste0("X", 1:4), data = df,
               method = "grf", preprocess = FALSE, random_state = 0L,
               eval_nuisance = "ridge")
  auto <- do.call(cf_fit, base)
  cont <- do.call(cf_fit, c(base, list(outcome_type = "continuous")))
  bin  <- do.call(cf_fit, c(base, list(outcome_type = "binary")))
  # A 0/1 outcome is detected as binary, and the explicit value is what is
  # recorded rather than the detected one.
  expect_identical(auto$meta$estimator$outcome_type, "binary")
  expect_identical(bin$meta$estimator$outcome_type, "binary")
  expect_identical(cont$meta$estimator$outcome_type, "continuous")
  # The override has to reach the scoring nuisance, which is fitted as a
  # classifier for a binary outcome and as a regression otherwise.
  expect_false(isTRUE(all.equal(bin$res$nuisance$mu0_oof,
                                cont$res$nuisance$mu0_oof)))
  expect_error(do.call(cf_fit, c(base, list(outcome_type = "nonsense"))))
})


# ---------------------------------------------------------------------------
# A value that matches nothing is an error wherever several are accepted
# ---------------------------------------------------------------------------

test_that("an unmatched value is refused, alone or mixed with good ones", {
  skip_if_not_installed("grf")
  set.seed(7)
  n <- 120L
  X <- matrix(rnorm(n * 4L), n, 4L); colnames(X) <- paste0("X", 1:4)
  T <- rbinom(n, 1L, 0.5)
  df <- data.frame(Y = X[, 1] + T * X[, 2] + rnorm(n), T = T, X)
  cv <- colnames(X)
  cmp <- suppressMessages(cf_compare(
    Y = "Y", T = "T", X_baseline = cv, data = df, methods = "grf",
    n_folds = 3L, verbose = FALSE, feature_method = "none",
    params = list(num.trees = 150L, tune.parameters = "none")))
  m <- cmp$models$grf

  # match.arg(several.ok = TRUE) keeps the good values and drops the bad one,
  # so the same misspelling errored alone and passed in a vector.
  expect_error(cf_evaluate(m, metrics = c("rate", "blp", "gates", "policy")),
               class = "cateomics_invalid_input")
  expect_error(cf_evaluate(m, metrics = "policy"),
               class = "cateomics_invalid_input")
  expect_error(suppressMessages(cf_interpret(m, surrogate = c("ridge", "nope"))),
               class = "cateomics_invalid_input")
  expect_error(cf_attach_results(df, m, what = c("tau", "nope")),
               class = "cateomics_invalid_input")
  expect_error(suppressMessages(cf_compare(
    Y = "Y", T = "T", X_baseline = cv, data = df, methods = "grf",
    n_folds = 3L, verbose = FALSE, feature_method = "none",
    metric = c("r_loss", "nope"))),
    class = "cateomics_invalid_input")

  # Partial matching, which match.arg allowed, still works.
  ev <- cf_evaluate(m, metrics = c("ra", "bl"))
  expect_setequal(setdiff(names(ev), "metadata"), c("rate", "blp"))
})
