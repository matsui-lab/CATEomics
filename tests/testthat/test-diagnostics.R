# Smoke tests for cf_diagnostics(), which returns the overlap, balance,
# per-variable and E-value views in one call.

make_diag_demo <- function(n = 200L, p = 5L, seed = 7L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(0.8 * X[, 1]))
  Y <- 0.5 * T + 0.4 * X[, 1] + 0.3 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X), cv = colnames(X))
}


test_that("cf_diagnostics returns the four-view contract on a fitted model", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 200L,
                tune.parameters = "none"))
  ))
  dg <- cf_diagnostics(fit, variables = c("X1", "X2"), top_k = 3L)
  expect_s3_class(dg, "cateomics_diagnostics")
  expect_true(all(c("overlap", "smd", "variable_overlap", "evalue",
                    "metadata") %in% names(dg)))
  # Overlap view: the model's own (here held-out) propensity, and how many
  # observations the estimator had to clip.
  expect_true(is.numeric(dg$overlap$propensity_scores))
  expect_length(dg$overlap$overlap_region, 2L)
  expect_identical(dg$overlap$ps_type, "oof")
  expect_identical(dg$overlap$propensity_scores, fit$res$ps)
  expect_true(is.numeric(dg$overlap$n_clipped))
  expect_equal(dg$overlap$clipped_pct,
               dg$overlap$n_clipped / length(dg$overlap$propensity_scores))
  expect_false("n_extreme" %in% names(dg$overlap))
  # SMD Love view: top_k rows, ranked by |SMD|.
  expect_s3_class(dg$smd, "data.frame")
  expect_lte(nrow(dg$smd), 3L)
  expect_true(all(diff(dg$smd$abs_smd) <= 1e-9))  # descending
  # Requested-variable overlap.
  expect_true(is.data.frame(dg$variable_overlap))
  expect_setequal(unique(dg$variable_overlap$variable), c("X1", "X2"))
  # E-value.
  expect_true(isTRUE(dg$evalue$available))
  # The inputs of the formula are pinned to the fit, not read back out of
  # the diagnostic: standardising by the spread of tau-hat instead of the
  # spread of Y would otherwise agree with itself.
  ev <- dg$evalue
  expect_equal(ev$ate, fit$res$ate, tolerance = 1e-12)
  # The denominator is the pooled within-arm SD, recomputed here from the fit.
  # It is not the SD of the whole outcome: variance splits into a within-arm
  # and a between-arm part, and the second is made by the effect itself.
  y_fit <- as.numeric(fit$data$Y); w_fit <- as.integer(fit$data$T)
  n1 <- sum(w_fit == 1L); n0 <- sum(w_fit == 0L)
  pooled_expected <- sqrt(((n1 - 1) * stats::var(y_fit[w_fit == 1L]) +
                           (n0 - 1) * stats::var(y_fit[w_fit == 0L])) /
                          (n1 + n0 - 2))
  expect_equal(ev$sd_y, pooled_expected, tolerance = 1e-12)
  expect_identical(ev$sd_basis, "pooled within-arm SD")
  expect_lt(ev$sd_y, stats::sd(y_fit))
  rr_expected <- exp(0.91 * abs(ev$ate) / ev$sd_y)
  ev_expected <- if (rr_expected <= 1) 1 else
    rr_expected + sqrt(rr_expected * (rr_expected - 1))
  expect_equal(ev$evalue_point, ev_expected, tolerance = 1e-8)
})


test_that("cf_diagnostics evalue_alpha rescales the CI-based E-value", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 200L,
                tune.parameters = "none"))
  ))

  dg_05 <- cf_diagnostics(fit, evalue_alpha = 0.05)
  dg_10 <- cf_diagnostics(fit, evalue_alpha = 0.10)
  ev05 <- dg_05$evalue
  ev10 <- dg_10$evalue

  # (a) Both views are available and the point E-value is unaffected by alpha.
  expect_true(isTRUE(ev05$available))
  expect_true(isTRUE(ev10$available))
  expect_equal(ev05$evalue_point, ev10$evalue_point, tolerance = 1e-10)

  # (b) evalue_alpha reaches behavior: the CI-based E-value differs across
  # alpha, and a wider alpha (0.10) gives a narrower CI, moving the CI-bound
  # E-value toward the point E-value (i.e. no smaller than at alpha = 0.05).
  expect_false(isTRUE(all.equal(ev05$evalue_ci, ev10$evalue_ci)))
  expect_gte(ev10$evalue_ci, ev05$evalue_ci)
  expect_lte(abs(ev10$evalue_ci - ev10$evalue_point),
             abs(ev05$evalue_ci - ev05$evalue_point))

  # (c) Invariant: alpha = 0.05 reproduces the model's stored 95% ATE CI, so
  # the CI-based E-value equals the value derived from the stored bounds.
  cil <- unname(fit$res$summary$ate_ci_lower)
  ciu <- unname(fit$res$summary$ate_ci_upper)
  expect_equal(unname(ev05$ci[["lower"]]), cil, tolerance = 1e-10)
  expect_equal(unname(ev05$ci[["upper"]]), ciu, tolerance = 1e-10)
})


test_that("cf_diagnostics E-value uses the 0.91 continuous coefficient", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 200L,
                tune.parameters = "none"))
  ))
  dg <- cf_diagnostics(fit)
  ev <- dg$evalue
  # rr_point should equal exp(0.91 * d), with the effect and the spread taken
  # from the fit rather than from the diagnostic's own record of them.
  y_fit <- as.numeric(fit$data$Y); w_fit <- as.integer(fit$data$T)
  n1 <- sum(w_fit == 1L); n0 <- sum(w_fit == 0L)
  pooled_expected <- sqrt(((n1 - 1) * stats::var(y_fit[w_fit == 1L]) +
                           (n0 - 1) * stats::var(y_fit[w_fit == 0L])) /
                          (n1 + n0 - 2))
  expect_equal(ev$sd_y, pooled_expected, tolerance = 1e-12)
  expect_equal(ev$rr_point,
               exp(0.91 * abs(fit$res$ate) / pooled_expected),
               tolerance = 1e-8)
})


test_that("cf_diagnostics says what the interval contributed, or why nothing", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 200L,
                tune.parameters = "none"))
  ))

  # An interval that excludes the null contributes the limit nearer to it.
  ev <- cf_diagnostics(fit)$evalue
  if (isTRUE(fit$res$summary$ate_ci_lower > 0) ||
      isTRUE(fit$res$summary$ate_ci_upper < 0)) {
    expect_identical(ev$ci_basis, "confidence limit nearer the null")
    expect_true(is.finite(ev$evalue_ci))
  }

  # An interval that holds the null gives 1, and says so.
  crossing <- fit
  crossing$res$summary$ate_ci_lower <- -abs(fit$res$ate) - 1
  crossing$res$summary$ate_ci_upper <-  abs(fit$res$ate) + 1
  ev_cross <- cf_diagnostics(crossing)$evalue
  expect_equal(ev_cross$evalue_ci, 1)
  expect_identical(ev_cross$ci_basis, "the interval crosses the null")

  # No interval at all is NA, not 1: "the evidence does not survive the
  # interval" and "there is no interval" are different statements.
  none <- fit
  none$res$summary$ate_ci_lower <- NA_real_
  none$res$summary$ate_ci_upper <- NA_real_
  none$meta$estimator$bootstrap <- FALSE
  ev_none <- cf_diagnostics(none)$evalue
  expect_true(is.na(ev_none$evalue_ci))
  expect_true(is.na(ev_none$rr_ci))
  expect_match(ev_none$ci_basis, "no interval on the average effect")
  expect_match(ev_none$ci_basis, "bootstrap_ci = FALSE", fixed = TRUE)
  # The point estimate is unaffected by the missing interval.
  expect_equal(ev_none$evalue_point, ev$evalue_point, tolerance = 1e-12)
  # print survives the NA.
  expect_output(print(cf_diagnostics(none)), "not available", fixed = TRUE)

  # Requested but absent is a different sentence.
  none2 <- none
  none2$meta$estimator$bootstrap <- TRUE
  expect_identical(cf_diagnostics(none2)$evalue$ci_basis,
                   "no interval on the average effect")
})


test_that("cf_diagnostics names the spread it standardised by, and falls back", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 200L,
                tune.parameters = "none"))
  ))
  # Both arms present: the pooled within-arm SD.
  ev <- cf_diagnostics(fit)$evalue
  expect_identical(ev$sd_basis, "pooled within-arm SD")

  # One arm reduced to a single row: the pooled variance is undefined, so the
  # whole-sample SD is used and the record says so rather than keeping the
  # label of a quantity that was not computed.
  thin <- fit
  keep <- c(which(fit$data$T == 1L), which(fit$data$T == 0L)[1L])
  thin$data$Y <- fit$data$Y[keep]
  thin$data$T <- fit$data$T[keep]
  ev_thin <- cf_diagnostics(thin)$evalue
  expect_identical(ev_thin$sd_basis, "SD of the observed outcome")
  expect_equal(ev_thin$sd_y, stats::sd(as.numeric(thin$data$Y)),
               tolerance = 1e-12)
})


test_that("cf_diagnostics rejects non-cf_model input and unknown variables", {
  expect_error(cf_diagnostics(list(a = 1)),
               class = "cateomics_invalid_input")
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 100L,
                tune.parameters = "none"))
  ))
  expect_error(cf_diagnostics(fit, variables = "not_a_feature"),
               class = "cateomics_diagnostics_unknown_variables")
})


test_that("cf_diagnostics print and plots run", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "oof",
            random_state = 1L, params = list(num.trees = 100L,
                tune.parameters = "none"))
  ))
  dg <- cf_diagnostics(fit, variables = "X1")
  expect_output(print(dg), "<cateomics_diagnostics>")
  for (ty in c("overlap", "smd", "variable")) {
    p <- cf_plot_diagnostics(dg, type = ty)
    expect_s3_class(p, "ggplot")
  }
})


test_that("the clipped count follows the clip bounds, not a fixed threshold", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit_tight <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "none",
           clip = c(0.4, 0.6), random_state = 1L,
           params = list(num.trees = 200L, tune.parameters = "none"))))
  ov <- cf_diagnostics(fit_tight)$overlap
  expect_identical(ov$clip, c(0.4, 0.6))
  expect_gt(ov$n_clipped, 0L)
  expect_true(all(ov$propensity_scores >= 0.4 & ov$propensity_scores <= 0.6))
})


test_that("a method without a propensity loses only the overlap view", {
  skip_if_no_cf_python()
  d <- make_diag_demo(n = 60L, p = 3L)
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "slearner", partition = "none",
           random_state = 1L)))
  expect_message(dg <- cf_diagnostics(fit, variables = "X1"),
                 class = "cateomics_overlap_view_skipped")
  expect_null(dg$overlap)
  expect_s3_class(dg$smd, "data.frame")
  expect_true(is.data.frame(dg$variable_overlap))
  expect_true(isTRUE(dg$evalue$available))
  expect_output(print(dg), "not available")
})


test_that("an ensemble reports the shared evaluation nuisance", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()
  d <- make_diag_demo(n = 80L, p = 4L)
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf", "bcf"),
               n_folds = 3L, random_state = 1L,
               params = list(grf = list(num.trees = 100L,
                                        tune.parameters = "none"),
                             bcf = list(nburn = 10L, nsim = 10L,
                                        n_chains = 1L)))))
  ens <- cf_ensemble(cmp, method = "mean")
  ov  <- cf_diagnostics(ens)$overlap
  expect_identical(ov$ps_type, "evaluation nuisance")
  expect_identical(ov$propensity_scores, as.numeric(ens$res$nuisance$pi_oof))
  expect_identical(ov$clip, c(1e-4, 1 - 1e-4))
})


test_that("the imbalance count covers every covariate, not just the shown rows", {
  skip_if_not_installed("grf")
  set.seed(2)
  n <- 80L; p <- 40L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(0.9 * X[, 1]))
  Y <- T + X[, 1] + rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = colnames(X), method = "grf",
    partition = "none", feature_method = "none",
    params = list(num.trees = 150L, tune.parameters = "none")))
  dg <- cf_diagnostics(fit, top_k = 5L)
  expect_equal(nrow(dg$smd), 5L)
  expect_equal(attr(dg$smd, "n_total"), p)
  # The table reports the distribution over every covariate and calls none of
  # them balanced: at n << p an |SMD| of 0.1 is what randomisation alone
  # produces for most of them.
  expect_null(dg$smd$balanced)
  expect_null(attr(dg$smd, "n_imbalanced"))
  expect_true(attr(dg$smd, "max") >= attr(dg$smd, "p90"))
  expect_true(attr(dg$smd, "p90") >= attr(dg$smd, "median"))
  txt <- paste(capture.output(print(dg)), collapse = " ")
  expect_match(txt, "|SMD| median", fixed = TRUE)
  expect_false(grepl("covariates |SMD| >=", txt, fixed = TRUE))

  # A level outside (0, 1) used to be ignored, leaving the stored interval.
  expect_error(cf_diagnostics(fit, evalue_alpha = 0),
               class = "cateomics_invalid_input")

  # The stored flag for "the ranges do not meet" is gone; support_fraction
  # says it.
  expect_null(dg$overlap$no_overlap)
})

test_that("a known propensity is reported as the one value it is", {
  skip_if_not_installed("grf")
  set.seed(4)
  n <- 80L
  X <- matrix(rnorm(n * 5L), n, 5L); colnames(X) <- paste0("X", 1:5)
  T <- rbinom(n, 1, 0.5)
  Y <- T + X[, 1] + rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = colnames(X), method = "grf",
    partition = "none", feature_method = "none", propensity_score = 0.5,
    params = list(num.trees = 150L, tune.parameters = "none")))
  dg <- cf_diagnostics(fit)
  txt <- paste(capture.output(print(dg)), collapse = " ")
  # "100% in support" on a constant score is a tautology, not a finding.
  expect_match(txt, "the propensity takes one value", fixed = TRUE)
  expect_false(grepl("in support", txt, fixed = TRUE))
  # The figure says the same thing, with its reasons.
  expect_error(cf_plot_diagnostics(dg, type = "overlap"),
               class = "cateomics_plot_no_data")
})
