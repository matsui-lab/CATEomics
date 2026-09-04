# Smoke tests for the unified cf_diagnostics().
# Replaces the separate cf_balance_check / cf_overlap_check /
# cf_ps_diagnostics / cf_sensitivity / cf_diagnostics_report diagnostics.

make_diag_demo <- function(n = 200L, p = 5L, seed = 7L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(0.8 * X[, 1]))
  Y <- 0.5 * T + 0.4 * X[, 1] + 0.3 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X),
       formula = stats::as.formula(
         paste("Y ~ T |", paste(colnames(X), collapse = " + "))))
}


test_that("cf_diagnostics returns the four-view contract on a fitted model", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(d$formula, data = d$df, method = "grf",
           partition = "oof", random_state = 1L,
           num.trees = 200L, tune.parameters = "none")
  ))
  dg <- cf_diagnostics(fit, variables = c("X1", "X2"), top_k = 3L)
  expect_s3_class(dg, "cfomics_diagnostics")
  expect_true(all(c("overlap", "smd", "variable_overlap", "evalue",
                    "metadata") %in% names(dg)))
  # Overlap view.
  expect_true(is.numeric(dg$overlap$propensity_scores))
  expect_length(dg$overlap$overlap_region, 2L)
  # SMD Love view: top_k rows, ranked by |SMD|.
  expect_s3_class(dg$smd, "data.frame")
  expect_lte(nrow(dg$smd), 3L)
  expect_true(all(diff(dg$smd$abs_smd) <= 1e-9))  # descending
  # Requested-variable overlap.
  expect_true(is.data.frame(dg$variable_overlap))
  expect_setequal(unique(dg$variable_overlap$variable), c("X1", "X2"))
  # E-value.
  expect_true(isTRUE(dg$evalue$available))
  # evalue_point must equal the VanderWeele-Ding value of the RR derived from
  # the fit's ATE and outcome sd (not merely >= 1, which holds by construction).
  ev <- dg$evalue
  rr_expected <- exp(0.91 * abs(ev$ate) / ev$sd_y)
  ev_expected <- if (rr_expected <= 1) 1 else
    rr_expected + sqrt(rr_expected * (rr_expected - 1))
  expect_equal(ev$evalue_point, ev_expected, tolerance = 1e-8)
})


test_that("cf_diagnostics evalue_alpha rescales the CI-based E-value", {
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(d$formula, data = d$df, method = "grf",
           partition = "oof", random_state = 1L,
           num.trees = 200L, tune.parameters = "none")
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
    cf_fit(d$formula, data = d$df, method = "grf",
           partition = "oof", random_state = 1L,
           num.trees = 200L, tune.parameters = "none")
  ))
  dg <- cf_diagnostics(fit)
  ev <- dg$evalue
  # rr_point should equal exp(0.91 * |ATE| / sd_Y).
  expect_equal(ev$rr_point, exp(0.91 * abs(ev$ate) / ev$sd_y), tolerance = 1e-8)
})


test_that("cf_diagnostics rejects non-cf_model input and unknown variables", {
  expect_error(cf_diagnostics(list(a = 1)),
               class = "cfomics_invalid_input")
  skip_if_not_installed("grf")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(d$formula, data = d$df, method = "grf",
           partition = "oof", random_state = 1L,
           num.trees = 100L, tune.parameters = "none")
  ))
  expect_error(cf_diagnostics(fit, variables = "not_a_feature"),
               class = "cfomics_diagnostics_unknown_variables")
})


test_that("cf_diagnostics print and plots run", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_diag_demo()
  fit <- suppressMessages(suppressWarnings(
    cf_fit(d$formula, data = d$df, method = "grf",
           partition = "oof", random_state = 1L,
           num.trees = 100L, tune.parameters = "none")
  ))
  dg <- cf_diagnostics(fit, variables = "X1")
  expect_output(print(dg), "<cfomics_diagnostics>")
  for (ty in c("overlap", "smd", "variable", "evalue")) {
    p <- cf_plot_diagnostics(dg, type = ty)
    expect_s3_class(p, "ggplot")
  }
})
