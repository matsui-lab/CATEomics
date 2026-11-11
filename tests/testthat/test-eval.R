## Evaluation end-to-end smoke tests.
## Covers cf_compare / cf_select / cf_evaluate / cf_rule on an R-only
## path (grf).

make_demo <- function(n = 300, p = 4, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X), cv = colnames(X))
}

test_that("cf_compare returns cateomics_compare with scoreboard + models", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  expect_s3_class(cmp, "cateomics_compare")
  expect_true("r_loss" %in% colnames(cmp$scoreboard))
  expect_true("rank" %in% colnames(cmp$scoreboard))
  expect_named(cmp$models, "grf")
  expect_true(!is.null(cmp$models$grf$res$tau_hat_oof))
  expect_true(!is.null(cmp$models$grf$res$nuisance))
})

test_that("cf_compare requires methods argument", {
  d <- make_demo()
  expect_error(
    cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df),
    class = "cateomics_methods_required"
  )
})

test_that("cf_select picks the winner and reports the tie + between-fold CI", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  sel <- cf_select(cmp, metric = "r_loss")
  expect_s3_class(sel, "cateomics_select")
  expect_identical(sel$selected, "grf")
  expect_length(sel$selection_stats$between_fold_ci, 2L)
  # The old name is gone, asserted by name: a length check alone would pass
  # with both fields present.
  expect_null(sel$selection_stats[["naive_ci"]])
  expect_false(sel$selection_stats$tie_detected) # single method = no tie
  expect_identical(sel$selection_stats$tie_method, "report")
  # No honest_ci / honest_method / posi_multiplier fields.
  expect_null(sel$selection_stats$honest_ci)
  expect_null(sel$selection_stats$honest_method)
  expect_null(sel$selection_stats$posi_multiplier)
  # The config and the number of methods asked for are not kept a second time:
  # `model$meta$spec` and `length(selected)` already carry them.
  expect_null(sel[["spec"]])
  expect_null(sel[["top_n"]])
  expect_setequal(names(sel), c("selected", "model", "selection_stats", "metric"))

  # tie = "ignore" did not look, which is not the same as looked and separated.
  ign <- cf_select(cmp, metric = "r_loss", tie = "ignore")
  expect_true(is.na(ign$selection_stats$tie_detected))
  expect_identical(ign$selection_stats$tie_basis, "not assessed (tie = \'ignore\')")

  # The tie message states what was measured and steers nowhere.
  msg <- paste(capture.output(print(sel)), collapse = " ")
  expect_false(grepl("cf_ensemble", msg, fixed = TRUE))
})

test_that("cf_select rejects non-cateomics_compare input", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "none")
  expect_error(
    cf_select(fit),
    class = "cateomics_select_input_error"
  )
})

test_that("cf_evaluate returns RATE / BLP / GATES (merged metric)", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  ev <- cf_evaluate(cmp$models$grf,
                     metrics = c("rate", "blp", "gates"))
  expect_s3_class(ev, "cateomics_evaluate")
  # One ci_method for the object, in metadata; the settings live there too.
  expect_named(ev$rate, c("toc", "se", "ci_lower", "ci_upper",
                           "p_value", "p_value_one_sided", "autoc", "curve"))
  expect_named(ev$rate$autoc, c("value", "se", "ci_lower", "ci_upper",
                                 "p_value", "p_value_one_sided"))
  expect_null(ev$rate$ci_method)
  expect_null(ev$rate$target_quantile)
  expect_null(ev$rate$autoc$ci_method)
  expect_identical(ev$metadata$ci_method, "half-sample bootstrap SE")
  expect_identical(ev$metadata$target_quantile, 0.2)
  expect_null(ev$metadata$eval_sl_weights)
  # Calibration + BLP-HTE are merged into a single BLP metric that
  # reports the slope once and tests it against H0 = 1 (calibration) and
  # H0 = 0 (heterogeneity).
  expect_true("slope" %in% names(ev$blp))
  expect_true(all(c("calibration", "heterogeneity_loading") %in% names(ev$blp)))
  expect_true("p_value" %in% names(ev$blp$calibration))
  expect_true("p_value" %in% names(ev$blp$heterogeneity_loading))
  # beta = 0 has two readings, so the field is named after the coefficient.
  expect_null(ev$blp[["heterogeneity"]])
  # The calibration view draws `points`; the binned copy of what `gates`
  # already reports is gone.
  expect_null(ev$blp$plot_data)
  expect_null(ev$blp$slope$ci_method)
  expect_null(ev$calibration)   # removed standalone metric
  expect_null(ev$blp_hte)       # removed standalone metric
  expect_true("group_ate" %in% names(ev$gates))
  # The curve sits inside the quantity it is the curve of.
  expect_null(ev[["rate_curve"]])
  expect_s3_class(ev$rate$curve, "data.frame")
  # The curve is TOC(q), so AUTOC is the mean of the column that is plotted.
  expect_named(ev$rate$curve, c("quantile", "toc", "ci_lower", "ci_upper"))
  expect_equal(mean(ev$rate$curve$toc), ev$rate$autoc$value)
})

test_that("AUTOC is self-consistent and finite (cf_evaluate$rate$autoc)", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  ev <- cf_evaluate(cmp$models$grf, metrics = "rate")
  a <- ev$rate$autoc
  expect_true(is.finite(a$value) && is.finite(a$se) && a$se >= 0)
  expect_lte(a$ci_lower, a$value)
  expect_gte(a$ci_upper, a$value)
  expect_identical(ev$metadata$ci_method, "half-sample bootstrap SE")
})

test_that("AUTOC point matches grf::rank_average_treatment_effect (same DR scores)", {
  # Feed grf and .compute_autoc the SAME AIPW scores + priorities + q-grid;
  # the point estimates must agree up to grid discretisation.
  skip_if_not_installed("grf")
  set.seed(1)
  n <- 400L; p <- 6L
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- X[, 2] + W * (1 + 2 * X[, 1]) + rnorm(n)
  cf <- grf::causal_forest(X, Y, W, W.hat = 0.5, num.trees = 800, seed = 1)
  tau    <- predict(cf)$predictions
  scores <- grf::get_scores(cf)
  qg     <- seq(0.1, 1, by = 0.1)
  g    <- grf::rank_average_treatment_effect(cf, priorities = tau,
                                             target = "AUTOC", q = qg, R = 50)
  mine <- CATEomics:::.compute_autoc(tau, scores, qg)
  expect_lt(abs(mine$value - g$estimate) / abs(g$estimate), 0.10)  # within 10%
})

test_that("cf_rule produces threshold/topk/tree rules with DR policy value", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  model <- cmp$models$grf

  r1 <- cf_rule(model, type = "threshold", threshold = 0)
  expect_s3_class(r1, "cateomics_rule")
  expect_true(!is.null(r1$policy_value))
  expect_true(is.numeric(r1$policy_value$estimate))

  r2 <- cf_rule(model, type = "topk", topk_frac = 0.3)
  expect_identical(r2$type, "topk")
  expect_lte(r2$treated_frac, 0.35)

  r3 <- cf_rule(model, type = "tree", tree_depth = 2L)
  expect_identical(r3$type, "tree")
  expect_true(!is.null(r3$tree))
})

test_that("print.cateomics_rule shows all three contrasts and names a degenerate rule", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  model <- cmp$models$grf

  mid <- stats::median(model$res$tau_hat_oof, na.rm = TRUE)
  txt <- paste(capture.output(print(cf_rule(model, type = "threshold",
                                            threshold = mid))), collapse = "\n")
  # An improvement over treat-all can be the one that matters; it used not to
  # be printed at all.
  expect_match(txt, "improvement over same-budget random")
  expect_match(txt, "improvement over treat-all")
  expect_match(txt, "improvement over treat-none")
  expect_false(grepl("zero by construction", txt, fixed = TRUE))

  # A rule that treats everyone is one of the fixed policies it is compared
  # with, so two contrasts are zero by definition rather than by measurement.
  allt <- paste(capture.output(print(cf_rule(model, type = "threshold",
                                             threshold = -1e6))), collapse = "\n")
  expect_match(allt, "this rule is treat-all", fixed = TRUE)
  expect_match(allt, "zero by construction", fixed = TRUE)

  nont <- paste(capture.output(print(cf_rule(model, type = "threshold",
                                             threshold = 1e6))), collapse = "\n")
  expect_match(nont, "this rule is treat-none", fixed = TRUE)
})

test_that("cf_rule(grid=TRUE) returns a quantile-based default grid for threshold", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
          n_folds = 3L, verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  rg <- cf_rule(cmp$models$grf, type = "threshold", grid = TRUE)
  expect_s3_class(rg, "cateomics_rule_grid")
  expect_true(nrow(rg$results) > 5L)
  expect_named(rg$results,
                c("parameter", "n_treated", "policy_value",
                  "ci_lower", "ci_upper"))
})


test_that("the metric intervals are reproducible and follow n_bootstrap", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                    verbose = FALSE,
                    params = list(num.trees = 200L, tune.parameters = "none"))
  m <- cmp$models$grf
  a <- cf_evaluate(m, random_state = 7L, n_bootstrap = 100L)
  b <- cf_evaluate(m, random_state = 7L, n_bootstrap = 100L)
  expect_identical(a$rate$se, b$rate$se)
  expect_identical(a$gates$top_bottom_test$se, b$gates$top_bottom_test$se)
  expect_identical(a$blp$slope$se, b$blp$slope$se)

  # A different seed moves the standard error but not the point estimate.
  cc <- cf_evaluate(m, random_state = 8L, n_bootstrap = 100L)
  expect_identical(cc$rate$toc, a$rate$toc)
  expect_false(identical(cc$rate$se, a$rate$se))

  # Every interval brackets its point estimate, and the two BLP p-values come
  # from the bootstrap standard error.
  expect_lte(a$blp$slope$ci_lower, a$blp$slope$estimate)
  expect_gte(a$blp$slope$ci_upper, a$blp$slope$estimate)
  expect_equal(a$blp$heterogeneity_loading$z,
               a$blp$slope$estimate / a$blp$slope$se)
  expect_equal(a$blp$calibration$z,
               (a$blp$slope$estimate - 1) / a$blp$slope$se)
  expect_identical(a$metadata$ci_method, "half-sample bootstrap SE")
  expect_null(a$gates$top_bottom_test$ci_method)
  expect_true(all(is.finite(a$rate$curve$ci_lower)))
})


test_that("cf_rule resamples every component that depends on the rule", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                    verbose = FALSE,
                    params = list(num.trees = 200L, tune.parameters = "none"))
  m <- cmp$models$grf

  topk <- cf_rule(m, type = "topk", topk_frac = 0.3,
                  n_bootstrap = 100L, random_state = 5L)
  expect_identical(topk$policy_value$ci_method, "half-sample bootstrap SE")
  again <- cf_rule(m, type = "topk", topk_frac = 0.3,
                   n_bootstrap = 100L, random_state = 5L)
  expect_identical(topk$policy_value$ci_lower, again$policy_value$ci_lower)
  moved <- cf_rule(m, type = "topk", topk_frac = 0.3,
                   n_bootstrap = 100L, random_state = 6L)
  expect_identical(moved$policy_value$estimate, topk$policy_value$estimate)
  expect_false(identical(moved$policy_value$ci_lower,
                         topk$policy_value$ci_lower))

  # The replicates reflect the sample, not a changing rule, so a rule stated
  # as a fixed threshold is resampled too: the same rule is applied to the
  # drawn rows. The two fixed policies are plain means of per-row values and
  # keep the closed form whatever the rule is.
  thr <- cf_rule(m, type = "threshold", threshold = 0,
                 n_bootstrap = 100L, random_state = 5L)
  for (ty in list(topk, thr)) {
    expect_identical(ty$policy_value$ci_method, "half-sample bootstrap SE")
    expect_identical(ty$policy_value$baseline_treat_all$ci_method,
                     "influence function SE")
    expect_identical(ty$policy_value$baseline_treat_none$ci_method,
                     "influence function SE")
  }

  # A p-value on the three differences, none on the four levels, and one-sided
  # only where the alternative has a single natural direction.
  pv <- topk$policy_value
  for (nm in c("estimate", "se", "ci_lower", "ci_upper", "ci_method"))
    expect_true(nm %in% names(pv$baseline_random))
  expect_false("p_value" %in% names(pv$baseline_random))
  expect_false("p_value" %in% names(pv$baseline_treat_all))
  for (nm in c("improvement_over_random", "improvement_over_none")) {
    expect_true(all(c("p_value", "p_value_one_sided") %in% names(pv[[nm]])))
    # the p-value is the normal approximation on the standard error reported
    z <- pv[[nm]]$estimate / pv[[nm]]$se
    expect_equal(pv[[nm]]$p_value, 2 * stats::pnorm(-abs(z)), tolerance = 1e-12)
    expect_equal(pv[[nm]]$p_value_one_sided, stats::pnorm(-z),
                 tolerance = 1e-12)
  }
  expect_true("p_value" %in% names(pv$improvement_over_all))
  expect_false("p_value_one_sided" %in% names(pv$improvement_over_all))

  # The same-budget random baseline is defined through the budget the rule
  # actually spends, and that budget is realised rather than given, so it is
  # resampled even for a rule whose decisions are fixed. The decisive evidence
  # that the standard error came from the replicates and not from a closed form
  # is that it moves with the seed, which a closed form cannot do.
  mid <- stats::median(m$res$tau_hat_oof)
  half <- cf_rule(m, type = "threshold", threshold = mid,
                  n_bootstrap = 100L, random_state = 5L)
  expect_gt(half$treated_frac, 0)
  expect_lt(half$treated_frac, 1)
  for (nm in c("baseline_random", "improvement_over_random"))
    expect_identical(half$policy_value[[nm]]$ci_method,
                     "half-sample bootstrap SE")
  half2 <- cf_rule(m, type = "threshold", threshold = mid,
                   n_bootstrap = 100L, random_state = 6L)
  expect_identical(half2$policy_value$improvement_over_random$estimate,
                   half$policy_value$improvement_over_random$estimate)
  expect_false(identical(half2$policy_value$improvement_over_random$se,
                         half$policy_value$improvement_over_random$se))

  # A rule that treats everybody makes the same-budget random policy the
  # treat-all policy, so the two budget components are exact and keep the
  # closed form instead of reporting a bootstrap that every replicate refused.
  all_in <- cf_rule(m, type = "threshold", threshold = -Inf,
                    n_bootstrap = 100L, random_state = 5L)
  expect_identical(all_in$treated_frac, 1)
  expect_identical(all_in$policy_value$baseline_random$ci_method,
                   "influence function SE")
  expect_identical(all_in$policy_value$improvement_over_random$ci_method,
                   "influence function SE")
  expect_equal(all_in$policy_value$improvement_over_random$estimate, 0)
  all_in2 <- cf_rule(m, type = "threshold", threshold = -Inf,
                     n_bootstrap = 100L, random_state = 6L)
  expect_identical(all_in2$policy_value$baseline_random$se,
                   all_in$policy_value$baseline_random$se)

  # The label names the standard error that was actually used. Below 20 rows
  # the replicates are unusable, so a topk rule keeps the influence-function
  # errors and must say so rather than claim a bootstrap it did not get.
  small <- make_demo(n = 18L)
  cmp_s <- suppressWarnings(suppressMessages(
    cf_compare(Y = "Y", T = "T", X_baseline = small$cv, data = small$df, methods = "grf",
               n_folds = 3L, verbose = FALSE,
               eval_nuisance = "rf",
               params = list(num.trees = 100L,
                             tune.parameters = "none"))))
  tiny <- suppressWarnings(
    cf_rule(cmp_s$models$grf, type = "topk", topk_frac = 0.3,
            n_bootstrap = 50L, random_state = 5L))
  expect_identical(tiny$policy_value$ci_method, "influence function SE")
  expect_true(is.finite(tiny$policy_value$ci_lower))

  for (nm in c("improvement_over_random", "improvement_over_all",
               "improvement_over_none")) {
    b <- topk$policy_value[[nm]]
    expect_lte(b$ci_lower, b$estimate)
    expect_gte(b$ci_upper, b$estimate)
  }
})


test_that("the tie test differences the two configs fold by fold", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf",
                    base_learners = NULL, n_folds = 5L, verbose = FALSE,
                    params = list(num.trees = 150L, tune.parameters = "none"))
  # Two configs that share folds and nuisance: duplicate the scoreboard row
  # is not possible, so compare grf against a second grf-like config via the
  # base-learner grid on a pluggable learner is Python-bound; instead check
  # the helper directly on the stored per-fold matrix.
  pf <- cmp$metadata$per_fold_loss$r_loss
  expect_true(is.matrix(pf))
  expect_equal(ncol(pf), 5L)
  expect_setequal(rownames(pf), cmp$scoreboard$method)

  # A config differenced against itself has a zero difference and is a tie.
  pf2 <- rbind(a = pf[1, ], b = pf[1, ])
  md  <- list(per_fold_loss = list(r_loss = pf2))
  out <- CATEomics:::.cf_select_paired_diff(md, "r_loss", "a", "b")
  expect_equal(out$estimate, 0)
  expect_equal(out$k, 5L)

  # A constant shift moves the estimate but leaves the spread at zero, so the
  # interval excludes zero however small the shift is.
  pf3 <- rbind(a = pf[1, ], b = pf[1, ] + 0.01)
  out3 <- CATEomics:::.cf_select_paired_diff(
    list(per_fold_loss = list(r_loss = pf3)), "r_loss", "a", "b")
  expect_equal(out3$estimate, 0.01)
  expect_gt(out3$ci[1], 0)

  # Missing per-fold losses (e.g. an older object) disable the test.
  expect_null(CATEomics:::.cf_select_paired_diff(list(), "r_loss", "a", "b"))
})


test_that("cf_select reports the paired basis and a tighter interval", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  skip_if_skip_bcf()
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf", "bcf"),
                    n_folds = 5L, verbose = FALSE,
                    params = list(grf = list(num.trees = 150L,
                                             tune.parameters = "none"),
                                  bcf = list(nburn = 20L, nsim = 40L,
                                             n_chains = 1L)))
  ss <- suppressMessages(cf_select(cmp))$selection_stats
  expect_identical(ss$tie_basis, "fold-paired difference")
  expect_true(all(is.finite(ss$runner_up_diff_ci)))

  # The paired interval is narrower than adding the two between-fold
  # intervals, which is the point of pairing.
  sb <- cmp$scoreboard
  half <- (sb$r_loss_ci_upper - sb$r_loss_ci_lower) / 2
  unpaired_half <- sqrt(sum(half[order(sb$r_loss)][1:2]^2))
  paired_half <- diff(ss$runner_up_diff_ci) / 2
  expect_lt(paired_half, unpaired_half)
})


test_that("a constant tau-hat scores zero targeting", {
  # ecdf(tau)(tau) is 1 for every row when tau is constant, so selecting on
  # `>= 1 - q` used to put the whole sample in the "top q" and report a large
  # positive RATE for a method that found no heterogeneity at all.
  set.seed(1)
  n <- 200L
  Y_DR <- rnorm(n, mean = 1)
  tau  <- rep(0.7, n)
  expect_equal(CATEomics:::.toc_point(tau, Y_DR, 0.2), 0, tolerance = 1e-12)
  expect_equal(CATEomics:::.autoc_point(tau, Y_DR, seq(0.1, 1, by = 0.1)),
               0, tolerance = 1e-12)

  # The weight of the targeted set is exactly q * n, ties or no ties.
  expect_equal(sum(CATEomics:::.cf_top_q_weight(tau, 0.2)), 0.2 * n)
  set.seed(2)
  tau2 <- rnorm(n)
  expect_equal(sum(CATEomics:::.cf_top_q_weight(tau2, 0.35)), 0.35 * n)
  # The weighted set is the TOP of tau-hat, not the bottom: CATEomics scores the
  # units with the largest tau-hat and has no argument that flips that.
  w <- CATEomics:::.cf_top_q_weight(tau2, 0.25)
  expect_gt(mean(tau2[w > 0]), mean(tau2[w == 0]))
  expect_equal(min(tau2[w > 0]), sort(tau2, decreasing = TRUE)[0.25 * n])
})

test_that("cf_evaluate and cf_rule take no direction argument", {
  # The tail is not selectable: the convention is that a larger Y is the better
  # outcome, so a larger tau-hat is the better effect, and every function that
  # reads tau-hat follows it. Naming the removed argument is an error rather
  # than a silently ignored request for the other tail.
  expect_false("direction" %in% names(formals(cf_evaluate)))
  expect_false("direction" %in% names(formals(cf_rule)))
  expect_false("target_low" %in% names(formals(CATEomics:::.cf_top_q_weight)))
  expect_false("target_low" %in% names(formals(CATEomics:::.toc_point)))
  expect_false("target_low" %in% names(formals(CATEomics:::.gates_point)))

  d   <- make_demo(n = 120L)
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:4), data = d$df,
                method = "grf", partition = "oof", n_folds = 3L,
                eval_nuisance = "ridge", random_state = 1L,
                params = list(num.trees = 200L, tune.parameters = "none"))
  # cf_evaluate takes no dots, so R itself refuses the name; pin the reason so
  # the assertion cannot pass on an unrelated failure.
  expect_error(cf_evaluate(fit, direction = "negative"),
               regexp = "unused argument")
  expect_error(cf_rule(fit, type = "topk", direction = "negative"),
               class = "cateomics_unknown_argument")
  # and no object reports a direction any more
  ev <- cf_evaluate(fit, metrics = "rate", random_state = 1L)
  expect_false("direction" %in% names(ev$metadata))
  rl <- suppressMessages(cf_rule(fit, type = "topk", topk_frac = 0.3,
                                 compute_policy_value = FALSE,
                                 random_state = 1L))
  expect_false("direction" %in% names(rl))
  gr <- suppressMessages(cf_rule(fit, type = "topk", grid = c(0.2, 0.4),
                                 compute_policy_value = FALSE,
                                 random_state = 1L))
  expect_false("direction" %in% names(gr))
})


# cf_rule: the policy value is measured, not just bracketed ----------------

.rule_fixture <- function(n = 300L, seed = 3L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  T_vec <- rbinom(n, 1L, 0.5)
  # tau = 2 * x1: who benefits is known, so a rule that treats the top of
  # tau-hat must beat a coin flip and its mirror image must lose to one.
  Y <- x2 + T_vec * (2 * x1) + rnorm(n, sd = 0.3)
  df <- data.frame(Y = Y, T = T_vec, x1 = x1, x2 = x2)
  suppressMessages(cf_fit(data = df, Y = "Y", T = "T",
                          X_baseline = c("x1", "x2"), method = "grf",
                          partition = "oof", random_state = 1L,
                          propensity_score = 0.5))
}

test_that("treating the top of tau-hat beats random, and the mirror rule loses", {
  fit <- .rule_fixture()
  good <- suppressMessages(cf_rule(fit, type = "topk", topk_frac = 0.3,
                                   random_state = 1L))
  # The mirror rule treats the same number of units from the other end. There
  # is no argument for that any more, so it is built here from the model's own
  # out-of-fold effects and scored through the same policy-value helper.
  tau  <- fit$res$tau_hat_oof
  k    <- max(1L, floor(length(tau) * 0.3))
  d_bad <- as.integer(tau <= sort(tau)[k])
  bad_pv <- CATEomics:::.compute_policy_value(fit, d_bad, n_bootstrap = 200L,
                                            random_state = 1L)

  # An implementation that treats the complement, or ignores the decision,
  # would still return an interval that brackets its own estimate.
  expect_gt(good$policy_value$improvement_over_random$estimate, 0.2)
  expect_lt(bad_pv$improvement_over_random$estimate, -0.2)
  expect_gt(good$policy_value$estimate, bad_pv$estimate)
  # Same budget on both sides, so the gap is the ranking, not the count.
  expect_equal(sum(good$decisions), sum(d_bad))
})

test_that("the treat-all and treat-none baselines are the AIPW arm means", {
  fit <- .rule_fixture()
  rule <- suppressMessages(cf_rule(fit, type = "topk", topk_frac = 0.3,
                                   random_state = 1L))
  nu <- fit$res$nuisance
  Y  <- fit$data$Y
  T_vec <- fit$data$T
  # Written out in plain arithmetic rather than through the package helper,
  # so a wrong pseudo-outcome cannot agree with itself.
  ey1 <- mean(nu$mu1_oof + T_vec / nu$pi_oof * (Y - nu$mu1_oof))
  ey0 <- mean(nu$mu0_oof + (1 - T_vec) / (1 - nu$pi_oof) * (Y - nu$mu0_oof))
  expect_equal(rule$policy_value$baseline_treat_all$estimate, ey1,
               tolerance = 1e-10)
  expect_equal(rule$policy_value$baseline_treat_none$estimate, ey0,
               tolerance = 1e-10)
})


# The top-q weights, including how they split a tie ----------------------

test_that("the top-q weights select by rank and share a tie across the cut", {
  # Four of ten units are targeted; three of them are tied on the value at
  # the cut, so the two remaining slots are split between the three.
  tau <- c(5, 4, 3, 3, 3, 2, 1, 0, -1, -2)
  w <- CATEomics:::.cf_top_q_weight(tau, q = 0.4)
  expect_equal(sum(w), 4)
  expect_equal(w[1:2], c(1, 1))
  expect_equal(w[3:5], rep(2 / 3, 3))
  expect_true(all(w[6:10] == 0))

  # Taking `tau >= cut` with weight 1, which is what the grf cross-check
  # cannot see, would target five units instead of four.
  expect_false(sum(w) == 5)

  # No ties: a plain top-k.
  expect_equal(CATEomics:::.cf_top_q_weight(c(5, 4, 3, 2, 1), q = 0.4),
               c(1, 1, 0, 0, 0))
  # Constant tau: everyone is equally targeted, which is what makes RATE 0.
  expect_equal(CATEomics:::.cf_top_q_weight(rep(2, 10), q = 0.3), rep(0.3, 10))
  # Always the top: the order of the input does not move the selection.
  expect_equal(CATEomics:::.cf_top_q_weight(c(1, 2, 3, 4, 5), q = 0.4),
               c(0, 0, 0, 1, 1))
})

# ---------------------------------------------------------------------------
# One-sided p-values beside the two-sided ones. The reported intervals stay
# two-sided, so both statements are available and neither is reconstructed by
# the reader.
# ---------------------------------------------------------------------------

test_that(".cf_normal_p is one construction for both sides", {
  pv <- CATEomics:::.cf_normal_p(estimate = 0.5, se = 0.2)
  expect_equal(pv$z, 2.5)
  expect_equal(pv$two_sided, 2 * stats::pnorm(-2.5))
  expect_equal(pv$one_sided, stats::pnorm(-2.5))
  # the two sides of the one-sided test partition the probability
  expect_equal(pv$one_sided + stats::pnorm(2.5), 1)
  # a null other than zero shifts the estimate, not the construction
  pc <- CATEomics:::.cf_normal_p(estimate = 0.5, se = 0.2, h0 = 1)
  expect_equal(pc$z, -2.5)
  expect_equal(pc$two_sided, 2 * stats::pnorm(-2.5))
  # not estimable in, not estimable out
  pn <- CATEomics:::.cf_normal_p(estimate = 0.5, se = NA_real_)
  expect_true(is.na(pn$two_sided))
  expect_true(is.na(pn$one_sided))
})

test_that("cf_evaluate reports a one-sided p beside the two-sided one", {
  skip_if_not_installed("grf")
  d   <- make_demo(n = 200L)
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:4), data = d$df,
                method = "grf", partition = "oof", n_folds = 3L,
                eval_nuisance = "ridge", random_state = 1L,
                params = list(num.trees = 200L, tune.parameters = "none"))
  ev <- cf_evaluate(fit, metrics = c("rate", "blp", "gates"),
                    random_state = 1L)

  # the four quantities that have a signed direction of interest
  quads <- list(
    toc   = list(est = ev$rate$toc,                   obj = ev$rate),
    autoc = list(est = ev$rate$autoc$value,           obj = ev$rate$autoc),
    gates = list(est = ev$gates$top_bottom_test$diff, obj = ev$gates$top_bottom_test)
  )
  for (nm in names(quads)) {
    o <- quads[[nm]]$obj; est <- quads[[nm]]$est
    expect_true("p_value" %in% names(o), info = nm)
    expect_true("p_value_one_sided" %in% names(o), info = nm)
    z <- est / o$se
    expect_equal(o$p_value, 2 * stats::pnorm(-abs(z)), info = nm)
    expect_equal(o$p_value_one_sided, stats::pnorm(-z), info = nm)
  }
  # the BLP heterogeneity test reads the slope and the slope's SE
  z_het <- ev$blp$slope$estimate / ev$blp$slope$se
  expect_equal(ev$blp$heterogeneity$p_value, 2 * stats::pnorm(-abs(z_het)))
  expect_equal(ev$blp$heterogeneity$p_value_one_sided, stats::pnorm(-z_het))

  # the calibration null has no one-sided form: a slope either side of 1 is a
  # failure, so only the two-sided p-value is reported
  expect_false("p_value_one_sided" %in% names(ev$blp$calibration))
  expect_equal(ev$blp$calibration$p_value,
               2 * stats::pnorm(-abs((ev$blp$slope$estimate - 1) /
                                       ev$blp$slope$se)))

  # intervals stay two-sided
  expect_true(all(c("ci_lower", "ci_upper") %in% names(ev$rate)))
  expect_equal(ev$rate$ci_lower, ev$rate$toc - 1.96 * ev$rate$se)

  # and print shows both
  expect_output(print(ev), "one-sided")
})

test_that("print does not drop a line when a reported field is absent", {
  # sprintf() returns character(0) when any argument is NULL, so an object
  # written before a field existed used to lose its whole line -- TOC, AUTOC,
  # the heterogeneity test and GATES all vanished at once, with nothing saying
  # so. The line must still appear, with NA where the field is missing.
  skip_if_not_installed("grf")
  d   <- make_demo(n = 120L)
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:4), data = d$df,
                method = "grf", partition = "oof", n_folds = 3L,
                eval_nuisance = "ridge", random_state = 1L,
                params = list(num.trees = 200L, tune.parameters = "none"))
  ev  <- cf_evaluate(fit, random_state = 1L)

  older <- ev
  older$rate$p_value <- NULL
  older$rate$p_value_one_sided <- NULL
  older$rate$autoc$p_value <- NULL
  older$rate$autoc$p_value_one_sided <- NULL
  older$blp$heterogeneity$p_value_one_sided <- NULL
  older$gates$top_bottom_test$p_value_one_sided <- NULL

  out <- capture.output(print(older))
  expect_true(any(grepl("TOC @", out)))
  expect_true(any(grepl("AUTOC", out)))
  expect_true(any(grepl("heterogeneity", out)))
  expect_true(any(grepl("GATES", out)))
  expect_true(any(grepl("one-sided NA", out)))
  # the current object prints a number in the same slot
  expect_output(print(ev), "one-sided [0-9]")
})


# A tree that never split --------------------------------------------------

test_that("print says the tree did not split, whichever action the stump took", {
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  mk <- function(sign) {
    set.seed(5)
    n <- 120L
    df <- data.frame(T = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
    # One sign for everyone, so the stump assigns the same action to all rows.
    df$Y <- df$x2 + sign * df$T * (2 + abs(df$x1)) + rnorm(n, sd = 0.3)
    fit <- suppressMessages(cf_fit(
      data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "grf",
      partition = "oof", n_folds = 3L, random_state = 1L,
      propensity_score = 0.5,
      params = list(num.trees = 200L, tune.parameters = "none")))
    suppressMessages(cf_rule(fit, type = "tree", tree_depth = 1L,
                             min_node_size = 1000L))
  }
  treat_all  <- mk(+1)
  treat_none <- mk(-1)
  expect_length(treat_all$tree$nodes, 1L)
  expect_length(treat_none$tree$nodes, 1L)
  expect_equal(treat_all$treated_frac, 1)
  expect_equal(treat_none$treated_frac, 0)
  # The note exists to say the tree never split; a stump that assigns control
  # needs it just as much as one that treats.
  expect_output(print(treat_all),
                "the tree did not split: this rule is treat-all", fixed = TRUE)
  expect_output(print(treat_none),
                "the tree did not split: this rule is treat-none", fixed = TRUE)
})
