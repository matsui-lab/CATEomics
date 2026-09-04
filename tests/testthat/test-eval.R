## Evaluation end-to-end smoke tests.
## Covers cf_compare / cf_select / cf_evaluate / cf_rule on an R-only
## path (grf).

make_demo <- function(n = 300, p = 4, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X),
       formula = stats::as.formula(
         paste("Y ~ T |", paste(colnames(X), collapse = " + "))
       ))
}

test_that("cf_compare returns cfomics_compare with scoreboard + models", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  expect_s3_class(cmp, "cfomics_compare")
  expect_true("r_loss" %in% colnames(cmp$scoreboard))
  expect_true("rank" %in% colnames(cmp$scoreboard))
  expect_named(cmp$models, "grf")
  expect_true(!is.null(cmp$models$grf$res$tau_hat_oof))
  expect_true(!is.null(cmp$models$grf$res$nuisance))
})

test_that("cf_compare requires methods argument", {
  d <- make_demo()
  expect_error(
    cf_compare(d$formula, data = d$df),
    class = "cfomics_methods_required"
  )
})

test_that("cf_select picks the winner and reports tie + naive CI", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  sel <- cf_select(cmp, metric = "r_loss")
  expect_s3_class(sel, "cfomics_select")
  expect_identical(sel$selected, "grf")
  expect_length(sel$selection_stats$naive_ci, 2L)
  expect_false(sel$selection_stats$tie_detected) # single method = no tie
  expect_identical(sel$selection_stats$tie_method, "report")
  # No honest_ci / honest_method / posi_multiplier fields.
  expect_null(sel$selection_stats$honest_ci)
  expect_null(sel$selection_stats$honest_method)
  expect_null(sel$selection_stats$posi_multiplier)
})

test_that("cf_select rejects non-cfomics_compare input", {
  skip_if_not_installed("grf")
  d <- make_demo()
  fit <- cf_fit(d$formula, data = d$df, method = "grf", partition = "none")
  expect_error(
    cf_select(fit),
    class = "cfomics_select_input_error"
  )
})

test_that("cf_evaluate returns RATE / BLP / GATES (merged metric)", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  ev <- cf_evaluate(cmp$models$grf,
                     metrics = c("rate", "blp", "gates"))
  expect_s3_class(ev, "cfomics_evaluate")
  expect_named(ev$rate, c("value", "ci_lower", "ci_upper", "target_quantile",
                           "ci_method", "autoc"))
  expect_named(ev$rate$autoc, c("value", "se", "ci_lower", "ci_upper",
                                 "ci_method"))
  # Calibration + BLP-HTE are merged into a single BLP metric that
  # reports the slope once and tests it against H0 = 1 (calibration) and
  # H0 = 0 (heterogeneity).
  expect_true("slope" %in% names(ev$blp))
  expect_true(all(c("calibration", "heterogeneity") %in% names(ev$blp)))
  expect_true("p_value" %in% names(ev$blp$calibration))
  expect_true("p_value" %in% names(ev$blp$heterogeneity))
  expect_null(ev$calibration)   # removed standalone metric
  expect_null(ev$blp_hte)       # removed standalone metric
  expect_true("group_ate" %in% names(ev$gates))
  expect_s3_class(ev$rate_curve, "data.frame")
  expect_named(ev$rate_curve, c("quantile", "rate", "ci_lower", "ci_upper"))
})

test_that("AUTOC is self-consistent and finite (cf_evaluate$rate$autoc)", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  ev <- cf_evaluate(cmp$models$grf, metrics = "rate")
  a <- ev$rate$autoc
  expect_true(is.finite(a$value) && is.finite(a$se) && a$se >= 0)
  expect_lte(a$ci_lower, a$value)
  expect_gte(a$ci_upper, a$value)
  expect_match(a$ci_method, "influence-function")
})

test_that("AUTOC point matches grf::rank_average_treatment_effect (same DR scores)", {
  # Feed grf and .compute_autoc the SAME AIPW scores + priorities + q-grid;
  # the point estimates must agree up to grid discretisation. grf's SE is
  # bootstrap-based (propagates ranking uncertainty); cfomics' IF SE conditions
  # on the ranking, so only the POINT is compared here.
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
  mine <- cfomics:::.compute_autoc(tau, scores, qg)
  expect_lt(abs(mine$value - g$estimate) / abs(g$estimate), 0.10)  # within 10%
})

test_that("cf_rule produces threshold/topk/tree rules with DR policy value", {
  skip_if_not_installed("grf")
  skip_if_not_installed("rpart")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  model <- cmp$models$grf

  r1 <- cf_rule(model, type = "threshold", threshold = 0)
  expect_s3_class(r1, "cfomics_rule")
  expect_true(!is.null(r1$policy_value))
  expect_true(is.numeric(r1$policy_value$estimate))

  r2 <- cf_rule(model, type = "topk", topk_frac = 0.3)
  expect_identical(r2$type, "topk")
  expect_lte(r2$treated_frac, 0.35)

  r3 <- cf_rule(model, type = "tree", tree_depth = 2L)
  expect_identical(r3$type, "tree")
  expect_true(!is.null(r3$tree))
})

test_that("cf_rule(grid=TRUE) returns a quantile-based default grid for threshold", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = c("grf"),
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  rg <- cf_rule(cmp$models$grf, type = "threshold", grid = TRUE)
  expect_s3_class(rg, "cfomics_rule_grid")
  expect_true(nrow(rg$results) > 5L)
  expect_named(rg$results,
                c("parameter", "n_treated", "policy_value",
                  "ci_lower", "ci_upper"))
})
