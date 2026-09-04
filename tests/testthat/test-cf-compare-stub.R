## Thin contract regression test for cf_compare / cf_select: the
## `cfomics_not_yet_implemented` placeholder does not appear, and input-shape
## errors fire with their dedicated classes instead.

test_that("cf_compare requires methods argument (replaces 'not_yet_implemented' stub)", {
  set.seed(1)
  data <- data.frame(
    Y  = rnorm(20),
    T  = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_compare(Y ~ T | X1, data),
    class = "cfomics_methods_required"
  )
})

test_that("cf_select rejects non-cfomics_compare input (replaces stub)", {
  expect_error(
    cf_select(list(scoreboard = data.frame(method = "grf", r_loss = 1)),
              metric = "r_loss"),
    class = "cfomics_select_input_error"
  )
})

## cf_select reports tie_detected only (no formal selection-inference correction).


## ---------------------------------------------------------------------------
## cf_select regression tests
## ---------------------------------------------------------------------------

# Local demo helper (mirrors test-eval.R, scoped to this file).
.make_demo_compare <- function(n = 120L, seed = 42L) {
  set.seed(seed)
  X1 <- rnorm(n); X2 <- rnorm(n)
  T_vec <- rbinom(n, 1L, 0.5)
  Y <- 0.2 * X1 + 0.5 * X2 * T_vec + rnorm(n, sd = 0.5)
  list(
    df = data.frame(Y = Y, T = T_vec, X1 = X1, X2 = X2),
    formula = Y ~ T | X1 + X2
  )
}

test_that("cf_compare default metric is the c('r_loss','dr_loss','factual_loss') triplet", {
  # No data needed: just inspect the default value via formals().
  default_metric <- eval(formals(cf_compare)$metric)
  expect_setequal(default_metric, c("r_loss", "dr_loss", "factual_loss"))
})

test_that("cf_compare scoreboard exposes rank_<metric> when multi-metric", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare()
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(d$formula, data = d$df, methods = c("grf"),
                metric = c("r_loss", "dr_loss"), n_folds = 3L, verbose = FALSE,
                num.trees = 100L, tune.parameters = "none")
  ))
  sb <- cmp$scoreboard
  expect_true(all(c("rank", "rank_r_loss", "rank_dr_loss") %in% names(sb)))
})

test_that("bcf formals filter strips unknown args (R8B-follow regression)", {
  skip_if_not_installed("bcf")
  # Pass a base-learner unification arg that bcf::bcf() does not accept.
  # The formals filter should strip it before forwarding.
  d <- .make_demo_compare()
  expect_no_error(
    suppressMessages(suppressWarnings(
      cf_fit(d$formula, data = d$df, method = "bcf",
              partition = "none",
              model_propensity = "logistic",  # unknown to bcf::bcf
              random_state = 42L,
              n_burn = 100L, n_sim = 100L)
    ))
  )
})

test_that("cf_compare forwards device to per-method dispatch", {
  # Smoke: verify the formal arg now lives downstream by inspecting the
  # device match.arg vector. The actual GPU smoke is covered by helper-gpu.R.
  expect_true("device" %in% names(formals(cf_compare)))
  expect_setequal(eval(formals(cf_compare)$device), c("auto", "cpu", "gpu"))
})


## eval_nuisance separation + known propensity ---------------------

test_that(".cfomics_resolve_eval_nuisance broadcasts strings and validates lists", {
  s <- cfomics:::.cfomics_resolve_eval_nuisance("rf")
  expect_equal(s, list(outcome = "rf", propensity = "rf"))

  s2 <- cfomics:::.cfomics_resolve_eval_nuisance(list(outcome = "xgboost",
                                                      propensity = "known"))
  expect_equal(s2, list(outcome = "xgboost", propensity = "known"))

  # Missing partner fills from sibling.
  s3 <- cfomics:::.cfomics_resolve_eval_nuisance(list(outcome = "ridge"))
  expect_equal(s3$propensity, "ridge")

  expect_error(
    cfomics:::.cfomics_resolve_eval_nuisance("not_a_learner"),
    class = "cfomics_eval_nuisance_unknown"
  )
})

test_that("cf_compare(eval_nuisance = list(propensity = 'known'), eval_propensity_score = 0.5) works", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 80L)
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(d$formula, data = d$df, methods = c("grf"),
                eval_nuisance = list(outcome = "rf", propensity = "known"),
                eval_propensity_score = 0.5,
                n_folds = 3L, verbose = FALSE,
                num.trees = 100L, tune.parameters = "none")
  ))
  expect_s3_class(cmp, "cfomics_compare")
  # All OOF propensity values should equal 0.5 (within numeric tolerance).
  expect_true(all(abs(cmp$metadata$nuisance_oof$pi_oof - 0.5) < 1e-8))
})

test_that("known propensity without eval_propensity_score errors", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  expect_error(
    suppressMessages(suppressWarnings(
      cf_compare(d$formula, data = d$df, methods = c("grf"),
                  eval_nuisance = list(propensity = "known"),
                  n_folds = 3L, verbose = FALSE,
                  num.trees = 100L, tune.parameters = "none")
    )),
    class = "cfomics_eval_known_ps_missing"
  )
})


## propensity surface (consumed by cf_diagnostics overlap) ----

test_that("grf surfaces W.hat as fit$res$propensity", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  fit <- suppressMessages(suppressWarnings(
    cf_fit(d$formula, data = d$df, method = "grf",
            partition = "none",
            num.trees = 100L, tune.parameters = "none")
  ))
  ps <- fit$res$propensity
  expect_true(is.numeric(ps))
  expect_equal(length(ps), 60L)
  expect_true(all(ps > 0 & ps < 1))
})


# factual_loss is method-specific (per-method OOF mu0/mu1) ------------------

test_that(".build_compare_scoreboard factual_loss uses each method's own mu", {
  set.seed(11)
  n  <- 40L
  Y  <- rbinom(n, 1L, 0.5)
  T  <- rbinom(n, 1L, 0.5)
  fold_id <- rep(1:4, length.out = n)
  nuis <- list(
    mu0_oof = runif(n, 0.2, 0.5),
    mu1_oof = runif(n, 0.4, 0.7),
    pi_oof  = rep(0.5, n),
    m_oof   = runif(n, 0.3, 0.6)
  )
  ids     <- c("mA", "mB")
  oof_tau <- list(mA = rnorm(n, 0.1, 0.05), mB = rnorm(n, 0.2, 0.05))
  # Two methods with DIFFERENT held-out potential outcomes.
  oof_mu0 <- list(mA = rep(0.3, n),  mB = rep(0.6, n))
  oof_mu1 <- list(mA = rep(0.5, n),  mB = rep(0.2, n))

  sb <- cfomics:::.build_compare_scoreboard(
    ids = ids, metric = "factual_loss", oof_tau = oof_tau,
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    oof_mu0 = oof_mu0, oof_mu1 = oof_mu1
  )
  fl <- sb$factual_loss
  expect_length(fl, 2L)
  expect_true(all(is.finite(fl)))
  # Different mu per method must yield different factual loss (the old bug made
  # these identical because a shared nuisance mu was used for every method).
  expect_false(isTRUE(all.equal(fl[1], fl[2])))
})

test_that(".build_compare_scoreboard factual_loss is NA when a method has no mu", {
  set.seed(12)
  n  <- 30L
  Y  <- rbinom(n, 1L, 0.5)
  T  <- rbinom(n, 1L, 0.5)
  fold_id <- rep(1:3, length.out = n)
  nuis <- list(mu0_oof = runif(n), mu1_oof = runif(n),
               pi_oof = rep(0.5, n), m_oof = runif(n))
  sb <- cfomics:::.build_compare_scoreboard(
    ids = "mA", metric = "factual_loss", oof_tau = list(mA = rnorm(n)),
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    oof_mu0 = list(mA = rep(NA_real_, n)), oof_mu1 = list(mA = rep(NA_real_, n))
  )
  expect_true(is.na(sb$factual_loss[1]))
})
