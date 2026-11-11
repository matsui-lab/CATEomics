## Thin contract regression test for cf_compare / cf_select: the
## `cateomics_not_yet_implemented` placeholder does not appear, and input-shape
## errors fire with their dedicated classes instead.

test_that("cf_compare requires methods argument (replaces 'not_yet_implemented' stub)", {
  set.seed(1)
  data <- data.frame(
    Y  = rnorm(20),
    T  = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_compare(Y = "Y", T = "T", X_baseline = "X1", data = data),
    class = "cateomics_methods_required"
  )
})

test_that("cf_select rejects non-cateomics_compare input (replaces stub)", {
  expect_error(
    cf_select(list(scoreboard = data.frame(method = "grf", r_loss = 1)),
              metric = "r_loss"),
    class = "cateomics_select_input_error"
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
    cv = c("X1", "X2")
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
    cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = c("grf"),
            metric = c("r_loss", "dr_loss"), n_folds = 3L, verbose = FALSE,
            params = list(num.trees = 100L, tune.parameters = "none"))
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
      cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "bcf", partition = "none",
              model_propensity = "logistic", random_state = 42L,
              params = list(nburn = 100L, nsim = 100L))
    ))
  )
})

test_that("cf_compare hands the per-method fits a resolved device", {
  # A literal "auto" must never reach a worker: the workers compare the
  # string, so an unresolved value pins every Python method to the CPU no
  # matter what the caller asked for. Asserting that `.resolve_device()`
  # does not return "auto" would prove nothing -- it has four returns and
  # all of them are "cpu" or "gpu". Capture what cf_compare hands the
  # dispatcher instead, which is where the value was actually lost.
  seen <- NULL
  testthat::local_mocked_bindings(
    .kfold_oof_predict_matrix = function(..., device) {
      seen <<- device
      stop("the forwarded device is all this test needs")
    },
    .package = "CATEomics"
  )
  set.seed(8)
  df0 <- data.frame(Y = rnorm(40), T = rbinom(40, 1L, 0.5),
                    x1 = rnorm(40), x2 = rnorm(40))
  suppressMessages(suppressWarnings(try(cf_compare(
    data = df0, Y = "Y", T = "T", X_baseline = c("x1", "x2"),
    methods = "grf", device = "auto", propensity_score = 0.5,
    random_state = 1L, verbose = FALSE), silent = TRUE)))
  expect_false(is.null(seen))
  expect_false(identical(seen, "auto"))
  expect_true(seen %in% c("cpu", "gpu"))
})

test_that("a per-method fit runs on the device cf_compare was asked for", {
  skip_if_no_cf_python()
  set.seed(8)
  n <- 60L
  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1L, 0.5),
                   x1 = rnorm(n), x2 = rnorm(n))
  fit_dev <- function(dev) suppressMessages(suppressWarnings(cf_compare(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"),
    methods = "tarnet", device = dev, propensity_score = 0.5,
    n_folds = 2L, random_state = 1L, verbose = FALSE)))

  # Reading `formals(cf_compare)`, as this test first did, passes while the
  # value never reaches a backend. So does reading a meta field cf_compare
  # never writes: cf_compare dispatches `.kfold_oof_predict_matrix`, not
  # `cf_fit`, so `meta$provenance$device` stays NULL and any expectation
  # against it holds vacuously. Read what the backend itself reported.
  expect_identical(fit_dev("cpu")$models[["tarnet"]]$res$device_effective,
                   "cpu")

  # The helper's own default is "cpu", so a CPU-only host cannot tell a
  # forwarded "cpu" from no forwarding at all. Only a host with a GPU can.
  skip_if_not(identical(
    CATEomics:::.resolve_device("auto", method = "python",
                              r_native_methods = character(0)), "gpu"))
  # expect_false(identical(x, "cpu")) would also hold for NULL or NA, which
  # is what a backend that stopped reporting would give. Require the string.
  expect_match(fit_dev("gpu")$models[["tarnet"]]$res$device_effective, "^cuda")
  expect_match(fit_dev("auto")$models[["tarnet"]]$res$device_effective, "^cuda")
})


## eval_nuisance separation + known propensity ---------------------

test_that(".cateomics_resolve_eval_nuisance broadcasts strings and validates lists", {
  s <- CATEomics:::.cateomics_resolve_eval_nuisance("rf")
  expect_equal(s, list(outcome = "rf", propensity = "rf"))

  s2 <- CATEomics:::.cateomics_resolve_eval_nuisance(list(outcome = "xgboost",
                                                      propensity = "known"))
  expect_equal(s2, list(outcome = "xgboost", propensity = "known"))

  # Missing partner fills from sibling.
  s3 <- CATEomics:::.cateomics_resolve_eval_nuisance(list(outcome = "ridge"))
  expect_equal(s3$propensity, "ridge")

  expect_error(
    CATEomics:::.cateomics_resolve_eval_nuisance("not_a_learner"),
    class = "cateomics_eval_nuisance_unknown"
  )
})

test_that("a known evaluation propensity is used by both entries", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 80L)
  prm <- list(num.trees = 100L, tune.parameters = "none")
  # The score is the only public way to say "this propensity is known", and
  # it has to reach the evaluation nuisance from either entry. If it does not,
  # the losses and the metrics score an estimated propensity in a randomised
  # trial, with nothing in the output saying so.
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf",
               eval_propensity_score = 0.5, n_folds = 3L, verbose = FALSE,
               params = prm)
  ))
  expect_s3_class(cmp, "cateomics_compare")
  expect_true(all(abs(cmp$metadata$nuisance_oof$pi_oof - 0.5) < 1e-8))

  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf",
           eval_propensity_score = 0.5, n_folds = 3L, params = prm)
  ))
  expect_true(all(abs(fit$res$nuisance$pi_oof - 0.5) < 1e-8))
})

test_that("a vector known evaluation propensity lands on its own rows", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 80L)
  pv <- seq(0.2, 0.8, length.out = 80L)
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf",
               eval_propensity_score = pv, n_folds = 3L, verbose = FALSE,
               params = list(num.trees = 100L, tune.parameters = "none"))
  ))
  # Each fold writes the held-out rows, so a permuted vector would show up as a
  # mismatch here rather than as a plausible-looking propensity.
  expect_equal(cmp$metadata$nuisance_oof$pi_oof, pv, tolerance = 1e-12)
})

test_that("an out-of-range or mis-sized known propensity is refused", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  prm <- list(num.trees = 100L, tune.parameters = "none")
  expect_error(
    suppressMessages(suppressWarnings(
      cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                 verbose = FALSE, eval_propensity_score = 1.2, params = prm)
    )),
    class = "cateomics_eval_known_ps_out_of_range"
  )
  expect_error(
    suppressMessages(suppressWarnings(
      cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                 verbose = FALSE, eval_propensity_score = c(0.4, 0.6),
                 params = prm)
    )),
    class = "cateomics_eval_known_ps_length_mismatch"
  )
})

test_that("eval_nuisance takes a learner per slot at either entry", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  prm <- list(num.trees = 100L, tune.parameters = "none")
  # The estimator side names the outcome and the propensity learner
  # separately, and `eval_library` / `eval_params` already take the slot form,
  # so `eval_nuisance` takes it too. Both entries accept the same shapes.
  for (ok in list(list(outcome = "rf", propensity = "ridge"),
                  list(outcome = "rf"),
                  list(propensity = "ridge"))) {
    cmp <- suppressMessages(suppressWarnings(
      cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                 verbose = FALSE, eval_nuisance = ok, params = prm)))
    fit <- suppressMessages(suppressWarnings(
      cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", n_folds = 3L,
             eval_nuisance = ok, params = prm)))
    expect_identical(cmp$models[[1L]]$meta$estimator$eval_nuisance_resolved,
                     fit$meta$estimator$eval_nuisance_resolved)
  }
  # "known" carries no score, so it is not a learner a caller can name; the
  # known propensity goes to `eval_propensity_score` as a value.
  for (bad in list(list(outcome = "rf", propensity = "known"),
                   list(outcome = "nope"), list(foo = "rf"))) {
    expect_error(
      suppressMessages(suppressWarnings(
        cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
                   verbose = FALSE, eval_nuisance = bad, params = prm)
      )),
      class = "cateomics_invalid_eval_nuisance"
    )
    expect_error(
      suppressMessages(suppressWarnings(
        cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", n_folds = 3L,
               eval_nuisance = bad, params = prm)
      )),
      class = "cateomics_invalid_eval_nuisance"
    )
  }
})

test_that("a known propensity on one side only is reported", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  prm <- list(num.trees = 100L, tune.parameters = "none")
  expect_message(
    suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf",
                            n_folds = 3L, propensity_score = 0.5,
                            params = prm)),
    class = "cateomics_known_ps_one_side"
  )
  # Both sides set: nothing to report.
  expect_no_message(
    suppressWarnings(cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf",
                            n_folds = 3L, propensity_score = 0.5,
                            eval_propensity_score = 0.5, params = prm))
  )
})


## propensity surface (consumed by cf_diagnostics overlap) ----

test_that("grf surfaces W.hat as fit$res$ps", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 60L)
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", partition = "none",
            params = list(num.trees = 100L, tune.parameters = "none"))
  ))
  ps <- fit$res$ps
  expect_true(is.numeric(ps))
  expect_equal(length(ps), 60L)
  expect_true(all(ps > 0 & ps < 1))
  expect_identical(fit$res$ps_type, "oob")
  expect_identical(predict(fit, type = "ps"), ps)
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

  sb <- CATEomics:::.build_compare_scoreboard(
    ids = ids, metric = "factual_loss", oof_tau = oof_tau,
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    oof_mu0 = oof_mu0, oof_mu1 = oof_mu1,
    config_df = data.frame(label = ids, method = ids,
                           base = NA_character_, stringsAsFactors = FALSE)
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
  sb <- CATEomics:::.build_compare_scoreboard(
    ids = "mA", metric = "factual_loss", oof_tau = list(mA = rnorm(n)),
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    oof_mu0 = list(mA = rep(NA_real_, n)), oof_mu1 = list(mA = rep(NA_real_, n)),
    config_df = data.frame(label = "mA", method = "mA",
                           base = NA_character_, stringsAsFactors = FALSE)
  )
  expect_true(is.na(sb$factual_loss[1]))
})

# ---------------------------------------------------------------------------
# Fold-variation columns: how the spread of the out-of-fold tau-hat divides
# between folds and within them.
# ---------------------------------------------------------------------------

test_that(".build_compare_scoreboard splits tau_sd between and within folds", {
  set.seed(21)
  n  <- 60L
  Y  <- rnorm(n)
  T  <- rbinom(n, 1L, 0.5)
  fold_id <- rep(1:5, length.out = n)
  nuis <- list(mu0_oof = rnorm(n), mu1_oof = rnorm(n),
               pi_oof  = rep(0.5, n), m_oof = rnorm(n))
  # three deliberately different shapes of the same spread:
  #   purely between folds (constant inside each fold),
  #   purely within folds  (every fold has the same mean),
  #   and a mixture.
  fold_means  <- c(-2, -1, 0, 1, 2)[fold_id]
  within_only <- rep(c(-1, 1), length.out = n)
  oof_tau <- list(
    between = fold_means,
    within  = within_only,
    mixed   = fold_means + within_only
  )
  ids <- names(oof_tau)
  sb <- CATEomics:::.build_compare_scoreboard(
    ids = ids, metric = "r_loss", oof_tau = oof_tau,
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    config_df = data.frame(label = ids, method = ids,
                           base = NA_character_, stringsAsFactors = FALSE)
  )
  expect_true(all(c("within_fold_sd", "between_fold_sd",
                    "between_fold_share") %in% colnames(sb)))
  rownames(sb) <- sb$label
  # a score that is constant inside each fold is all between-fold
  expect_equal(sb["between", "between_fold_share"], 1, tolerance = 1e-12)
  expect_equal(sb["between", "within_fold_sd"], 0, tolerance = 1e-12)
  # a score whose fold means are all equal is all within-fold
  expect_equal(sb["within", "between_fold_share"], 0, tolerance = 1e-12)
  expect_equal(sb["within", "between_fold_sd"], 0, tolerance = 1e-12)
  # the mixture sits strictly between the two
  expect_gt(sb["mixed", "between_fold_share"], 0)
  expect_lt(sb["mixed", "between_fold_share"], 1)
  # and the split is exact for every row
  expect_equal(sb$within_fold_sd^2 + sb$between_fold_sd^2, sb$tau_sd^2,
               tolerance = 1e-10)
})

test_that(".build_compare_scoreboard takes both fold variances over the same rows", {
  # A method whose out-of-fold effect is missing on some rows: the share must
  # still be the analysis-of-variance decomposition of the rows it does have,
  # not a ratio of two different subsets.
  set.seed(23)
  n  <- 60L
  Y  <- rnorm(n)
  T  <- rbinom(n, 1L, 0.5)
  fold_id <- rep(1:4, length.out = n)
  nuis <- list(mu0_oof = rnorm(n), mu1_oof = rnorm(n),
               pi_oof  = rep(0.5, n), m_oof = rnorm(n))
  tau <- rnorm(n)
  tau[c(3L, 17L, 41L)] <- NA_real_
  sb <- CATEomics:::.build_compare_scoreboard(
    ids = "gappy", metric = "r_loss", oof_tau = list(gappy = tau),
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    config_df = data.frame(label = "gappy", method = "gappy",
                           base = NA_character_, stringsAsFactors = FALSE))
  ok <- !is.na(tau)
  t_ok <- tau[ok]; g_ok <- fold_id[ok]
  expected <- sum((stats::ave(t_ok, g_ok, FUN = mean) - mean(t_ok))^2) /
    sum((t_ok - mean(t_ok))^2)
  expect_equal(sb$between_fold_share, expected, tolerance = 1e-12)
  # a sum of squares over a subset of itself: always a proportion
  expect_gte(sb$between_fold_share, 0)
  expect_lte(sb$between_fold_share, 1)
  # and the split still reproduces tau_sd, which is also taken over those rows
  expect_equal(sb$within_fold_sd^2 + sb$between_fold_sd^2, sb$tau_sd^2,
               tolerance = 1e-10)
  expect_equal(sb$tau_sd, stats::sd(t_ok))
})

test_that(".build_compare_scoreboard reports NA fold variation for a constant tau-hat", {
  set.seed(22)
  n  <- 40L
  Y  <- rnorm(n)
  T  <- rbinom(n, 1L, 0.5)
  fold_id <- rep(1:4, length.out = n)
  nuis <- list(mu0_oof = rnorm(n), mu1_oof = rnorm(n),
               pi_oof  = rep(0.5, n), m_oof = rnorm(n))
  sb <- CATEomics:::.build_compare_scoreboard(
    ids = "flat", metric = "r_loss", oof_tau = list(flat = rep(0.5, n)),
    Y = Y, T = T, nuis = nuis, fold_id = fold_id, failed = list(),
    config_df = data.frame(label = "flat", method = "flat",
                           base = NA_character_, stringsAsFactors = FALSE)
  )
  # There is no spread to divide, so the share is undefined rather than 0 or 1;
  # tau_sd = 0 sits beside it and says why.
  expect_equal(sb$tau_sd, 0)
  expect_true(is.na(sb$between_fold_share))
  expect_true(is.na(sb$within_fold_sd))
  expect_true(is.na(sb$between_fold_sd))
})

test_that("cf_compare carries the fold-variation columns and prints the share", {
  se  <- cf_example_data(n = 60L, p = 100L, n_prog = 8L, n_pred = 4L, seed = 3L)
  cmp <- cf_compare(Y = "outcome", T = "treatment", data = se, omics = "all",
                    X_baseline = c("age", "sex"), methods = "grf",
                    n_folds = 3L, n_features = 40L, eval_nuisance = "ridge",
                    random_state = 0L, BPPARAM = BiocParallel::SerialParam())
  sb <- cmp$scoreboard
  expect_true(all(c("within_fold_sd", "between_fold_sd",
                    "between_fold_share") %in% colnames(sb)))
  # recomputable from what the object stores, so the column is not a free number
  tau <- cmp$metadata$oof_tau_matrix[, sb$label[1]]
  fid <- cmp$metadata$fold_id
  expect_equal(sb$between_fold_share[1],
               stats::var(stats::ave(tau, fid)) / stats::var(tau),
               tolerance = 1e-12)
  # the compact print shows the share; the full print shows every column
  expect_output(print(cmp), "between_fold_share")
  expect_output(print(cmp, full = TRUE), "within_fold_sd")
})

test_that("the three displays of a comparison say different things", {
  skip_if_not_installed("grf")
  d <- .make_demo_compare(n = 80L)
  cmp <- suppressMessages(suppressWarnings(
    cf_compare(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
               methods = "grf", n_folds = 3L, verbose = FALSE,
               params = list(num.trees = 100L, tune.parameters = "none"))))
  # summary() answers what a selection rests on: which candidate comes first,
  # under which metric, and in which direction that metric is read.
  out <- capture.output(summary(cmp))
  expect_true(any(grepl("^  best: ", out)))
  expect_true(any(grepl("selection metric: r_loss \\(lower is better\\)", out)))
  expect_true(any(grepl("every column: print", out)))
  # and it is not the scoreboard: the columns only the full print shows must
  # not appear, or summary() would be print(full = TRUE) under another name.
  expect_false(any(grepl("within_fold_sd|factual_squared_score", out)))
  # it fits on a screen
  expect_lt(length(out), 20L)
  # the full print still shows every column
  expect_output(print(cmp, full = TRUE), "within_fold_sd")
})
