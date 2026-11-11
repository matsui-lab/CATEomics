# The Super Learner evaluation nuisance.
#
# The assertions are measurements in situations whose answer is known: a linear
# outcome makes a penalised linear member the right regression, and a noise
# outcome with a level makes the intercept the right one. Replacing the solved
# weights with equal weights fails those two and the restricted-library one.
# The randomised-treatment case asserts the prediction rather than the weights,
# because with no signal the penalised members collapse onto the intercept and
# the solution over the resulting collinear columns is not unique.

sl_data <- function(n = 120L, p = 8L, seed = 7L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  T <- stats::rbinom(n, 1L, 0.5)
  list(X = X, T = T,
       y_linear = 4 + 2 * X[, 1] - 1.5 * X[, 2] + stats::rnorm(n, sd = 0.3),
       y_noise  = 4 + stats::rnorm(n),
       folds = .create_stratified_folds(T, n_folds = 5L, random_state = 0L))
}

lib_all <- function() .CF_SL_LIBRARY_DEFAULT
no_params <- function(members) stats::setNames(rep(list(list()), length(members)), members)

test_that("the two entries that fit a nuisance default to the Super Learner", {
  expect_identical(formals(cf_fit)$eval_nuisance, "super_learner")
  expect_identical(formals(cf_compare)$eval_nuisance, "super_learner")
  expect_true("eval_library" %in% names(formals(cf_fit)))
  expect_true("eval_library" %in% names(formals(cf_compare)))
  expect_identical(.CF_SL_LIBRARY_DEFAULT,
                   c("mean", "ridge", "elastic_net", "lasso", "rf",
                     "xgboost_d2", "xgboost_d4"))
})

test_that("cf_evaluate takes no nuisance arguments and needs one on the model", {
  # The nuisance is fitted where the model is, so cf_evaluate reads the one the
  # model was scored with. Configuring a second one here would let the
  # comparison and the evaluation disagree.
  for (a in c("eval_nuisance", "eval_library", "eval_params",
              "eval_propensity_score")) {
    expect_false(a %in% names(formals(cf_evaluate)))
  }
  d <- sl_data()
  df <- data.frame(Y = d$y_linear, T = d$T, d$X)
  fit <- cf_fit(data = df, Y = "Y", T = "T",
                X_baseline = colnames(d$X), method = "grf",
                preprocess = FALSE, random_state = 0L)
  expect_false(is.null(fit$res$nuisance))
  stripped <- fit
  stripped$res$nuisance <- NULL
  expect_error(cf_evaluate(stripped, metrics = "blp"),
               class = "cateomics_evaluate_missing_nuisance")
})

test_that("a linear outcome puts the weight on a penalised linear member", {
  d <- sl_data()
  r <- .cf_sl_oof_quantity(d$X, d$y_linear, d$folds, lib_all(),
                           type = "regression", random_state = 0L,
                           member_params = no_params(lib_all()))
  expect_named(r, c("pred", "weights"))
  expect_equal(sum(r$weights), 1, tolerance = 1e-8)
  expect_true(all(r$weights >= 0))
  linear_members <- c("ridge", "elastic_net", "lasso")
  # Equal weights would give these three 3/7 = 0.43 together, and the fit
  # would be far worse than this.
  expect_gt(sum(r$weights[linear_members]), 0.9)
  r2 <- 1 - mean((d$y_linear - r$pred)^2) / stats::var(d$y_linear)
  expect_gt(r2, 0.9)
})

test_that("a randomised treatment keeps the propensity near the marginal rate", {
  d <- sl_data()
  r <- .cf_sl_oof_quantity(d$X, d$T, d$folds, lib_all(),
                           type = "classification", random_state = 0L,
                           member_params = no_params(lib_all()))
  # T is independent of X here, so the honest answer is the marginal rate and
  # the combination has to stay there even though two members (measured: the
  # boosted ones, up to 0.37 away) wander.
  expect_true(all(r$pred > 0 & r$pred < 1))
  expect_lt(max(abs(r$pred - mean(d$T))), 0.2)
  # The weight allocation is deliberately not asserted: with no signal the
  # penalised classifiers collapse onto the intercept, so their columns agree
  # with the `mean` column to three decimals and the non-negative least
  # squares solution over collinear columns is not unique.
})

test_that("a noise outcome with a level puts the weight on the intercept", {
  d <- sl_data()
  r <- .cf_sl_oof_quantity(d$X, d$y_noise, d$folds, lib_all(),
                           type = "regression", random_state = 0L,
                           member_params = no_params(lib_all()))
  expect_gt(r$weights[["mean"]], 0.5)
})

test_that("restricting the library restricts the weights", {
  d <- sl_data()
  small <- c("mean", "ridge")
  r <- .cf_sl_oof_quantity(d$X, d$y_linear, d$folds, small,
                           type = "regression", random_state = 0L,
                           member_params = no_params(small))
  expect_named(r$weights, small)
  expect_equal(sum(r$weights), 1, tolerance = 1e-8)
  expect_gt(r$weights[["ridge"]], 0.9)
})

test_that("the default seed reproduces the weights and the predictions", {
  d <- sl_data()
  # The ambient RNG is moved between the two calls: CATEomics restores
  # `.Random.seed`, so two back-to-back calls would agree even if the seed were
  # not honoured.
  set.seed(111)
  a <- .cf_sl_oof_quantity(d$X, d$y_linear, d$folds, lib_all(),
                           type = "regression", random_state = 0L,
                           member_params = no_params(lib_all()))
  set.seed(222)
  b <- .cf_sl_oof_quantity(d$X, d$y_linear, d$folds, lib_all(),
                           type = "regression", random_state = 0L,
                           member_params = no_params(lib_all()))
  expect_identical(a$weights, b$weights)
  expect_equal(a$pred, b$pred, tolerance = 0)
})

test_that("cf_fit records the library and the weights, and keeps the nuisance contract", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  fit <- cf_fit(data = df, Y = "y", T = "t",
                X_baseline = colnames(d$X), method = "grf",
                preprocess = FALSE, random_state = 0L)
  est <- fit$meta$estimator
  # The request is not recorded beside the resolved pair: the pair carries it
  # (both slots name the learner that ran) and says more, because a supplied
  # propensity turns that slot into "known". Read it with [[ ]]: `$` would
  # partial-match `eval_nuisance_resolved` and return the pair instead of
  # failing.
  expect_null(est[["eval_nuisance"]])
  expect_identical(est$eval_nuisance_resolved$outcome, "super_learner")
  expect_identical(est$eval_nuisance_resolved$propensity, "super_learner")
  expect_identical(est$eval_library$outcome, .CF_SL_LIBRARY_DEFAULT)
  expect_setequal(names(est$eval_sl_weights), c("pi", "mu0", "mu1", "m"))
  for (q in names(est$eval_sl_weights)) {
    expect_equal(sum(est$eval_sl_weights[[q]]), 1, tolerance = 1e-8)
  }
  # The weights describe the nuisance but are not one of its vectors.
  expect_false("sl_weights" %in% names(fit$res$nuisance))
  expect_setequal(names(fit$res$nuisance),
                  c("mu0_oof", "mu1_oof", "pi_oof", "m_oof"))
})

test_that("partition = 'none' still solves weights, from an internal split", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  fit <- cf_fit(data = df, Y = "y", T = "t",
                X_baseline = colnames(d$X), method = "grf",
                preprocess = FALSE, random_state = 0L, partition = "none")
  w <- fit$meta$estimator$eval_sl_weights
  expect_true(length(w) > 0L)
  expect_equal(sum(w$mu0), 1, tolerance = 1e-8)
  expect_false("sl_weights" %in% names(fit$res$nuisance))
})

test_that("a single learner is unaffected and records no library", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  fit <- cf_fit(data = df, Y = "y", T = "t",
                X_baseline = colnames(d$X), method = "grf",
                preprocess = FALSE, random_state = 0L,
                eval_nuisance = "rf", eval_library = c("mean", "ridge"))
  expect_identical(fit$meta$estimator$eval_nuisance_resolved$outcome, "rf")
  expect_identical(fit$meta$estimator$eval_nuisance_resolved$propensity, "rf")
  expect_null(fit$meta$estimator$eval_library)
  expect_null(fit$meta$estimator$eval_sl_weights)
})

test_that("eval_library is validated", {
  expect_identical(.cf_sl_resolve_library(NULL, "outcome"), .CF_SL_LIBRARY_DEFAULT)
  expect_identical(.cf_sl_resolve_library(c("mean", "rf"), "outcome"),
                   c("mean", "rf"))
  expect_identical(
    .cf_sl_resolve_library(list(outcome = c("mean", "rf")), "propensity"),
    c("mean", "rf"))   # the missing slot falls back to its partner
  expect_error(.cf_sl_resolve_library(c("mean", "nope"), "outcome"),
               class = "cateomics_invalid_eval_library")
  expect_error(.cf_sl_resolve_library(c("mean", "mean"), "outcome"),
               class = "cateomics_invalid_eval_library")
  expect_error(.cf_sl_resolve_library(character(0), "outcome"),
               class = "cateomics_invalid_eval_library")
})

test_that("under the Super Learner eval_params is keyed by member", {
  lib <- lib_all()
  mem <- function(ep) .cf_eval_params_members(ep, "outcome", lib, type = "regression")

  expect_identical(mem(list(rf = list(num.trees = 1000L)))$rf,
                   list(num.trees = 1000L))
  expect_identical(mem(list(outcome = list(rf = list(num.trees = 1000L))))$rf,
                   list(num.trees = 1000L))
  two <- mem(list(outcome = list(rf = list(num.trees = 500L),
                                 xgboost_d2 = list(min_child_weight = 5L))))
  expect_identical(two$rf, list(num.trees = 500L))
  expect_identical(two$xgboost_d2, list(min_child_weight = 5L))
  expect_identical(two$ridge, list())

  # A bare hyperparameter has no single learner to be checked against.
  expect_error(mem(list(outcome = list(max_depth = 3L))),
               class = "cateomics_sl_params_need_member")
  expect_error(mem(list(outcome = list(xgboost_d8 = list(nrounds = 10L)))),
               class = "cateomics_sl_params_need_member")
  expect_error(.cf_eval_params_members(list(outcome = list(rf = list())),
                                       "outcome", c("mean", "ridge"),
                                       type = "regression"),
               class = "cateomics_sl_params_need_member")
  # The member label is the name of the setting it carries.
  expect_error(mem(list(outcome = list(xgboost_d2 = list(max_depth = 5L)))),
               class = "cateomics_duplicate_param")
  expect_error(mem(list(outcome = list(mean = list(nrounds = 10L)))),
               class = "cateomics_unknown_param")
  expect_error(mem(list(outcome = list(rf = list(nope = 1L)))),
               class = "cateomics_unknown_param")
  expect_error(mem(list(outcome = list(rf = list(seed = 7L)))),
               class = "cateomics_rng_param")
})

test_that("a single learner refuses a learner block, which is the mirror mistake", {
  expect_error(
    .cf_eval_params_slot(list(rf = list(num.trees = 10L)), "outcome", "rf",
                         type = "regression"),
    class = "cateomics_eval_params_learner_block")
  # and still accepts the hyperparameters written directly
  expect_identical(
    .cf_eval_params_slot(list(num.trees = 10L), "outcome", "rf",
                         type = "regression"),
    list(num.trees = 10L))
})

test_that("member hyperparameters reach the fit", {
  d <- sl_data()
  big <- .cf_sl_fit_member(d$X, d$y_linear, type = "regression", member = "rf",
                           random_state = 0L, params = list(num.trees = 23L))
  expect_equal(big$model$num.trees, 23)   # ranger stores it as a double
  shallow <- .cf_sl_fit_member(d$X, d$y_linear, type = "regression",
                               member = "xgboost_d2", random_state = 0L)
  deep <- .cf_sl_fit_member(d$X, d$y_linear, type = "regression",
                            member = "xgboost_d4", random_state = 0L)
  # The two boosted members differ only in capacity, so their predictions must
  # not coincide.
  expect_gt(max(abs(.cf_sl_predict_member(shallow, d$X) -
                    .cf_sl_predict_member(deep, d$X))), 1e-6)
})

test_that("the intercept member predicts the training mean", {
  d <- sl_data()
  m <- .cf_sl_fit_member(d$X, d$y_linear, type = "regression", member = "mean")
  p <- .cf_sl_predict_member(m, d$X[1:5, , drop = FALSE])
  expect_equal(p, rep(mean(d$y_linear), 5L))
})

# --- a member that cannot fit, and the group check ---

sl_fail_params <- function(members, bad = "rf") {
  p <- no_params(members)
  p[[bad]] <- list(num.trees = 0L)   # ranger refuses this at fit time
  p
}

test_that("a member that cannot fit is dropped, warned about, and re-weighted", {
  d <- sl_data()
  expect_warning(
    r <- .cf_sl_oof_quantity(d$X, d$y_linear, d$folds, lib_all(),
                             type = "regression", random_state = 0L,
                             member_params = sl_fail_params(lib_all())),
    class = "cateomics_sl_member_failed")
  expect_equal(r$weights[["rf"]], 0)
  expect_equal(sum(r$weights), 1, tolerance = 1e-8)
  expect_true(all(is.finite(r$pred)))
})

test_that("every member failing is an error, not a silent zero", {
  d <- sl_data()
  p <- no_params(c("rf"))
  p$rf <- list(num.trees = 0L)
  expect_error(
    suppressWarnings(.cf_sl_oof_quantity(d$X, d$y_linear, d$folds, "rf",
                                         type = "regression", random_state = 0L,
                                         member_params = p)),
    class = "cateomics_sl_all_members_failed")
})

test_that("a propensity-only library reaches the propensity slot", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  fit <- cf_fit(data = df, Y = "y", T = "t", X_baseline = colnames(d$X),
                method = "grf", preprocess = FALSE, random_state = 0L,
                eval_library = list(outcome = c("mean", "ridge"),
                                    propensity = c("mean", "lasso")))
  lib <- fit$meta$estimator$eval_library
  expect_identical(lib$outcome, c("mean", "ridge"))
  expect_identical(lib$propensity, c("mean", "lasso"))
  expect_setequal(names(fit$meta$estimator$eval_sl_weights$pi),
                  c("mean", "lasso"))
  expect_setequal(names(fit$meta$estimator$eval_sl_weights$mu0),
                  c("mean", "ridge"))
})

test_that("a binary outcome makes the outcome members probability learners", {
  d <- sl_data()
  yb <- as.integer(d$y_linear > stats::median(d$y_linear))
  df <- data.frame(y = yb, t = d$T, d$X)
  fit <- suppressMessages(
    cf_fit(data = df, Y = "y", T = "t", X_baseline = colnames(d$X),
           method = "grf", preprocess = FALSE, random_state = 0L))
  expect_identical(fit$meta$estimator$outcome_type, "binary")
  nu <- fit$res$nuisance
  expect_true(all(nu$mu0_oof >= 0 & nu$mu0_oof <= 1))
  expect_true(all(nu$mu1_oof >= 0 & nu$mu1_oof <= 1))
  expect_equal(sum(fit$meta$estimator$eval_sl_weights$mu1), 1, tolerance = 1e-8)
})

test_that("cf_compare records the weights and leaves the nuisance contract alone", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  cmp <- suppressWarnings(
    cf_compare(data = df, Y = "y", T = "t", X_baseline = colnames(d$X),
               methods = "grf", preprocess = FALSE, random_state = 0L,
               verbose = FALSE))
  expect_setequal(names(cmp$metadata$eval_sl_weights), c("pi", "mu0", "mu1", "m"))
  expect_identical(cmp$metadata$eval_library$outcome, .CF_SL_LIBRARY_DEFAULT)
  expect_false("sl_weights" %in% names(cmp$metadata$nuisance_oof))
  expect_false("sl_weights" %in% names(cmp$models[[1]]$res$nuisance))
})

test_that("eval_library beside a single learner is refused rather than ignored", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  expect_warning(
    cf_fit(data = df, Y = "y", T = "t", X_baseline = colnames(d$X),
           method = "grf", preprocess = FALSE, random_state = 0L,
           eval_nuisance = "rf", eval_library = c("mean", "ridge")),
    class = "cateomics_eval_library_unused")
})

test_that("the resolved spec is what is recorded, not the request", {
  d <- sl_data()
  df <- data.frame(y = d$y_linear, t = d$T, d$X)
  fit <- cf_fit(data = df, Y = "y", T = "t", X_baseline = colnames(d$X),
                method = "grf", preprocess = FALSE, random_state = 0L,
                eval_propensity_score = 0.5)
  est <- fit$meta$estimator
  expect_null(est[["eval_nuisance"]])
  # the propensity slot actually ran as the supplied score
  expect_identical(est$eval_nuisance_resolved$propensity, "known")
  expect_identical(est$eval_nuisance_resolved$outcome, "super_learner")
  expect_null(est$eval_library$propensity)
  expect_null(est$eval_sl_weights$pi)
})


# Per-slot learners for the evaluation nuisance -----------------------------
#
# The estimator side lets the outcome and the propensity be fitted by
# different learners (`model_outcome` / `model_propensity`). The evaluation
# side has to offer the same, which it does through
# `eval_nuisance = list(outcome = , propensity = )` -- the shape
# `eval_library` and `eval_params` already take. The assertion that matters is
# not that the call is accepted but that the two slots really used different
# learners, so each test reads the fitted nuisance, not only the metadata.

ps_data <- function(n = 90L, p = 6L, seed = 5L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  T <- stats::rbinom(n, 1L, 0.5)
  Y <- X[, 2] + T * (1 + X[, 1]) + stats::rnorm(n)
  list(df = data.frame(Y = Y, T = T, X), cv = colnames(X))
}

test_that("cf_fit fits the two nuisance slots with the learners named per slot", {
  d <- ps_data()
  f <- function(...) cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
                            method = "grf", preprocess = FALSE, ...)
  one  <- f(eval_nuisance = "ridge")
  split <- f(eval_nuisance = list(outcome = "ridge", propensity = "rf"))
  expect_identical(unlist(one$meta$estimator$eval_nuisance_resolved,
                          use.names = FALSE), c("ridge", "ridge"))
  expect_identical(unlist(split$meta$estimator$eval_nuisance_resolved,
                          use.names = FALSE), c("ridge", "rf"))
  # Naming a different propensity learner has to move the propensity and leave
  # the outcome alone; metadata alone would not show that.
  expect_false(isTRUE(all.equal(one$res$nuisance$pi_oof,
                                split$res$nuisance$pi_oof)))
  expect_equal(one$res$nuisance$mu0_oof, split$res$nuisance$mu0_oof)
  expect_equal(one$res$nuisance$mu1_oof, split$res$nuisance$mu1_oof)
  # One slot named applies to both, which is the resolver's broadcast.
  expect_identical(unlist(f(eval_nuisance = list(outcome = "lasso"))$meta$
                            estimator$eval_nuisance_resolved,
                          use.names = FALSE), c("lasso", "lasso"))
})

test_that("cf_compare takes the per-slot form the same way cf_fit does", {
  d <- ps_data()
  g <- function(...) cf_compare(data = d$df, Y = "Y", T = "T",
                                X_baseline = d$cv, methods = "grf",
                                preprocess = FALSE, n_folds = 3L, ...)
  one   <- g(eval_nuisance = "ridge")$models[[1L]]
  split <- g(eval_nuisance = list(outcome = "ridge",
                                  propensity = "rf"))$models[[1L]]
  expect_identical(unlist(split$meta$estimator$eval_nuisance_resolved,
                          use.names = FALSE), c("ridge", "rf"))
  expect_false(isTRUE(all.equal(one$res$nuisance$pi_oof,
                                split$res$nuisance$pi_oof)))
})

test_that("a known evaluation propensity and a named propensity learner conflict", {
  d <- ps_data()
  # Both entries refuse it, with the same condition, the way the estimator
  # side refuses `propensity_score` together with `model_propensity`.
  expect_error(
    cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv, method = "grf",
           preprocess = FALSE, eval_nuisance = list(propensity = "rf"),
           eval_propensity_score = 0.5),
    class = "cateomics_eval_propensity_args_conflict")
  expect_error(
    cf_compare(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
               methods = "grf", preprocess = FALSE, n_folds = 3L,
               eval_nuisance = list(propensity = "rf"),
               eval_propensity_score = 0.5),
    class = "cateomics_eval_propensity_args_conflict")
  # Naming only the outcome leaves the known score free to take the other slot.
  f <- cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv, method = "grf",
              preprocess = FALSE, eval_nuisance = list(outcome = "ridge"),
              eval_propensity_score = 0.5)
  expect_identical(unlist(f$meta$estimator$eval_nuisance_resolved,
                          use.names = FALSE), c("ridge", "known"))
})

test_that("the per-slot form refuses what the single name refuses", {
  d <- ps_data()
  f <- function(v) cf_fit(data = d$df, Y = "Y", T = "T", X_baseline = d$cv,
                          method = "grf", preprocess = FALSE,
                          eval_nuisance = v)
  expect_error(f(list(outcome = "nope")), class = "cateomics_invalid_eval_nuisance")
  expect_error(f(list(propensity = "nope")), class = "cateomics_invalid_eval_nuisance")
  # "known" carries no score, so it is not a learner a caller can name.
  expect_error(f(list(propensity = "known")), class = "cateomics_invalid_eval_nuisance")
  expect_error(f(list(foo = "rf")), class = "cateomics_invalid_eval_nuisance")
  expect_error(f(list()), class = "cateomics_invalid_eval_nuisance")
  expect_error(f(c("ridge", "rf")), class = "cateomics_invalid_eval_nuisance")
})
