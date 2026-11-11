# Smoke tests for the cf_interpret() surrogate-SHAP path.
# Goals (smoke only): both surrogates run, return signed SHAP features +
# fidelity, accept ensemble inputs. Fidelity is held out (see below).

set.seed(99)
.make_cate_fixture <- function(n = 80L, p = 6L) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  # Build a tau hat that responds to X[, 1] and X[, 2] mostly.
  tau <- 0.5 * X[, 1] - 0.4 * X[, 2] + 0.1 * rnorm(n)
  list(X = X, tau = tau)
}

.make_mock_cf_model <- function(method = "grf", n = 80L, p = 6L) {
  fx <- .make_cate_fixture(n, p)
  structure(
    list(
      method  = method,
      res     = list(tau_hat_oof = fx$tau, tau_hat = fx$tau),
      data    = list(X_preprocessed = fx$X, Y = rnorm(n),
                     T = rbinom(n, 1L, 0.5)),
      meta    = list(partition = "oof")
    ),
    class = "cf_model"
  )
}


# Ridge + XGB run + output contract --------------------------------------

test_that("cf_interpret returns vi$ridge and vi$xgb with the standard fields", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("xgboost")
  m <- .make_mock_cf_model()
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, top_n = 5L, params = list(xgb = list(nrounds = 30L)))
  ))
  expect_s3_class(out, "cateomics_interpret")
  expect_setequal(names(out$vi), c("ridge", "xgb"))

  for (s in c("ridge", "xgb")) {
    feats <- out$vi[[s]]$features
    expect_s3_class(feats, "data.frame")
    # features table is feature_name / importance / direction / rank; there is
    # no signed-mean shap_value column.
    expect_true(all(c("feature_name", "importance", "direction", "rank") %in%
                      names(feats)))
    expect_false("shap_value" %in% names(feats))
    expect_true(all(feats$importance >= 0))          # mean|SHAP| is unsigned
    expect_true(all(feats$direction %in% c(-1, 0, 1)))
    # The full per-observation SHAP matrix + feature values are retained for
    # as_shapviz().
    expect_true(is.matrix(out$vi[[s]]$shap))
    # One design matrix for the object, not one per surrogate.
    expect_equal(nrow(out$vi[[s]]$shap), nrow(out$X))
    expect_null(out$vi[[s]][["X"]])
    expect_null(out$vi[[s]][["coef"]])
    expect_null(out$vi[[s]][["intercept"]])
    # Held-out R^2: bounded above by 1, unbounded below (a surrogate can do
    # worse than the mean on rows it did not see).
    expect_true(is.numeric(out$vi[[s]]$fidelity))
    expect_true(is.finite(out$vi[[s]]$fidelity))
    expect_lte(out$vi[[s]]$fidelity, 1.01)
  }
})

test_that("as_shapviz errors informatively without shapviz, works with it", {
  skip_if_not_installed("glmnet")
  m <- .make_mock_cf_model()
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, surrogate = "ridge", top_n = 5L)))
  if (requireNamespace("shapviz", quietly = TRUE)) {
    sv <- as_shapviz(out, surrogate = "ridge")
    expect_s3_class(sv, "shapviz")
  } else {
    expect_error(as_shapviz(out, surrogate = "ridge"),
                 class = "cateomics_optional_dependency_missing")
  }
})


test_that("cf_interpret picks up the strongest tau driver via Ridge sign", {
  skip_if_not_installed("glmnet")
  m <- .make_mock_cf_model(n = 200L, p = 5L)
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, surrogate = "ridge", top_n = 5L)
  ))
  # X1 has the largest positive driver; expect X1 in top-3 by abs SHAP.
  top3 <- out$vi$ridge$features$feature_name[1:3]
  expect_true("X1" %in% top3)
})


# Ensemble input ---------------------------------------------------------

test_that("cf_interpret accepts cateomics_ensemble input", {
  skip_if_not_installed("glmnet")
  fx <- .make_cate_fixture(n = 80L, p = 5L)
  ens <- structure(
    list(
      method  = "ensemble",
      res     = list(tau_hat_oof = fx$tau, tau_hat = fx$tau),
      data    = list(X = fx$X, Y = rnorm(80), T = rbinom(80, 1L, 0.5)),
      meta    = list(partition = "oof")
    ),
    class = c("cateomics_ensemble", "cf_model")
  )
  out <- suppressMessages(suppressWarnings(
    cf_interpret(ens, surrogate = "ridge", top_n = 5L)
  ))
  expect_equal(out$metadata$method, "ensemble")
})


# Missing X / tau --------------------------------------------------------

test_that("cf_interpret errors when the cached design matrix is absent", {
  m <- .make_mock_cf_model()
  m$data$X_preprocessed <- NULL
  expect_error(
    cf_interpret(m),
    class = "cateomics_interpret_no_X"
  )
})


test_that("cf_interpret errors when no usable tau exists", {
  m <- .make_mock_cf_model()
  m$res$tau_hat_oof <- rep(NA_real_, length(m$res$tau_hat_oof))
  m$res$tau_hat <- rep(NA_real_, length(m$res$tau_hat_oof))
  expect_error(
    cf_interpret(m),
    class = "cateomics_interpret_no_tau"
  )
})


# print / summary smoke --------------------------------------------------

test_that("print and summary methods run on the interpret object", {
  skip_if_not_installed("glmnet")
  m <- .make_mock_cf_model()
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, surrogate = "ridge", top_n = 5L)
  ))
  expect_output(print(out), "cateomics_interpret")
  s <- summary(out)
  expect_s3_class(s, "summary.cateomics_interpret")
  expect_output(print(s), "cf_interpret summary")
})


test_that("fidelity is held out, so a surrogate cannot memorise its way to 1", {
  skip_if_not_installed("xgboost")
  # tau hat is pure noise: nothing in X predicts it, so an honest fidelity is
  # around zero even though a deep tree can interpolate the training rows.
  set.seed(4)
  n <- 80L; p <- 6L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  m <- structure(
    list(method = "grf", res = list(tau_hat_oof = rnorm(n)),
         data = list(X = X, Y = rnorm(n), T = rbinom(n, 1L, 0.5)),
         meta = list(partition = "oof")),
    class = "cf_model")
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, surrogate = "xgb", top_n = 3L, random_state = 1L,
                 params = list(xgb = list(nrounds = 200L, max_depth = 6L)))))
  expect_lt(out$vi$xgb$fidelity, 0.3)

  # In-sample R^2 of the same surrogate on the same rows, for contrast.
  bst <- xgboost::xgb.train(
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6L,
                  verbosity = 0L),
    data = xgboost::xgb.DMatrix(data = out$X,
                                label = m$res$tau_hat_oof),
    nrounds = 200L, verbose = 0L)
  in_sample <- CATEomics:::.cateomics_r_squared(
    m$res$tau_hat_oof,
    as.numeric(stats::predict(bst, newdata = out$X)))
  expect_gt(in_sample, 0.9)
})


# Fidelity is bounded from below as well as from above --------------------

test_that("a surrogate that genuinely fits scores high, one that cannot scores low", {
  skip_if_not_installed("glmnet")
  set.seed(31)
  n <- 200L; p <- 5L
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("g", seq_len(p))))
  mk <- function(tau) structure(
    list(method = "grf", res = list(tau_hat_oof = tau, tau_hat = tau),
         data = list(X_preprocessed = X, Y = rnorm(n),
                     T = rbinom(n, 1L, 0.5)),
         meta = list(partition = "oof")), class = "cf_model")

  # tau is an exact linear function of X, so a ridge surrogate must recover
  # it out of fold. Only bounding fidelity from above would let a surrogate
  # that always reported zero pass.
  linear <- suppressMessages(suppressWarnings(
    cf_interpret(mk(as.numeric(X %*% c(2, -1, 0.5, 0, 0))),
                 surrogate = "ridge", top_n = 3L)))
  expect_gt(linear$vi$ridge$fidelity, 0.9)

  noise <- suppressMessages(suppressWarnings(
    cf_interpret(mk(rnorm(n)), surrogate = "ridge", top_n = 3L)))
  expect_lt(noise$vi$ridge$fidelity, 0.3)
})

test_that("surrogate hyperparameters reach the fit and are checked", {
  skip_if_not_installed("grf")
  skip_if_not_installed("xgboost")
  set.seed(42)
  n <- 150L; p <- 10L
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("g", 1:p)))
  T_vec <- rbinom(n, 1L, 0.5)
  Y <- X[, 1] + T_vec * (1 + X[, 2]) + rnorm(n)
  df <- data.frame(Y = Y, T = T_vec, X)
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = colnames(X), method = "grf",
    partition = "oof", random_state = 1L, preprocess = FALSE))

  base <- suppressMessages(cf_interpret(fit, surrogate = c("ridge", "xgb"),
                                        random_state = 11L))
  # The resolved values are recorded, so a figure or a number can be traced
  # back to the setting that produced it.
  expect_equal(base$vi$ridge$params$nfolds, 5L)
  expect_equal(base$vi$xgb$params$nrounds, 100L)
  expect_equal(base$vi$xgb$params$max_depth, 4L)

  # Reaching the fit, not just being stored: a much smaller forest has to
  # change the held-out fidelity.
  small <- suppressMessages(cf_interpret(
    fit, surrogate = "xgb", random_state = 11L,
    params = list(xgb = list(nrounds = 20L, max_depth = 2L, eta = 0.3))))
  expect_equal(small$vi$xgb$params$nrounds, 20L)
  expect_false(isTRUE(all.equal(small$vi$xgb$fidelity, base$vi$xgb$fidelity)))

  # Names come from the R function that runs the surrogate, derived from its
  # formals rather than a list kept here, so a forwarded argument works...
  expect_no_error(suppressMessages(cf_interpret(
    fit, surrogate = "ridge", random_state = 11L,
    params = list(ridge = list(nfolds = 10L, penalty.factor = rep(1, p))))))
  # ...a typo does not, and neither do the values CATEomics fixes itself.
  expect_error(cf_interpret(fit, surrogate = "xgb",
                            params = list(xgb = list(not_a_param = 1))),
               class = "cateomics_unknown_param")
  expect_error(cf_interpret(fit, surrogate = "xgb",
                            params = list(xgb = list(objective = "x"))),
               class = "cateomics_unknown_param")
  expect_error(cf_interpret(fit, surrogate = "ridge",
                            params = list(ridge = list(alpha = 1))),
               class = "cateomics_unknown_param")
  expect_error(cf_interpret(fit, surrogate = "ridge",
                            params = list(tree = list(a = 1))),
               class = "cateomics_unknown_params_slot")
  expect_error(cf_interpret(fit, surrogate = "ridge",
                            params = list(xgb = list(nrounds = 5L))),
               class = "cateomics_unused_params")
})

test_that("the surrogates come back in the order they were asked for", {
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 15, random_state = 0L,
                propensity_score = 0.5, eval_propensity_score = 0.5)
  # cf_plot_interpret() and as_shapviz() take the first surrogate present when
  # none is named, so the order decides which panel a plot shows.
  both <- suppressWarnings(cf_interpret(fit, surrogate = c("xgb", "ridge"),
                                        random_state = 0L))
  expect_identical(names(both$vi), c("xgb", "ridge"))
  rev <- suppressWarnings(cf_interpret(fit, surrogate = c("ridge", "xgb"),
                                       random_state = 0L))
  expect_identical(names(rev$vi), c("ridge", "xgb"))
})
