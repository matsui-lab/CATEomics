# Smoke tests for the cf_interpret() surrogate-SHAP path.
# Goals (smoke only): both surrogates run, return signed SHAP features +
# fidelity, accept ensemble inputs.

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
      data    = list(X = fx$X, Y = rnorm(n), T = rbinom(n, 1L, 0.5)),
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
    cf_interpret(m, top_n = 5L, xgb_nrounds = 30L)
  ))
  expect_s3_class(out, "cfomics_interpret")
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
    expect_equal(nrow(out$vi[[s]]$shap), nrow(out$vi[[s]]$X))
    expect_true(is.numeric(out$vi[[s]]$fidelity))
    expect_true(out$vi[[s]]$fidelity >= -1 &&
                  out$vi[[s]]$fidelity <= 1.01)  # numeric slack
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
                 class = "cfomics_missing_dependency")
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

test_that("cf_interpret accepts cfomics_ensemble input", {
  skip_if_not_installed("glmnet")
  fx <- .make_cate_fixture(n = 80L, p = 5L)
  ens <- structure(
    list(
      method  = "ensemble",
      res     = list(tau_hat_oof = fx$tau, tau_hat = fx$tau),
      data    = list(X = fx$X, Y = rnorm(80), T = rbinom(80, 1L, 0.5)),
      meta    = list(partition = "oof")
    ),
    class = c("cfomics_ensemble", "cf_model")
  )
  out <- suppressMessages(suppressWarnings(
    cf_interpret(ens, surrogate = "ridge", top_n = 5L)
  ))
  expect_equal(out$metadata$method, "ensemble")
})


# Missing X / tau --------------------------------------------------------

test_that("cf_interpret errors when model$data$X is absent", {
  m <- .make_mock_cf_model()
  m$data$X <- NULL
  expect_error(
    cf_interpret(m),
    class = "cfomics_interpret_no_X"
  )
})


test_that("cf_interpret errors when no usable tau exists", {
  m <- .make_mock_cf_model()
  m$res$tau_hat_oof <- rep(NA_real_, length(m$res$tau_hat_oof))
  m$res$tau_hat <- rep(NA_real_, length(m$res$tau_hat_oof))
  expect_error(
    cf_interpret(m),
    class = "cfomics_interpret_no_tau"
  )
})


# print / summary smoke --------------------------------------------------

test_that("print and summary methods run on the interpret object", {
  skip_if_not_installed("glmnet")
  m <- .make_mock_cf_model()
  out <- suppressMessages(suppressWarnings(
    cf_interpret(m, surrogate = "ridge", top_n = 5L)
  ))
  expect_output(print(out), "cfomics_interpret")
  s <- summary(out)
  expect_s3_class(s, "summary.cfomics_interpret")
  expect_output(print(s), "cf_interpret summary")
})
