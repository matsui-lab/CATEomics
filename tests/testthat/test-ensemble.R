# Smoke tests for cf_ensemble.
# Goals (smoke only):
#   - 5 methods (mean / softmax / r_stacking / dr_stacking / cba) run
#   - output contract (class, oof_ite, weights, fit$res${ite,ate,...})
#   - input contract (cfomics_compare required, cf_model rejected)
#   - downstream connectivity (cf_evaluate / cf_benefit_groups accept it)
# No PEHE / accuracy checks.

set.seed(789)
.make_mock_compare <- function(n = 40L, K = 4L,
                                include_models = TRUE,
                                include_nuisance = TRUE,
                                include_oof = TRUE) {
  Y <- rnorm(n, sd = 1)
  T_vec <- rbinom(n, 1L, prob = 0.5)
  mu0 <- rnorm(n, sd = 0.5)
  mu1 <- mu0 + 0.5 + rnorm(n, sd = 0.2)
  pi  <- pmin(pmax(plogis(rnorm(n, sd = 0.3)), 1e-4), 1 - 1e-4)

  oof_tau_matrix <- matrix(rnorm(n * K, sd = 0.5), nrow = n, ncol = K)
  colnames(oof_tau_matrix) <- paste0("m", seq_len(K))

  scoreboard <- data.frame(
    method       = colnames(oof_tau_matrix),
    r_loss       = runif(K, 0.5, 1.5),
    dr_loss      = runif(K, 0.5, 1.5),
    factual_loss = runif(K, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  X_mat <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  colnames(X_mat) <- paste0("X", seq_len(5))

  models <- if (include_models) {
    out <- list()
    for (k in seq_len(K)) {
      mdl <- list(data = list(Y = Y, T = T_vec, X = X_mat),
                  meta = list(covariate_names = colnames(X_mat),
                              n = n, p = 5L))
      out[[paste0("m", k)]] <- mdl
    }
    out
  } else {
    list()
  }

  metadata <- list(
    n_folds        = 5L,
    nuisance_oof   = if (include_nuisance) list(mu0_oof = mu0, mu1_oof = mu1,
                                                 pi_oof = pi,
                                                 m_oof = rnorm(n)) else NULL,
    oof_tau_matrix = if (include_oof) oof_tau_matrix else NULL,
    parsed         = list(covariate_names = colnames(X_mat)),
    formula        = NULL,
    random_state   = NULL
  )

  structure(
    list(scoreboard = scoreboard,
         models     = models,
         metadata   = metadata),
    class = c("cfomics_compare", "list")
  )
}


# 5 method smoke runs --------------------------------------------------------

test_that("cf_ensemble runs for all 5 methods and returns proper class", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)

  for (m in c("mean", "softmax", "r_stacking", "dr_stacking", "cba")) {
    out <- expect_no_error(
      suppressMessages(cf_ensemble(cmp, method = m))
    )
    expect_s3_class(out, "cfomics_ensemble")
    expect_s3_class(out, "cf_model")
    expect_s3_class(out, "cf_model")
    expect_equal(out$method, "ensemble")
    expect_equal(out$meta$ensemble_method, m)
    expect_equal(length(out$res$tau_hat_oof), 40L)
    expect_equal(length(out$weights), 4L)
    # Weights are non-negative. r_stacking is NNLS and is NOT renormalised to
    # sum 1; the convex methods (mean/softmax/dr_stacking) sum to 1, and cba
    # surfaces consensus_members instead of a weight simplex.
    expect_true(all(out$weights >= -1e-8))
    if (!m %in% c("r_stacking")) {
      expect_true(abs(sum(out$weights) - 1) < 1e-6 ||
                  length(out$meta$consensus_members) > 0L)
    }
  }
})


# Output contract -----------------------------------------------------------

test_that("cf_ensemble output has the minimum contract for downstream", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  out <- suppressMessages(cf_ensemble(cmp, method = "mean"))

  # Required top-level fields (four-pillar).
  expect_true(all(c("method", "res", "model", "data", "meta",
                    "weights") %in% names(out)))
  # data triplet.
  expect_true(all(c("Y", "T", "X") %in% names(out$data)))
  # res with tau_hat / tau_hat_oof / mu0_hat / mu1_hat / ate / summary.
  res <- out$res
  expect_true(all(c("tau_hat", "tau_hat_oof", "ate", "mu0_hat", "mu1_hat",
                    "summary") %in% names(res)))
  expect_equal(length(res$tau_hat), 30L)
  expect_equal(length(res$mu0_hat), 30L)
  expect_equal(length(res$mu1_hat), 30L)
  expect_true(is.numeric(res$ate))
  # Summary fields (IQR form).
  expect_true(all(c("tau_mean", "tau_sd", "tau_iqr") %in%
                    names(res$summary)))
})


# Input contract ------------------------------------------------------------

test_that("cf_ensemble rejects non-cfomics_compare inputs", {
  expect_error(
    cf_ensemble(list(scoreboard = NULL), method = "mean"),
    class = "cfomics_invalid_compare_result"
  )
  cf_model_like <- structure(list(method = "grf"), class = "cf_model")
  expect_error(
    cf_ensemble(cf_model_like, method = "mean"),
    class = "cfomics_ensemble_single_model"
  )
})

test_that("cf_ensemble errors when oof_tau_matrix is unusable", {
  cmp <- .make_mock_compare(include_oof = FALSE)
  expect_error(
    suppressMessages(cf_ensemble(cmp, method = "mean")),
    class = "cfomics_ensemble_empty_input"
  )
})


# softmax dispatch with eta -------------------------------------------------

test_that("softmax eta = 0 reduces to mean, large eta concentrates", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)

  out_eq <- suppressMessages(cf_ensemble(cmp, method = "softmax", eta = 0))
  expect_true(all(abs(out_eq$weights - 1 / 3) < 1e-8))

  out_hi <- suppressMessages(cf_ensemble(cmp, method = "softmax", eta = 1e6))
  # The single highest score (lowest loss) should dominate.
  w <- out_hi$weights
  expect_true(max(w) > 0.95)
})


# stacking weights satisfy the simplex constraint --------------------------

test_that("r_stacking weights satisfy w >= 0 (NNLS free scale)", {
  cmp <- .make_mock_compare(n = 50L, K = 3L)
  out <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))
  expect_true(all(out$weights >= -1e-8))
  # r_stacking is non-negative least squares; sum(w) is NOT constrained to 1.
  # Only require a non-degenerate positive total.
  expect_true(sum(out$weights) > 0)
})

test_that("dr_stacking weights satisfy w >= 0 and sum(w) = 1", {
  cmp <- .make_mock_compare(n = 50L, K = 3L)
  out <- suppressMessages(cf_ensemble(cmp, method = "dr_stacking"))
  expect_true(all(out$weights >= -1e-8))
  expect_true(abs(sum(out$weights) - 1) < 1e-6)
})


# CBA selects a consensus subset --------------------------------------------

test_that("cba returns a consensus_members subset and uniform weights over it", {
  cmp <- .make_mock_compare(n = 50L, K = 4L)
  out <- suppressMessages(cf_ensemble(cmp, method = "cba"))
  expect_true(length(out$meta$consensus_members) >= 2L)
  expect_true(length(out$meta$consensus_members) <= 4L)
  # Weights are uniform across selected members.
  picked <- out$weights[out$meta$consensus_members]
  expect_true(all(abs(picked - picked[1]) < 1e-12))
  # Non-selected methods have weight 0.
  others <- setdiff(names(out$weights), out$meta$consensus_members)
  if (length(others) > 0L) {
    expect_true(all(out$weights[others] == 0))
  }
})


# Failed methods are dropped silently with an info note ---------------------

test_that("cf_ensemble drops all-NA columns via OOF-matrix cleaning", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  # Mark m1 as a "failed" method (all NA).
  cmp$metadata$oof_tau_matrix[, "m1"] <- NA_real_

  out <- suppressMessages(cf_ensemble(cmp, method = "mean"))
  expect_true("m1" %in% out$meta$methods_dropped)
  expect_equal(length(out$meta$methods_used), 3L)
  expect_equal(out$weights[["m1"]], 0)
})


# print / summary smoke -----------------------------------------------------

test_that("print and summary methods do not error", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  out <- suppressMessages(cf_ensemble(cmp, method = "softmax", eta = 1.0))
  expect_output(print(out), "<cfomics_ensemble>")
  s <- summary(out)
  expect_s3_class(s, "summary.cfomics_ensemble")
  expect_output(print(s), "cf_ensemble summary")
})


# predict() on an ensemble --------------------------------------------------

test_that("predict() returns the stored ensemble tau (no newdata)", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))

  tau <- predict(ens, type = "tau")
  expect_type(tau, "double")
  expect_length(tau, 40L)
  expect_identical(tau, ens$res$tau_hat)

  all <- predict(ens, type = "all")
  expect_true(is.list(all))
  expect_identical(all$tau_hat, ens$res$tau_hat)
  expect_identical(all$ate, ens$res$ate)
})

test_that("predict() on an ensemble rejects newdata", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "mean"))
  nd  <- matrix(rnorm(40L * 5L), nrow = 40L, ncol = 5L)
  expect_error(
    predict(ens, newdata = nd),
    class = "cfomics_ensemble_no_newdata"
  )
})

test_that("cf_plot_cate() renders an ensemble CATE distribution", {
  skip_if_not_installed("ggplot2")
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))
  p   <- cf_plot_cate(ens)
  expect_s3_class(p, "ggplot")
})
