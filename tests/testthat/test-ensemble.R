# Smoke tests for cf_ensemble.
# Goals (smoke only):
#   - 5 methods (mean / softmax / r_stacking / dr_stacking / cba) run
#   - output contract (class, weights, res$tau_hat_oof, res$ate)
#   - input contract (cateomics_compare required, cf_model rejected)
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
    class = c("cateomics_compare", "list")
  )
}


# 5 method smoke runs --------------------------------------------------------

test_that("cf_ensemble runs for all 5 methods and returns proper class", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)

  for (m in c("mean", "softmax", "r_stacking", "dr_stacking", "cba")) {
    out <- expect_no_error(
      suppressMessages(cf_ensemble(cmp, method = m))
    )
    expect_s3_class(out, "cateomics_ensemble")
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
  expect_true(all(c("Y", "T", "X_preprocessed") %in% names(out$data)))
  # res with tau_hat_oof / mu0_hat / mu1_hat / ate / summary; an ensemble has
  # no backend to re-predict with, so tau_hat is NULL.
  res <- out$res
  expect_true(all(c("tau_hat_oof", "ate", "mu0_hat", "mu1_hat",
                    "summary") %in% names(res)))
  expect_null(res$tau_hat)
  expect_equal(length(res$tau_hat_oof), 30L)
  expect_equal(length(res$mu0_hat), 30L)
  expect_equal(length(res$mu1_hat), 30L)
  expect_true(is.numeric(res$ate))
  # Summary fields (IQR form).
  expect_true(all(c("tau_mean", "tau_sd", "tau_iqr") %in%
                    names(res$summary)))
})


# Input contract ------------------------------------------------------------

test_that("cf_ensemble rejects non-cateomics_compare inputs", {
  expect_error(
    cf_ensemble(list(scoreboard = NULL), method = "mean"),
    class = "cateomics_invalid_compare_result"
  )
  cf_model_like <- structure(list(method = "grf"), class = "cf_model")
  expect_error(
    cf_ensemble(cf_model_like, method = "mean"),
    class = "cateomics_ensemble_single_model"
  )
})

test_that("cf_ensemble errors when oof_tau_matrix is unusable", {
  cmp <- .make_mock_compare(include_oof = FALSE)
  expect_error(
    suppressMessages(cf_ensemble(cmp, method = "mean")),
    class = "cateomics_ensemble_empty_input"
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
  expect_output(print(out), "<cateomics_ensemble>")
  s <- summary(out)
  expect_s3_class(s, "summary.cateomics_ensemble")
  expect_output(print(s), "cf_ensemble summary")
})


# predict() on an ensemble --------------------------------------------------

test_that("predict() returns the stored ensemble tau (no newdata)", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))

  tau <- predict(ens, type = "tau")
  expect_type(tau, "double")
  expect_length(tau, 40L)
  expect_identical(tau, ens$res$tau_hat_oof)

  all <- predict(ens, type = "all")
  expect_true(is.list(all))
  expect_identical(all$tau_hat, ens$res$tau_hat_oof)
  expect_identical(all$ate, ens$res$ate)
})

test_that("predict() on an ensemble stops when no member can score newdata", {
  # The mock's "models" are bare lists with no backend, so every member fails
  # to score. One member surviving is enough to return an effect, so the abort
  # has to require that none did.
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "mean"))
  nd  <- matrix(rnorm(40L * 5L), nrow = 40L, ncol = 5L)
  expect_error(
    suppressWarnings(predict(ens, newdata = nd)),
    class = "cateomics_ensemble_newdata_failed"
  )
})

test_that("an ensemble carries the members that hold non-zero weight", {
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "mean"))
  expect_true(is.list(ens$model$members))
  expect_identical(ens$model$kind, "ensemble")
  # mean weights every method, so every model is carried.
  expect_setequal(names(ens$model$members), names(ens$weights))
  expect_true(all(ens$weights[names(ens$model$members)] > 0))
})

test_that("stacking refuses a library it cannot identify", {
  # The weights are a regression of a target on K candidate columns, so K rows
  # over K methods does not determine them. Same condition the Super Learner
  # applies to its own library.
  cmp <- .make_mock_compare(n = 4L, K = 4L)
  expect_error(
    suppressMessages(cf_ensemble(cmp, method = "r_stacking")),
    class = "cateomics_ensemble_library_too_large"
  )
  # `mean` needs no weights, so it still runs on the same input.
  expect_s3_class(suppressMessages(cf_ensemble(cmp, method = "mean")),
                  "cateomics_ensemble")
})

test_that("cf_plot_cate() renders an ensemble CATE distribution", {
  skip_if_not_installed("ggplot2")
  cmp <- .make_mock_compare(n = 40L, K = 4L)
  ens <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))
  p   <- cf_plot_cate(ens)
  expect_s3_class(p, "ggplot")
})


# Weights: the solvers are measured against a known answer -----------------
#
# The mock is built so the right weight vector is known in advance. The
# nuisance helpers are used to CONSTRUCT the inputs, not to form the
# expectation: the expectation is "the solver puts its weight on the column
# that solves its own objective exactly", which does not depend on those
# helpers being right.

.make_solvable_compare <- function(n = 200L, seed = 5L) {
  set.seed(seed)
  T_vec <- rbinom(n, 1L, 0.5)
  pi    <- rep(0.5, n)                 # keeps T - pi bounded away from zero
  m_oof <- rnorm(n)
  Yres  <- rnorm(n)
  Y     <- m_oof + Yres
  Tres  <- T_vec - pi
  mu0   <- rnorm(n, sd = 0.3)
  mu1   <- mu0 + 0.4
  Y_DR  <- CATEomics:::.compute_dr_pseudo_outcome(Y, T_vec, mu0, mu1, pi)

  # r-stacking minimises (Yres - sum_k w_k Tres tau_k)^2, so Yres / Tres is
  # the exact answer; dr-stacking minimises (sum_k w_k tau_k - Y_DR)^2.
  M <- cbind(signal_r  = Yres / Tres,
             signal_dr = Y_DR,
             noise1    = rnorm(n),
             noise2    = rnorm(n))
  X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("X", 1:3)))

  structure(list(
    scoreboard = data.frame(method = colnames(M), r_loss = runif(4),
                            dr_loss = runif(4), factual_loss = runif(4),
                            stringsAsFactors = FALSE),
    models = stats::setNames(lapply(colnames(M), function(k)
      list(data = list(Y = Y, T = T_vec, X_preprocessed = X),
           meta = list(covariate_names = colnames(X), n = n, p = 3L))),
      colnames(M)),
    metadata = list(
      n_folds = 5L,
      nuisance_oof = list(mu0_oof = mu0, mu1_oof = mu1, pi_oof = pi,
                          m_oof = m_oof),
      oof_tau_matrix = M,
      parsed = list(covariate_names = colnames(X)),
      formula = NULL, random_state = NULL)),
    class = c("cateomics_compare", "list"))
}

test_that("r_stacking puts its weight on the column that solves the R-loss", {
  cmp <- .make_solvable_compare()
  w <- suppressMessages(cf_ensemble(cmp, method = "r_stacking"))$weights
  expect_gt(w[["signal_r"]], 0.95)
  expect_lt(max(w[c("signal_dr", "noise1", "noise2")]), 0.05)
  # Equal weights would satisfy "w >= 0" just as well, so pin the contrast.
  expect_gt(w[["signal_r"]], 10 * max(w[c("noise1", "noise2")]) + 0.1)
})

test_that("dr_stacking puts its weight on the column that solves the DR loss", {
  cmp <- .make_solvable_compare()
  w <- suppressMessages(cf_ensemble(cmp, method = "dr_stacking"))$weights
  expect_gt(w[["signal_dr"]], 0.95)
  expect_lt(max(w[c("signal_r", "noise1", "noise2")]), 0.05)
  expect_equal(sum(w), 1, tolerance = 1e-6)
})

test_that("cba drops a method that agrees with none of the others", {
  set.seed(9)
  n <- 200L
  base <- rnorm(n)
  M <- cbind(a   = base + rnorm(n, sd = 0.05),
             b   = base + rnorm(n, sd = 0.05),
             c   = base + rnorm(n, sd = 0.05),
             odd = rnorm(n))
  T_vec <- rbinom(n, 1L, 0.5); pi <- rep(0.5, n)
  m_oof <- rnorm(n); Y <- m_oof + rnorm(n)
  mu0 <- rnorm(n, sd = 0.3); mu1 <- mu0 + 0.4
  X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("X", 1:3)))
  cmp <- structure(list(
    scoreboard = data.frame(method = colnames(M), r_loss = runif(4),
                            dr_loss = runif(4), factual_loss = runif(4),
                            stringsAsFactors = FALSE),
    models = stats::setNames(lapply(colnames(M), function(k)
      list(data = list(Y = Y, T = T_vec, X_preprocessed = X),
           meta = list(covariate_names = colnames(X), n = n, p = 3L))),
      colnames(M)),
    metadata = list(n_folds = 5L,
      nuisance_oof = list(mu0_oof = mu0, mu1_oof = mu1, pi_oof = pi,
                          m_oof = m_oof),
      oof_tau_matrix = M, parsed = list(covariate_names = colnames(X)),
      formula = NULL, random_state = NULL)),
    class = c("cateomics_compare", "list"))

  ens <- suppressMessages(cf_ensemble(cmp, method = "cba"))
  # Keeping every method and weighting them equally would pass a
  # "weights are non-negative and the members are a subset" check.
  expect_false("odd" %in% ens$meta$consensus_members)
  expect_setequal(ens$meta$consensus_members, c("a", "b", "c"))
  expect_equal(unname(ens$weights[["odd"]]), 0)
})

test_that("the ensemble records how many rows stood in for fewer members", {
  # A row whose out-of-fold column is missing for some member is combined from
  # the members that are there. That substitution is not recoverable from the
  # returned tau, so the count is part of the result.
  cmp <- .make_mock_compare(n = 60L, K = 3L)
  ens <- suppressWarnings(suppressMessages(cf_ensemble(cmp, method = "mean")))
  expect_identical(ens$meta$n_rows_rescaled, 0L)

  # Punch holes into one member's column: those rows are then combined from the
  # other two, at the intended total weight, and the member is still in use.
  holes <- c(2L, 7L, 19L)
  cmp2 <- cmp
  cmp2$metadata$oof_tau_matrix[holes, 2L] <- NA_real_
  ens2 <- suppressWarnings(suppressMessages(cf_ensemble(cmp2, method = "mean")))
  expect_identical(ens2$meta$n_rows_rescaled, length(holes))
  expect_setequal(ens2$meta$methods_used, c("m1", "m2", "m3"))
  expect_true(all(is.finite(ens2$res$tau_hat_oof)))
  # The stored estimate carries no attributes of its own.
  expect_null(attributes(ens2$res$tau_hat_oof))

  # A row where no member is available is not tolerated: the summary quantiles
  # refuse it, so such a row stops the call rather than being carried as NA.
  cmp3 <- cmp
  cmp3$metadata$oof_tau_matrix[5L, ] <- NA_real_
  expect_error(suppressWarnings(suppressMessages(
    cf_ensemble(cmp3, method = "mean"))))
})

test_that("a weighting that found nothing is reported, not replaced", {
  # The NNLS branch keeps the answer it was given. An exact zero used to be
  # replaced with equal weights, which both invented an answer and hid it:
  # the report below only speaks when no weight clears the floor.
  wm <- .cf_ensemble_weighted_members
  expect_length(wm(c(a = 0, b = 0)), 0L)
  expect_length(wm(c(a = 1.11e-16, b = 1.11e-16)), 0L)
  expect_length(wm(c(a = 0.3, b = 0)), 1L)
  ens <- structure(
    list(method = "ensemble",
         res = list(tau_hat_oof = rep(0, 5)),
         model = list(members = list()),
         weights = c(a = 0, b = 0), meta = list()),
    class = c("cateomics_ensemble", "cf_model"))
  # And predict() names that cause rather than telling the caller to rebuild
  # an object that was built correctly.
  expect_error(predict(ens, newdata = data.frame(x = 1:5)),
               class = "cateomics_ensemble_zero_weights_newdata")
})
