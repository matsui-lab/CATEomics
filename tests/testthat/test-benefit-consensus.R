# Smoke tests for cf_benefit_groups(consensus = TRUE).
# Goals: consensus stats shape + ensemble hard-error + select-route works.

set.seed(31415)
.make_mock_compare_consensus <- function(n = 50L, K = 4L) {
  tau_mat <- matrix(rnorm(n * K, sd = 0.5), nrow = n, ncol = K)
  colnames(tau_mat) <- paste0("m", seq_len(K))
  list(
    compare = structure(
      list(
        scoreboard = data.frame(method = colnames(tau_mat),
                                 stringsAsFactors = FALSE),
        models     = list(),
        metadata   = list(oof_tau_matrix = tau_mat,
                          nuisance_oof   = list(mu0_oof = rnorm(n),
                                                 mu1_oof = rnorm(n),
                                                 pi_oof  = rep(0.5, n)))
      ),
      class = c("cfomics_compare", "list")
    ),
    tau_select = tau_mat[, 1L]
  )
}

.make_select_model <- function(tau, model_method = "grf") {
  structure(
    list(
      method  = model_method,
      res     = list(tau_hat_oof = tau, tau_hat = tau),
      meta    = list(formula = NULL)
    ),
    class = "cf_model"
  )
}


# ensemble rejection -------------------------------------------------------

test_that("cf_benefit_groups(consensus=TRUE) errors on ensemble input", {
  ens <- structure(
    list(
      method  = "ensemble",
      res     = list(tau_hat_oof = rnorm(40), tau_hat = rnorm(40)),
      meta    = list(formula = NULL)
    ),
    class = c("cfomics_ensemble", "cf_model")
  )
  mock <- .make_mock_compare_consensus(n = 40L, K = 3L)
  expect_error(
    cf_benefit_groups(ens, consensus = TRUE,
                         consensus_source = mock$compare),
    class = "cfomics_consensus_not_applicable_to_ensemble"
  )
})


# Consensus stats shape ----------------------------------------------------

test_that("consensus stats are returned in select-route median split", {
  mock <- .make_mock_compare_consensus(n = 60L, K = 4L)
  model <- .make_select_model(mock$tau_select)
  out <- cf_benefit_groups(model, method = "median",
                              consensus = TRUE,
                              consensus_source = mock$compare,
                              consensus_q      = 0.25)
  expect_s3_class(out, "cfomics_groups")
  cn <- out$consensus
  expect_true(!is.null(cn))
  expect_equal(length(cn$consensus_fraction), 60L)
  expect_equal(length(cn$sign_agreement),     60L)
  expect_equal(length(cn$model_sd),           60L)
  expect_equal(length(cn$unstable_flag),      60L)
  expect_setequal(cn$methods_used, c("m1", "m2", "m3", "m4"))
  expect_equal(cn$consensus_q, 0.25)
  # Consensus fractions are in [0, 1].
  expect_true(all(cn$consensus_fraction >= 0 &
                    cn$consensus_fraction <= 1))
  expect_true(all(cn$sign_agreement >= 0 & cn$sign_agreement <= 1,
                    na.rm = TRUE))
})


test_that("consensus stats work with extreme_quartile method", {
  mock <- .make_mock_compare_consensus(n = 50L, K = 3L)
  model <- .make_select_model(mock$tau_select)
  out <- cf_benefit_groups(model, method = "extreme_quartile",
                              consensus = TRUE,
                              consensus_source = mock$compare)
  expect_s3_class(out, "cfomics_groups")
  expect_true(!is.null(out$consensus))
  expect_equal(length(out$consensus$consensus_fraction), 50L)
})


# Method filter -----------------------------------------------------------

test_that("consensus_methods filter restricts the cohort", {
  mock <- .make_mock_compare_consensus(n = 50L, K = 4L)
  model <- .make_select_model(mock$tau_select)
  out <- cf_benefit_groups(model,
                              consensus = TRUE,
                              consensus_source  = mock$compare,
                              consensus_methods = c("m1", "m2"))
  expect_setequal(out$consensus$methods_used, c("m1", "m2"))
})

test_that("unknown consensus_methods raise an error", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  expect_error(
    cf_benefit_groups(model, consensus = TRUE,
                         consensus_source  = mock$compare,
                         consensus_methods = c("mghost")),
    class = "cfomics_consensus_unknown_methods"
  )
})


# Invalid q --------------------------------------------------------------

test_that("invalid consensus_q raises", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  expect_error(
    cf_benefit_groups(model, consensus = TRUE,
                         consensus_source = mock$compare,
                         consensus_q      = 0.8),
    class = "cfomics_consensus_invalid_q"
  )
})


# Missing source --------------------------------------------------------

test_that("missing consensus_source raises", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  expect_error(
    cf_benefit_groups(model, consensus = TRUE,
                         consensus_source = NULL),
    class = "cfomics_consensus_source_missing"
  )
})


# print smoke ------------------------------------------------------------

test_that("print method handles consensus output", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  out <- cf_benefit_groups(model, consensus = TRUE,
                              consensus_source = mock$compare)
  expect_output(print(out), "consensus")
})


# stability_filter ------------------------------------------------

test_that("stability_filter writes unstable individuals back as NA in groups", {
  # Make m4 a noisy column so some individuals have low sign_agreement.
  set.seed(2026)
  n <- 60L
  base <- rnorm(n)
  tau_mat <- cbind(
    m1 = base + rnorm(n, sd = 0.2),
    m2 = base + rnorm(n, sd = 0.2),
    m3 = base + rnorm(n, sd = 0.2),
    m4 = rnorm(n, sd = 1)
  )
  cmp <- structure(
    list(
      scoreboard = data.frame(method = colnames(tau_mat),
                               stringsAsFactors = FALSE),
      models     = list(),
      metadata   = list(oof_tau_matrix = tau_mat,
                        nuisance_oof   = list(mu0_oof = rnorm(n),
                                               mu1_oof = rnorm(n),
                                               pi_oof  = rep(0.5, n)))
    ),
    class = c("cfomics_compare", "list")
  )
  model <- .make_select_model(tau_mat[, 1L])

  base_out <- cf_benefit_groups(model, method = "median",
                                    consensus = TRUE,
                                    consensus_source = cmp)
  stab_out <- cf_benefit_groups(model, method = "median",
                                    consensus = TRUE,
                                    consensus_source = cmp,
                                    stability_filter = TRUE,
                                    stability_min_agreement = 0.6)

  expect_true(stab_out$metadata$stability_filter)
  expect_equal(stab_out$metadata$stability_min_agreement, 0.6)
  # The stable run has at least the same number of NAs as the base run.
  expect_true(sum(is.na(stab_out$groups)) >= sum(is.na(base_out$groups)))
})


test_that("stability_filter without consensus errors", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  expect_error(
    cf_benefit_groups(model, stability_filter = TRUE),
    class = "cfomics_stability_requires_consensus"
  )
})


test_that("stability_min_agreement bounds are enforced", {
  mock <- .make_mock_compare_consensus()
  model <- .make_select_model(mock$tau_select)
  expect_error(
    cf_benefit_groups(model, consensus = TRUE,
                         consensus_source = mock$compare,
                         stability_filter = TRUE,
                         stability_min_agreement = 1.5),
    class = "cfomics_stability_invalid_threshold"
  )
})
