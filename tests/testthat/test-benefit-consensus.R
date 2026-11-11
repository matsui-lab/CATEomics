# cf_benefit_groups: the split itself, and how much the other candidates of a
# comparison support it.
#
# The support answers "how stable is this grouping across methods". It is
# counted per method as same group or not, the method the groups came from is
# not one of the voters, and a method that leaves a sample ungrouped does not
# support it.

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
      class = c("cateomics_compare", "list")
    ),
    tau_select = tau_mat[, 1L]
  )
}

.mock_compare_from <- function(tau_mat) {
  structure(
    list(scoreboard = data.frame(method = colnames(tau_mat),
                                 stringsAsFactors = FALSE),
         models   = list(),
         metadata = list(oof_tau_matrix = tau_mat)),
    class = c("cateomics_compare", "list"))
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


# the split on its own -----------------------------------------------------

test_that("the object reports the split once, with the source of the effects", {
  mock <- .make_mock_compare_consensus(n = 60L, K = 3L)
  g <- cf_benefit_groups(.make_select_model(mock$tau_select), method = "median")
  expect_s3_class(g, "cateomics_groups")
  expect_setequal(names(g), c("groups", "labels", "thresholds", "n_per_group",
                              "tau", "consensus", "metadata"))
  # The split is named once, in metadata, next to where the effects came from.
  expect_null(g[["method"]])
  expect_identical(g$metadata$method, "median")
  expect_identical(g$metadata$tau_source, "oof")
  expect_identical(g$metadata$source_method, "grf")
  expect_equal(g$n_per_group, unname(tabulate(g$groups, nbins = 2L)))
  # Without a source there is nothing to count agreement over.
  expect_null(g$consensus)
})

test_that("in-sample effects are used with a warning and recorded as such", {
  m <- structure(list(method = "grf",
                      res = list(tau_hat_oof = NULL, tau_hat = rnorm(40))),
                 class = "cf_model")
  expect_warning(g <- cf_benefit_groups(m),
                 class = "cateomics_partition_none_warning")
  expect_identical(g$metadata$tau_source, "in-sample")
})

test_that("thresholds are refused unless the split is the custom one", {
  mock <- .make_mock_compare_consensus(n = 40L, K = 2L)
  m <- .make_select_model(mock$tau_select)
  # The quantile splits would silently ignore them.
  expect_error(cf_benefit_groups(m, method = "median", thresholds = 0.1),
               class = "cateomics_invalid_input")
  expect_error(cf_benefit_groups(m, method = "custom"),
               class = "cateomics_invalid_input")
  g <- cf_benefit_groups(m, method = "custom", thresholds = c(-0.2, 0.2))
  expect_length(g$labels, 3L)
  expect_equal(g$thresholds, c(-0.2, 0.2))
})


# the support --------------------------------------------------------------

test_that("the support is a share over the other methods, one row per sample", {
  mock <- .make_mock_compare_consensus(n = 60L, K = 4L)
  g <- cf_benefit_groups(.make_select_model(mock$tau_select),
                         method = "median", consensus_source = mock$compare)
  cn <- g$consensus
  expect_setequal(names(cn), c("group_support", "methods_used",
                               "K_used", "filter"))
  # There is no `methods_dropped`: the only two reasons a method does not vote
  # are that it is the method the groups came from and that its estimates are
  # unusable, and one field holding both would be read as the second.
  expect_null(cn[["methods_dropped"]])
  expect_length(cn$group_support, 60L)
  expect_true(all(cn$group_support >= 0 & cn$group_support <= 1))
  # The method the groups came from does not vote on them.
  expect_identical(cn$K_used, 3L)
  expect_false("m1" %in% cn$methods_used)
  # The support takes K + 1 values and nothing else.
  expect_true(all(round(cn$group_support * cn$K_used) %in% 0:3))
})

test_that("the method that produced the groups is dropped by its estimates, not its name", {
  n <- 50L
  tau_mat <- matrix(rnorm(n * 3L), nrow = n,
                    dimnames = list(NULL, c("grf_big", "grf_sml", "bcf")))
  cmp <- .mock_compare_from(tau_mat)
  # A grid comparison labels the same method once per base learner, so the
  # model's own `method` ("grf") matches no column name.
  m <- .make_select_model(tau_mat[, "grf_sml"], model_method = "grf")
  cn <- cf_benefit_groups(m, consensus_source = cmp)$consensus
  expect_identical(cn$K_used, 2L)
  expect_setequal(cn$methods_used, c("grf_big", "bcf"))
})

test_that("support is 1 when every other method agrees and (K-1)/K when one dissents", {
  n <- 40L
  base <- seq(-1, 1, length.out = n)
  # Three columns ordering the samples the same way: every median split agrees.
  tau_mat <- cbind(a = base, b = base * 2, c = base + 5)
  cmp <- .mock_compare_from(tau_mat)
  cn <- cf_benefit_groups(.make_select_model(base), consensus_source = cmp)$consensus
  expect_identical(cn$K_used, 2L)
  expect_true(all(cn$group_support == 1))

  # One column reversed: it puts every sample in the other group.
  tau_mat2 <- cbind(a = base, b = base * 2, d = -base)
  cn2 <- cf_benefit_groups(.make_select_model(base),
                           consensus_source = .mock_compare_from(tau_mat2))$consensus
  expect_identical(cn2$K_used, 2L)
  expect_true(all(cn2$group_support == 0.5))
})

test_that("a sample no method grouped is reported as ungrouped, not as agreement", {
  mock <- .make_mock_compare_consensus(n = 60L, K = 3L)
  g <- cf_benefit_groups(.make_select_model(mock$tau_select),
                         method = "extreme_quartile",
                         consensus_source = mock$compare)
  # The middle half is in no group on the reported split.
  expect_true(anyNA(g$groups))
  expect_true(all(is.na(g$consensus$group_support[is.na(g$groups)])))
  # Where the reported split does give a group, the share is defined and the
  # denominator is still the number of voters: a voter whose own middle half
  # covers the sample simply does not support it.
  ok <- !is.na(g$groups)
  expect_false(anyNA(g$consensus$group_support[ok]))
  expect_true(all(g$consensus$group_support[ok] <= 1))
})

test_that("consensus_methods restricts the voters and an unknown name is an error", {
  mock <- .make_mock_compare_consensus(n = 50L, K = 4L)
  m <- .make_select_model(mock$tau_select)
  cn <- cf_benefit_groups(m, consensus_source = mock$compare,
                          consensus_methods = c("m1", "m2", "m3"))$consensus
  expect_setequal(cn$methods_used, c("m2", "m3"))
  expect_identical(cn$K_used, 2L)
  expect_error(
    cf_benefit_groups(m, consensus_source = mock$compare,
                      consensus_methods = c("m2", "nope")),
    class = "cateomics_consensus_unknown_methods")
})

test_that("a comparison holding only the reporting method warns and returns NA", {
  n <- 30L
  tau <- rnorm(n)
  cmp <- .mock_compare_from(cbind(only = tau))
  expect_warning(
    cn <- cf_benefit_groups(.make_select_model(tau),
                            consensus_source = cmp)$consensus,
    class = "cateomics_consensus_no_clean_methods")
  expect_identical(cn$K_used, 0L)
  expect_true(all(is.na(cn$group_support)))
})


# the filter ---------------------------------------------------------------

test_that("the threshold leaves weakly supported samples ungrouped and says how many", {
  mock <- .make_mock_compare_consensus(n = 80L, K = 4L)
  m <- mock$compare$metadata$oof_tau_matrix[, 1L]
  g_all <- cf_benefit_groups(.make_select_model(m), method = "tertile",
                             consensus_source = mock$compare)
  g_f <- cf_benefit_groups(.make_select_model(m), method = "tertile",
                           consensus_source = mock$compare,
                           consensus_min_agreement = 0.5)
  drop <- g_all$consensus$group_support < 0.5
  expect_identical(sum(is.na(g_f$groups)), sum(drop))
  expect_identical(g_f$consensus$filter$n_dropped, sum(drop))
  expect_identical(g_f$consensus$filter$min_agreement, 0.5)
  # The support itself is untouched by the filter; only `groups` changes.
  expect_identical(g_f$consensus$group_support, g_all$consensus$group_support)
  expect_equal(sum(g_f$n_per_group), nrow(mock$compare$metadata$oof_tau_matrix) - sum(drop))
})

test_that("the filter needs a source and a threshold in [0, 1]", {
  mock <- .make_mock_compare_consensus(n = 40L, K = 3L)
  m <- .make_select_model(mock$tau_select)
  expect_error(cf_benefit_groups(m, consensus_min_agreement = 0.5),
               class = "cateomics_consensus_source_missing")
  expect_error(cf_benefit_groups(m, consensus_source = mock$compare,
                                 consensus_min_agreement = 1.5),
               class = "cateomics_consensus_invalid_threshold")
})


# ensemble -----------------------------------------------------------------

test_that("an ensemble can be grouped, with or without a source to count over", {
  mock <- .make_mock_compare_consensus(n = 40L, K = 3L)
  combined <- rowMeans(mock$compare$metadata$oof_tau_matrix)
  ens <- structure(
    list(method  = "ensemble",
         res     = list(tau_hat_oof = combined, tau_hat = NULL),
         meta    = list(formula = NULL)),
    class = c("cateomics_ensemble", "cf_model"))
  expect_s3_class(suppressWarnings(cf_benefit_groups(ens)), "cateomics_groups")
  # Every candidate votes: an ensemble is not one of them, so nothing is
  # excluded. Its estimates are built from theirs, which is why the support
  # runs high; that is arithmetic, not agreement between methods.
  cn <- cf_benefit_groups(ens, consensus_source = mock$compare)$consensus
  expect_identical(cn$K_used, 3L)
  expect_length(cn$group_support, 40L)
})


# print --------------------------------------------------------------------

test_that("print reports the counts rather than a summary of them", {
  mock <- .make_mock_compare_consensus(n = 60L, K = 4L)
  g <- cf_benefit_groups(.make_select_model(mock$tau_select),
                         consensus_source = mock$compare,
                         consensus_min_agreement = 0.5)
  txt <- paste(capture.output(print(g)), collapse = "\n")
  expect_match(txt, "consensus (3 other methods)", fixed = TRUE)
  expect_match(txt, "0/3 agree", fixed = TRUE)
  expect_match(txt, "3/3 agree", fixed = TRUE)
  expect_match(txt, "filtered at 0.5", fixed = TRUE)
  # A median or a mean of a K + 1 valued share would hide the shape.
  expect_false(grepl("median =", txt, fixed = TRUE))
})


# row alignment --------------------------------------------------------------

test_that("the support is matched to the comparison by sample name, not by position", {
  skip_if_not_installed("grf")
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 60L, p = 80L, seed = 1L)
  g <- function(nt, mn = 5L) list(
    method = "grf",
    params = list(num.trees = nt, tune.parameters = "none", min.node.size = mn))
  cmp <- suppressMessages(cf_compare(
    data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
    n_features = 40L, n_folds = 3L, verbose = FALSE, random_state = 1L,
    configs = list(a = g(300L), b = g(120L), c = g(200L, 20L))))
  ids <- CATEomics:::.cf_consensus_row_ids(cmp, ncol(se))
  expect_length(ids, ncol(se))
  tau <- cmp$models$a$res$tau_hat_oof
  grp <- CATEomics:::.cf_benefit_cut(tau, "median")$groups
  f <- CATEomics:::.cf_benefit_group_support
  plain <- suppressMessages(f(cmp, NULL, grp, tau, "median", NULL, ids))
  set.seed(3)
  perm <- sample(ncol(se))
  # The same samples and the same estimates, handed over in another order: the
  # fit is not redone, so only the alignment is under test.
  shuffled <- suppressMessages(
    f(cmp, NULL, grp[perm], tau[perm], "median", NULL, ids[perm]))
  expect_identical(shuffled$K_used, plain$K_used)
  expect_equal(shuffled$group_support[order(perm)], plain$group_support)
  # Without names there is nothing to align on, and the numbers do move --
  # which is what equal lengths alone used to allow through in silence.
  blind <- suppressMessages(
    f(cmp, NULL, grp[perm], tau[perm], "median", NULL, NULL))
  expect_false(isTRUE(all.equal(blind$group_support[order(perm)],
                                plain$group_support)))
  # A comparison over other samples is refused rather than lined up by row.
  expect_error(
    suppressMessages(f(cmp, NULL, grp, tau, "median", NULL,
                       paste0("other", seq_len(ncol(se))))),
    class = "cateomics_consensus_sample_mismatch")
})
