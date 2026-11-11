## Bioconductor core smoke tests for the SE / MAE round-trip:
##   cf_benefit_groups, cf_aggregate_vi, cf_attach_results.
## (Downstream biology is delegated to standard Bioconductor tools via the
## vignette recipes.)

make_demo <- function(n = 200, p = 6, G = 50, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  df <- data.frame(Y = Y, T = T, X)
  expr <- matrix(rpois(G * n, lambda = 20), G, n)
  rownames(expr) <- paste0("gene", seq_len(G))
  list(df = df, cv = colnames(X), expr = expr)
}

test_that("cf_benefit_groups partitions samples by tau-hat quantile", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  grp_med <- cf_benefit_groups(cmp$models$grf, method = "median")
  expect_s3_class(grp_med, "cateomics_groups")
  expect_identical(grp_med$labels, c("low_benefit", "high_benefit"))
  expect_length(grp_med$groups, length(cmp$models$grf$res$tau_hat_oof))

  grp_q <- cf_benefit_groups(cmp$models$grf, method = "quartile")
  expect_length(grp_q$labels, 4L)
})

## Differential expression and enrichment are delegated to standard
## Bioconductor tools (limma / fgsea) via the omics-workflow vignette; the
## upstream honest ingredient is cf_benefit_groups.

test_that("cf_aggregate_vi runs the RRA consensus", {
  skip_if_not_installed("grf")
  skip_if_not_installed("RobustRankAggreg")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  # RRA is the single cross-method VI consensus.
  agg <- suppressWarnings(cf_aggregate_vi(cmp))
  expect_s3_class(agg, "cateomics_aggregate_vi")
  expect_s3_class(agg$results, "data.frame")
  # What was aggregated, and how, is recorded once, with the other settings.
  expect_equal(agg$metadata$method, "rra")
  expect_null(agg[["method"]])
  expect_null(agg[["methods"]])
  expect_null(agg$metadata$n_methods_used)
  expect_equal(length(agg$metadata$methods), 1L)
  expect_true(all(c("feature_name", "score", "rank") %in% colnames(agg$results)))
  # The seed the surrogates were fitted with is the one this call used.
  expect_equal(agg$metadata$random_state, 0L)
  # The dropped signed_mean / rank-only variants are no longer accepted.
  expect_error(suppressWarnings(cf_aggregate_vi(cmp, method = "signed_mean")))
})

test_that("cf_attach_results writes cf_tau / cf_group to data.frame", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  grp <- cf_benefit_groups(cmp$models$grf, method = "median")
  df2 <- cf_attach_results(d$df, cmp$models$grf, what = "all",
                            groups = grp)
  expect_true(all(c("cf_tau", "cf_group") %in% colnames(df2)))
  expect_true(is.factor(df2$cf_group))
  # The average effect is one number for the whole fit, so it is not written
  # as n copies of itself.
  expect_false("cf_ate" %in% colnames(df2))

  # A field named explicitly whose input is missing is an error, not silence.
  expect_error(cf_attach_results(d$df, cmp$models$grf, what = "group"),
               class = "cateomics_attach_missing_input")
  expect_error(cf_attach_results(d$df, cmp$models$grf, what = "vi"),
               class = "cateomics_attach_missing_input")
  # Through "all" the same gap is reported and skipped.
  expect_message(df3 <- cf_attach_results(d$df, cmp$models$grf, what = "all"),
                 class = "cateomics_attach_all_partial")
  expect_true("cf_tau" %in% colnames(df3))
  expect_false("cf_group" %in% colnames(df3))
})

## cf_aggregate_vi returns a cateomics_aggregate_vi; pathway analysis of the
## ranked features is delegated to Bioconductor tools (fgsea) in the
## omics-workflow vignette.
