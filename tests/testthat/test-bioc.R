## Bioconductor core smoke tests for the SE / MAE round-trip:
##   cf_benefit_groups, cf_aggregate_vi, cf_attach_results, as_cf_data.
## (Downstream biology is delegated to standard Bioconductor tools via the
## vignette recipes.)

make_demo <- function(n = 200, p = 6, G = 50, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  df <- data.frame(Y = Y, T = T, X)
  formula <- stats::as.formula(
    paste("Y ~ T |", paste(colnames(X), collapse = " + "))
  )

  expr <- matrix(rpois(G * n, lambda = 20), G, n)
  rownames(expr) <- paste0("gene", seq_len(G))
  list(df = df, formula = formula, expr = expr)
}

test_that("cf_benefit_groups partitions samples by tau-hat quantile", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  grp_med <- cf_benefit_groups(cmp$models$grf, method = "median")
  expect_s3_class(grp_med, "cfomics_groups")
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
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  # RRA is the single cross-method VI consensus.
  agg <- suppressWarnings(cf_aggregate_vi(cmp))
  expect_s3_class(agg, "cfomics_aggregate_vi")
  expect_s3_class(agg$results, "data.frame")
  expect_equal(agg$method, "rra")
  expect_true(all(c("feature_name", "score", "rank") %in% colnames(agg$results)))
  # The dropped signed_mean / rank-only variants are no longer accepted.
  expect_error(suppressWarnings(cf_aggregate_vi(cmp, method = "signed_mean")))
})

test_that("cf_attach_results writes cf_tau / cf_ate / cf_group to data.frame", {
  skip_if_not_installed("grf")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  grp <- cf_benefit_groups(cmp$models$grf, method = "median")
  df2 <- cf_attach_results(d$df, cmp$models$grf, what = "all",
                            groups = grp)
  expect_true(all(c("cf_tau", "cf_ate", "cf_group") %in% colnames(df2)))
  expect_true(is.factor(df2$cf_group))
})

## cf_aggregate_vi returns a cfomics_aggregate_vi; pathway analysis of the
## ranked features is delegated to Bioconductor tools (fgsea) in the
## omics-workflow vignette.
