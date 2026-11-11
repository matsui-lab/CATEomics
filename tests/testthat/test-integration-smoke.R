## End-to-end integration smoke test:
##   Exercises the full CATEomics workflow from cf_fit -> cf_compare ->
##   cf_select -> cf_ensemble -> cf_evaluate / cf_interpret / cf_rule ->
##   cf_benefit_groups -> cf_aggregate_vi, and writing the results back onto
##   a SummarizedExperiment with cf_attach_results.

test_that("End-to-end CATEomics workflow (case-B smoke)", {
  skip_if_not_installed("grf")
  skip_if_not_installed("limma")

  set.seed(42)
  n <- 200; p <- 6; G <- 30
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  df <- data.frame(Y = Y, T = T, X)
  # Step 1: cf_compare (single CATE method to keep smoke fast).
  cmp <- cf_compare(data = df, Y = "Y", T = "T", X_baseline = colnames(X),
          methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  expect_s3_class(cmp, "cateomics_compare")

  # Step 2: cf_select.
  sel <- cf_select(cmp, metric = "r_loss")
  expect_s3_class(sel, "cateomics_select")
  expect_identical(sel$selected, "grf")
  model <- sel$model

  # Step 3 (Route A): cf_evaluate + cf_interpret + cf_rule.
  ev <- cf_evaluate(model)
  expect_s3_class(ev, "cateomics_evaluate")

  itp <- cf_interpret(model, surrogate = "ridge", top_n = 5L)
  expect_s3_class(itp, "cateomics_interpret")

  rule <- cf_rule(model, type = "topk", topk_frac = 0.3)
  expect_s3_class(rule, "cateomics_rule")
  expect_true(is.numeric(rule$policy_value$estimate))

  # Step 4 (Route B): cf_benefit_groups. Downstream differential expression /
  # enrichment are delegated to Bioconductor recipes in the omics-workflow
  # vignette.
  grp <- cf_benefit_groups(model, method = "median")
  expect_s3_class(grp, "cateomics_groups")
  expect_true(all(grp$groups %in% c(1L, 2L)))

  # Step 6: cf_aggregate_vi smoke (single-method cmp -> RRA aggregation works).
  if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    agg <- suppressWarnings(cf_aggregate_vi(cmp))
    expect_s3_class(agg, "cateomics_aggregate_vi")
  }

  # Step 7: cf_attach_results (data.frame backend).
  df2 <- cf_attach_results(df, model, what = "all", groups = grp)
  expect_true(all(c("cf_tau", "cf_group") %in% colnames(df2)))
  expect_false("cf_ate" %in% colnames(df2))
  # The documented default call. `what` defaults to the mode "all", not to the
  # list of every item: a default holding "group" and "vi" was read as the
  # caller naming both, so this call aborted asking for inputs nobody
  # mentioned, and no shipped caller used the bare default to notice.
  expect_message(df3 <- cf_attach_results(df, model),
                 class = "cateomics_attach_all_partial")
  expect_true("cf_tau" %in% colnames(df3))
  expect_false("cf_group" %in% colnames(df3))
})
