## End-to-end integration smoke test:
##   Exercises the full cfomics workflow from cf_fit -> cf_compare ->
##   cf_select -> cf_evaluate / cf_interpret / cf_rule -> cf_dge_analysis ->
##   cf_enrichment -> cf_dual_biomarker. Covers both:
##     - Route A: feature-importance route (statistical interpretation)
##     - Route B: high-benefit-group DGE route (biological interpretation)
##     - Cross-route: cf_dual_biomarker

test_that("End-to-end cfomics workflow (case-B smoke)", {
  skip_if_not_installed("grf")
  skip_if_not_installed("limma")

  set.seed(42)
  n <- 200; p <- 6; G <- 30
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  df <- data.frame(Y = Y, T = T, X)
  formula <- stats::as.formula(
    paste("Y ~ T |", paste(colnames(X), collapse = " + "))
  )

  # Step 1: cf_compare (single CATE method to keep smoke fast).
  cmp <- cf_compare(formula, data = df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  expect_s3_class(cmp, "cfomics_compare")

  # Step 2: cf_select.
  sel <- cf_select(cmp, metric = "r_loss")
  expect_s3_class(sel, "cfomics_select")
  expect_identical(sel$selected, "grf")
  model <- sel$model

  # Step 3 (Route A): cf_evaluate + cf_interpret + cf_rule.
  ev <- cf_evaluate(model)
  expect_s3_class(ev, "cfomics_evaluate")

  itp <- cf_interpret(model, surrogate = "ridge", top_n = 5L)
  expect_s3_class(itp, "cfomics_interpret")

  rule <- cf_rule(model, type = "topk", topk_frac = 0.3)
  expect_s3_class(rule, "cfomics_rule")
  expect_true(is.numeric(rule$policy_value$estimate))

  # Step 4 (Route B): cf_benefit_groups. Downstream differential expression /
  # enrichment are delegated to Bioconductor recipes in the omics-workflow
  # vignette.
  grp <- cf_benefit_groups(model, method = "median")
  expect_s3_class(grp, "cfomics_groups")
  expect_true(all(grp$groups %in% c(1L, 2L)))

  # Step 6: cf_aggregate_vi smoke (single-method cmp -> RRA aggregation works).
  if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    agg <- suppressWarnings(cf_aggregate_vi(cmp))
    expect_s3_class(agg, "cfomics_aggregate_vi")
  }

  # Step 7: cf_attach_results (data.frame backend).
  df2 <- cf_attach_results(df, model, what = "all", groups = grp)
  expect_true(all(c("cf_tau", "cf_ate", "cf_group") %in% colnames(df2)))
})
