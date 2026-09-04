## VIZ topic smoke tests: cf_plot_compare / cf_plot_evaluate(type=) /
## cf_plot_policy_curve / cf_plot_interpret / cf_plot_stability /
## cf_plot_rule / cf_plot_cate / cf_plot_diagnostics.

make_demo <- function(n = 300, p = 6, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X),
       formula = stats::as.formula(
         paste("Y ~ T |", paste(colnames(X), collapse = " + "))
       ))
}

test_that("cf_plot_rule renders rpart tree for cf_rule(type='tree')", {
  skip_if_not_installed("grf")
  skip_if_not_installed("rpart")
  d <- make_demo()
  m <- cf_fit(d$formula, data = d$df, method = "grf",
              random_state = 42L, num.trees = 200L,
              tune.parameters = "none")
  rule <- cf_rule(m, type = "tree", tree_depth = 3L)
  pdf(tempfile(fileext = ".pdf"))
  expect_invisible(cf_plot_rule(rule))
  dev.off()
  # base backend works without rpart.plot
  pdf(tempfile(fileext = ".pdf"))
  expect_invisible(cf_plot_rule(rule, backend = "base"))
  dev.off()
})

test_that("cf_plot_rule errors on non-tree rule", {
  skip_if_not_installed("grf")
  d <- make_demo()
  m <- cf_fit(d$formula, data = d$df, method = "grf",
              random_state = 42L, num.trees = 200L,
              tune.parameters = "none")
  rule_th <- cf_rule(m, type = "threshold", threshold = 0)
  expect_error(cf_plot_rule(rule_th),
               class = "cfomics_invalid_input")
})

test_that("cf_plot_compare returns a ranked-bar ggplot", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  p <- cf_plot_compare(cmp)
  expect_s3_class(p, "ggplot")
  # The returned plot must not capture the (large) source object.
  expect_lt(as.numeric(utils::object.size(p)),
            as.numeric(utils::object.size(cmp)))
})

test_that("cf_plot_evaluate (rate/gates/blp) from cf_evaluate", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  ev <- cf_evaluate(cmp$models$grf)
  expect_s3_class(cf_plot_evaluate(ev, type = "blp"), "ggplot")
  expect_s3_class(cf_plot_evaluate(ev, type = "rate"), "ggplot")
  expect_s3_class(cf_plot_evaluate(ev, type = "gates"), "ggplot")
})

test_that("cf_plot_policy_curve from cf_rule(grid=TRUE)", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  rg <- cf_rule(cmp$models$grf, type = "topk", grid = TRUE)
  expect_s3_class(cf_plot_policy_curve(rg), "ggplot")
})

test_that("cf_plot_interpret returns a signed-SHAP bar ggplot", {
  skip_if_not_installed("grf")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("shapviz")   # cf_plot_interpret delegates rendering to shapviz
  d <- make_demo()
  cmp <- cf_compare(d$formula, data = d$df, methods = "grf",
                     n_folds = 3L, verbose = FALSE,
                     num.trees = 200L, tune.parameters = "none")
  itp <- suppressMessages(suppressWarnings(
    cf_interpret(cmp$models$grf, surrogate = "ridge", top_n = 5L)))
  p <- cf_plot_interpret(itp)
  expect_s3_class(p, "ggplot")
})

# cf_plot_volcano and cf_plot_heatmap are not part of the plot family;
# cf_plot_interpret covers the interpretation views.
