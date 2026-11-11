## VIZ topic smoke tests: cf_plot_compare / cf_plot_evaluate(type=) /
## cf_plot_policy_curve / cf_plot_interpret / cf_plot_stability /
## cf_plot_rule / cf_plot_cate / cf_plot_diagnostics.

make_demo <- function(n = 300, p = 6, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  T <- rbinom(n, 1, plogis(X[, 1]))
  Y <- 0.4 * T + 0.3 * X[, 1] + 0.5 * X[, 1] * T + rnorm(n, sd = 0.3)
  list(df = data.frame(Y = Y, T = T, X), cv = colnames(X))
}

test_that("cf_plot_rule renders the policy tree for cf_rule(type='tree')", {
  skip_if_not_installed("grf")
  d <- make_demo()
  m <- cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", random_state = 42L,
          params = list(num.trees = 200L, tune.parameters = "none"))
  rule <- cf_rule(m, type = "tree", tree_depth = 3L)
  pdf(tempfile(fileext = ".pdf"))
  # With DiagrammeR the widget is handed back so it can be rendered outside an
  # auto-printing console; without it the text summary is printed instead and
  # NULL comes back invisibly. Both paths are the shipped behaviour.
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    expect_s3_class(cf_plot_rule(rule), "htmlwidget")
  } else {
    expect_invisible(cf_plot_rule(rule))
  }
  dev.off()
  # a misspelled argument is rejected rather than ignored
  expect_error(cf_plot_rule(rule, backend = "base"),
               class = "cateomics_unknown_argument")
})

test_that("cf_plot_rule errors on non-tree rule", {
  skip_if_not_installed("grf")
  d <- make_demo()
  m <- cf_fit(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, method = "grf", random_state = 42L,
          params = list(num.trees = 200L, tune.parameters = "none"))
  rule_th <- cf_rule(m, type = "threshold", threshold = 0)
  expect_error(cf_plot_rule(rule_th),
               class = "cateomics_invalid_input")
})

test_that("cf_plot_compare returns a ranked-bar ggplot", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  p <- cf_plot_compare(cmp)
  expect_s3_class(p, "ggplot")
  # The returned plot must not capture the source object: printing one that
  # does made RStudio segfault while introspecting it. Compare what the plot
  # actually holds, not object.size(), which walks the whole environment
  # chain into the package namespace and so measures nothing useful here.
  expect_lt(nrow(p$data), 100L)
  expect_false(any(vapply(ls(p$plot_env),
                          function(nm) inherits(get(nm, envir = p$plot_env),
                                                "cateomics_compare"),
                          logical(1))))
})

test_that("cf_plot_evaluate (rate/gates/blp) from cf_evaluate", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  ev <- cf_evaluate(cmp$models$grf)
  expect_s3_class(cf_plot_evaluate(ev, type = "blp"), "ggplot")
  expect_s3_class(cf_plot_evaluate(ev, type = "rate"), "ggplot")
  expect_s3_class(cf_plot_evaluate(ev, type = "gates"), "ggplot")
})

test_that("cf_plot_policy_curve from cf_rule(grid=TRUE)", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  rg <- cf_rule(cmp$models$grf, type = "topk", grid = TRUE)
  expect_s3_class(cf_plot_policy_curve(rg), "ggplot")
})

test_that("cf_plot_interpret returns a signed-SHAP bar ggplot", {
  skip_if_not_installed("grf")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("shapviz")   # cf_plot_interpret delegates rendering to shapviz
  d <- make_demo()
  cmp <- cf_compare(Y = "Y", T = "T", X_baseline = d$cv, data = d$df, methods = "grf", n_folds = 3L,
          verbose = FALSE, params = list(num.trees = 200L,
              tune.parameters = "none"))
  itp <- suppressMessages(suppressWarnings(
    cf_interpret(cmp$models$grf, surrogate = "ridge", top_n = 5L)))
  p <- cf_plot_interpret(itp)
  expect_s3_class(p, "ggplot")
  # shapviz titles nothing, so the title is the package's; every view the
  # family returns as a ggplot carries one.
  expect_identical(p$labels$title, "SHAP importance (surrogate model)")
  pb <- cf_plot_interpret(itp, type = "beeswarm")
  expect_identical(pb$labels$title, "SHAP importance (surrogate model)")
  pdep <- cf_plot_interpret(
    itp, type = "dependence",
    feature = itp$vi[[1L]]$features$feature_name[1L])
  expect_identical(pdep$labels$title, "SHAP dependence")
})

# cf_plot_volcano and cf_plot_heatmap are not part of the plot family;
# cf_plot_interpret covers the interpretation views.


# cf_plot_evaluate honours `type` ----------------------------------------

test_that("each evaluate view is the one that was asked for", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("grf")
  set.seed(3)
  n <- 200L
  df <- data.frame(T = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  df$Y <- df$x2 + df$T * (2 * df$x1) + rnorm(n, sd = 0.3)
  fit <- suppressMessages(cf_fit(data = df, Y = "Y", T = "T",
                                 X_baseline = c("x1", "x2"), method = "grf",
                                 partition = "oof", random_state = 1L,
                                 propensity_score = 0.5))
  ev <- suppressMessages(cf_evaluate(fit, random_state = 1L))

  plots <- lapply(c("rate", "gates", "blp"),
                  function(ty) suppressMessages(cf_plot_evaluate(ev, type = ty)))
  # Returning the same view for every `type` would pass an
  # expect_s3_class(p, "ggplot") check three times over.
  expect_match(plots[[1]]$labels$y, "TOC")
  expect_match(plots[[2]]$labels$y, "Group ATE")
  expect_match(plots[[3]]$labels$x, "Predicted CATE")
  expect_false(identical(plots[[1]]$data, plots[[2]]$data))
  expect_false(identical(plots[[2]]$data, plots[[3]]$data))
})


# The calibration view is per-sample, not per-bin -------------------------

test_that("the BLP view draws one point per sample and names its lines", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("grf")
  set.seed(5)
  n <- 150L
  df <- data.frame(T = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  df$Y <- df$x2 + df$T * (2 * df$x1) + rnorm(n, sd = 0.3)
  fit <- suppressMessages(cf_fit(data = df, Y = "Y", T = "T",
                                 X_baseline = c("x1", "x2"), method = "grf",
                                 partition = "oof", random_state = 1L,
                                 propensity_score = 0.5))
  ev <- suppressMessages(cf_evaluate(fit, random_state = 1L))

  # The evaluation carries the pairs the regression was fitted on.
  expect_equal(nrow(ev$blp$points), n)
  expect_named(ev$blp$points, c("predicted", "observed"))
  # The stored coefficients are that regression: re-fitting it on the raw
  # scale has to give back the slope, and the intercept once it is shifted
  # off the centred regressor.
  m <- stats::lm(observed ~ predicted, data = ev$blp$points)
  expect_equal(unname(stats::coef(m)[2L]), ev$blp$slope$estimate)
  expect_equal(unname(stats::coef(m)[1L]),
               ev$blp$intercept -
                 ev$blp$slope$estimate * mean(ev$blp$points$predicted))
  # The bins are still reported, and they are the group effects GATES draws,
  # which is why the calibration view no longer draws them too.
  # The binned copy is gone; the view reads the per-sample pairs instead and
  # the group effects are reported once, by `gates`.
  expect_null(ev$blp$plot_data)
  expect_s3_class(ev$blp$points, "data.frame")

  p <- suppressMessages(cf_plot_evaluate(ev, type = "blp"))
  expect_equal(nrow(p$data), n)
  expect_match(p$labels$y, "pseudo-outcome")
  expect_match(p$labels$subtitle, "dashed = slope 1", fixed = TRUE)
  # No fixed aspect ratio: the two axes share a unit but not a scale.
  expect_false(inherits(p$coordinates, "CoordFixed"))
  # Two reference lines (slope 1 and the fit) plus the points.
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_equal(sum(geoms == "GeomAbline"), 2L)
  expect_equal(sum(geoms == "GeomPoint"), 1L)

  pg <- suppressMessages(cf_plot_evaluate(ev, type = "gates"))
  expect_match(pg$labels$subtitle, "top - bottom", fixed = TRUE)
})


# What each view draws ------------------------------------------------------
# Each block checks the built plot rather than the picture.

make_viz_fixture <- function() {
  set.seed(7)
  n <- 150L
  df <- data.frame(T = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  df$Y <- df$x2 + df$T * (2 * df$x1) + rnorm(n, sd = 0.3)
  prm <- list(num.trees = 200L, tune.parameters = "none")
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "grf",
    partition = "oof", n_folds = 3L, random_state = 1L, params = prm))
  list(df = df, fit = fit)
}

geoms_of <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], character(1))

test_that("the overlap view shows the scores themselves, in the arm palette", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  fx <- make_viz_fixture()
  # The requested-variable view reads the set the diagnostics were asked for.
  dg <- suppressMessages(cf_diagnostics(fx$fit, variables = "x1"))
  p  <- suppressMessages(cf_plot_diagnostics(dg, type = "overlap"))
  # A density draws each arm past its own data, so the picture alone cannot say
  # where the observations are; the rug can.
  expect_true("GeomRug" %in% geoms_of(p))
  built <- ggplot2::ggplot_build(p)
  cols <- unlist(lapply(built$data, function(d)
    c(d$fill, d$colour)[nzchar(c(d$fill, d$colour))]))
  expect_true(all(c("#009E73", "#CC79A7") %in% cols))
  # The three diagnostics views use the family's theme, not a bare one.
  for (ty in c("overlap", "smd", "variable")) {
    pp <- suppressMessages(cf_plot_diagnostics(dg, type = ty))
    expect_identical(pp$theme$plot.title$face, "bold", info = ty)
    expect_identical(pp$theme$legend.position, "bottom", info = ty)
  }
})

test_that("the calibration view shades the slopes inside the interval", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  fx <- make_viz_fixture()
  ev <- suppressMessages(cf_evaluate(fx$fit, random_state = 1L))
  p  <- suppressMessages(cf_plot_evaluate(ev, type = "blp"))
  expect_true("GeomRibbon" %in% geoms_of(p))
  # The band is a pencil of lines through the centroid, so at the centre it
  # has no width and it widens outwards.
  band <- ggplot2::ggplot_build(p)$data[[which(geoms_of(p) == "GeomRibbon")]]
  w <- band$ymax - band$ymin
  expect_equal(min(w), 0, tolerance = 0.02 * max(w))
  expect_gt(max(w), 0)
  # The subtitle names both reference lines and the band, and wraps so that
  # the clause naming the band is not cut at the default figure width.
  expect_match(p$labels$subtitle, "dashed = slope 1", fixed = TRUE)
  expect_match(p$labels$subtitle, "shaded = slopes inside the interval",
               fixed = TRUE)
  expect_match(p$labels$subtitle, "\n", fixed = TRUE)
})

test_that("the CATE view marks the average effect", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  fx <- make_viz_fixture()
  p <- suppressMessages(cf_plot_cate(fx$fit))
  i <- which(geoms_of(p) == "GeomVline")
  expect_length(i, 1L)
  expect_equal(ggplot2::ggplot_build(p)$data[[i]]$xintercept,
               unname(fx$fit$res$ate))
  expect_match(p$labels$subtitle, "average effect", fixed = TRUE)
})

test_that("the stability heatmap is titled for what it shows", {
  skip_if_not_installed("grf")
  skip_if_not_installed("ggplot2")
  fx <- make_viz_fixture()
  g <- function(nt) list(method = "grf",
                         params = list(num.trees = nt, tune.parameters = "none"))
  cmp <- suppressMessages(cf_compare(
    data = fx$df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), n_folds = 3L,
    verbose = FALSE, random_state = 1L,
    configs = list(a = g(300L), b = g(120L), c = g(200L))))
  stab <- suppressMessages(cf_model_stability(cmp))
  p <- cf_plot_stability(stab)
  expect_identical(p$labels$title, "Pairwise tau correlation across methods")
  expect_match(p$labels$subtitle, "^Spearman correlation")
})

test_that("the policy tree says what its leaves do and on what scale", {
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  skip_if_not_installed("DiagrammeR")
  fx <- make_viz_fixture()
  rule <- suppressMessages(cf_rule(fx$fit, type = "tree", tree_depth = 2L))
  w <- cf_plot_rule(rule)
  dot <- w$x$diagram
  # The action numbers alone do not say which one treats.
  expect_match(dot, "do not treat", fixed = TRUE)
  expect_match(dot, "Policy tree: who to treat", fixed = TRUE)
  expect_match(dot, "preprocessed covariate scale", fixed = TRUE)
  expect_false(grepl("action = ", dot, fixed = TRUE))
})
