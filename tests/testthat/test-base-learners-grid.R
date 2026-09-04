# Tests for cf_compare(base_learners=): grid over (method x base learner).
# The pluggable meta/DR learners fan out to one config per base learner
# (label = method_base); grf / bcf / NN stay single. `label` is the uniform
# external key; each config keeps its bare `method` field. base_learners = NULL
# reproduces the legacy single-config behaviour.

grid_demo <- function(n = 70L, p = 45L, seed = 3L) {
  cf_example_data(n = n, p = p, n_prog = 8L, n_pred = 4L, seed = seed)
}

pluggable <- c("slearner", "tlearner", "xlearner", "drlearner")


test_that(".cfomics_build_configs expands only pluggable methods", {
  # NULL = one row per method, base NA, label == method.
  cfg0 <- cfomics:::.cfomics_build_configs(c("grf", "xlearner"), NULL, pluggable)
  expect_identical(cfg0$label, c("grf", "xlearner"))
  expect_identical(cfg0$method, c("grf", "xlearner"))
  expect_true(all(is.na(cfg0$base)))

  # Grid = pluggable methods fan out; non-pluggable stay single (base NA).
  cfg <- cfomics:::.cfomics_build_configs(
    c("grf", "xlearner", "drlearner"), c("ridge", "rf"), pluggable)
  expect_identical(
    cfg$label,
    c("grf", "xlearner_ridge", "xlearner_rf", "drlearner_ridge", "drlearner_rf"))
  expect_identical(cfg$method,
                   c("grf", "xlearner", "xlearner", "drlearner", "drlearner"))
  expect_identical(cfg$base, c(NA, "ridge", "rf", "ridge", "rf"))
})


test_that("base_learners = NULL reproduces the legacy scoreboard layout", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  se <- grid_demo()
  cmp <- suppressWarnings(suppressMessages(cf_compare(
    data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
    methods = c("grf", "bcf"), n_features = 20L,
    n_burn = 5L, n_iter = 5L, n_chains = 1L,
    propensity_score = 0.5, random_state = 1, verbose = FALSE)))
  # Legacy: identity column is `method`; no label / base columns.
  expect_true("method" %in% colnames(cmp$scoreboard))
  expect_false("label" %in% colnames(cmp$scoreboard))
  expect_false("base" %in% colnames(cmp$scoreboard))
  expect_setequal(names(cmp$models), c("grf", "bcf"))
  expect_setequal(cmp$metadata$methods, c("grf", "bcf"))
  expect_null(cmp$metadata$base_learners)
})


test_that("grid mode adds label/method/base columns and keys by label", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  se <- grid_demo()
  # grf + bcf are both non-pluggable, so base_learners has no effect but the
  # grid plumbing (label/method/base, label keying) is exercised without Python.
  cmp <- suppressWarnings(suppressMessages(cf_compare(
    data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
    methods = c("grf", "bcf"), base_learners = c("ridge", "rf"),
    n_features = 20L, n_burn = 5L, n_iter = 5L, n_chains = 1L,
    propensity_score = 0.5, random_state = 1, verbose = FALSE)))
  expect_true(all(c("label", "method", "base") %in% colnames(cmp$scoreboard)))
  expect_setequal(cmp$scoreboard$label, c("grf", "bcf"))
  expect_true(all(is.na(cmp$scoreboard$base)))          # non-pluggable
  expect_setequal(names(cmp$models), c("grf", "bcf"))
  expect_setequal(cmp$metadata$methods, c("grf", "bcf"))
  expect_identical(cmp$metadata$base_learners, c("ridge", "rf"))
  # Each config keeps its BARE method (predict dispatch key), not the label.
  expect_identical(cmp$models[["grf"]]$method, "grf")
})


test_that("grid-singleton fits match the legacy fits (base ignored, byte parity)", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  se <- grid_demo()
  common <- list(data = se, Y = "outcome", T = "treatment",
                 X_baseline = c("age", "sex"), methods = c("grf", "bcf"),
                 n_features = 20L, n_burn = 5L, n_iter = 5L, n_chains = 1L,
                 propensity_score = 0.5, random_state = 1, verbose = FALSE)
  cmp0 <- suppressWarnings(suppressMessages(do.call(cf_compare, common)))
  cmp1 <- suppressWarnings(suppressMessages(
    do.call(cf_compare, c(common, list(base_learners = c("ridge", "rf"))))))
  # grf / bcf ignore base learners -> identical r_loss under both calls.
  expect_equal(cmp0$scoreboard$r_loss[cmp0$scoreboard$method == "grf"],
               cmp1$scoreboard$r_loss[cmp1$scoreboard$label == "grf"],
               tolerance = 1e-10)
  expect_equal(cmp0$scoreboard$r_loss[cmp0$scoreboard$method == "bcf"],
               cmp1$scoreboard$r_loss[cmp1$scoreboard$label == "bcf"],
               tolerance = 1e-10)
})


test_that("base_learners is mutually exclusive with model_outcome/model_final", {
  se <- grid_demo(n = 40L, p = 20L)
  expect_error(
    suppressWarnings(suppressMessages(cf_compare(
      data = se, Y = "outcome", T = "treatment", X_baseline = "age",
      methods = c("grf", "xlearner"), base_learners = "rf",
      model_outcome = "ridge", verbose = FALSE))),
    class = "cfomics_base_learners_conflict")
})


test_that("unknown base_learners abort; all-singleton set informs no effect", {
  se <- grid_demo(n = 40L, p = 20L)
  expect_error(
    suppressWarnings(suppressMessages(cf_compare(
      data = se, Y = "outcome", T = "treatment", X_baseline = "age",
      methods = "grf", base_learners = "not_a_learner", verbose = FALSE))),
    class = "cfomics_base_learners_unknown")
  skip_if_not_installed("grf")
  expect_message(
    suppressWarnings(cf_compare(
      data = se, Y = "outcome", T = "treatment", X_baseline = "age",
      methods = "grf", base_learners = "rf", n_features = 15L,
      propensity_score = 0.5, random_state = 1, verbose = FALSE)),
    class = "cfomics_base_learners_no_pluggable")
})


test_that("Python pluggable method actually varies with the base learner", {
  skip_if_no_cf_python()
  se <- grid_demo(n = 80L, p = 50L, seed = 2L)
  cmp <- suppressWarnings(suppressMessages(cf_compare(
    data = se, Y = "outcome", T = "treatment", X_baseline = c("age", "sex"),
    methods = c("grf", "xlearner"), base_learners = c("ridge", "rf"),
    n_features = 25L, propensity_score = 0.5, random_state = 1, verbose = FALSE)))
  expect_setequal(names(cmp$models),
                  c("grf", "xlearner_ridge", "xlearner_rf"))
  xr <- cmp$models[["xlearner_ridge"]]
  xf <- cmp$models[["xlearner_rf"]]
  # bare method kept; spec records the base learner.
  expect_identical(xr$method, "xlearner")
  expect_identical(xr$meta$spec$model_outcome, "ridge")
  expect_identical(xf$meta$spec$model_outcome, "rf")
  # Different base learners -> genuinely different OOF tau.
  expect_gt(max(abs(xr$res$tau_hat_oof - xf$res$tau_hat_oof), na.rm = TRUE), 1e-3)
  # cf_select returns the winning label + a reproducible spec.
  sel <- suppressMessages(cf_select(cmp, metric = "r_loss"))
  expect_true(sel$selected %in% names(cmp$models))
  expect_identical(sel$spec$method,
                   cmp$models[[sel$selected]]$method)
  # cf_ensemble spans the full grid.
  ens <- suppressWarnings(suppressMessages(cf_ensemble(cmp, method = "mean")))
  expect_setequal(ens$meta$methods_used,
                  c("grf", "xlearner_ridge", "xlearner_rf"))
  # predict dispatches on the bare method of a grid model.
  pr <- suppressWarnings(suppressMessages(predict(xf, type = "tau")))
  expect_length(pr, ncol(se))
  expect_true(all(is.finite(pr)))
})
