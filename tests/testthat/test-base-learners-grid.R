# Tests for cf_compare(base_learners=): grid over (method x base learner).
# The pluggable meta/DR learners fan out to one config per base learner
# (label = method_base); grf / bcf / NN stay single. `label` is the uniform
# external key; each config keeps its bare `method` field. base_learners = NULL
# reproduces the legacy single-config behaviour.

grid_demo <- function(n = 70L, p = 45L, seed = 3L) {
  cf_example_data(n = n, p = p, n_prog = 8L, n_pred = 4L, seed = seed)
}

pluggable <- c("slearner", "tlearner", "xlearner", "drlearner")


test_that(".cateomics_build_configs expands only pluggable methods", {
  # NULL = one row per method, base NA, label == method.
  cfg0 <- CATEomics:::.cateomics_build_configs(c("grf", "xlearner"), NULL, pluggable)
  expect_identical(cfg0$label, c("grf", "xlearner"))
  expect_identical(cfg0$method, c("grf", "xlearner"))
  expect_true(all(is.na(cfg0$base)))

  # Grid = pluggable methods fan out; non-pluggable stay single (base NA).
  cfg <- CATEomics:::.cateomics_build_configs(
    c("grf", "xlearner", "drlearner"), c("ridge", "rf"), pluggable)
  expect_identical(
    cfg$label,
    c("grf", "xlearner_ridge", "xlearner_rf", "drlearner_ridge", "drlearner_rf"))
  expect_identical(cfg$method,
                   c("grf", "xlearner", "xlearner", "drlearner", "drlearner"))
  expect_identical(cfg$base, c(NA, "ridge", "rf", "ridge", "rf"))
})


test_that("base_learners = NULL keys by method and leaves base empty", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  se <- grid_demo()
  cmp <- suppressWarnings(suppressMessages(cf_compare(data = se, Y = "outcome", T = "treatment",
          X_baseline = c("age", "sex"), methods = c("grf",
              "bcf"), n_features = 20L, propensity_score = 0.5,
          random_state = 1, verbose = FALSE, params = list(bcf = list(nburn = 5L,
              nsim = 5L, n_chains = 1L)))))
  # The decomposition is always present: `label` is the key, `method` the bare
  # method, `base` the base learner that actually ran. grf and bcf have no
  # pluggable outcome slot, so `base` is NA for both whether or not a grid was
  # asked for.
  expect_true(all(c("label", "method", "base") %in% colnames(cmp$scoreboard)))
  expect_identical(cmp$scoreboard$label, cmp$scoreboard$method)
  expect_true(all(is.na(cmp$scoreboard$base)))
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
  cmp <- suppressWarnings(suppressMessages(cf_compare(data = se, Y = "outcome", T = "treatment",
          X_baseline = c("age", "sex"), methods = c("grf",
              "bcf"), base_learners = c("ridge", "rf"), n_features = 20L,
          propensity_score = 0.5, random_state = 1, verbose = FALSE,
          params = list(bcf = list(nburn = 5L, nsim = 5L,
              n_chains = 1L)))))
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
                 n_features = 20L, propensity_score = 0.5, random_state = 1,
                 verbose = FALSE,
                 params = list(bcf = list(nburn = 5L, nsim = 5L,
                                          n_chains = 1L)))
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
    class = "cateomics_base_learners_conflict")
})


test_that("unknown base_learners abort; all-singleton set informs no effect", {
  se <- grid_demo(n = 40L, p = 20L)
  expect_error(
    suppressWarnings(suppressMessages(cf_compare(
      data = se, Y = "outcome", T = "treatment", X_baseline = "age",
      methods = "grf", base_learners = "not_a_learner", verbose = FALSE))),
    class = "cateomics_base_learners_unknown")
  skip_if_not_installed("grf")
  expect_message(
    suppressWarnings(cf_compare(
      data = se, Y = "outcome", T = "treatment", X_baseline = "age",
      methods = "grf", base_learners = "rf", n_features = 15L,
      propensity_score = 0.5, random_state = 1, verbose = FALSE)),
    class = "cateomics_base_learners_no_pluggable")
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
  # The config behind the winner is reachable through the model it returns;
  # `cf_select` does not keep a second copy of it.
  expect_identical(sel$model$meta$spec$method,
                   cmp$models[[sel$selected]]$method)
  expect_null(sel[["spec"]])
  # cf_ensemble spans the full grid.
  ens <- suppressWarnings(suppressMessages(cf_ensemble(cmp, method = "mean")))
  expect_setequal(ens$meta$methods_used,
                  c("grf", "xlearner_ridge", "xlearner_rf"))
  # predict dispatches on the bare method of a grid model.
  pr <- suppressWarnings(suppressMessages(predict(xf, type = "tau")))
  expect_length(pr, ncol(se))
  expect_true(all(is.finite(pr)))
})


# ---------------------------------------------------------------------------
# cf_compare(configs = ): naming the candidate set outright
# ---------------------------------------------------------------------------

test_that("the shortcut arguments and configs resolve to the same table", {
  # One internal representation: whatever the caller wrote, everything
  # downstream reads the table this returns.
  sugar <- CATEomics:::.cf_resolve_comparison(
    methods = c("grf", "xlearner"), base_learners = NULL, configs = NULL,
    params = NULL, model_outcome = "rf", model_final = "rf",
    pluggable_methods = pluggable,
    cate_methods_all = c("grf", "bcf", "slearner", "tlearner", "xlearner",
                         "drlearner", "cevae", "ganite", "tarnet", "dragonnet"))
  explicit <- CATEomics:::.cf_resolve_comparison(
    methods = NULL, base_learners = NULL,
    configs = list(grf = list(method = "grf"),
                   xlearner = list(method = "xlearner", model_outcome = "rf",
                                   model_final = "rf")),
    params = NULL, model_outcome = NULL, model_final = NULL,
    pluggable_methods = pluggable,
    cate_methods_all = c("grf", "bcf", "slearner", "tlearner", "xlearner",
                         "drlearner", "cevae", "ganite", "tarnet", "dragonnet"))
  expect_identical(sugar$table, explicit$table)
  # grf has no pluggable outcome slot, so naming one for it is unreachable;
  # the shortcut form applies it only where it lands.
  expect_identical(sugar$specs[[2]]$model_outcome, "rf")
  expect_identical(explicit$specs[[2]]$model_outcome, "rf")
})

test_that("the resolved table records the learner that actually runs", {
  res <- CATEomics:::.cf_resolve_comparison(
    methods = c("grf", "xlearner", "drlearner"), base_learners = NULL,
    configs = NULL, params = NULL, model_outcome = NULL, model_final = NULL,
    pluggable_methods = pluggable,
    cate_methods_all = c("grf", "bcf", "slearner", "tlearner", "xlearner",
                         "drlearner", "cevae", "ganite", "tarnet", "dragonnet"))
  # Nothing was asked for, so `base` carries each fitter's own default rather
  # than NA -- the request and the effective value are not the same thing.
  expect_true(is.na(res$table$base[res$table$method == "grf"]))
  expect_identical(res$table$base[res$table$method == "xlearner"], "ridge")
  expect_identical(res$table$base[res$table$method == "drlearner"], "ridge")
})

test_that("configs is mutually exclusive with the arguments that build one", {
  d <- cf_example_data(n = 40L, p = 12L, n_prog = 4L, n_pred = 2L, seed = 5L)
  cfgs <- list(grf = list(method = "grf"))
  for (extra in list(list(methods = "grf"),
                     list(base_learners = "ridge"),
                     list(model_outcome = "ridge"),
                     list(model_final = "ridge"),
                     list(params = list(grf = list(num.trees = 10L))))) {
    expect_error(
      do.call(cf_compare, c(list(data = d, Y = "outcome", T = "treatment",
                                 X_baseline = "age", configs = cfgs,
                                 verbose = FALSE), extra)),
      class = "cateomics_configs_conflict"
    )
  }
  expect_error(
    cf_compare(data = d, Y = "outcome", T = "treatment",
               X_baseline = "age", verbose = FALSE),
    class = "cateomics_methods_required"
  )
})

test_that("configs validates its labels, methods and entries", {
  d <- cf_example_data(n = 40L, p = 12L, n_prog = 4L, n_pred = 2L, seed = 5L)
  bad <- function(cfgs) {
    expect_error(
      cf_compare(data = d, Y = "outcome", T = "treatment",
                 X_baseline = "age", configs = cfgs, verbose = FALSE),
      class = "cateomics_configs_invalid"
    )
  }
  bad(list(list(method = "grf")))                     # unnamed
  bad(stats::setNames(list(list(method = "grf")), "")) # empty name
  bad(list(a = list(method = "grf"), a = list(method = "bcf")))  # duplicate
  bad(list(a = "grf"))                                 # not a list
  bad(list(a = list(method = "nope")))                 # unknown method
  bad(list(a = list(method = "grf", n_folds = 3L)))    # unknown entry
  expect_error(
    cf_compare(data = d, Y = "outcome", T = "treatment",
               X_baseline = "age", verbose = FALSE,
               configs = list(a = list(method = "xlearner",
                                        model_outcome = "nope"))),
    class = "cateomics_base_learner_error"
  )
})

test_that("the entry guards run on the resolved configs, not just on methods", {
  d <- cf_example_data(n = 40L, p = 12L, n_prog = 4L, n_pred = 2L, seed = 5L)
  # grf has no `bootstrap_ci` formal, so the request is refused at the entry.
  # The guards read the resolved config table, so naming the candidate set
  # either way has to be refused the same: before the table was built first,
  # the configs path reached the fits with the request silently dropped.
  expect_error(
    cf_compare(data = d, Y = "outcome", T = "treatment", X_baseline = "age",
               methods = "grf", bootstrap_ci = TRUE, verbose = FALSE),
    class = "cateomics_bootstrap_unsupported")
  expect_error(
    cf_compare(data = d, Y = "outcome", T = "treatment", X_baseline = "age",
               configs = list(g = list(method = "grf")), bootstrap_ci = TRUE,
               verbose = FALSE),
    class = "cateomics_bootstrap_unsupported")
})

test_that("the recorded spec names only the learners that actually run", {
  # One resolver serves both entry points. A slot that no learner fills stays
  # NA: writing the fitter's default there would name a learner that never ran.
  slots <- function(method, model_outcome = NULL, model_final = NULL,
                    model_propensity = NULL, propensity_score = NULL,
                    clip = c(0.025, 0.975)) {
    CATEomics:::.cf_spec_slots(method, model_outcome, model_final,
                             model_propensity, propensity_score, clip)
  }

  # grf centres with its own forest: nothing is pluggable, but it does clip.
  g <- slots("grf")
  expect_true(is.na(g$model_outcome))
  expect_true(is.na(g$model_final))
  expect_true(is.na(g$model_propensity))
  expect_identical(g$clip, c(0.025, 0.975))

  # Defaults live on each fitter, not on the public formals, so an unnamed
  # slot still records the value that runs.
  x <- slots("xlearner")
  expect_identical(x$model_outcome, "ridge")
  expect_identical(x$model_final, "ridge")
  expect_identical(x$model_propensity, "logistic")
  expect_identical(slots("bcf")$model_propensity, "logistic")

  # A known propensity score means no propensity learner is fitted at all.
  expect_true(is.na(slots("xlearner", propensity_score = 0.5)$model_propensity))

  # `clip` bounds a propensity score, so a method without that formal must not
  # record bounds it never applied.
  for (m in c("slearner", "tlearner", "tarnet", "dragonnet")) {
    expect_true(is.na(slots(m)$clip), info = m)
  }
})
