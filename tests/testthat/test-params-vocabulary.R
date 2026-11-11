# The vocabulary of `params`: which names are accepted, which CATEomics keeps
# for itself, and that the seed is not one of the accepted ones on any route.

test_that("a random seed cannot be set through params on any route", {
  # grf and bcf expose their own seed argument, the scikit-learn estimators
  # behind the meta-learners expose `random_state`, and the neural methods
  # forward whatever they are given to their Python constructor. All of them
  # would have taken over from `random_state` for that one fitter.
  expect_error(.cf_params_check("grf", list(seed = 1L)),
               class = "cateomics_rng_param")
  expect_error(.cf_params_check("bcf", list(random_seed = 1L)),
               class = "cateomics_rng_param")
  expect_error(.cf_params_check("tlearner", list(outcome = list(random_state = 1L))),
               class = "cateomics_rng_param")
  for (m in c("cevae", "ganite", "tarnet", "dragonnet")) {
    expect_error(.cf_params_check(m, list(seed = 1L)),
                 class = "cateomics_rng_param")
  }
  expect_error(.cf_learner_params_check("rf", list(seed = 1L), "eval_params"),
               class = "cateomics_rng_param")
  expect_error(.cf_learner_params_check("xgboost", list(seed = 1L), "params$xgb"),
               class = "cateomics_rng_param")
})

test_that("params keys carry the upstream package's spelling", {
  skip_if_not_installed("grf")
  skip_if_not_installed("bcf")
  # grf's own dotted names, and bcf's own names rather than renamed ones.
  expect_no_error(.cf_params_check("grf", list(num.trees = 300L,
                                               min.node.size = 10L)))
  expect_no_error(.cf_params_check("bcf", list(nburn = 20L, nsim = 20L)))
  # `n_chains` has no upstream counterpart: CATEomics runs the chains itself.
  expect_no_error(.cf_params_check("bcf", list(n_chains = 1L)))
  # The names CATEomics used to rename them to are gone.
  expect_error(.cf_params_check("bcf", list(n_burn = 20L)),
               class = "cateomics_unknown_param")
  expect_error(.cf_params_check("bcf", list(n_iter = 20L)),
               class = "cateomics_unknown_param")
  # The orthogonalization predictions stay closed: CATEomics fits them and
  # passes them positionally, and a known propensity has its own argument.
  expect_error(.cf_params_check("grf", list(Y.hat = rep(0, 4))),
               class = "cateomics_unknown_param")
  expect_error(.cf_params_check("grf", list(W.hat = rep(0.5, 4))),
               class = "cateomics_unknown_param")
})

test_that("the linear learner refuses the one name that would break it", {
  # `qr` is the only name lm takes and glm does not, and it breaks both
  # branches of `.fit_linear()`: glm forwards it to glm.control(), which has
  # no `...`, and a regression fitted with qr = FALSE fails in predict.lm().
  expect_error(.cf_learner_params_check("linear", list(qr = TRUE),
                                        "eval_params"),
               class = "cateomics_unknown_param")
  expect_no_error(.cf_learner_params_check("linear", list(singular.ok = TRUE),
                                           "eval_params"))
  expect_no_error(.cf_learner_params_check("logistic", list(epsilon = 1e-9),
                                           "eval_params"))
  # The data and the family stay CATEomics's.
  expect_error(.cf_learner_params_check("linear", list(family = "x"),
                                        "eval_params"),
               class = "cateomics_unknown_param")
})

test_that("xgboost's two spellings of the learning rate do not collide", {
  skip_if_not_installed("xgboost")
  set.seed(8)
  X <- matrix(rnorm(80 * 6), 80, 6)
  y <- X[, 1] + rnorm(80)
  rate <- function(m) {
    as.numeric(xgboost::xgb.config(m$model)$learner$gradient_booster$tree_train_param$eta)
  }
  # The native spelling used to be accepted by the validator and then
  # overridden by the default, because both names reached the same list.
  expect_equal(rate(.fit_base_learner(X, y, learner = "xgboost",
                                      random_state = 1L, eta = 0.9)),
               0.9, tolerance = 1e-6)
  expect_equal(rate(.fit_base_learner(X, y, learner = "xgboost",
                                      random_state = 1L, learning_rate = 0.9)),
               0.9, tolerance = 1e-6)
  expect_error(.fit_base_learner(X, y, learner = "xgboost", random_state = 1L,
                                 eta = 0.9, learning_rate = 0.1),
               class = "cateomics_duplicate_param")
})

test_that("the xgboost base learner carries the SuperLearner wrapper's shape", {
  skip_if_not_installed("xgboost")
  set.seed(4)
  X <- matrix(rnorm(80 * 8), 80, 8)
  y <- X[, 1] + rnorm(80)
  tp <- function(m) {
    xgboost::xgb.config(m$model)$learner$gradient_booster$tree_train_param
  }
  d <- tp(.fit_base_learner(X, y, learner = "xgboost", random_state = 1L))
  expect_equal(as.integer(d$max_depth), 4L)
  expect_equal(as.numeric(d$eta), 0.1, tolerance = 1e-6)
  expect_equal(as.numeric(d$min_child_weight), 10)
  expect_equal(as.numeric(d$colsample_bytree), 0.5, tolerance = 1e-6)
  # And every one of them is still the caller's to change.
  o <- tp(.fit_base_learner(X, y, learner = "xgboost", random_state = 1L,
                            max_depth = 2L, colsample_bytree = 1,
                            min_child_weight = 1L, learning_rate = 0.3))
  expect_equal(as.integer(o$max_depth), 2L)
  expect_equal(as.numeric(o$colsample_bytree), 1, tolerance = 1e-6)
  expect_equal(as.numeric(o$min_child_weight), 1)
  expect_equal(as.numeric(o$eta), 0.3, tolerance = 1e-6)
})

test_that("the Python factory resolves xgboost's native parameter names", {
  skip_if_no_cf_python()
  # The R side derives its accepted names from xgb.params(), which lists the
  # native spelling; `get_params()` on the scikit-learn wrapper lists only the
  # wrapper spelling, so the same key used to be accepted on one side of the
  # package and refused on the other.
  got <- .cf_py_run(function(pypath) {
    import_sys <- reticulate::import("sys", convert = FALSE)
    import_sys$path$insert(0L, pypath)
    lf <- reticulate::import("learner_factory", convert = TRUE)
    xgb <- reticulate::import("xgboost", convert = TRUE)
    out <- list()
    for (pair in list(c("eta", "learning_rate"), c("min_split_loss", "gamma"),
                      c("lambda", "reg_lambda"), c("alpha", "reg_alpha"))) {
      est <- lf$`_apply_params`(xgb$XGBRegressor(), "xgboost",
                                 stats::setNames(list(0.25), pair[1L]))
      out[[pair[1L]]] <- est$get_params()[[pair[2L]]]
    }
    out
  }, pypath = system.file("python", package = "CATEomics"))
  expect_equal(unlist(got), c(eta = 0.25, min_split_loss = 0.25,
                              lambda = 0.25, alpha = 0.25))
})

test_that("the default seed reproduces every R-side learner", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("xgboost")
  # `random_state = 0L` is the package default, and ranger reads `seed = 0` as
  # "no seed": passing the value through left the default evaluation nuisance
  # following whatever RNG state the caller was in. The ambient state is moved
  # between the two fits, which is what a second session looks like.
  set.seed(1)
  n <- 150L; p <- 8L
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", seq_len(p))))
  y <- X[, 1] + rnorm(n)
  yb <- as.integer(y > 0)
  twice <- function(learner, type, target) {
    f <- function() {
      m <- .fit_base_learner(X, target, type = type, learner = learner,
                             random_state = 0L)
      as.numeric(.predict_base_learner(m, X))
    }
    set.seed(111); a <- f()
    set.seed(222); b <- f()
    max(abs(a - b))
  }
  for (l in c("rf", "xgboost", "ridge", "lasso")) {
    expect_equal(twice(l, "regression", y), 0,
                 info = paste(l, "regression"))
  }
  expect_equal(twice("rf", "classification", yb), 0)
  # And the seed still does something: a different one has to give a
  # different forest, or the first assertion would hold trivially.
  set.seed(111)
  m0 <- .fit_base_learner(X, y, learner = "rf", random_state = 0L)
  set.seed(111)
  m1 <- .fit_base_learner(X, y, learner = "rf", random_state = 1L)
  expect_gt(max(abs(.predict_base_learner(m0, X) - .predict_base_learner(m1, X))), 1e-8)
})


# ---- catenets validation split: the R formal and the Python signature agree --

test_that("val_split_prop is an argument of both catenets fitters", {
  # The split catenets takes changes how many rows train the weights, so the
  # default is an argument rather than something reachable only through
  # `params`. 0.3 is catenets' own default, so the fits are unchanged.
  for (fitter in list(CATEomics:::cf_fit_tarnet, CATEomics:::cf_fit_dragonnet)) {
    fo <- formals(fitter)
    expect_true("val_split_prop" %in% names(fo))
    expect_equal(eval(fo$val_split_prop), 0.3)
  }
})

test_that("the catenets wrapper declares val_split_prop and forwards it", {
  # The R side names the argument explicitly, so the Python side must take it
  # rather than swallow it into **kwargs: this reads the file instead of
  # trusting the R call.
  f <- system.file("python", "catenets_wrapper.py", package = "CATEomics")
  if (!nzchar(f)) f <- "../../inst/python/catenets_wrapper.py"
  skip_if_not(file.exists(f))
  src <- readLines(f, warn = FALSE)
  # One declaration per fitter (fit_tarnet / fit_dragonnet) ...
  expect_length(grep("^\\s*val_split_prop=0\\.3,\\s*$", src), 2L)
  # ... and one hand-off into the constructor arguments for each.
  expect_length(grep("val_split_prop=float\\(val_split_prop\\)", src), 2L)
})

# A hyperparameter the Python side requires to be whole has to arrive as an
# integer. R writes `10` as a double, so the natural literal is refused, and
# the refusal reads as an upstream type error or, through torch, as an error
# about a tensor. Neither names the fix.

test_that("the scan finds the whole-number doubles and nothing else", {
  hits <- .cf_whole_doubles(list(
    # Not lists, so never scanned: these are the data, not hyperparameters.
    X = matrix(1, 3, 3), Y = c(1, 2, 3), T = c(0, 1, 0),
    clip = c(0.025, 0.975),
    params = list(outcome = list(cv = 10, eps = 0.01,
                                 alphas = c(0.1, 1, 10),
                                 max_iter = 2000L)),
    extra = list(n_units_out = 8)))
  expect_setequal(hits, c("outcome$cv = 10", "n_units_out = 8"))
  # A block with nothing to report gives an empty vector, not NULL.
  expect_identical(.cf_whole_doubles(list(params = list(eps = 0.01))),
                   character(0))
  expect_identical(.cf_whole_doubles(list()), character(0))
})

test_that("a double where Python wants an integer is named as such", {
  skip_if_no_cf_python()
  se <- cf_example_data(n = 60, p = 80, n_prog = 4, n_pred = 4, seed = 1)
  fit_with <- function(v) {
    suppressWarnings(suppressMessages(
      cf_fit(data = se, Y = "outcome", T = "treatment",
             X_baseline = c("age", "sex"), method = "tlearner",
             n_features = 40, model_outcome = "lasso",
             params = list(tlearner = list(outcome = list(cv = v))))))
  }
  err <- expect_error(fit_with(10), class = "cateomics_python_param_type")
  expect_match(conditionMessage(err), "10L", fixed = TRUE)
  expect_match(conditionMessage(err), "cv = 10", fixed = TRUE)
  # The same call with an integer runs, so the message is about the type and
  # not about the name.
  expect_s3_class(fit_with(10L), "cf_model")
})

test_that("an unrelated Python failure keeps its own condition", {
  skip_if_no_cf_python()
  se <- cf_example_data(n = 60, p = 80, n_prog = 4, n_pred = 4, seed = 1)
  err <- expect_error(suppressWarnings(suppressMessages(
    cf_fit(data = se, Y = "outcome", T = "treatment",
           X_baseline = c("age", "sex"), method = "tlearner",
           n_features = 40, model_outcome = "lasso",
           params = list(tlearner = list(outcome = list(nope = 1L)))))))
  expect_false(inherits(err, "cateomics_python_param_type"))
  expect_match(conditionMessage(err), "nope")
})
