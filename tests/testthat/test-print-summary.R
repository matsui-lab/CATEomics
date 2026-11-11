# print.cf_model / summary.cf_model / print.summary.cf_model

make_fitted <- function() {
  skip_if_not_installed("grf")
  set.seed(2026)
  data <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50)
  )
  cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
}

test_that("summary.cf_model returns a structured list", {
  fit <- make_fitted()
  s <- summary(fit)
  expect_s3_class(s, "summary.cf_model")
  # No estimand: every method here estimates a CATE, so the field recorded a
  # constant and the line it fed printed something the object did not hold.
  expect_named(s, c("method", "n", "p", "ate", "ate_ci",
                    "ci_method", "tau_summary"))
  expect_length(s$ate_ci, 2)
  expect_named(s$tau_summary, c("mean", "sd", "iqr"))
})

test_that("print.summary.cf_model produces formatted output", {
  fit <- make_fitted()
  s <- summary(fit)
  expect_output(print(s), "Causal Inference Model Summary")
  expect_output(print(s), "ATE 95% CI")
})

test_that("print.cf_model returns the object invisibly", {
  fit <- make_fitted()
  expect_output(p <- print(fit), "<cf_model:")
  expect_identical(p, fit)
})

test_that("summary.cf_model aborts when fit$res$summary is NULL", {
  fit <- make_fitted()
  fit$res$summary <- NULL
  expect_error(
    summary(fit),
    class = "cateomics_invalid_result"
  )
})


test_that("a cached fold fit prints its shape and can be summarised", {
  skip_if_not_installed("grf")
  # keep_kfold_models hands these objects to the caller, so they go through the
  # same print and summary. They carry no `data` pillar, so the shape has to be
  # in their own meta: without it sprintf() returned character(0) and the whole
  # shape line vanished from print, while summary stopped with "unsupported
  # type".
  set.seed(2)
  n <- 60L
  df <- data.frame(T = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  df$Y <- df$x2 + df$T * df$x1 + rnorm(n, sd = 0.3)
  fit <- suppressMessages(cf_fit(
    data = df, Y = "Y", T = "T", X_baseline = c("x1", "x2"), method = "grf",
    partition = "oof", n_folds = 3L, random_state = 1L,
    keep_kfold_models = TRUE,
    params = list(num.trees = 150L, tune.parameters = "none")))
  fold <- fit$fold_cache[[1L]]
  expect_true(all(c("n", "p") %in% names(fold$meta)))
  expect_output(print(fold), "N = [0-9]+, p = 2")
  expect_s3_class(summary(fold), "summary.cf_model")
  expect_output(print(summary(fold)), "Method")
})

test_that("a print survives a field an older object does not carry", {
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 15, random_state = 0L,
                propensity_score = 0.5, eval_propensity_score = 0.5)
  # `sprintf` fed a NULL returns character(0) and `cat()` prints nothing, so
  # the line used to disappear; `$` also partial-matches a longer sibling.
  older <- fit
  older$meta[["n"]] <- NULL
  older$meta[["p"]] <- NULL
  out <- capture.output(print(older))
  expect_true(any(grepl("N = NA, p = NA", out, fixed = TRUE)))
  sm <- capture.output(print(summary(older)))
  expect_true(any(grepl("N samples: NA", sm, fixed = TRUE)))

  ev <- suppressWarnings(cf_evaluate(fit, n_bootstrap = 25L, random_state = 0L))
  thin <- ev
  thin$rate$se <- NULL
  thin$blp$calibration$p_value <- NULL
  thin$metadata$n_bootstrap <- NULL
  eo <- capture.output(print(thin))
  expect_true(any(grepl("TOC @", eo, fixed = TRUE)))
  expect_true(any(grepl("calibration", eo, fixed = TRUE)))
  expect_true(any(grepl("intervals:", eo, fixed = TRUE)))
})

test_that("no print stops part way when a field is absent", {
  # A print says what it cannot find rather than dying: an object written by an
  # earlier version, or one a caller has edited, still has to be readable
  # enough to show what is wrong with it. This pins the shapes that used to
  # take the whole print down.
  se <- cf_example_data(n = 60, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cdd <- SummarizedExperiment::colData(se)
  df <- data.frame(Y = cdd$outcome, T = cdd$treatment, age = cdd$age,
                   sex = cdd$sex)
  q <- function(e) suppressWarnings(suppressMessages(e))
  fit <- q(cf_fit(data = df, Y = "Y", T = "T", X_baseline = c("age", "sex"),
                  method = "grf", random_state = 0L, propensity_score = 0.5,
                  eval_propensity_score = 0.5))
  cmp <- q(cf_compare(data = df, Y = "Y", T = "T",
                      X_baseline = c("age", "sex"),
                      configs = list(a = list(method = "grf"),
                                     b = list(method = "grf")),
                      n_folds = 3L, random_state = 0L, propensity_score = 0.5,
                      eval_propensity_score = 0.5))
  drop1 <- function(obj, path) {
    if (length(path) == 1L) obj[[path]] <- NULL
    else obj[[path[1L]]] <- drop1(obj[[path[1L]]], path[-1L])
    obj
  }
  cases <- list(
    list(summary(cmp), "scoreboard"),
    list(summary(cmp), c("metadata", "metric")),
    list(q(cf_select(cmp)), "selection_stats"),
    list(q(cf_select(cmp)), c("selection_stats", "runner_up")),
    list(q(cf_ensemble(cmp, method = "mean")), "meta"),
    list(q(cf_ensemble(cmp, method = "mean")), c("meta", "ensemble_method")),
    list(q(cf_rule(fit, type = "topk", n_bootstrap = 25L)),
         c("policy_value", "improvement_over_random", "estimate")),
    list(q(cf_rule(fit, type = "topk", grid = c(0.2, 0.4),
                   n_bootstrap = 25L)), "results"),
    list(q(cf_diagnostics(fit)), "smd"),
    list(q(cf_diagnostics(fit)), c("evalue", "ate"))
  )
  for (cs in cases) {
    broken <- drop1(cs[[1L]], cs[[2L]])
    expect_output(print(broken))
  }

  # And the interpretation print reads the surrogates it holds, not the index
  # of them, so a damaged index cannot make it show less than the object does.
  itp <- q(cf_interpret(fit, surrogate = "ridge", random_state = 0L))
  bent <- itp
  bent$metadata$surrogates <- NA_character_
  expect_equal(length(capture.output(print(bent))),
               length(capture.output(print(itp))))
})
