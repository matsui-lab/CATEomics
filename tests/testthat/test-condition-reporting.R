## How often a condition is signalled. The rule is one per call: a fold loop
## raises a fold's warning once per fold, but the caller asked for one fit.
## Nothing is capped by the clock or by the session, so a second call with the
## same problem says so again.

library(testthat)

.cr_tally <- function(expr) {
  # Ask rlang for every occurrence, so the count is the code's own and not
  # whatever an earlier test in this session already spent.
  op <- options(rlib_warning_verbosity = "verbose",
                rlib_message_verbosity = "verbose")
  on.exit(options(op), add = TRUE)
  hits <- character(0)
  keep <- function(cnd) {
    cl <- setdiff(class(cnd), c("warning", "message", "condition",
                                "simpleWarning", "simpleMessage",
                                "rlang_warning", "rlang_message"))
    if (length(cl)) hits <<- c(hits, cl[[1L]])
  }
  withCallingHandlers(
    force(expr),
    warning = function(c) { keep(c); invokeRestart("muffleWarning") },
    message = function(c) { keep(c); invokeRestart("muffleMessage") })
  table(hits)
}
.cr_n <- function(tab, cls) if (cls %in% names(tab)) as.integer(tab[[cls]]) else 0L

.cr_se <- function() cf_example_data(n = 60, p = 60, n_prog = 3, n_pred = 3,
                                     seed = 1)

test_that("a fold loop says a fold's caveat once, not once per fold", {
  # `model_final` is not an argument grf takes. The dispatcher sees that on
  # every fold and again on the full fit.
  tab <- .cr_tally(cf_fit(data = .cr_se(), Y = "outcome", T = "treatment",
                          X_baseline = c("age", "sex"), method = "grf",
                          n_features = 20, model_final = "ridge"))
  expect_equal(.cr_n(tab, "cateomics_control_arg_unreachable"), 1L)
})

test_that("a grid sweep says a caveat about its model once", {
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = .cr_se(), Y = "outcome", T = "treatment",
           X_baseline = c("age", "sex"), method = "grf", n_features = 20,
           partition = "none")))
  one  <- .cr_tally(cf_rule(fit, type = "threshold"))
  grid <- .cr_tally(cf_rule(fit, type = "threshold", grid = TRUE))
  expect_equal(.cr_n(one,  "cateomics_partition_none_warning"), 1L)
  # The sweep builds a rule per grid point from the same model.
  expect_equal(.cr_n(grid, "cateomics_partition_none_warning"), 1L)
})

test_that("the same problem in a second call is reported again", {
  call_once <- function() .cr_tally(
    cf_fit(data = .cr_se(), Y = "outcome", T = "treatment",
           X_baseline = c("age", "sex"), method = "grf", n_features = 20,
           model_final = "ridge"))
  expect_equal(.cr_n(call_once(), "cateomics_control_arg_unreachable"), 1L)
  expect_equal(.cr_n(call_once(), "cateomics_control_arg_unreachable"), 1L)
})

test_that("a different wording of one class is not collapsed away", {
  # The key is the class and the message together: one class covers several
  # messages, and keying on the class alone would keep only the last.
  mk <- function(cls, msg) structure(
    class = c(cls, "warning", "condition"),
    list(message = msg, call = NULL))
  a <- mk("cateomics_sl_zero_weights", "slot A")
  b <- mk("cateomics_sl_zero_weights", "slot B")
  expect_false(identical(.cf_condition_key(a), .cf_condition_key(b)))
  expect_identical(.cf_condition_key(a), .cf_condition_key(mk(
    "cateomics_sl_zero_weights", "slot A")))
  expect_identical(.cf_condition_class(a), "cateomics_sl_zero_weights")
  # A condition with no class of its own keys on its text.
  bare <- simpleWarning("plain")
  expect_identical(.cf_condition_class(bare), "")
})

test_that("cf_compare reports its candidates' conditions once, by candidate", {
  # No clinical covariate, so the propensity falls back and says so -- inside
  # the shared preprocessing, not inside a candidate.
  tab <- .cr_tally(cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
                              omics = "all", methods = "grf", n_features = 20))
  expect_equal(.cr_n(tab, "cateomics_propensity_X_fallback"), 1L)
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, model_final = "ridge")))
  # Raised on every fold of that candidate; kept once, named by candidate.
  expect_true("grf" %in% names(cmp$metadata$config_warnings))
  got <- cmp$metadata$config_warnings$grf
  expect_true("cateomics_control_arg_unreachable" %in% names(got))
  expect_equal(sum(names(got) == "cateomics_control_arg_unreachable"), 1L)
})

test_that("a comparison with nothing to report carries no warnings record", {
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_null(cmp$metadata$config_warnings)
})

test_that("a candidate's messages are reported apart from its warnings", {
  # A backend's progress output is a message. Calling it a warning tells the
  # reader a caveat was raised when none was.
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = c("grf", "bcf"),
               n_features = 20, n_folds = 3L, random_state = 0L,
               params = list(bcf = list(nburn = 20L, nsim = 20L,
                                        n_threads = 1L)))))
  msg <- cmp$metadata$config_messages
  expect_true("bcf" %in% names(msg))
  expect_true(any(grepl("BCF", msg$bcf, fixed = TRUE)))
  # The same strings must not also be filed as warnings.
  expect_false(any(grepl("Initializing BCF", cmp$metadata$config_warnings$bcf,
                         fixed = TRUE)))
})

test_that("cf_evaluate names every quantity it could not resample, once", {
  se  <- cf_example_data(n = 30, p = 25, n_prog = 2, n_pred = 2, seed = 1)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 15, random_state = 0L, n_folds = 2L,
                propensity_score = 0.5, eval_propensity_score = 0.5)
  tab <- .cr_tally(cf_evaluate(fit, n_bootstrap = 5L, random_state = 0L))
  expect_equal(.cr_n(tab, "cateomics_bootstrap_se_unavailable"), 1L)
  # And the one line says which quantities, which a warning raised from inside
  # the resampler could not do.
  w <- tryCatch(cf_evaluate(fit, n_bootstrap = 5L, random_state = 0L),
                warning = function(w) w)
  expect_s3_class(w, "cateomics_bootstrap_se_unavailable")
  expect_match(conditionMessage(w), "quantit")
})

test_that("cf_compare says once that a GPU request cannot be served", {
  tab <- .cr_tally(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, n_folds = 3L, random_state = 0L,
               device = "gpu", propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_equal(.cr_n(tab, "cateomics_gpu_unsupported_warning"), 1L)
  tab2 <- .cr_tally(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, n_folds = 3L, random_state = 0L,
               device = "cpu", propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_equal(.cr_n(tab2, "cateomics_gpu_unsupported_warning"), 0L)
})

test_that("the set-level control report is not repeated by the candidate line", {
  tab <- .cr_tally(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, n_folds = 3L, random_state = 0L,
               model_final = "ridge", propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_equal(.cr_n(tab, "cateomics_control_arg_unreachable"), 1L)
  # The per-candidate record is still kept, because candidates can differ.
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = .cr_se(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 20, n_folds = 3L, random_state = 0L,
               model_final = "ridge", propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_true("cateomics_control_arg_unreachable" %in%
                names(cmp$metadata$config_warnings$grf))
  expect_equal(.cr_n(tab, "cateomics_compare_config_warnings"), 0L)
})

test_that("a thin inner split is reported once for the call, learners named", {
  # One call fits a penalised learner once per fold, once per Super Learner
  # member and once per candidate; all of them are the same fact about the
  # same data, so the line names the learners and the smallest count.
  se <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  cd <- SummarizedExperiment::colData(se)
  df <- data.frame(Y = cd$outcome, treatment = cd$treatment, age = cd$age,
                   sex = cd$sex)
  tab <- .cr_tally(suppressMessages(
    cf_compare(data = df, Y = "Y", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf", n_folds = 3L,
               random_state = 0L, propensity_score = 0.5,
               eval_propensity_score = 0.5)))
  expect_equal(.cr_n(tab, "cateomics_cv_split_thin"), 1L)
  w <- NULL
  suppressMessages(withCallingHandlers(
    cf_compare(data = df, Y = "Y", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf", n_folds = 3L,
               random_state = 0L, propensity_score = 0.5,
               eval_propensity_score = 0.5),
    warning = function(cnd) {
      if (inherits(cnd, "cateomics_cv_split_thin")) w <<- cnd
      invokeRestart("muffleWarning")
    }))
  expect_match(conditionMessage(w), "penalised fit")
  expect_match(conditionMessage(w), "ridge")
  expect_match(conditionMessage(w), "lasso")

  # Called on its own, outside any public entry, it still says it at once and
  # in the singular form.
  one <- tryCatch(.cf_check_cv_split(5L, "lasso"), warning = function(x) x)
  expect_s3_class(one, "cateomics_cv_split_thin")
  expect_match(conditionMessage(one), "The lasso selects its penalty")
})
