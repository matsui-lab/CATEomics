# Tests for input validation and error handling

test_that("cf_fit rejects non-numeric outcome", {
  data <- data.frame(
    Y = letters[1:20],
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "numeric"
  )
})

test_that("cf_fit rejects non-binary treatment", {
  data <- data.frame(
    Y = rnorm(20),
    T = 1:20,
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "binary"
  )
})

test_that("cf_fit rejects treatment with values other than 0/1", {
  data <- data.frame(
    Y = rnorm(20),
    T = sample(c(1, 2), 20, replace = TRUE),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "binary"
  )
})

test_that("cf_fit accepts logical treatment and converts to 0/1", {
  set.seed(42)
  data <- data.frame(
    Y = rnorm(50),
    T = sample(c(TRUE, FALSE), 50, replace = TRUE),
    X1 = rnorm(50)
  )

  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
})

test_that("cf_fit auto-expands factor / character covariates", {
  # factor / character covariates are auto one-hot encoded via
  # model.matrix(~ . - 1). `.expand_factor_columns` records the train levels in
  # `meta$factor_levels`, and cf_predict.R re-applies them to newdata. The
  # fixture uses a 3-level factor (a 20-level factor would be rank-deficient at
  # n = 20).
  set.seed(42)
  data <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    X1 = sample(c("a", "b", "c"), 40, replace = TRUE),
    stringsAsFactors = FALSE
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
  expect_true(!is.null(fit$meta$factor_levels))
  expect_true("X1" %in% names(fit$meta$factor_levels))
  expect_setequal(fit$meta$factor_levels$X1, c("a", "b", "c"))
})

test_that("cf_fit handles NA values in outcome", {
  data <- data.frame(
    Y = c(NA, rnorm(19)),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  # NA in outcome is detected by parse_cf_formula and aborts.
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_na_in_variable"
  )
})

test_that("cf_fit handles NA values in treatment", {
  data <- data.frame(
    Y = rnorm(20),
    T = c(NA, rbinom(19, 1, 0.5)),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_na_in_variable"
  )
})

test_that("cf_fit handles NA values in covariates", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = c(NA, rnorm(19))
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_na_in_variable"
  )
})

test_that("cf_fit handles empty data", {
  data <- data.frame(
    Y = numeric(0),
    T = integer(0),
    X1 = numeric(0)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_no_observations"
  )
})

test_that("cf_fit handles single observation", {
  data <- data.frame(
    Y = 1.0,
    T = 1L,
    X1 = 0.5
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_insufficient_n"
  )
})

test_that("cf_fit handles single treatment group (all treated)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(1L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_group_empty"
  )
})

test_that("cf_fit handles single treatment group (all control)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(0L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    class = "cateomics_group_empty"
  )
})

test_that("cf_fit rejects unknown method", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "unknown_method"),
    "arg"
  )
})


test_that("parse_cf_formula extracts correct components", {
  set.seed(123)
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = rnorm(50),
    X2 = rnorm(50)
  )

  result <- parse_cf_formula(Y ~ T | X1 + X2, data)

  expect_equal(result$outcome_name, "Y")
  expect_equal(result$treatment_name, "T")
  expect_equal(result$covariate_names, c("X1", "X2"))
  expect_equal(length(result$Y), 50)
  expect_equal(length(result$T), 50)
  expect_equal(ncol(result$X), 2)
  expect_true(is.numeric(result$Y))
  expect_true(is.integer(result$T))
  expect_true(is.matrix(result$X))
})


test_that("cf_fit handles data with many covariates (p >> n)", {
  set.seed(42)
  # n=150 chosen so GRF honest splitting (B9 default) has enough rows per leaf
  # after tune_forest's internal subsampling. Smaller n triggers grf's
  # "honesty fraction too close to 1 or 0" error.
  n <- 150
  p <- 100
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("X", 1:p)

  data <- as.data.frame(X)
  data$Y <- rnorm(n)
  data$T <- rbinom(n, 1, 0.5)

  skip_if_not_installed("grf")

  fit <- cf_fit(data = data, Y = "Y", T = "T", X_baseline = paste0("X", 1:p),
                method = "grf")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$meta$p, p)
})

test_that("cf_fit rejects very small sample (n=5) at formula parser", {
  # parse_cf_formula enforces a minimum of 10 observations at the API boundary
  # rather than delegating to downstream methods with tiny n.
  set.seed(123)
  data <- data.frame(
    Y = rnorm(5),
    T = c(0L, 0L, 1L, 1L, 1L),
    X1 = rnorm(5)
  )

  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", partition = "none"),
    "Insufficient sample size"
  )
})

test_that("cf_fit handles zero-variance covariate", {
  set.seed(42)
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = rep(1, 50),  # constant
    X2 = rnorm(50)
  )

  # Should still fit (method may or may not use the constant variable)
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
})

test_that("`method` has no automatic choice", {
  # Resolving the estimator from what happens to be installed made the same
  # script run a different method on a different machine, and choosing one
  # silently is the opposite of what this package is for.
  expect_false(exists("select_best_method",
                      envir = asNamespace("CATEomics"), inherits = FALSE))
  expect_false("auto" %in% eval(formals(cf_fit)$method))
})

# Allow dot-containing variable names (`gene.A`, `snp.rs123`, ...).
test_that("cf_fit accepts dot-containing covariate names", {
  set.seed(42)
  data <- data.frame(
    Y       = rnorm(40),
    T       = rbinom(40, 1, 0.5),
    gene.A  = rnorm(40),
    gene.B  = rnorm(40)
  )
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("gene.A", "gene.B"),
                data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
  expect_setequal(fit$meta$covariate_names, c("gene.A", "gene.B"))
})

# Y type advisory: binary 0/1 Y is a supported mode (risk difference), so it
# is not warned about; count-like Y is flagged as treated-as-continuous.
test_that("parse_cf_formula does not warn for binary 0/1 Y", {
  set.seed(7)
  data <- data.frame(
    Y  = rbinom(40, 1, 0.5),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40)
  )
  expect_no_warning(
    parse_cf_formula(Y ~ T | X1, data),
    class = "cateomics_y_type_warning"
  )
})

test_that("parse_cf_formula warns when Y is count-like", {
  set.seed(7)
  data <- data.frame(
    Y  = rpois(40, lambda = 3),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40)
  )
  expect_warning(
    parse_cf_formula(Y ~ T | X1, data),
    class = "cateomics_y_type_warning"
  )
})

test_that("parse_cf_formula does not warn for continuous Y", {
  set.seed(7)
  data <- data.frame(
    Y  = rnorm(40),
    T  = rbinom(40, 1, 0.5),
    X1 = rnorm(40)
  )
  expect_no_warning(
    parse_cf_formula(Y ~ T | X1, data)
  )
})

# Validator failures abort (no tryCatch demotion).
test_that("validate_cateomics_result aborts on missing summary fields", {
  bad <- structure(
    list(
      method = "grf",
      res    = list(
        ate     = 0.5,
        tau_hat = rnorm(10),
        mu0_hat  = rnorm(10),
        mu1_hat  = rnorm(10),
        summary = list(ate = 0.5, ate_ci_lower = 0, ate_ci_upper = 1)
      ),
      meta   = list(
        n = 10L, p = 1L,
        covariate_names  = "X1",
        named_args       = list(Y = "Y", T = "T", X_baseline = "X1"),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cateomics_result(bad),
    class = "cateomics_invalid_result"
  )
})

test_that("validate_cateomics_result aborts on malformed ite_quantiles", {
  bad <- structure(
    list(
      method = "grf",
      res    = list(
        ate     = 0.5,
        tau_hat = rnorm(10),
        mu0_hat  = rnorm(10),
        mu1_hat  = rnorm(10),
        summary = list(
          ate = 0.5, ate_ci_lower = 0, ate_ci_upper = 1,
          tau_mean = 0.5, tau_sd = 1.0,
          tau_iqr = list(q0 = 0, q50 = 0.5, q100 = 1)
        )
      ),
      meta = list(
        n = 10L, p = 1L,
        covariate_names  = "X1",
        named_args       = list(Y = "Y", T = "T", X_baseline = "X1"),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cateomics_result(bad),
    class = "cateomics_invalid_result"
  )
})

# ate_ci / ite_ci unified location validator.
test_that("validate_cateomics_result accepts NULL ate_ci / ite_ci (additive, optional)", {
  ## The validator requires `n`, `p`, `covariate_names`, `named_args` and the
  ## two buckets (`provenance`, `estimator`), so the synthetic `ok` object
  ## carries them. This still verifies that NULL `ate_ci` / `ite_ci` paths
  ## remain accepted.
  ok <- structure(
    list(
      method = "grf",
      res    = list(
        ate     = 0.5,
        tau_hat = rnorm(10),
        mu0_hat  = rnorm(10),
        mu1_hat  = rnorm(10),
        summary = list(
          ate = 0.5, ate_ci_lower = 0, ate_ci_upper = 1,
          tau_mean = 0.5, tau_sd = 1.0,
          tau_iqr = list(q0 = 0, q25 = 0.25, q50 = 0.5, q75 = 0.75, q100 = 1)
        )
      ),
      meta = list(
        n = 10L, p = 1L,
        covariate_names = "X1",
        named_args       = list(Y = "Y", T = "T", X_baseline = "X1"),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_true(validate_cateomics_result(ok))
})

test_that("validate_cateomics_result rejects malformed ate_ci", {
  bad <- structure(
    list(
      method = "grf",
      res    = list(
        ate     = 0.5,
        ate_ci  = c(1.0, 0.0),
        tau_hat = rnorm(10),
        mu0_hat  = rnorm(10),
        mu1_hat  = rnorm(10),
        summary = list(
          ate = 0.5, ate_ci_lower = 0, ate_ci_upper = 1,
          tau_mean = 0.5, tau_sd = 1.0,
          tau_iqr = list(q0 = 0, q25 = 0.25, q50 = 0.5, q75 = 0.75, q100 = 1)
        )
      ),
      meta = list(
        n = 10L, p = 1L,
        covariate_names  = "X1",
        named_args       = list(Y = "Y", T = "T", X_baseline = "X1"),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cateomics_result(bad),
    class = "cateomics_invalid_result"
  )
})

test_that("validate_cateomics_result rejects malformed ite_ci nrow", {
  bad <- structure(
    list(
      method = "grf",
      res    = list(
        ate     = 0.5,
        tau_hat = rnorm(10),
        tau_ci  = data.frame(lower = rnorm(5), upper = rnorm(5)),
        mu0_hat  = rnorm(10),
        mu1_hat  = rnorm(10),
        summary = list(
          ate = 0.5, ate_ci_lower = 0, ate_ci_upper = 1,
          tau_mean = 0.5, tau_sd = 1.0,
          tau_iqr = list(q0 = 0, q25 = 0.25, q50 = 0.5, q75 = 0.75, q100 = 1)
        )
      ),
      meta = list(
        n = 10L, p = 1L,
        covariate_names  = "X1",
        named_args       = list(Y = "Y", T = "T", X_baseline = "X1"),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cateomics_result(bad),
    class = "cateomics_invalid_result"
  )
})

# --- treatment groups too small for the fold structure (2026-10-08) --------
# Left unchecked this produced three different outcomes: the evaluation
# nuisance returned a hard zero for the affected rows and drove the
# doubly-robust pseudo-outcome to a mean of -1036 where the true effect was 0,
# grf returned a plausible-looking number, and bcf stopped with an unclassed
# glmnet message.

small_arm_df <- function(m_control, n = 40L, p = 4L, seed = 77L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  data.frame(y = 5 + X[, 1] + stats::rnorm(n),
             t = c(rep(0L, m_control), rep(1L, n - m_control)), X)
}
fit_small <- function(df, ...) {
  suppressWarnings(cf_fit(data = df, Y = "y", T = "t",
                          X_baseline = paste0("g", 1:4), method = "grf",
                          preprocess = FALSE, random_state = 0L, ...))
}

test_that("cross-fitting refuses a group that leaves under two rows to train on", {
  # 2 controls: a fold holds out one, leaving one. No n_folds can fix that.
  expect_error(fit_small(small_arm_df(2L)),
               class = "cateomics_group_too_small_for_folds")
  expect_error(fit_small(small_arm_df(2L), n_folds = 20L),
               class = "cateomics_group_too_small_for_folds")
  # 3 controls at 5 folds: one held out, two to train on. Checked with a single
  # learner, because the Super Learner has its own requirement on the number of
  # rows its weights are solved on (below).
  expect_no_error(fit_small(small_arm_df(3L), eval_nuisance = "rf"))
  # The same 3 at 2 folds holds out two and leaves one, so fewer folds is
  # stricter here, not looser.
  expect_error(fit_small(small_arm_df(3L), n_folds = 2L, eval_nuisance = "rf"),
               class = "cateomics_group_too_small_for_folds")
  expect_no_error(fit_small(small_arm_df(4L), n_folds = 2L, eval_nuisance = "rf"))
})

test_that("the Super Learner needs more rows than members to weight them", {
  # The fold structure is satisfied here and no row is predicted by a model
  # that saw it, so only the identification condition catches this.
  expect_error(fit_small(small_arm_df(3L)),
               class = "cateomics_sl_library_too_large")
  msg <- tryCatch(fit_small(small_arm_df(3L)), error = conditionMessage)
  expect_match(msg, "at most 2 member", fixed = TRUE)
  # A library the group can support goes through.
  expect_no_error(fit_small(small_arm_df(3L),
                            eval_library = c("mean", "ridge")))
})

test_that("the message names a reachable remedy", {
  msg2 <- tryCatch(fit_small(small_arm_df(2L)), error = conditionMessage)
  expect_match(msg2, "No number of folds can satisfy this", fixed = TRUE)
  msg3 <- tryCatch(fit_small(small_arm_df(3L), n_folds = 2L),
                   error = conditionMessage)
  expect_match(msg3, "3 or more would do", fixed = TRUE)
})

test_that("without folds a group still needs two rows to fit at all", {
  # partition = "none" fits on all of the group's rows, so two is enough for a
  # single learner and the Super Learner still cannot weight on held-out rows.
  expect_no_error(fit_small(small_arm_df(2L), partition = "none",
                            eval_nuisance = "rf"))
  expect_error(fit_small(small_arm_df(2L), partition = "none"),
               class = "cateomics_group_too_small_for_folds")
  expect_error(fit_small(small_arm_df(1L), partition = "none",
                         eval_nuisance = "rf"),
               class = "cateomics_group_too_small")
})

test_that("a glmnet penalty that cannot be cross-validated says so", {
  set.seed(77); n <- 40L; p <- 4L
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  Tv <- c(0L, 0L, rep(1L, n - 2L))
  # The shared rule now reaches this first: two rows of a class cannot fill
  # the folds, so the penalty cannot be selected.
  expect_error(
    .fit_base_learner(X, Tv, type = "classification", learner = "lasso",
                      random_state = 0L),
    class = "cateomics_cv_split_too_small")
  # The wrapper around glmnet's own refusal is a backstop for failures this
  # package does not model, and it stays reachable: a degenerate input is
  # refused by glmnet for a reason the row and class counts say nothing about.
  # Measured: raising `nfolds` or thinning the minority class no longer gets
  # here, because the shared rule reaches those first.
  expect_error(
    suppressWarnings(
      .fit_base_learner(X, rep(1, n), type = "regression", learner = "lasso",
                        random_state = 0L)),
    class = "cateomics_glmnet_cv_failed")
  Xz <- X; Xz[, ] <- 0
  expect_error(
    suppressWarnings(
      .fit_base_learner(Xz, stats::rnorm(n), type = "regression",
                        learner = "lasso", random_state = 0L)),
    class = "cateomics_glmnet_cv_failed")
  # the same data through a learner that does not split internally
  expect_no_error(suppressWarnings(
    .fit_base_learner(X, Tv, type = "classification", learner = "logistic",
                      random_state = 0L)))
})

# --- splits CATEomics parameterises on the Python side (2026-10-09) ----------
# Two different things, kept apart: the internal cross-fit estimates the
# nuisances out of fold so the final stage is orthogonal, and a penalised
# learner's own cross-validation selects a penalty. Only the DR learner has
# the former (verified: `cv` appears in no other fitter's formals, and
# econml's metalearners take no `cv`).

py_arm_df <- function(m_control, n = 60L, p = 4L, seed = 3L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  data.frame(y = 2 + X[, 1] + stats::rnorm(n),
             t = c(rep(0L, m_control), rep(1L, n - m_control)), X)
}
fit_py <- function(df, method, ...) {
  suppressMessages(suppressWarnings(
    cf_fit(data = df, Y = "y", T = "t", X_baseline = paste0("g", 1:4),
           method = method, preprocess = FALSE, random_state = 0L,
           eval_nuisance = "rf", ...)))
}

test_that("the internal cross-fit refuses a group it cannot split", {
  skip_if_no_cf_python()
  # m = 3 at 5 outer folds leaves 2, and cv = 5 inside leaves 1 of those.
  expect_error(fit_py(py_arm_df(3L), "drlearner"),
               class = "cateomics_group_too_small_for_cv")
  # Lowering cv does not help: a smaller inner split holds out a bigger share.
  expect_error(fit_py(py_arm_df(3L), "drlearner",
                      params = list(drlearner = list(cv = 2L))),
               class = "cateomics_group_too_small_for_cv")
  # m = 4 leaves 3 after the outer split and 2 after the inner one.
  expect_no_error(fit_py(py_arm_df(4L), "drlearner", model_outcome = "rf",
                         model_final = "rf", model_propensity = "logistic"))
  # The metalearners do not cross-fit internally, so the same data goes
  # through: `cv` appears in no fitter's formals but the DR learner's.
  # `model_outcome = "rf"` isolates the cross-fit: the default outcome learner
  # is a penalty, which the shared rule refuses at this arm size on its own.
  expect_no_error(fit_py(py_arm_df(3L), "tlearner", model_outcome = "rf"))
  expect_null(formals(CATEomics:::.cf_method_fitter("tlearner"))$cv)
  expect_false(is.null(formals(CATEomics:::.cf_method_fitter("drlearner"))$cv))
})

test_that("one fold count and one rule govern every penalty selection", {
  # Every split that selects a penalty uses the same count, so one pair of
  # thresholds describes all of them: below it the split cannot be formed at
  # all, below twice it the folds hold fewer than two rows each.
  cv <- CATEomics:::.CF_CV_FOLDS
  expect_identical(cv, 5L)
  expect_identical(CATEomics:::.CF_CV_TOKENS, c("lasso", "ridge", "elastic_net"))

  # The decision, stated directly on the count reaching the learner.
  expect_error(CATEomics:::.cf_check_cv_split(cv - 1L, "lasso"),
               class = "cateomics_cv_split_too_small")
  expect_warning(CATEomics:::.cf_check_cv_split(cv, "lasso"),
                 class = "cateomics_cv_split_thin")
  expect_warning(CATEomics:::.cf_check_cv_split(2L * cv - 1L, "lasso"),
                 class = "cateomics_cv_split_thin")
  expect_silent(CATEomics:::.cf_check_cv_split(2L * cv, "lasso"))
  # Learners that do not cross-validate are not subject to it at any count.
  for (tok in c("linear", "logistic", "rf", "xgboost")) {
    expect_silent(CATEomics:::.cf_check_cv_split(1L, tok))
  }
})

test_that("the count a learner is checked against has the splits above it removed", {
  # Folds nest, so the rows a penalised nuisance actually receives are what is
  # left after the outer fold and, where the method has one, after the
  # internal cross-fit.
  expect_identical(CATEomics:::.cf_cv_rows_left(20L), 20L)
  expect_identical(CATEomics:::.cf_cv_rows_left(20L, n_folds = 5L), 16L)
  expect_identical(CATEomics:::.cf_cv_rows_left(20L, n_folds = 5L, cv = 5L), 12L)
  # A count below the floor at one level is not rescued by skipping the next.
  expect_identical(CATEomics:::.cf_cv_rows_left(8L, n_folds = 5L, cv = 5L), 4L)
  # An absent or degenerate split count subtracts nothing.
  expect_identical(CATEomics:::.cf_cv_rows_left(9L, n_folds = NULL, cv = 1L), 9L)
})

test_that("a penalised token reaching the Python factory is checked, in either role", {
  # The propensity's target is always a class label, so the split is
  # stratified; a continuous outcome within an arm is split on rows. Both
  # thresholds are the same number now, and ridge is no longer exempt because
  # the factory gives RidgeCV the same fold count as the others.
  expect_error(fit_py(py_arm_df(4L), "tlearner", partition = "none",
                      model_propensity = "ridge"),
               class = "cateomics_cv_split_too_small")
  expect_no_error(fit_py(py_arm_df(4L), "tlearner", partition = "none",
                         model_propensity = "logistic", model_outcome = "rf"))
  expect_no_error(fit_py(py_arm_df(20L), "tlearner", partition = "none",
                         model_propensity = "ridge"))

  fit_c <- function(df, ...) suppressMessages(suppressWarnings(
    cf_fit(data = df, Y = "y", T = "t", X_baseline = paste0("g", 1:4),
           method = "tlearner", preprocess = FALSE, random_state = 0L,
           eval_nuisance = "rf", partition = "none", ...)))
  for (tok in c("lasso", "elastic_net")) {
    expect_error(fit_c(py_arm_df(4L), model_outcome = tok),
                 class = "cateomics_cv_split_too_small")
  }
  # The one exemption, and the reason it is one: on a continuous target in
  # Python, ridge is RidgeCV at scikit-learn's default, the closed-form
  # leave-one-out, which forms no partition at all.
  expect_no_error(fit_c(py_arm_df(4L), model_outcome = "ridge"))
  expect_false("ridge" %in% setdiff(CATEomics:::.CF_CV_TOKENS,
                                    CATEomics:::.CF_CV_NO_SPLIT_PY_ROWS))
  expect_true("ridge" %in% CATEomics:::.CF_CV_TOKENS)
  expect_silent(CATEomics:::.cf_check_cv_rows(1L, NULL, NULL, "ridge", "model_outcome"))
  # but it is not exempt where it does split: stratified in Python, and either
  # role in R.
  expect_error(CATEomics:::.cf_check_cv_classes(c(rep(0L, 40L), rep(1L, 3L)),
                                              NULL, NULL, "ridge", "model_propensity"),
               class = "cateomics_cv_split_too_small")
  set.seed(8)
  Xr <- matrix(stats::rnorm(16L), 4L, 4L)
  expect_error(.fit_base_learner(Xr, stats::rnorm(4L), type = "regression",
                                 learner = "ridge", random_state = 0L),
               class = "cateomics_cv_split_too_small")
  expect_no_error(fit_c(py_arm_df(4L), model_outcome = "rf"))
  expect_no_error(fit_c(py_arm_df(10L), model_outcome = "lasso"))
})

test_that("the internal cross-fit nests inside the penalty check", {
  # The DR learner splits again before its nuisance learner is fitted, so the
  # same arm that passes for a metalearner can fail for it. m = 10 leaves 8
  # after the outer split, which the metalearners hand straight to the
  # learner, and 6 after the DR learner's own cv.
  expect_identical(CATEomics:::.cf_cv_rows_left(10L, 5L, NULL), 8L)
  expect_identical(CATEomics:::.cf_cv_rows_left(10L, 5L, 5L), 6L)
  # m = 8 is the boundary: 6 for a metalearner, 4 for the DR learner.
  expect_identical(CATEomics:::.cf_cv_rows_left(8L, 5L, NULL), 6L)
  expect_identical(CATEomics:::.cf_cv_rows_left(8L, 5L, 5L), 4L)
  expect_warning(CATEomics:::.cf_check_cv_rows(8L, 5L, NULL, "lasso", "model_outcome"),
                 class = "cateomics_cv_split_thin")
  expect_error(CATEomics:::.cf_check_cv_rows(8L, 5L, 5L, "lasso", "model_outcome"),
               class = "cateomics_cv_split_too_small")
  # The classification form reads the smallest class, not the total.
  strata <- c(rep(0L, 40L), rep(1L, 8L))
  expect_error(CATEomics:::.cf_check_cv_classes(strata, 5L, 5L, "ridge", "model_propensity"),
               class = "cateomics_cv_split_too_small")
  expect_silent(CATEomics:::.cf_check_cv_classes(c(rep(0L, 40L), rep(1L, 40L)),
                                               5L, 5L, "ridge", "model_propensity"))
})

test_that("the fold count in the Python factory has not drifted from the R copy", {
  # The factory writes the count in Python and R cannot read it across the
  # process boundary, so it is restated. This reads the file rather than
  # trusting the restatement, and checks that every penalised branch uses the
  # one constant rather than a literal of its own.
  f <- system.file("python", "learner_factory.py", package = "CATEomics")
  if (!nzchar(f)) f <- "../../inst/python/learner_factory.py"
  skip_if_not(file.exists(f))
  src <- readLines(f, warn = FALSE)
  decl <- grep("^CV_FOLDS\\s*=\\s*[0-9]+", src, value = TRUE)
  expect_length(decl, 1L)
  expect_identical(as.integer(sub("[^0-9]*([0-9]+).*", "\\1", decl)),
                   as.integer(CATEomics:::.CF_CV_FOLDS))
  # No penalised constructor may carry its own numeric fold count.
  ctor <- grep("LassoCV\\(|RidgeCV\\(|ElasticNetCV\\(|LogisticRegressionCV\\(", src)
  block <- unlist(lapply(ctor, function(i) src[i:min(i + 4L, length(src))]))
  expect_false(any(grepl("cv\\s*=\\s*[0-9]", block)))
  # Five of the six reach the constant. The sixth is RidgeCV on a continuous
  # target, left at scikit-learn's closed-form leave-one-out.
  expect_length(grep("cv=CV_FOLDS", src), 5L)
  expect_length(grep("return RidgeCV\\(\\)", src), 1L)
})

test_that("the R penalty selection takes the same count and the same rule", {
  # cv.glmnet defaults to 10 folds of its own and, on a continuous target,
  # returns a fit without stopping or warning when the rows cannot fill them
  # (measured at 3, 5 and 8 rows). The package passes its own count and
  # applies the shared rule, on the rows actually in hand.
  cv <- CATEomics:::.CF_CV_FOLDS
  set.seed(1)
  thin <- function(n) {
    X <- matrix(stats::rnorm(n * 4L), n, 4L)
    list(X = X, y = stats::rnorm(n))
  }
  d <- thin(cv - 1L)
  expect_error(.fit_base_learner(d$X, d$y, type = "regression",
                                 learner = "lasso", random_state = 0L),
               class = "cateomics_cv_split_too_small")
  d <- thin(cv)
  expect_warning(.fit_base_learner(d$X, d$y, type = "regression",
                                   learner = "lasso", random_state = 0L),
                 class = "cateomics_cv_split_thin")
  # still a fit in the thin band
  f <- suppressWarnings(.fit_base_learner(d$X, d$y, type = "regression",
                                          learner = "lasso", random_state = 0L))
  expect_s3_class(f, "cf_base_learner")
  # above twice the count there is nothing to report
  d <- thin(2L * cv)
  expect_no_condition(
    suppressWarnings(.fit_base_learner(d$X, d$y, type = "regression",
                                       learner = "lasso", random_state = 0L)),
    class = "cateomics_cv_split_thin")
  # The fold count actually handed to glmnet is the package's, not glmnet's
  # default of 10, and a caller's own `nfolds` still wins. Read off the
  # recorded call rather than inferred from the selected lambda, which can
  # coincide across fold counts.
  d <- thin(40L)
  f <- .fit_base_learner(d$X, d$y, type = "regression", learner = "lasso",
                         random_state = 0L)
  expect_identical(as.integer(f$model$call$nfolds), cv)
  f3 <- .fit_base_learner(d$X, d$y, type = "regression", learner = "lasso",
                          random_state = 0L, .params = list(nfolds = 3L))
  expect_identical(as.integer(f3$model$call$nfolds), 3L)
})

test_that("the binary R penalty selection reads the smallest class", {
  # The same rule, on the class that has to fill the folds rather than on the
  # total rows.
  set.seed(2)
  X <- matrix(stats::rnorm(40L * 4L), 40L, 4L)
  y <- c(rep(0L, 36L), rep(1L, 4L))
  expect_error(.fit_base_learner(X, y, type = "classification",
                                 learner = "ridge", random_state = 0L),
               class = "cateomics_cv_split_too_small")
  y <- c(rep(0L, 33L), rep(1L, 7L))
  expect_warning(.fit_base_learner(X, y, type = "classification",
                                   learner = "ridge", random_state = 0L),
                 class = "cateomics_cv_split_thin")
  y <- c(rep(0L, 20L), rep(1L, 20L))
  expect_no_condition(
    suppressWarnings(.fit_base_learner(X, y, type = "classification",
                                       learner = "ridge", random_state = 0L)),
    class = "cateomics_cv_split_thin")
})

test_that("the rule follows the fold count actually in force, not the default", {
  # The count is this package's default and a caller can lower it, so a check
  # pinned to the constant would refuse a fit the learner can perform. Measured:
  # cv.glmnet returns a fit for four rows over three folds.
  cv <- CATEomics:::.CF_CV_FOLDS
  expect_error(CATEomics:::.cf_check_cv_split(3L, "lasso", cv_folds = cv),
               class = "cateomics_cv_split_too_small")
  expect_warning(CATEomics:::.cf_check_cv_split(3L, "lasso", cv_folds = 3L),
                 class = "cateomics_cv_split_thin")
  expect_silent(CATEomics:::.cf_check_cv_split(6L, "lasso", cv_folds = 3L))
  # A degenerate count decides nothing rather than stopping on a comparison.
  expect_silent(CATEomics:::.cf_check_cv_split(1L, "lasso", cv_folds = 1L))

  set.seed(4)
  X <- matrix(stats::rnorm(16L), 4L, 4L); y <- stats::rnorm(4L)
  expect_error(.fit_base_learner(X, y, type = "regression", learner = "lasso",
                                 random_state = 0L),
               class = "cateomics_cv_split_too_small")
  expect_warning(.fit_base_learner(X, y, type = "regression", learner = "lasso",
                                   random_state = 0L, .params = list(nfolds = 3L)),
                 class = "cateomics_cv_split_thin")
})

test_that("a Python slot's own fold count is read, not the internal cross-fit's", {
  # `.cf_apply_params()` splits the two by the fitter's formals: a flat `cv` is
  # the DR learner's internal cross-fit, a slot's `cv` is that learner's own.
  skip_if_no_cf_python()
  d <- py_arm_df(6L)
  fit <- function(...) suppressMessages(suppressWarnings(
    cf_fit(data = d, Y = "y", T = "t", X_baseline = paste0("g", 1:4),
           method = "tlearner", preprocess = FALSE, random_state = 0L,
           eval_nuisance = "rf", model_outcome = "lasso", ...)))
  expect_error(fit(), class = "cateomics_cv_split_too_small")
  expect_no_error(fit(params = list(outcome = list(cv = 3L))))
})


# The treatment column's type ------------------------------------------------
#
# `as.integer()` on a factor returns level codes, which start at 1, so a
# two-level factor used to pass the 0/1 check by character coercion and then
# become 1/2. Every arm-wise fit then looked for a control group that was not
# there, and the caller saw a Super Learner failure rather than a word about
# the treatment. The types whose arms are unambiguous are converted; the ones
# that would need a guess are refused.

test_that("unambiguous treatment types give the same fit as numeric 0/1", {
  set.seed(7); n <- 80L; p <- 5L
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("g", seq_len(p))))
  Tn <- stats::rbinom(n, 1L, 0.5)
  Y <- X[, 2] + Tn * (1 + X[, 1]) + stats::rnorm(n)
  fit_with <- function(tv) cf_fit(data = data.frame(Y = Y, T = tv, X),
                                  Y = "Y", T = "T", X_baseline = colnames(X),
                                  method = "grf", preprocess = FALSE)$res$ate
  base <- fit_with(Tn)
  expect_equal(fit_with(Tn == 1L), base)                      # logical
  expect_equal(fit_with(factor(Tn, levels = c(0, 1))), base)  # factor "0"/"1"
  expect_equal(fit_with(as.character(Tn)), base)              # character "0"/"1"
})

test_that("a treatment whose arms are not named 0/1 is refused, not guessed", {
  set.seed(7); n <- 60L
  X <- matrix(stats::rnorm(n * 4L), n, 4L,
              dimnames = list(NULL, paste0("g", seq_len(4L))))
  Tn <- stats::rbinom(n, 1L, 0.5)
  Y <- X[, 2] + Tn + stats::rnorm(n)
  bad <- function(tv) cf_fit(data = data.frame(Y = Y, T = tv, X), Y = "Y",
                             T = "T", X_baseline = colnames(X),
                             method = "grf", preprocess = FALSE)
  # Which level is treated is the caller's to state: level order is theirs to
  # set, so taking the second one would silently invert some callers' arms.
  expect_error(bad(factor(ifelse(Tn == 1L, "treat", "ctrl"))),
               class = "cateomics_treatment_not_binary")
  expect_error(bad(ifelse(Tn == 1L, "treat", "ctrl")),
               class = "cateomics_treatment_not_binary")
  expect_error(bad(Tn + 1L), class = "cateomics_treatment_not_binary")
  # The message has to name the treatment rather than whatever fails next.
  expect_error(bad(factor(ifelse(Tn == 1L, "treat", "ctrl"))), "numeric 0/1")
})

test_that("cf_evaluate refuses a treated fraction outside (0, 1]", {
  skip_if_not_installed("grf")
  se  <- cf_example_data(n = 40, p = 40, n_prog = 6, n_pred = 3, seed = 2)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 30, partition = "oof",
                params = list(num.trees = 100L, tune.parameters = "none"))
  # 0 divides by zero in the 1/q factor and used to return NaN silently;
  # anything above 1 asks for more units than there are and used to surface a
  # raw "missing value where TRUE/FALSE needed" from the weighting helper.
  expect_error(cf_evaluate(fit, metrics = "rate", target_quantile = 0),
               class = "cateomics_invalid_quantile")
  expect_error(cf_evaluate(fit, metrics = "rate", target_quantile = 1.5),
               class = "cateomics_invalid_quantile")
  expect_error(cf_evaluate(fit, metrics = "rate", target_quantile = NA_real_),
               class = "cateomics_invalid_quantile")
  expect_error(cf_evaluate(fit, metrics = "rate",
                           rate_quantiles = c(0, 0.5, 1)),
               class = "cateomics_invalid_quantile")
  # q = 1 is the whole sample, which is a defined (and zero) TOC.
  ev <- cf_evaluate(fit, metrics = "rate", target_quantile = 1,
                    n_bootstrap = 20L)
  expect_equal(ev$rate$toc, 0, tolerance = 1e-10)
})

test_that("cf_compare reports a control argument no compared method takes", {
  skip_if_not_installed("grf")
  d <- data.frame(Y = rnorm(60), T = rbinom(60, 1, 0.5),
                  g1 = rnorm(60), g2 = rnorm(60), g3 = rnorm(60))
  nm <- c("g1", "g2", "g3")
  # The dispatcher's per-fit report is swallowed by cf_compare's fold helper,
  # so the set-level report has to be raised at the entry.
  expect_warning(
    cf_compare(data = d, Y = "Y", T = "T", X_baseline = nm,
               methods = c("slearner", "tlearner"), n_folds = 2L,
               clip = c(0.1, 0.9), verbose = FALSE),
    class = "cateomics_control_arg_unreachable")
  # A value that reaches at least one compared method was used, so nothing is
  # said: a comparison legitimately mixes methods that differ in what they take.
  # Scoped to the class under test: a small split also draws an unrelated
  # glmnet warning from the Super Learner members.
  expect_no_warning(
    cf_compare(data = d, Y = "Y", T = "T", X_baseline = nm,
               methods = "grf", n_folds = 2L, clip = c(0.1, 0.9),
               verbose = FALSE,
               params = list(num.trees = 100L, tune.parameters = "none")),
    class = "cateomics_control_arg_unreachable")
  # bcf takes a known propensity through its native `pihat`, so it is not
  # reported even though it has no `propensity_score` formal.
  expect_no_warning(.cf_check_unreachable_control("bcf", "propensity_score"),
                    class = "cateomics_control_arg_unreachable")
  expect_warning(
    .cf_check_unreachable_control(c("slearner", "tlearner"), "model_final"),
    class = "cateomics_control_arg_unreachable")
  expect_no_warning(
    .cf_check_unreachable_control(c("slearner", "xlearner"), "model_final"),
    class = "cateomics_control_arg_unreachable")
})

test_that("topk_frac is checked as a quantile, in both directions", {
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 15, random_state = 0L,
                propensity_score = 0.5, eval_propensity_score = 0.5)
  # Unchecked, 0 treated one unit and 1.5 made every decision NA.
  expect_error(cf_rule(fit, type = "topk", topk_frac = 0),
               class = "cateomics_invalid_quantile")
  expect_error(cf_rule(fit, type = "topk", topk_frac = 1.5),
               class = "cateomics_invalid_quantile")
  expect_error(cf_rule(fit, type = "topk", grid = c(0.2, 1.5)),
               class = "cateomics_invalid_quantile")
  # A valid fraction still goes through, and so does the boundary.
  expect_no_error(suppressMessages(suppressWarnings(
    cf_rule(fit, type = "topk", topk_frac = 1, compute_policy_value = FALSE))))
})

test_that("cf_plot_policy_curve refuses an unknown argument, as documented", {
  skip_if_not_installed("ggplot2")
  se  <- cf_example_data(n = 40, p = 30, n_prog = 3, n_pred = 3, seed = 1)
  fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
                X_baseline = c("age", "sex"), method = "grf",
                n_features = 15, random_state = 0L,
                propensity_score = 0.5, eval_propensity_score = 0.5)
  rg <- suppressMessages(suppressWarnings(
    cf_rule(fit, type = "topk", grid = c(0.2, 0.4))))
  expect_error(cf_plot_policy_curve(rg, nonsense_argument = 1),
               class = "cateomics_unknown_argument")
  expect_no_error(cf_plot_policy_curve(rg))
})

test_that("a mixed bootstrap request is served where it can be", {
  chk <- .cf_check_ci_policy
  # One method, and it cannot resample: refused, because a single fit has
  # nothing to fall back on.
  expect_error(chk("grf", TRUE, NULL, allow_partial = FALSE),
               class = "cateomics_bootstrap_unsupported")
  # A set in which some can: served where it can be, said once.
  expect_warning(chk(c("grf", "drlearner"), TRUE, list(NULL, NULL),
                     allow_partial = TRUE),
                 class = "cateomics_bootstrap_partial")
  # A set in which none can: still refused.
  expect_error(chk(c("grf", "bcf"), TRUE, list(NULL, NULL),
                   allow_partial = TRUE),
               class = "cateomics_bootstrap_unsupported")
  # One warning for the set, not one per candidate.
  n <- 0L
  withCallingHandlers(
    chk(c("grf", "bcf", "drlearner"), TRUE, list(NULL, NULL, NULL),
        allow_partial = TRUE),
    warning = function(w) {
      if (inherits(w, "cateomics_bootstrap_partial")) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
  expect_equal(n, 1L)
})
