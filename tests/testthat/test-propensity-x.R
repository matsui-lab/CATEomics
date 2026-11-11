## Tests for propensity_X across the PS-using CATE methods and the ignored
## warning for PS-not-used / PS-architecturally-fixed methods.

library(testthat)

.make_df <- function(n = 120L, p = 6L, seed = 42L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  T <- stats::rbinom(n, 1, 0.5)
  Y <- as.numeric(T) * 1.5 + rowSums(X) + 0.1 * stats::rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  colnames(df) <- c("Y", "T", paste0("X", seq_len(p)))
  df
}

.fit_pX <- function(method, pX = c("X1", "X3")) {
  df <- .make_df()
  # `tune.parameters` is grf-specific, so it only goes into `params` when the
  # method is grf (another method would reject the unknown name).
  args <- list(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
    propensity_X = pX,
    data = df, method = method,
    random_state = 42L, partition = "none"
  )
  if (identical(method, "grf")) args$params <- list(tune.parameters = "none")
  suppressWarnings(suppressMessages(do.call(cf_fit, args)))
}

test_that("grf accepts propensity_X (forest.W on subset)", {
  fit <- .fit_pX("grf")
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
  expect_equal(fit$meta$propensity_X, c("X1", "X3"))
})

test_that("bcf accepts propensity_X (ps_model on subset)", {
  skip_if_not_installed("bcf")
  fit <- .fit_pX("bcf")
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
})

test_that("drlearner accepts propensity_X (ColumnSubsetPropensityModel)", {
  skip_if_no_cf_python()
  fit <- .fit_pX("drlearner")
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
})

test_that("xlearner accepts propensity_X (ColumnSubsetPropensityModel)", {
  skip_if_no_cf_python()
  fit <- .fit_pX("xlearner")
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
})

test_that("7-method warning fires for PS-ignored methods", {
  df <- .make_df()
  ignored_methods <- c("slearner", "tlearner",
                       "dragonnet", "cevae", "ganite", "tarnet")
  for (m in ignored_methods) {
    # Reset warning verbosity so each method's once-per-session warning fires.
    tryCatch(rlang::reset_warning_verbosity(paste0("propensity_X_", m)),
             error = function(e) NULL)
    warned <- FALSE
    .res <- tryCatch(
      withCallingHandlers(
        suppressMessages(cf_fit(
          Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
          propensity_X = c("X1", "X3"),
          data = df, method = m,
          random_state = 42L, partition = "none"
        )),
        cateomics_propensity_X_ignored = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    expect_true(warned, info = sprintf("method='%s'", m))
  }
})

test_that("propensity_X missing columns are rejected", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
            propensity_X = c("X1", "X_does_not_exist"), data = df,
            method = "grf", random_state = 42L, partition = "none",
            params = list(tune.parameters = "none")),
    class = "cateomics_propensity_X_missing"
  )
})

test_that("propensity_X invalid type is rejected", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
            propensity_X = 1L, data = df, method = "grf", random_state = 42L,
            partition = "none", params = list(tune.parameters = "none")),
    class = "cateomics_invalid_propensity_X"
  )
})


# propensity_X reaches the propensity model ------------------------------
#
# The tests above end at expect_s3_class(fit, "cf_model") and a finite ATE,
# which a method that discarded the subset would also satisfy. These fit the
# same data twice with different subsets and require the stored propensity
# to differ.

.px_confounded <- function(n = 120L, seed = 21L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  # x3 drives assignment; x1 and x2 do not.
  T_vec <- rbinom(n, 1L, plogis(1.8 * x3))
  list(df = data.frame(Y = x2 + T_vec * (1 + x1) + rnorm(n), T = T_vec,
                       x1 = x1, x2 = x2, x3 = x3),
       true_ps = plogis(1.8 * x3))
}

test_that("grf's propensity changes with the propensity_X subset", {
  skip_if_not_installed("grf")
  d <- .px_confounded()
  without <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3"),
    method = "grf", propensity_X = "x1", partition = "none",
    random_state = 1L))
  with_it <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3"),
    method = "grf", propensity_X = c("x1", "x3"), partition = "none",
    random_state = 1L))

  expect_gt(max(abs(without$res$ps - with_it$res$ps)), 0.1)
  # Adding the variable that actually drives assignment has to help.
  expect_gt(stats::cor(with_it$res$ps, d$true_ps), 0.8)
  expect_gt(stats::cor(with_it$res$ps, d$true_ps),
            stats::cor(without$res$ps, d$true_ps))
})

test_that("bcf's propensity changes with the propensity_X subset", {
  skip_if_not_installed("bcf")
  d <- .px_confounded()
  short <- list(bcf = list(nburn = 10L, nsim = 10L, n_chains = 1L))
  # bcf's propensity model needs at least two columns, so both subsets have
  # two or more and they differ by the one that matters.
  without <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3"),
    method = "bcf", params = short, propensity_X = c("x1", "x2"),
    partition = "none", random_state = 1L))
  with_it <- suppressMessages(cf_fit(
    data = d$df, Y = "Y", T = "T", X_baseline = c("x1", "x2", "x3"),
    method = "bcf", params = short, propensity_X = c("x1", "x2", "x3"),
    partition = "none", random_state = 1L))

  expect_gt(max(abs(without$res$ps - with_it$res$ps)), 0.1)
  expect_gt(stats::cor(with_it$res$ps, d$true_ps), 0.8)
})

# --------------------------------------------------------------------------
# The default: with `propensity_X` unset the propensity goes on the clinical
# covariates alone, and the omics block is left to the outcome models.
# --------------------------------------------------------------------------

.se_small <- function() {
  cf_example_data(n = 80, p = 100, n_prog = 4, n_pred = 4, seed = 1)
}
.fit_se <- function(...) {
  suppressWarnings(suppressMessages(
    cf_fit(data = .se_small(), Y = "outcome", T = "treatment",
           method = "grf", n_features = 50, ...)))
}

test_that("unset, the propensity is fitted on the clinical covariates", {
  fit <- .fit_se(X_baseline = c("age", "sex"))
  # `sex` is a factor, so the expansion names are what reaches the model.
  expect_setequal(fit$meta$propensity_X, c("age", "sexF", "sexM"))
  expect_gt(ncol(fit$data$X_preprocessed), length(fit$meta$propensity_X))
  # Naming every column reproduces what earlier versions did by default, and
  # it is a different propensity, so the default is doing something.
  all_cols <- .fit_se(X_baseline = c("age", "sex"),
                      propensity_X = colnames(fit$data$X_preprocessed))
  expect_gt(max(abs(fit$res$ps - all_cols$res$ps)), 0.01)
})

test_that("cf_fit and cf_compare resolve the same columns", {
  fit <- .fit_se(X_baseline = c("age", "sex"))
  cmp <- suppressWarnings(suppressMessages(
    cf_compare(data = .se_small(), Y = "outcome", T = "treatment",
               X_baseline = c("age", "sex"), methods = "grf",
               n_features = 50)))
  expect_identical(fit$meta$propensity_X, cmp$models$grf$meta$propensity_X)
})

test_that("the evaluation nuisance propensity follows the same default", {
  se <- .se_small()
  common <- list(data = se, Y = "outcome", T = "treatment",
                 X_baseline = c("age", "sex"), methods = "grf",
                 n_features = 50)
  by_default <- suppressWarnings(suppressMessages(do.call(cf_compare, common)))
  all_cols <- suppressWarnings(suppressMessages(do.call(cf_compare,
    c(common, list(propensity_X = colnames(
      by_default$models$grf$data$X_preprocessed))))))
  expect_gt(max(abs(by_default$metadata$nuisance_oof$pi_oof -
                      all_cols$metadata$nuisance_oof$pi_oof)), 0.01)
})

test_that("with no clinical covariate the fit says so and uses every column", {
  w <- NULL
  fit <- withCallingHandlers(
    suppressMessages(cf_fit(data = .se_small(), Y = "outcome", T = "treatment",
                            omics = "all", method = "grf", n_features = 50)),
    cateomics_propensity_X_fallback = function(c) {
      w <<- conditionMessage(c); invokeRestart("muffleWarning")
    })
  expect_true(!is.null(w))
  expect_match(w, "propensity_X", fixed = TRUE)
  expect_match(w, "propensity_score", fixed = TRUE)
  # NULL is the record that means every column, which is what it fell back to.
  expect_null(fit$meta$propensity_X)
})

test_that("a known score on both sides restricts nothing and says nothing", {
  w <- FALSE
  fit <- withCallingHandlers(
    suppressMessages(cf_fit(data = .se_small(), Y = "outcome", T = "treatment",
                            omics = "all", method = "grf", n_features = 50,
                            propensity_score = 0.5,
                            eval_propensity_score = 0.5)),
    cateomics_propensity_X_fallback = function(c) {
      w <<- TRUE; invokeRestart("muffleWarning")
    })
  expect_false(w)
  expect_null(fit$meta$propensity_X)
})

test_that("the ignored-argument warning is about the argument, not the default", {
  # The warning says `propensity_X` did not reach the estimator. Firing it on
  # a call that never named it would make it meaningless.
  fires <- function(...) {
    hit <- FALSE
    # The warning is throttled per method, and an earlier block in this file
    # has already spent slearner's turn, so ask rlang for every occurrence.
    op <- options(rlib_warning_verbosity = "verbose")
    on.exit(options(op), add = TRUE)
    withCallingHandlers(
      suppressMessages(cf_fit(data = .se_small(), Y = "outcome", T = "treatment",
                              X_baseline = c("age", "sex"), method = "slearner",
                              n_features = 50, ...)),
      cateomics_propensity_X_ignored = function(c) {
        hit <<- TRUE; invokeRestart("muffleWarning")
      })
    hit
  }
  skip_if_no_cf_python()
  expect_false(fires())
  expect_true(fires(propensity_X = "age"))
})
