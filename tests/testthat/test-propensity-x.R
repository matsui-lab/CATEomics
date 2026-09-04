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
  # `tune.parameters = "none"` is grf-specific; passing it to other backends
  # via `...` triggers an "unused argument" error (e.g. bcf). Apply it
  # only when method = "grf".
  args <- list(
    Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
    propensity_X = pX,
    data = df, method = method,
    random_state = 42L, partition = "none"
  )
  if (identical(method, "grf")) args$tune.parameters <- "none"
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
        cfomics_propensity_X_ignored = function(w) {
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
           propensity_X = c("X1", "X_does_not_exist"),
           data = df, method = "grf",
           random_state = 42L, partition = "none",
           tune.parameters = "none"),
    class = "cfomics_propensity_X_missing"
  )
})

test_that("propensity_X invalid type is rejected", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = paste0("X", 1:6),
           propensity_X = 1L,  # numeric, not character
           data = df, method = "grf",
           random_state = 42L, partition = "none",
           tune.parameters = "none"),
    class = "cfomics_invalid_propensity_X"
  )
})
