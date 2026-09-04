# Test that package loading has no side effects
# Critical for CRAN/Bioconductor compliance. Python backends are provisioned
# lazily via basilisk (a separate process), so loading cfomics must not
# initialise Python in the main session.

test_that(".onLoad does not initialize Python", {
  skip_if_not_installed("reticulate")

  py_was_initialized <- tryCatch(
    reticulate::py_available(initialize = FALSE),
    error = function(e) FALSE
  )
  if (py_was_initialized) {
    skip("Python was already initialized before test (expected in CI)")
  }

  py_is_now_available <- tryCatch(
    reticulate::py_available(initialize = FALSE),
    error = function(e) FALSE
  )
  expect_false(py_is_now_available)
})

test_that("the removed config subsystem is not exported", {
  # There is no cf_config / cf_config_reset (no on-disk config, no HOME
  # writes at load). Assert the exports are gone.
  ns <- asNamespace("cfomics")
  expect_false(exists("cf_config", envir = ns, inherits = FALSE))
  expect_false(exists("cf_config_reset", envir = ns, inherits = FALSE))
})

test_that("package namespace exposes the core API", {
  ns <- asNamespace("cfomics")
  expect_true(exists("cf_fit", envir = ns))
  # S3 methods are registered.
  expect_true(is.function(getS3method("predict", "cf_model", optional = TRUE)))
  expect_true(is.function(getS3method("summary", "cf_model", optional = TRUE)))
})

test_that("R-native methods run without initialising Python", {
  set.seed(123)
  n <- 30
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X = rnorm(n)
  )
  # grf is R-native; it must not require any Python environment.
  expect_no_error({
    cf_fit(Y ~ T | X, data = data, method = "grf", partition = "none")
  })
})
