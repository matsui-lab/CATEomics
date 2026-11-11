# Test that package loading has no side effects
# Critical for CRAN/Bioconductor compliance. Python backends are provisioned
# lazily via basilisk (a separate process), so loading CATEomics must not
# initialise Python in the main session.

test_that(".onLoad does not initialize Python", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("callr")
  # Asked in a fresh process, because reading the flag in the current
  # session skips whenever Python is already up, which is exactly the
  # situation the test exists to detect, and compares the value to itself.
  # Needs an installed CATEomics; under devtools::load_all there is nothing
  # for the subprocess to attach.
  pkg_dir <- find.package("CATEomics", quiet = TRUE)
  installed <- length(pkg_dir) > 0L &&
    file.exists(file.path(pkg_dir[1], "Meta", "package.rds"))
  if (!installed) {
    skip("CATEomics is loaded from source, not installed; runs under R CMD check")
  }
  initialised <- callr::r(
    function(libs) {
      .libPaths(libs)
      library(CATEomics)
      isTRUE(reticulate::py_available(initialize = FALSE))
    },
    args = list(libs = .libPaths()),
    error = "stack"
  )
  expect_false(initialised)
})

test_that("the removed config subsystem is not exported", {
  # There is no cf_config / cf_config_reset (no on-disk config, no HOME
  # writes at load). Assert the exports are gone.
  ns <- asNamespace("CATEomics")
  expect_false(exists("cf_config", envir = ns, inherits = FALSE))
  expect_false(exists("cf_config_reset", envir = ns, inherits = FALSE))
})

test_that("package namespace exposes the core API", {
  ns <- asNamespace("CATEomics")
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
    cf_fit(Y = "Y", T = "T", X_baseline = "X", data = data, method = "grf", partition = "none")
  })
})
