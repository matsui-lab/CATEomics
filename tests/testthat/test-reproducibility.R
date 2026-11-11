test_that("Results are reproducible with same random_state (grf)", {
  skip_if_not_installed("grf")

  set.seed(999)
  data <- data.frame(
    Y = rnorm(100),
    T = rbinom(100, 1, 0.5),
    X1 = rnorm(100)
  )

  fit1 <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", random_state = 42)
  fit2 <- cf_fit(Y = "Y", T = "T", X_baseline = "X1", data = data, method = "grf", random_state = 42)

  expect_equal(fit1$res$ate, fit2$res$ate)
  expect_equal(fit1$res$tau_hat, fit2$res$tau_hat)
})

## The 4-bucket metadata layout stores `software_versions` under
## `meta$provenance$software_versions`; these tests verify that nested path
## end-to-end.

test_that("meta contains software versions (grf)", {
  skip_if_not_installed("grf")

  set.seed(123)
  ## Use large enough n and p so grf::regression_forest tune_forest does not
  ## abort with "honesty fraction is too close to 1 or 0" (grf small-n +
  ## single-feature edge case).
  data <- data.frame(Y = rnorm(300), T = rbinom(300, 1, 0.5),
                     X1 = rnorm(300), X2 = rnorm(300), X3 = rnorm(300))
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3"), data = data, method = "grf",
                random_state = 1L)

  sv <- fit$meta$provenance$software_versions
  expect_true("software_versions" %in% names(fit$meta$provenance))
  expect_true("R" %in% names(sv))
  # Every package that computed something is in `packages`, CATEomics included.
  expect_true("CATEomics" %in% names(sv$packages))
  expect_true("timestamp" %in% names(sv))
  expect_match(sv$R, "^[0-9]+\\.[0-9]+")
  expect_s3_class(sv$timestamp, "POSIXct")
})

test_that("software_versions includes method-specific packages", {
  skip_if_not_installed("grf")

  set.seed(123)
  ## Large enough n and p (same grf small-n edge case as above).
  data <- data.frame(Y = rnorm(300), T = rbinom(300, 1, 0.5),
                     X1 = rnorm(300), X2 = rnorm(300), X3 = rnorm(300))
  fit <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2", "X3"), data = data, method = "grf",
                random_state = 1L)

  expect_true("packages" %in% names(fit$meta$provenance$software_versions))
  expect_true("grf" %in% names(fit$meta$provenance$software_versions$packages))
})

# ---------------------------------------------------------------------------
# RNG-pollution guard.
#
# For each method the contract is:
#   (a) `.Random.seed` is unchanged after `cf_fit(..., random_state = s)`
#   (b) Running cf_fit twice with the same seed yields identical tau/ate
#
# grf / cevae are the in-scope methods here. GANITE pollutes tf.random by
# design (see methods_ganite.R roxygen caveat), so it intentionally does not
# appear here.
# ---------------------------------------------------------------------------

.make_reproducibility_data <- function(n = 100L, seed = 7L) {
  set.seed(seed)
  data.frame(
    Y  = stats::rnorm(n),
    T  = stats::rbinom(n, 1L, 0.5),
    X1 = stats::rnorm(n),
    X2 = stats::rnorm(n)
  )
}

.snapshot_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
}

test_that("cf_fit(grf) does not pollute global RNG and is reproducible", {
  skip_if_not_installed("grf")
  dat <- .make_reproducibility_data()

  set.seed(11L)
  before <- .snapshot_seed()
  fit1 <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = dat, method = "grf",
                 random_state = 42L)
  expect_identical(.snapshot_seed(), before)

  set.seed(11L)
  before2 <- .snapshot_seed()
  fit2 <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = dat, method = "grf",
                 random_state = 42L)
  expect_identical(.snapshot_seed(), before2)

  expect_equal(fit1$res$tau_hat, fit2$res$tau_hat)
  expect_equal(fit1$res$ate, fit2$res$ate)
})

test_that("cf_fit(cevae) does not pollute global RNG and is reproducible", {
  skip_if_no_cf_python()
  dat <- .make_reproducibility_data(n = 80L)

  set.seed(11L)
  before <- .snapshot_seed()
  fit1 <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = dat, method = "cevae",
          random_state = 42L, params = list(num_epochs = 3L,
              latent_dim = 4L, hidden_dim = 16L, num_samples = 4L,
              batch_size = 32L))
  expect_identical(.snapshot_seed(), before)

  set.seed(11L)
  before2 <- .snapshot_seed()
  fit2 <- cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = dat, method = "cevae",
          random_state = 42L, params = list(num_epochs = 3L,
              latent_dim = 4L, hidden_dim = 16L, num_samples = 4L,
              batch_size = 32L))
  expect_identical(.snapshot_seed(), before2)

  expect_equal(fit1$res$ate, fit2$res$ate)
  expect_equal(fit1$res$tau_hat, fit2$res$tau_hat)
})
