# =============================================================================
# test-sup-6-meta-buckets.R
# =============================================================================
#
# Metadata bucket structure (provenance / estimator). The top-level
# identifiers kept are n / p / covariate_names; the variables the caller named
# are in named_args. The validator enforces the bucket field presence.

test_that("create_cf_meta returns 4-bucket structure with required fields", {
  meta <- create_cf_meta(
    n              = 100L,
    p              = 1L,
    covariate_names = "X1",
    method         = "grf",
    random_state   = 42L
  )

  # Top-level identifiers preserved. The method name is not among them: the
  # result carries it at the top level (`fit$method`), and the argument is read
  # here only for the backend name and the package version.
  expect_false("method" %in% names(meta))
  # The outcome and the treatment are named in `named_args`, not copied to the
  # top level: the same string under two names is one more thing to keep in step.
  expect_null(meta$outcome_name)
  expect_null(meta$treatment_name)
  expect_identical(meta$covariate_names, "X1")
  expect_identical(meta$provenance$backend, "grf (R)")

  # Required buckets are provenance / estimator (diagnostics are on-demand via
  # cf_diagnostics).
  expect_true(all(c("provenance", "estimator") %in% names(meta)))

  # The estimand bucket is gone: it held one constant, "CATE", which every
  # method reports, so it recorded nothing. `target_population` went the same
  # way earlier, for the same reason.
  expect_null(meta$estimand)

  # provenance bucket structure.
  expect_true(is.list(meta$provenance))
  expect_true(!is.null(meta$provenance$backend))
  expect_identical(meta$provenance$rng_seed, 42L)
  expect_true(is.list(meta$provenance$software_versions))
  # Every package that computed something is in `packages`, CATEomics included,
  # so there is one place to look rather than a field beside the list.
  expect_null(meta$provenance$software_versions$CATEomics)
  expect_true(!is.null(meta$provenance$software_versions$packages$CATEomics))

  # estimator bucket (method-specific config, partial populated).
  expect_true(is.list(meta$estimator))
})

test_that(".infer_backend dispatches across the CATE methods", {
  expect_match(CATEomics:::.infer_backend("grf"),       "grf")
  expect_match(CATEomics:::.infer_backend("bcf"),       "bcf")
  expect_match(CATEomics:::.infer_backend("drlearner"), "econml.dr")
  # xlearner (like s/t-learner) maps to the econml.metalearners backend, which
  # is distinct from drlearner's econml.dr label.
  expect_match(CATEomics:::.infer_backend("xlearner"),  "metalearners")
  # ganite is implemented in this package, so its label has to name that
  # module rather than any upstream package.
  expect_match(CATEomics:::.infer_backend("ganite"), "ganite_impl")
  expect_match(CATEomics:::.infer_backend("cevae"),     "pyro")
  expect_match(CATEomics:::.infer_backend("tarnet"),    "catenets")
  expect_match(CATEomics:::.infer_backend("dragonnet"), "catenets")
  expect_identical(CATEomics:::.infer_backend("unknown_method"), "unknown")
})

test_that("validate_cateomics_result enforces the bucket field presence", {
  # Four-pillar layout (res / model / meta / data), IQR summary; the meta set
  # does not require an assumption_check bucket.
  base <- list(
    method = "grf",
    res = list(
      tau_hat = numeric(5), ate = 0.1, mu0_hat = numeric(5), mu1_hat = numeric(5),
      summary = list(ate = 0.1, ate_ci_lower = 0, ate_ci_upper = 0.2,
                      tau_mean = 0, tau_sd = 1,
                      tau_iqr = list(q0 = -1, q25 = -0.5, q50 = 0,
                                     q75 = 0.5, q100 = 1))
    ),
    model = list(),
    meta = list(
      n = 5L, p = 1L, covariate_names = "X1",
      method = "grf",
      named_args = list(Y = "Y", T = "T", X_baseline = "X1"),
      provenance = list(),
      estimator  = list()
    ),
    data = list()
  )
  class(base) <- "cf_model"

  # Valid result: should pass.
  expect_invisible(validate_cateomics_result(base))

  # Missing required bucket: should fail.
  bad2 <- base
  bad2$meta$provenance <- NULL
  expect_error(validate_cateomics_result(bad2),
                 class = "cateomics_invalid_result")

  bad3 <- base
  bad3$meta$estimator <- NULL
  expect_error(validate_cateomics_result(bad3),
                 class = "cateomics_invalid_result")

  # A result with no per-unit estimate at all: should fail. This is what the
  # estimand field used to gate.
  bad4 <- base
  bad4$res$tau_hat <- rep(NA_real_, length(bad4$res$tau_hat))
  expect_error(validate_cateomics_result(bad4),
                 class = "cateomics_invalid_result")
})

test_that("cf_fit produces the meta buckets end-to-end (grf)", {
  skip_if_not_installed("grf")
  set.seed(7)
  df <- data.frame(
    Y = rnorm(40), T = rbinom(40, 1, 0.5),
    X1 = rnorm(40), X2 = rnorm(40)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y = "Y", T = "T", X_baseline = c("X1", "X2"), data = df, method = "grf",
            partition = "none", random_state = 42L)
  ))
  expect_true(all(c("provenance", "estimator") %in% names(fit$meta)))
  expect_null(fit$meta$estimand)
  expect_match(fit$meta$provenance$backend, "grf")
})
