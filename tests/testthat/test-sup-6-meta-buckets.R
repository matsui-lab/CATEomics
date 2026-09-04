# =============================================================================
# test-sup-6-meta-buckets.R
# =============================================================================
#
# Metadata 4-bucket structure (estimand / assumption_check / provenance /
# estimator). Top-level identifiers (formula / n / p / outcome_name /
# treatment_name / covariate_names / method) are kept. The validator enforces
# the 4-bucket field presence and the estimand$type value range.

test_that("create_cf_meta returns 4-bucket structure with required fields", {
  meta <- create_cf_meta(
    formula        = Y ~ T | X1,
    n              = 100L,
    p              = 1L,
    outcome_name   = "Y",
    treatment_name = "T",
    covariate_names = "X1",
    method         = "grf",
    estimand       = "ATE",
    random_state   = 42L
  )

  # Top-level identifiers preserved.
  expect_identical(meta$method, "grf")
  expect_identical(meta$outcome_name, "Y")

  # Required buckets are estimand / provenance / estimator (diagnostics are
  # on-demand via cf_diagnostics).
  expect_true(all(c("estimand", "provenance",
                     "estimator") %in% names(meta)))

  # estimand bucket structure.
  expect_true(is.list(meta$estimand))
  expect_identical(meta$estimand$type, "ATE")
  expect_identical(meta$estimand$target_population, "all")

  # provenance bucket structure.
  expect_true(is.list(meta$provenance))
  expect_true(!is.null(meta$provenance$backend))
  expect_identical(meta$provenance$rng_seed, 42L)
  expect_true(is.list(meta$provenance$software_versions))
  expect_true(!is.null(meta$provenance$software_versions$cfomics))

  # estimator bucket (method-specific config, partial populated).
  expect_true(is.list(meta$estimator))
})

test_that(".infer_backend dispatches across the CATE methods", {
  expect_match(cfomics:::.infer_backend("grf"),       "grf")
  expect_match(cfomics:::.infer_backend("bcf"),       "bcf")
  expect_match(cfomics:::.infer_backend("drlearner"), "econml.dr")
  # xlearner (like s/t-learner) maps to the econml.metalearners backend, which
  # is distinct from drlearner's econml.dr label.
  expect_match(cfomics:::.infer_backend("xlearner"),  "metalearners")
  # ganite uses the PyPI ganite (PyTorch) backend.
  expect_match(cfomics:::.infer_backend("ganite"), "ganite|PyTorch")
  expect_match(cfomics:::.infer_backend("cevae"),     "pyro")
  expect_match(cfomics:::.infer_backend("tarnet"),    "catenets")
  expect_match(cfomics:::.infer_backend("dragonnet"), "catenets")
  expect_identical(cfomics:::.infer_backend("unknown_method"), "unknown")
})

test_that("validate_cfomics_result enforces 4-bucket field presence", {
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
      formula = Y ~ T | X1, n = 5L, p = 1L,
      outcome_name = "Y", treatment_name = "T", covariate_names = "X1",
      method = "grf",
      estimand = list(type = "ATE", target_population = "all"),
      provenance = list(),
      estimator  = list()
    ),
    data = list()
  )
  class(base) <- "cf_model"

  # Valid result: should pass.
  expect_invisible(validate_cfomics_result(base))

  # Missing required bucket: should fail.
  bad2 <- base
  bad2$meta$provenance <- NULL
  expect_error(validate_cfomics_result(bad2),
                 class = "cfomics_invalid_result")

  bad3 <- base
  bad3$meta$estimator <- NULL
  expect_error(validate_cfomics_result(bad3),
                 class = "cfomics_invalid_result")

  # estimand as a flat string (old structure): should fail.
  bad4 <- base
  bad4$meta$estimand <- "ATE"
  expect_error(validate_cfomics_result(bad4),
                 class = "cfomics_invalid_result")

  # estimand$type invalid: should fail.
  bad5 <- base
  bad5$meta$estimand$type <- "NOTAVALIDESTIMAND"
  expect_error(validate_cfomics_result(bad5),
                 class = "cfomics_invalid_result")
})

test_that("cf_fit produces 4-bucket meta end-to-end (grf)", {
  skip_if_not_installed("grf")
  set.seed(7)
  df <- data.frame(
    Y = rnorm(40), T = rbinom(40, 1, 0.5),
    X1 = rnorm(40), X2 = rnorm(40)
  )
  fit <- suppressMessages(suppressWarnings(
    cf_fit(Y ~ T | X1 + X2, data = df, method = "grf",
            partition = "none", random_state = 42L)
  ))
  expect_true(all(c("estimand",
                     "provenance", "estimator") %in% names(fit$meta)))
  expect_identical(fit$meta$estimand$type, "CATE")
  expect_match(fit$meta$provenance$backend, "grf")
})
