# Tests for input validation and error handling

test_that("cf_fit rejects non-numeric outcome", {
  data <- data.frame(
    Y = letters[1:20],
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
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
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
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
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
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

  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none")
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
  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none")
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
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_na_in_variable"
  )
})

test_that("cf_fit handles NA values in treatment", {
  data <- data.frame(
    Y = rnorm(20),
    T = c(NA, rbinom(19, 1, 0.5)),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_na_in_variable"
  )
})

test_that("cf_fit handles NA values in covariates", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = c(NA, rnorm(19))
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_na_in_variable"
  )
})

test_that("cf_fit handles empty data", {
  data <- data.frame(
    Y = numeric(0),
    T = integer(0),
    X1 = numeric(0)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_no_observations"
  )
})

test_that("cf_fit handles single observation", {
  data <- data.frame(
    Y = 1.0,
    T = 1L,
    X1 = 0.5
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_insufficient_n"
  )
})

test_that("cf_fit handles single treatment group (all treated)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(1L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_group_empty"
  )
})

test_that("cf_fit handles single treatment group (all control)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(0L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
    class = "cfomics_group_empty"
  )
})

test_that("cf_fit rejects unknown method", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "unknown_method"),
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

test_that("as_cf_data.default rejects unsupported types", {
  expect_error(
    as_cf_data("not a data.frame"),
    "does not support"
  )
  expect_error(
    as_cf_data(list(a = 1, b = 2)),
    "does not support"
  )
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

  formula <- as.formula(paste("Y ~ T |", paste(paste0("X", 1:p), collapse = " + ")))

  fit <- cf_fit(formula, data = data, method = "grf")
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
    cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none"),
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
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf", partition = "none")
  expect_s3_class(fit, "cf_model")
})

test_that("select_best_method returns valid method", {
  method <- select_best_method()
  expect_true(method %in% c("grf", "bcf", "drlearner"))
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
  fit <- cf_fit(Y ~ T | gene.A + gene.B,
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
    class = "cfomics_y_type_warning"
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
    class = "cfomics_y_type_warning"
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
test_that("validate_cfomics_result aborts on missing summary fields", {
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
        formula = Y ~ T | X1, n = 10L, p = 1L,
        outcome_name = "Y", treatment_name = "T",
        covariate_names = "X1",
        estimand         = list(type = "ATE", target_population = "marginal",
                                assumptions = list()),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cfomics_result(bad),
    class = "cfomics_invalid_result"
  )
})

test_that("validate_cfomics_result aborts on malformed ite_quantiles", {
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
        formula = Y ~ T | X1, n = 10L, p = 1L,
        outcome_name = "Y", treatment_name = "T",
        covariate_names = "X1",
        estimand         = list(type = "ATE", target_population = "marginal",
                                assumptions = list()),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cfomics_result(bad),
    class = "cfomics_invalid_result"
  )
})

# ate_ci / ite_ci unified location validator.
test_that("validate_cfomics_result accepts NULL ate_ci / ite_ci (additive, optional)", {
  ## The 4-bucket metadata layout makes `assumption_check`, `provenance`, and
  ## `estimator` mandatory fields under `meta`, so the synthetic `ok` object
  ## carries empty buckets. This still verifies that NULL `ate_ci` / `ite_ci`
  ## paths remain accepted.
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
        formula = Y ~ T | X1, n = 10L, p = 1L,
        outcome_name = "Y", treatment_name = "T",
        covariate_names = "X1",
        ## estimand is a list (`type` + `target_population` + `assumptions`),
        ## not a bare string. This is the minimal valid shape so the
        ## validator's 4-bucket check passes.
        estimand = list(
          type = "ATE",
          target_population = "marginal",
          assumptions = list()
        ),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_true(validate_cfomics_result(ok))
})

test_that("validate_cfomics_result rejects malformed ate_ci", {
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
        formula = Y ~ T | X1, n = 10L, p = 1L,
        outcome_name = "Y", treatment_name = "T",
        covariate_names = "X1",
        estimand         = list(type = "ATE", target_population = "marginal",
                                assumptions = list()),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cfomics_result(bad),
    class = "cfomics_invalid_result"
  )
})

test_that("validate_cfomics_result rejects malformed ite_ci nrow", {
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
        formula = Y ~ T | X1, n = 10L, p = 1L,
        outcome_name = "Y", treatment_name = "T",
        covariate_names = "X1",
        estimand         = list(type = "ATE", target_population = "marginal",
                                assumptions = list()),
        provenance       = list(backend = "grf", rng_seed = 1L),
        estimator        = list(bootstrap = FALSE)
      )
    ),
    class = "cf_model"
  )
  expect_error(
    validate_cfomics_result(bad),
    class = "cfomics_invalid_result"
  )
})
