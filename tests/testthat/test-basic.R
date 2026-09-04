test_that("parse_cf_formula validation", {
  # factor / character covariates are auto one-hot encoded via
  # `.expand_factor_columns`. C is a 3-level factor here (a 10-level factor
  # would give a degenerate 1-level-per-sample factor at n = 10).
  set.seed(7)
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X = rnorm(20),
    C = sample(c("a", "b", "c"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  fit_c <- cf_fit(Y ~ T | C, data = data, method = "grf", partition = "none")
  expect_s3_class(fit_c, "cf_model")
  expect_setequal(fit_c$meta$factor_levels$C, c("a", "b", "c"))
})
