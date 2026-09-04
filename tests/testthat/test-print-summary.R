# print.cf_model / summary.cf_model / print.summary.cf_model

make_fitted <- function() {
  skip_if_not_installed("grf")
  set.seed(2026)
  data <- data.frame(
    Y  = rnorm(50),
    T  = rbinom(50, 1, 0.5),
    X1 = rnorm(50)
  )
  cf_fit(Y ~ T | X1, data = data, method = "grf", partition = "none")
}

test_that("summary.cf_model returns a structured list", {
  fit <- make_fitted()
  s <- summary(fit)
  expect_s3_class(s, "summary.cf_model")
  expect_named(s, c("method", "n", "p", "estimand", "ate", "ate_ci",
                    "ci_method", "tau_summary"))
  expect_length(s$ate_ci, 2)
  expect_named(s$tau_summary, c("mean", "sd", "iqr"))
})

test_that("print.summary.cf_model produces formatted output", {
  fit <- make_fitted()
  s <- summary(fit)
  expect_output(print(s), "Causal Inference Model Summary")
  expect_output(print(s), "ATE 95% CI")
})

test_that("print.cf_model returns the object invisibly", {
  fit <- make_fitted()
  expect_output(p <- print(fit), "<cf_model:")
  expect_identical(p, fit)
})

test_that("summary.cf_model aborts when fit$res$summary is NULL", {
  fit <- make_fitted()
  fit$res$summary <- NULL
  expect_error(
    summary(fit),
    class = "cfomics_invalid_result"
  )
})
