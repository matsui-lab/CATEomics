# Formula parser strictness.
#
# `parse_cf_formula` uses `Formula::Formula` + `model.part`, so the
# `Y ~ T | X1 + X2` syntax is enforced at parse time rather than silently
# mis-interpreted.

make_dat <- function(n = 30L, seed = 1L) {
  set.seed(seed)
  data.frame(
    Y  = stats::rnorm(n),
    T  = stats::rbinom(n, 1L, 0.5),
    X1 = stats::rnorm(n),
    X2 = stats::rnorm(n),
    Z  = stats::rnorm(n)
  )
}

test_that("Y ~ T | X1 + X2 parses T and X correctly", {
  dat <- make_dat()
  p <- parse_cf_formula(Y ~ T | X1 + X2, dat)
  expect_equal(p$outcome_name, "Y")
  expect_equal(p$treatment_name, "T")
  expect_equal(p$covariate_names, c("X1", "X2"))
  expect_equal(ncol(p$X), 2L)
})

test_that("Y ~ T | X1 (single covariate) parses correctly", {
  dat <- make_dat()
  p <- parse_cf_formula(Y ~ T | X1, dat)
  expect_equal(p$treatment_name, "T")
  expect_equal(p$covariate_names, "X1")
  expect_equal(ncol(p$X), 1L)
})

test_that("formulas without `|` are rejected", {
  dat <- make_dat()
  # A one-part formula has no place to say which variable is the treatment, so
  # it is refused. Assert the class, not the
  # wording: the two-part formula is built internally now, so this guard
  # protects an invariant rather than instructing a caller, and its text is
  # free to say so.
  expect_error(
    parse_cf_formula(Y ~ T + X1 + X2, dat),
    class = "cateomics_formula_rhs_missing"
  )
})

test_that("treatment variable position is respected (edge case)", {
  # `Y ~ X1 | T + X2` means outcome=Y, treatment=X1, covariates=(T, X2).
  # X1 in this dataset is continuous, so it fails the binary-T check —
  # which is the correct behavior: the parser trusts the `|` placement.
  dat <- make_dat()
  expect_error(
    parse_cf_formula(Y ~ X1 | T + X2, dat),
    "binary"
  )

  # If we swap so that the LHS of `|` is actually binary, parsing succeeds
  # and identifies treatment from the `|` placement, not from variable
  # appearance order.
  dat2 <- dat
  dat2$Bin <- stats::rbinom(nrow(dat2), 1L, 0.5)
  p <- parse_cf_formula(Y ~ Bin | X1 + X2, dat2)
  expect_equal(p$treatment_name, "Bin")
  expect_equal(p$covariate_names, c("X1", "X2"))
})

test_that("three-part formula `Y ~ T | X | Z` is rejected", {
  dat <- make_dat()
  expect_error(
    parse_cf_formula(Y ~ T | X1 | Z, dat),
    class = "cateomics_formula_rhs_extra"
  )
})
