test_that("cf_example_data returns a valid SummarizedExperiment with ground truth", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  suppressPackageStartupMessages(library(SummarizedExperiment))

  se <- cf_example_data(n = 60L, p = 200L, n_prog = 20L, n_pred = 10L, seed = 1L)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(200L, 60L))

  # Counts: integer, non-negative (valid for lognorm / edgeR).
  cts <- SummarizedExperiment::assay(se, "counts")
  expect_true(is.integer(cts))
  expect_true(all(cts >= 0))

  # Real gene symbols as row names, no duplicates.
  expect_type(rownames(se), "character")
  expect_false(anyDuplicated(rownames(se)) > 0)

  # colData carries treatment / outcome / clinical / truth.
  cd <- SummarizedExperiment::colData(se)
  expect_true(all(c("treatment", "outcome", "age", "sex", "cf_tau_true") %in%
                    colnames(cd)))
  expect_true(all(cd$treatment %in% c(0L, 1L)))
  expect_true(is.factor(cd$sex))
  expect_setequal(levels(cd$sex), c("F", "M"))
  expect_length(cd$cf_tau_true, 60L)

  # rowData truth flags: exactly n_pred predictive, n_prog prognostic, disjoint.
  rd <- SummarizedExperiment::rowData(se)
  expect_true(all(c("is_prognostic", "is_predictive") %in% colnames(rd)))
  expect_equal(sum(rd$is_predictive), 10L)
  expect_equal(sum(rd$is_prognostic), 20L)
  expect_equal(sum(rd$is_predictive & rd$is_prognostic), 0L)

  # metadata: embedded pathway + truth.
  md <- S4Vectors::metadata(se)
  expect_true(all(c("cf_pathway", "cf_predictive_genes", "cf_prognostic_genes",
                    "cf_true_ate", "cf_snr") %in% names(md)))
  expect_equal(md$cf_pathway, "HALLMARK_APOPTOSIS")
  expect_length(md$cf_predictive_genes, 10L)
  # predictive genes are exactly the is_predictive rows.
  expect_setequal(md$cf_predictive_genes, rownames(se)[rd$is_predictive])
})

test_that("cf_example_data moments follow the documented formulas", {
  skip_if_not_installed("SummarizedExperiment")
  suppressPackageStartupMessages(library(SummarizedExperiment))

  # Nothing is rescaled after the fact, so sd(tau) = effect_size * sqrt(n_pred
  # + 0.25) and sd(mu0) = effect_size * sqrt(n_prog + 1) up to sampling noise
  # in the z-scored expression.
  se <- cf_example_data(n = 300L, p = 300L, n_prog = 5L, n_pred = 5L,
                        effect_size = 1.5, ate = 0, noise_sd = 1, seed = 7L)
  cd <- SummarizedExperiment::colData(se)
  tau <- cd$cf_tau_true
  mu0_plus_eps <- cd$outcome - cd$treatment * tau
  expect_equal(stats::sd(tau), 1.5 * sqrt(5.25), tolerance = 0.2)
  expect_equal(stats::sd(mu0_plus_eps), sqrt((1.5 * sqrt(6))^2 + 1),
               tolerance = 0.2)
  expect_equal(mean(tau), 0, tolerance = 0.5)

  # More predictive genes means more heterogeneity, not the same amount spread
  # thinner -- the property the rescaled version did not have.
  wide <- cf_example_data(n = 300L, p = 300L, n_pred = 20L, seed = 7L)
  expect_gt(stats::sd(SummarizedExperiment::colData(wide)$cf_tau_true),
            stats::sd(tau))

  # `ate` is the constant term of tau, and with ate = 0 about half are harmed.
  shifted <- cf_example_data(n = 300L, p = 300L, ate = 2, seed = 7L)
  expect_equal(mean(SummarizedExperiment::colData(shifted)$cf_tau_true), 2,
               tolerance = 0.5)
  expect_equal(mean(tau < 0), 0.5, tolerance = 0.15)

  # `noise_sd` moves the noise only, so the realised Var(tau)/Var(eps) falls.
  loud <- cf_example_data(n = 300L, p = 300L, noise_sd = 3, seed = 7L)
  expect_lt(S4Vectors::metadata(loud)$cf_snr, S4Vectors::metadata(se)$cf_snr)

  # true ATE is the mean individual effect.
  expect_equal(S4Vectors::metadata(se)$cf_true_ate, mean(tau))
})

test_that("cf_example_data is reproducible and leaves the caller's RNG alone", {
  skip_if_not_installed("SummarizedExperiment")
  se  <- cf_example_data(n = 60L, p = 200L, seed = 7L)
  se2 <- cf_example_data(n = 60L, p = 200L, seed = 7L)
  expect_identical(SummarizedExperiment::assay(se, "counts"),
                   SummarizedExperiment::assay(se2, "counts"))
  se3 <- cf_example_data(n = 60L, p = 200L, seed = 8L)
  expect_false(identical(SummarizedExperiment::assay(se, "counts"),
                         SummarizedExperiment::assay(se3, "counts")))

  set.seed(123L); a <- runif(1)
  set.seed(123L); invisible(cf_example_data(n = 40L, p = 60L)); b <- runif(1)
  expect_equal(a, b)
})

test_that("the signal genes stay inside their pathways", {
  skip_if_not_installed("SummarizedExperiment")
  pool <- CATEomics:::.cateomics_example_genes
  # Asking for more signal genes than the curated head used to top up from the
  # noise pool, which made `cf_pathway` only partly true; now the pools hold
  # the whole pathway and going past it is an error.
  se <- cf_example_data(n = 40L, p = 400L, n_pred = 40L, n_prog = 40L, seed = 1L)
  md <- S4Vectors::metadata(se)
  expect_true(all(md$cf_predictive_genes %in% pool$predictive))
  expect_true(all(md$cf_prognostic_genes %in% pool$prognostic))
  expect_length(md$cf_predictive_genes, 40L)
  # The ceilings are the pathway sizes the help names, and both of them work:
  # recomputing them here from the pools would mirror an off-by-one in the
  # implementation instead of catching it.
  expect_equal(length(pool$predictive), 161L)
  expect_equal(length(pool$prognostic), 194L)
  expect_s4_class(cf_example_data(n = 10L, p = 400L, n_pred = 161L), "SummarizedExperiment")
  expect_s4_class(cf_example_data(n = 10L, p = 400L, n_prog = 194L), "SummarizedExperiment")
  expect_error(cf_example_data(n_pred = 162L, p = 1000L), class = "cateomics_arg_error")
  expect_error(cf_example_data(n_prog = 195L, p = 1000L), class = "cateomics_arg_error")
})

test_that("cf_example_data validates arguments", {
  expect_error(cf_example_data(n = 2L), class = "cateomics_arg_error")
  expect_error(cf_example_data(p = 5L), class = "cateomics_arg_error")
  expect_error(cf_example_data(p = 25L, n_prog = 20L, n_pred = 10L),
               class = "cateomics_arg_error")
  expect_error(cf_example_data(effect_size = 0), class = "cateomics_arg_error")
  expect_error(cf_example_data(noise_sd = -1), class = "cateomics_arg_error")
  expect_error(cf_example_data(ate = NA), class = "cateomics_arg_error")
  # `p` beyond what the generator offers is named, not left to sample(). The
  # limit is a published constant, so it is written out here: 4000 has to work
  # for every signal count, including the smallest.
  expect_s4_class(cf_example_data(n = 10L, p = 4000L), "SummarizedExperiment")
  expect_s4_class(cf_example_data(n = 10L, p = 4000L, n_pred = 1L, n_prog = 1L),
                  "SummarizedExperiment")
  expect_error(cf_example_data(n = 10L, p = 4001L), class = "cateomics_arg_error")
  # A count given as a fraction used to be truncated in silence, and a seed
  # that cannot be coerced raised base R's error after the RNG-restoring
  # on.exit() had been registered, which added a second, unrelated warning.
  expect_error(cf_example_data(n_pred = 5.7), class = "cateomics_arg_error")
  expect_error(cf_example_data(n = 150.5), class = "cateomics_arg_error")
  expect_error(cf_example_data(seed = 1.9), class = "cateomics_arg_error")
  expect_error(cf_example_data(seed = NA), class = "cateomics_arg_error")
  expect_error(cf_example_data(seed = "a"), class = "cateomics_arg_error")
  # A vector where a single number belongs reached `||` and died there.
  expect_error(cf_example_data(n = c(20L, 30L)), class = "cateomics_arg_error")
  expect_error(cf_example_data(p = c(600L, 700L)), class = "cateomics_arg_error")
  # The recorded seed is the integer that was used, not the value as given.
  expect_identical(
    S4Vectors::metadata(cf_example_data(n = 20L, p = 50L, seed = 3))$cf_params$seed, 3L)
})

test_that("cf_example_data feeds cf_fit (R-native) end to end", {
  skip_if_not_installed("SummarizedExperiment")
  se <- cf_example_data(n = 80L, p = 150L, seed = 3L)
  fit <- suppressWarnings(suppressMessages(
    cf_fit(data = se, Y = "outcome", T = "treatment",
           X_baseline = c("age", "sex"), method = "grf",
           partition = "none", random_state = 1L)))
  expect_s3_class(fit, "cf_model")
  expect_true(is.finite(fit$res$ate))
  expect_length(fit$res$tau_hat, 80L)
})
