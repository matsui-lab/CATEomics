# The covariates may come from the omics block alone: an expression-only
# container has no clinical column to name, and a whole layer never enters the
# formula (it travels down the assay path), so the formula names no covariate
# at all. These tests pin that route on all three containers, pin that the
# design matrix is the one the clinical-column route already produced, and pin
# that a call with no omics block left is still refused.

make_se <- function(n = 60L, g = 8L, clinical = FALSE, seed = 11L) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)
  cnt <- matrix(stats::rpois(n * g, 50), nrow = g,
                dimnames = list(paste0("gene", seq_len(g)),
                                paste0("s", seq_len(n))))
  cd <- data.frame(resp = stats::rnorm(n) + 0.5 * rep(0:1, each = n / 2),
                   arm  = rep(0:1, each = n / 2),
                   row.names = colnames(cnt))
  if (clinical) cd$age <- stats::rnorm(n)
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cnt),
    colData = S4Vectors::DataFrame(cd))
}

make_df <- function(n = 60L, g = 8L, clinical = FALSE, seed = 11L) {
  set.seed(seed)
  genes <- as.data.frame(matrix(stats::rpois(n * g, 50), n, g,
                                dimnames = list(NULL, paste0("gene", seq_len(g)))))
  base <- data.frame(resp = stats::rnorm(n) + 0.5 * rep(0:1, each = n / 2),
                     arm  = rep(0:1, each = n / 2))
  if (clinical) base$age <- stats::rnorm(n)
  cbind(base, genes)
}

fit_cheap <- function(...) {
  cf_fit(..., method = "grf", preprocess = FALSE, random_state = 0L,
         eval_nuisance = "ridge")
}

test_that("a SummarizedExperiment with no clinical column fits on the assay", {
  se <- make_se()
  # omics is not named at all: its default for a SummarizedExperiment is "all".
  f1 <- fit_cheap(Y = "resp", T = "arm", data = se)
  expect_s3_class(f1, "cf_model")
  expect_identical(colnames(f1$data$X_preprocessed), paste0("gene", 1:8))
  # The same three ways of asking for the whole layer agree.
  f2 <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = se)
  f3 <- fit_cheap(Y = "resp", T = "arm", omics = "counts", data = se)
  expect_equal(f2$data$X_preprocessed, f1$data$X_preprocessed)
  expect_equal(f3$data$X_preprocessed, f1$data$X_preprocessed)
})

test_that("such a fit names no baseline covariate at all", {
  f <- fit_cheap(Y = "resp", T = "arm", data = make_se())
  # Nothing was named: the design comes entirely from the omics layer, which
  # travels down the assay path rather than being listed as covariates.
  expect_null(f$meta$named_args$X_baseline)
  expect_named(f$meta$named_args$omics, "counts")
  expect_length(f$meta$covariate_names, 8L)
  expect_identical(colnames(f$data$X_preprocessed), paste0("gene", 1:8))
})

test_that("the design matrix equals the clinical-column route on its genes", {
  bare <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = make_se())
  clin <- fit_cheap(Y = "resp", T = "arm", X_baseline = "age", omics = "all",
                    data = make_se(clinical = TRUE))
  genes <- paste0("gene", 1:8)
  expect_identical(colnames(bare$data$X_preprocessed), genes)
  expect_identical(setdiff(colnames(clin$data$X_preprocessed), genes), "age")
  # The central z-score is per column, so dropping one column cannot move the
  # others: the gene block must be bit-identical.
  expect_equal(bare$data$X_preprocessed[, genes],
               clin$data$X_preprocessed[, genes], tolerance = 0)
})

test_that("naming the features gives the same fit as taking the layer", {
  se <- make_se()
  whole <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = se)
  named <- fit_cheap(Y = "resp", T = "arm", omics = paste0("gene", 1:8),
                     data = se)
  expect_equal(named$data$X_preprocessed, whole$data$X_preprocessed,
               tolerance = 0)
})

test_that("a data.frame with no clinical column fits on the omics columns", {
  df <- make_df()
  # For a data.frame the omics columns are the assay matrix and never enter
  # the formula, so both ways of naming them leave the formula empty.
  f1 <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = df)
  f2 <- fit_cheap(Y = "resp", T = "arm", omics = paste0("gene", 1:8), data = df)
  expect_identical(colnames(f1$data$X_preprocessed), paste0("gene", 1:8))
  expect_equal(f2$data$X_preprocessed, f1$data$X_preprocessed, tolerance = 0)
  clin <- fit_cheap(Y = "resp", T = "arm", X_baseline = "age", omics = "all",
                    data = make_df(clinical = TRUE))
  genes <- paste0("gene", 1:8)
  expect_equal(f1$data$X_preprocessed[, genes],
               clin$data$X_preprocessed[, genes], tolerance = 0)
})

test_that("a MultiAssayExperiment with no clinical column fits on its layers", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  se <- make_se()
  m <- SummarizedExperiment::assay(se, "counts")
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(rnaseq = SummarizedExperiment::SummarizedExperiment(
                    assays = list(e = m))),
    colData = SummarizedExperiment::colData(se))
  f <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = mae)
  expect_s3_class(f, "cf_model")
  # MAE columns keep their layer prefix.
  expect_identical(colnames(f$data$X_preprocessed), paste0("rnaseq_gene", 1:8))
})

test_that("predict(newdata) works on the omics-only route", {
  se <- make_se()
  f <- fit_cheap(Y = "resp", T = "arm", omics = "all", data = se)
  p <- predict(f, newdata = se)
  expect_length(p$tau_hat, ncol(se))
  expect_true(all(is.finite(p$tau_hat)))
})

test_that("a call with no omics block left is still refused", {
  se <- make_se()
  # `omics = "none"` with nothing in X_baseline has no covariates anywhere.
  expect_error(fit_cheap(Y = "resp", T = "arm", omics = "none", data = se),
               class = "cateomics_no_X_supplied")
  # A data.frame defaults to `omics = "none"`, so the same holds there.
  expect_error(fit_cheap(Y = "resp", T = "arm", data = make_df()),
               class = "cateomics_no_X_supplied")
  # The message must not ask again for the argument the caller used: it names
  # the state (the omics block is out), not `omics`.
  err <- tryCatch(fit_cheap(Y = "resp", T = "arm", omics = "none", data = se),
                  error = function(e) e)
  expect_match(conditionMessage(err), "omics block", fixed = TRUE)
  expect_no_match(conditionMessage(err), "must be supplied via", fixed = TRUE)
})

test_that("the parser still refuses an empty covariate part by default", {
  dat <- data.frame(Y = stats::rnorm(20), T = rep(0:1, 10),
                    X1 = stats::rnorm(20))
  # The relaxation is opt-in: only the callers that build an assay matrix pass
  # it, so a bare formula with nothing after `|` is refused as before.
  expect_error(parse_cf_formula(Y ~ T | 1, dat), class = "cateomics_no_covariate")
  p <- parse_cf_formula(Y ~ T | 1, dat, allow_no_covariate = TRUE)
  expect_length(p$covariate_names, 0L)
  expect_identical(dim(p$X), c(20L, 0L))
})

test_that("cf_compare takes the same route", {
  se <- make_se()
  cmp <- cf_compare(Y = "resp", T = "arm", omics = "all", data = se,
                    methods = "grf", preprocess = FALSE, random_state = 0L,
                    eval_nuisance = "ridge")
  expect_s3_class(cmp, "cateomics_compare")
  expect_identical(cmp$metadata$n_covariates, 8L)
})
