## Tests for the data-class-aware X-slot resolution.
##
## `.resolve_X_columns()` enforces the slot-shape rules (df: any dedicated +
## no X_assays; SE: at most one dedicated, no X_assays; MAE: anything goes).
## Column-existence validation remains the job of the downstream
## `parse_cf_formula` / `as_cf_data` path and is not duplicated here.

library(testthat)

.make_df <- function(n = 120L, p = 4L, seed = 42L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  T <- stats::rbinom(n, 1, 0.5)
  Y <- as.numeric(T) * 1.5 + rowSums(X) + 0.1 * stats::rnorm(n)
  df <- data.frame(Y = Y, T = T, X)
  colnames(df) <- c("Y", "T", paste0("X", seq_len(p)))
  df
}

test_that("df: X_assays is rejected (no assay concept)", {
  df <- .make_df()
  expect_error(
    cf_fit(Y = "Y", T = "T", X_baseline = "X1",
           X_assays = list(atac = "X2"),
           data = df, method = "grf"),
    class = "cfomics_X_assays_unsupported_df"
  )
})

test_that("df: any combination of dedicated X_<omics> slots is accepted as plain column lookups", {
  df <- .make_df()
  # df has no real omics, but the resolver treats all X_<omics> values as
  # column names so passing them as a dedicated slot is fine. Use
  # `tune.parameters = "none"` so grf does not invoke `tune_forest`
  # (which is fragile on the synthetic 4-column toy data).
  fit <- cf_fit(Y = "Y", T = "T",
                X_baseline = c("X1", "X2"),
                X_rnaseq   = c("X3", "X4"),
                data = df, method = "grf",
                random_state = 42L, partition = "none",
                tune.parameters = "none")
  expect_s3_class(fit, "cf_model")
  expect_setequal(fit$meta$named_args$X_baseline, c("X1", "X2"))
  expect_setequal(fit$meta$named_args$X_rnaseq, c("X3", "X4"))
})

test_that("SE: at most one dedicated X_<omics> is allowed (multi-dedicated errors)", {
  skip_if_not_installed("SummarizedExperiment")
  # Build a tiny SE with colData + one assay, then probe the multi-dedicated
  # validation path. We do not need the fit to succeed; the validation runs
  # before parse_cf_formula / as_cf_data, so an early error is the goal.
  n  <- 40L
  ng <- 6L
  cd <- S4Vectors::DataFrame(
    Y  = stats::rnorm(n),
    T  = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n)
  )
  mat <- matrix(stats::rnorm(ng * n), nrow = ng, ncol = n,
                dimnames = list(paste0("G", seq_len(ng)), paste0("S", seq_len(n))))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = mat),
    colData = cd
  )

  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = "age",
           X_rnaseq   = c("G1", "G2"),
           X_mirna    = c("G3"),  # second dedicated → should fail
           data = se, method = "grf"),
    class = "cfomics_multiple_omics_unsupported_se"
  )
})

test_that("SE: X_assays is rejected (single-assay container)", {
  skip_if_not_installed("SummarizedExperiment")
  n  <- 40L
  ng <- 6L
  cd <- S4Vectors::DataFrame(
    Y  = stats::rnorm(n),
    T  = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n)
  )
  mat <- matrix(stats::rnorm(ng * n), nrow = ng, ncol = n,
                dimnames = list(paste0("G", seq_len(ng)), paste0("S", seq_len(n))))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = mat),
    colData = cd
  )

  expect_error(
    cf_fit(Y = "Y", T = "T",
           X_baseline = "age",
           X_assays   = list(atac = c("G1")),
           data = se, method = "grf"),
    class = "cfomics_X_assays_unsupported_se"
  )
})

test_that("MAE: multiple dedicated X_<omics> + X_assays are accepted (validation passes)", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  # Build a tiny MAE; we only need the validation path to clear without
  # error. The downstream parse / fit may or may not succeed depending on
  # bioc_integration handling — the point is that itself does not
  # complain.
  n  <- 40L
  cd_outer <- S4Vectors::DataFrame(
    Y   = stats::rnorm(n),
    T   = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n),
    row.names = paste0("S", seq_len(n))
  )
  make_se <- function(prefix, ng) {
    m <- matrix(stats::rnorm(ng * n), nrow = ng, ncol = n,
                dimnames = list(paste0(prefix, seq_len(ng)),
                                paste0("S", seq_len(n))))
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  }
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rnaseq = make_se("R", 5L),
                       mirna  = make_se("M", 4L)),
    colData     = cd_outer
  )
  # Validation should succeed (no abort) — full pipeline may still error
  # downstream but itself is permissive on MAE.
  expect_silent(
    .resolve_X_columns(
      data       = mae,
      X_rnaseq   = c("R1", "R2"),
      X_mirna    = c("M1"),
      X_assays   = list(extra = c("E1"))
    )
  )
})

test_that("unsupported data class is rejected", {
  expect_error(
    .resolve_X_columns(data = list(foo = 1)),
    class = "cfomics_invalid_data_type"
  )
})
