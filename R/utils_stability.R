# Shared helpers for stability diagnostics and ensemble.
#
# This file groups two helpers:
#
# - .cateomics_clean_oof_tau_matrix: drop all-NA / zero-sd columns +
#   complete.cases rows of the K-method OOF tau matrix returned by
#   cf_compare(...)$metadata$oof_tau_matrix.
# - .cateomics_tau_correlation_knee: pairwise correlation K x K + mean
#   off-diagonal ordering + knee detection. Used by cf_ensemble(method = "cba").
#
# The cleaner is called by cf_ensemble, cf_model_stability and
# cf_benefit_groups; both helpers operate purely on numeric matrices.


#' Clean the OOF tau matrix (drop failed methods + incomplete rows).
#'
#' Shared by `cf_ensemble`, `cf_model_stability` and `cf_benefit_groups` so
#' they use the same cleaning contract:
#'
#' 1. Drop columns whose values are all NA or whose standard deviation is 0
#'    (failed methods from cf_compare are stored as all-NA columns).
#' 2. Drop rows containing any remaining NA (`complete.cases`) so downstream
#'    correlation / QP / knee operations see a dense numeric matrix.
#'
#' The function preserves column names where available so callers can map
#' the cleaned matrix back to the original method labels.
#'
#' @param oof_tau_matrix Numeric `n x K` matrix from
#'   `cf_compare(...)$metadata$oof_tau_matrix`.
#' @return A list with elements:
#'   * `mat`              : the cleaned matrix (or `NULL` when the input is
#'                          unusable);
#'   * `kept_methods`     : character vector of column labels retained;
#'   * `dropped_methods`  : character vector of column labels dropped
#'                          (all-NA or zero-sd);
#'   * `n_rows_kept`      : integer, number of rows retained;
#'   * `n_rows_dropped`   : integer, number of rows removed by
#'                          `complete.cases`;
#'   * `status`           : character, "ok" / "no_input" / "no_columns" /
#'                          "no_rows" / "single_column".
#'
#' @keywords internal
#' @noRd
.cateomics_clean_oof_tau_matrix <- function(oof_tau_matrix) {
  if (is.null(oof_tau_matrix) || !is.matrix(oof_tau_matrix) ||
      !is.numeric(oof_tau_matrix)) {
    return(list(
      mat            = NULL,
      kept_methods   = character(0),
      dropped_methods = character(0),
      n_rows_kept    = 0L,
      n_rows_dropped = 0L,
      status         = "no_input"
    ))
  }

  all_methods <- colnames(oof_tau_matrix)
  if (is.null(all_methods)) {
    all_methods <- as.character(seq_len(ncol(oof_tau_matrix)))
  }

  # Step 1: drop columns with all NA or sd == 0.
  keep <- vapply(seq_len(ncol(oof_tau_matrix)), function(j) {
    col <- oof_tau_matrix[, j]
    any(!is.na(col)) && stats::sd(col, na.rm = TRUE) > 0
  }, logical(1))

  kept_methods    <- all_methods[keep]
  dropped_methods <- all_methods[!keep]

  if (!any(keep)) {
    return(list(
      mat            = NULL,
      kept_methods   = character(0),
      dropped_methods = dropped_methods,
      n_rows_kept    = 0L,
      n_rows_dropped = 0L,
      status         = "no_columns"
    ))
  }

  mat <- oof_tau_matrix[, keep, drop = FALSE]

  # Step 2: drop incomplete rows.
  cc <- stats::complete.cases(mat)
  n_rows_dropped <- sum(!cc)
  mat <- mat[cc, , drop = FALSE]
  n_rows_kept <- nrow(mat)

  if (n_rows_kept == 0L) {
    return(list(
      mat            = NULL,
      kept_methods   = kept_methods,
      dropped_methods = dropped_methods,
      n_rows_kept    = 0L,
      n_rows_dropped = n_rows_dropped,
      status         = "no_rows"
    ))
  }

  status <- if (ncol(mat) < 2L) "single_column" else "ok"

  list(
    mat            = mat,
    kept_methods   = kept_methods,
    dropped_methods = dropped_methods,
    n_rows_kept    = n_rows_kept,
    n_rows_dropped = n_rows_dropped,
    status         = status
  )
}


#' Pairwise correlation matrix + kneedle ordering for tau hat columns.
#'
#' Used by `cf_ensemble(method = "cba")` (Consensus-Based Averaging, Machluf
#' et al. 2024) to rank methods by their mean agreement with the others and
#' locate a knee in the sorted mean correlation curve. `cf_model_stability`
#' computes its own correlations and does not call this.
#'
#' Implementation notes:
#'
#' * Kendall is the default; Spearman and Pearson can be chosen instead.
#'   Kendall's pairwise complexity is O(n^2 K^2) against Spearman's
#'   O(n log n K^2).
#' * Knee detection uses the largest first difference of the sorted mean
#'   off-diagonal correlation curve.
#'
#' @param mat Cleaned numeric `n x K` matrix (output of
#'   `.cateomics_clean_oof_tau_matrix(...)$mat`).
#' @param cor_method Character; one of `"kendall"` (default), `"spearman"`,
#'   or `"pearson"`.
#' @param min_knee Integer; minimum number of methods to keep when a knee is
#'   detected (defaults to `2L`; keeps the implementation defensible when K
#'   is small).
#' @return A list with elements:
#'   * `cor_matrix`        : `K x K` correlation matrix (with 1 on the
#'                           diagonal);
#'   * `ordering`          : integer permutation of `seq_len(K)` sorting
#'                           methods by descending mean off-diagonal
#'                           correlation;
#'   * `sorted_mean_cor`   : numeric vector of the sorted mean correlations
#'                           (same length / order as `ordering`);
#'   * `knee_k`            : integer; number of top-ranked methods on the
#'                           "consensus" side of the knee (clamped to
#'                           `>= min_knee` and `<= K`);
#'   * `cor_method`        : the correlation method used.
#'
#' When `K < 2`, returns the same list with `cor_matrix = NULL`,
#' `ordering = NULL`, `sorted_mean_cor = NULL`, `knee_k = NA_integer_`.
#'
#' @keywords internal
#' @noRd
.cateomics_tau_correlation_knee <- function(mat,
                                          cor_method = c("kendall",
                                                         "spearman",
                                                         "pearson"),
                                          min_knee = 2L) {
  cor_method <- match.arg(cor_method)
  min_knee   <- as.integer(min_knee)

  if (is.null(mat) || !is.matrix(mat) || !is.numeric(mat) || ncol(mat) < 2L) {
    return(list(
      cor_matrix      = NULL,
      ordering        = NULL,
      sorted_mean_cor = NULL,
      knee_k          = NA_integer_,
      cor_method      = cor_method
    ))
  }

  cor_matrix <- stats::cor(mat, method = cor_method,
                           use = "pairwise.complete.obs")

  # Mean off-diagonal correlation per method (agreement with others).
  diag(cor_matrix) <- NA_real_
  mean_cor <- rowMeans(cor_matrix, na.rm = TRUE)
  diag(cor_matrix) <- 1

  ordering        <- order(mean_cor, decreasing = TRUE)
  sorted_mean_cor <- mean_cor[ordering]

  knee_k <- .cateomics_locate_knee(sorted_mean_cor, min_knee = min_knee)

  list(
    cor_matrix      = cor_matrix,
    ordering        = ordering,
    sorted_mean_cor = sorted_mean_cor,
    knee_k          = knee_k,
    cor_method      = cor_method
  )
}


#' Locate the knee of a monotonically decreasing curve.
#'
#' Internal helper for `.cateomics_tau_correlation_knee`. Implements a simple
#' kneedle-style heuristic without external dependencies:
#'
#' 1. Compute first differences of the sorted (descending) curve.
#' 2. The knee index is the position of the largest drop, i.e. the number of
#'    methods before it.
#' 3. Clamp the index to `[min_knee, length(x)]`.
#'
#' Computed here rather than through an external knee-detection package.
#'
#' @param x Numeric vector of sorted mean correlations (descending).
#' @param min_knee Integer; minimum number of methods to keep.
#' @return Integer in `[min_knee, length(x)]`.
#' @keywords internal
#' @noRd
.cateomics_locate_knee <- function(x, min_knee = 2L) {
  n <- length(x)
  if (n <= min_knee) return(n)
  diffs <- diff(x)
  diffs[!is.finite(diffs)] <- 0
  if (length(diffs) < 1L) return(n)
  drop_idx <- which.min(diffs)
  knee_k <- as.integer(drop_idx)
  knee_k <- max(knee_k, min_knee)
  knee_k <- min(knee_k, n)
  knee_k
}
