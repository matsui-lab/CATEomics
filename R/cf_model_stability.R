#' Diagnose estimator dependence across the methods in a cf_compare result.
#'
#' \code{cf_model_stability()} summarises how much the per-method OOF tau
#' hats agree with each other. It is a \strong{descriptive} diagnostic, not
#' a validation tool: the output must not be used as a model-selection
#' criterion.
#'
#' Three families of agreement are reported:
#'
#' \itemize{
#'   \item Pairwise correlation matrices (Pearson, Spearman, Kendall) over
#'         the K cleaned methods. Spearman is the primary rank-based
#'         measure; Kendall mirrors what \code{cf_ensemble(method = "cba")}
#'         uses internally so the two stay consistent.
#'   \item Top-q high-benefit-set Jaccard index between every pair of methods.
#'   \item Mean off-diagonal correlation per method (primary correlation)
#'         and a heuristic flag for outlier methods (mean agreement below
#'         the cohort median by more than one MAD).
#' }
#'
#' @param compare_result A \code{cfomics_compare} object.
#' @param cor_methods Character vector of correlation methods to compute;
#'   any subset of \code{c("pearson", "spearman", "kendall")}.
#' @param primary_cor The correlation method whose K-vector of mean
#'   off-diagonals drives the outlier flag and the summary statistics.
#' @param benefit_top_q Numeric in (0, 0.5]; quantile threshold for the
#'   top-q high-benefit set definition (default 0.2 = top 20\% by tau hat).
#' @param outlier_mad Numeric; an outlier method is one whose mean
#'   off-diagonal correlation is below the median minus \code{outlier_mad}
#'   times the MAD of mean correlations.
#' @param ... Currently unused.
#'
#' @return An object of class \code{c("cfomics_model_stability", "list")}
#'   with:
#'   \itemize{
#'     \item \code{cor_matrices}: named list of K x K correlation matrices,
#'           one per requested correlation method;
#'     \item \code{benefit_group_jaccard}: K x K matrix of top-q high-benefit Jaccard
#'           values;
#'     \item \code{per_method_tau_correlation}: numeric K-vector of mean
#'           off-diagonal correlations (primary);
#'     \item \code{outlier_methods}: character vector flagged as outliers
#'           (median - outlier_mad x MAD);
#'     \item \code{summary}: list with \code{median_agreement},
#'           \code{min_agreement}, \code{mean_agreement},
#'           \code{n_outliers}, \code{n_methods}, \code{n_rows};
#'     \item \code{meta}: methods_used, methods_dropped, primary_cor,
#'           benefit_top_q, outlier_mad.
#'   }
#'
#' @section Use:
#' Treat the output as descriptive; do not use it for model selection.
#' The intended uses are (a) flagging an outlier
#' method for investigation and (b) gauging how sensitive a downstream
#' claim is to estimator choice.
#'
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf) so agreement is
#' # non-degenerate; known randomised propensity. Runs without Python.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L,  # 1 chain keeps bcf fast
#'                   propensity_score = 0.5, random_state = 1)
#' cf_model_stability(cmp)
#'
#' @export
cf_model_stability <- function(compare_result,
                               cor_methods     = c("pearson", "spearman",
                                                   "kendall"),
                               primary_cor     = "spearman",
                               benefit_top_q = 0.2,
                               outlier_mad     = 2.0,   # Hampel MAD threshold (1.0 over-flags)
                               ...) {
  if (!inherits(compare_result, "cfomics_compare")) {
    rlang::abort(
      "cf_model_stability requires a cfomics_compare object.",
      class = "cfomics_invalid_compare_result"
    )
  }

  cor_methods    <- intersect(cor_methods, c("pearson", "spearman", "kendall"))
  if (length(cor_methods) == 0L) {
    rlang::abort(
      "cor_methods must include at least one of pearson / spearman / kendall.",
      class = "cfomics_stability_no_cor_method"
    )
  }
  primary_cor <- match.arg(primary_cor, c("pearson", "spearman", "kendall"))
  if (!primary_cor %in% cor_methods) cor_methods <- c(cor_methods, primary_cor)
  benefit_top_q <- as.numeric(benefit_top_q)[1L]
  if (!is.finite(benefit_top_q) ||
      benefit_top_q <= 0 || benefit_top_q > 0.5) {
    rlang::abort(
      "benefit_top_q must be in (0, 0.5].",
      class = "cfomics_stability_invalid_top_q"
    )
  }
  outlier_mad <- as.numeric(outlier_mad)[1L]

  cleaned <- .cfomics_clean_oof_tau_matrix(
    compare_result$metadata$oof_tau_matrix
  )
  if (cleaned$status %in% c("no_input", "no_columns", "no_rows")) {
    rlang::abort(
      sprintf(paste0("cf_model_stability: oof_tau_matrix is unusable ",
                     "(status = '%s')."), cleaned$status),
      class = "cfomics_stability_empty_input"
    )
  }
  if (cleaned$status == "single_column") {
    rlang::warn(
      paste0("cf_model_stability: only one method survived cleaning ('",
             cleaned$kept_methods, "'); pairwise diagnostics are degenerate."),
      class = "cfomics_stability_single_method"
    )
  }
  if (length(cleaned$dropped_methods) > 0L) {
    rlang::inform(
      sprintf("cf_model_stability: dropped %d failed method(s): %s",
              length(cleaned$dropped_methods),
              paste(cleaned$dropped_methods, collapse = ", ")),
      class = "cfomics_stability_dropped_methods"
    )
  }

  mat <- cleaned$mat
  K   <- ncol(mat)
  kept <- cleaned$kept_methods

  # Pairwise correlation matrices.
  cor_matrices <- lapply(cor_methods, function(cm) {
    stats::cor(mat, method = cm, use = "pairwise.complete.obs")
  })
  names(cor_matrices) <- cor_methods

  # Top-q high-benefit Jaccard matrix.
  benefit_group_jaccard <- .cfomics_top_q_jaccard(mat, top_q = benefit_top_q)

  # Per-method mean off-diagonal agreement under the primary correlation.
  prim <- cor_matrices[[primary_cor]]
  prim_off <- prim
  diag(prim_off) <- NA_real_
  per_method_tau_correlation <- if (K > 1L) {
    rowMeans(prim_off, na.rm = TRUE)
  } else {
    setNames(NA_real_, kept)
  }

  # Outlier detection: mean agreement < median - outlier_mad * MAD.
  outlier_methods <- character(0)
  if (K > 2L) {
    med <- stats::median(per_method_tau_correlation, na.rm = TRUE)
    mad <- stats::mad(per_method_tau_correlation, na.rm = TRUE)
    if (is.finite(mad) && mad > 0) {
      threshold <- med - outlier_mad * mad
      flag <- per_method_tau_correlation < threshold &
              is.finite(per_method_tau_correlation)
      outlier_methods <- kept[flag]
    }
  }

  summary <- list(
    n_methods         = K,
    n_rows            = cleaned$n_rows_kept,
    n_rows_dropped    = cleaned$n_rows_dropped,
    median_agreement  = if (K > 1L) stats::median(per_method_tau_correlation,
                                                  na.rm = TRUE) else NA_real_,
    min_agreement     = if (K > 1L) min(per_method_tau_correlation,
                                        na.rm = TRUE) else NA_real_,
    mean_agreement    = if (K > 1L) mean(per_method_tau_correlation,
                                         na.rm = TRUE) else NA_real_,
    n_outliers        = length(outlier_methods)
  )

  out <- list(
    cor_matrices         = cor_matrices,
    benefit_group_jaccard       = benefit_group_jaccard,
    per_method_tau_correlation = per_method_tau_correlation,
    outlier_methods      = outlier_methods,
    summary              = summary,
    meta                 = list(
      methods_used    = kept,
      methods_dropped = cleaned$dropped_methods,
      primary_cor     = primary_cor,
      benefit_top_q = benefit_top_q,
      outlier_mad     = outlier_mad,
      cor_methods     = cor_methods,
      n_folds         = compare_result$metadata$n_folds
    )
  )
  class(out) <- c("cfomics_model_stability", "list")
  out
}


# -------------------------------------------------------------------------
# Top-q high-benefit Jaccard
# -------------------------------------------------------------------------

#' Top-q high-benefit Jaccard helper.
#'
#' For each pair of method columns, define their top-q "high-benefit" sets as
#' the indices with the highest tau hat values (ties broken by row order),
#' and return the K x K Jaccard matrix of these sets.
#'
#' @param mat Numeric n x K matrix (cleaned, all-finite).
#' @param top_q Quantile threshold (0 < top_q <= 0.5).
#' @keywords internal
#' @noRd
.cfomics_top_q_jaccard <- function(mat, top_q = 0.2) {
  K <- ncol(mat)
  n <- nrow(mat)
  if (K == 0L || n == 0L) return(matrix(NA_real_, nrow = 0L, ncol = 0L))
  top_n <- max(1L, ceiling(n * top_q))
  sets <- lapply(seq_len(K), function(j) {
    ord <- order(mat[, j], decreasing = TRUE)
    ord[seq_len(top_n)]
  })
  jac <- matrix(NA_real_, nrow = K, ncol = K,
                dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i == j) {
        jac[i, j] <- 1
      } else {
        a <- sets[[i]]; b <- sets[[j]]
        inter <- length(intersect(a, b))
        uni   <- length(union(a, b))
        jac[i, j] <- if (uni > 0L) inter / uni else NA_real_
      }
    }
  }
  jac
}


# -------------------------------------------------------------------------
# S3 print / summary
# -------------------------------------------------------------------------

#' @export
print.cfomics_model_stability <- function(x, ...) {
  cat("<cfomics_model_stability>: K = ", x$summary$n_methods,
      ", n = ", x$summary$n_rows,
      ", primary_cor = ", x$meta$primary_cor, "\n", sep = "")
  cat("  median agreement: ",
      format(x$summary$median_agreement, digits = 3),
      ", min: ", format(x$summary$min_agreement, digits = 3), "\n", sep = "")
  if (x$summary$n_outliers > 0L) {
    cat("  outlier_methods (", x$summary$n_outliers, "): ",
        paste(x$outlier_methods, collapse = ", "), "\n", sep = "")
  } else {
    cat("  outlier_methods: none\n")
  }
  if (length(x$meta$methods_dropped) > 0L) {
    cat("  methods_dropped: ", paste(x$meta$methods_dropped, collapse = ", "),
        "\n", sep = "")
  }
  cat("  cor_methods: ", paste(x$meta$cor_methods, collapse = ", "), "\n",
      sep = "")
  cat("  benefit_top_q: ", x$meta$benefit_top_q, "\n", sep = "")
  cat("  NOTE: this is a descriptive diagnostic. Do not use for model ",
      "selection.\n", sep = "")
  invisible(x)
}


#' @export
summary.cfomics_model_stability <- function(object, ...) {
  res <- list(
    methods         = object$meta$methods_used,
    methods_dropped = object$meta$methods_dropped,
    summary         = object$summary,
    per_method      = object$per_method_tau_correlation,
    outliers        = object$outlier_methods,
    primary_cor     = object$meta$primary_cor,
    benefit_top_q = object$meta$benefit_top_q
  )
  class(res) <- "summary.cfomics_model_stability"
  res
}


#' @export
print.summary.cfomics_model_stability <- function(x, ...) {
  cat("cf_model_stability summary\n")
  cat("  n_methods: ", x$summary$n_methods,
      ", n_rows: ", x$summary$n_rows, "\n", sep = "")
  cat("  primary_cor: ", x$primary_cor, "\n", sep = "")
  cat("  median / min / mean agreement: ",
      format(x$summary$median_agreement, digits = 3), " / ",
      format(x$summary$min_agreement, digits = 3), " / ",
      format(x$summary$mean_agreement, digits = 3), "\n", sep = "")
  cat("  per-method agreement:\n")
  for (nm in names(x$per_method)) {
    cat("    ", nm, ": ", format(x$per_method[nm], digits = 3), "\n",
        sep = "")
  }
  if (length(x$outliers) > 0L) {
    cat("  outliers: ", paste(x$outliers, collapse = ", "), "\n", sep = "")
  }
  cat("  NOTE: descriptive only; do not use for model selection.\n")
  invisible(x)
}
