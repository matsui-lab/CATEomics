#' Diagnose estimator dependence across the methods in a cf_compare result.
#'
#' \code{cf_model_stability()} asks how much the candidates of a comparison
#' order the samples the same way. It reads one thing, the K columns of
#' out-of-fold tau-hat, and reports the correlation between every pair of
#' methods -- so a high value says the methods rank the samples alike, not
#' that any of them is right.
#'
#' \itemize{
#'   \item Pairwise correlation matrices (Pearson, Spearman, Kendall) over
#'         the K cleaned methods. Spearman is the primary rank-based
#'         measure; \code{cf_ensemble(method = "cba")} defaults to Kendall.
#'   \item The mean off-diagonal of the primary matrix for each method: how
#'         much that method agrees with the rest of the cohort.
#' }
#'
#' Which samples each method would call high-benefit is the question
#' \code{\link{cf_benefit_groups}} answers, per sample and for the split that
#' was reported; this function stays at the level of methods.
#'
#' @param compare_result A \code{cateomics_compare} object.
#' @param cor_methods Character vector of correlation methods to compute;
#'   any subset of \code{c("pearson", "spearman", "kendall")}, all three by
#'   default. \code{primary_cor} is appended when it is not among them,
#'   since the summary reads it.
#' @param primary_cor The correlation method (default \code{"spearman"})
#'   whose K-vector of mean off-diagonals is reported per method.
#' @param ... Not used; naming an unknown argument is an error.
#'
#' @details
#' Every reported quantity is a correlation between two methods, so at least
#' two of the compared methods must carry usable out-of-fold tau-hat. A
#' comparison that leaves one is refused with
#' \code{cateomics_stability_single_method}, and one that leaves none with
#' \code{cateomics_stability_empty_input}; methods whose fit failed are listed
#' in \code{cf_compare()}'s \code{metadata$failed_methods}.
#'
#' @return An object of class \code{c("cateomics_model_stability", "list")}
#'   with:
#'   \itemize{
#'     \item \code{cor_matrices}: named list of K x K correlation matrices,
#'           one per requested correlation method;
#'     \item \code{per_method_tau_correlation}: numeric K-vector of mean
#'           off-diagonal correlations (primary);
#'     \item \code{meta}: \code{methods_used}, \code{methods_dropped},
#'           \code{primary_cor}, \code{n_rows}, \code{n_rows_dropped},
#'           \code{n_folds}. The number of methods is
#'           \code{length(methods_used)} and the correlation methods computed
#'           are \code{names(cor_matrices)}.
#'   }
#'
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf) so agreement is
#' # non-degenerate; known randomised propensity. Runs without Python.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   # 1 chain keeps bcf fast
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)),
#'                   propensity_score = 0.5, eval_propensity_score = 0.5,
#'                   random_state = 1)
#' cf_model_stability(cmp)
#'
#' @export
cf_model_stability <- function(compare_result,
                               cor_methods     = c("pearson", "spearman",
                                                   "kendall"),
                               primary_cor     = "spearman",
                               ...) {
  .cf_reject_dots(..., .fn = "cf_model_stability()")
  if (!inherits(compare_result, "cateomics_compare")) {
    rlang::abort(
      "cf_model_stability requires a cateomics_compare object.",
      class = "cateomics_invalid_compare_result"
    )
  }

  # A name that matches nothing is refused rather than dropped: dropping it
  # leaves that matrix uncomputed, and the complaint then surfaces from
  # cf_plot_stability, about an object that does not hold it.
  cor_methods <- .cf_match_several(cor_methods,
                                   c("pearson", "spearman", "kendall"),
                                   "cor_methods", "cf_model_stability()")
  primary_cor <- match.arg(primary_cor, c("pearson", "spearman", "kendall"))
  if (!primary_cor %in% cor_methods) cor_methods <- c(cor_methods, primary_cor)
  cleaned <- .cateomics_clean_oof_tau_matrix(
    compare_result$metadata$oof_tau_matrix
  )
  if (cleaned$status %in% c("no_input", "no_columns", "no_rows")) {
    rlang::abort(
      sprintf(paste0("cf_model_stability: oof_tau_matrix is unusable ",
                     "(status = '%s')."), cleaned$status),
      class = "cateomics_stability_empty_input"
    )
  }
  # Every quantity this function reports is an agreement between methods, so
  # one surviving column leaves nothing to estimate: a correlation matrix is 1
  # by definition and the per-method correlation is NA. Returning that object
  # invites it to be read as "agreement unknown" when the right reading is
  # "there is no second method to agree with", so the request is refused.
  if (cleaned$status == "single_column") {
    rlang::abort(
      c(sprintf(paste0("cf_model_stability: only one method ('%s') has usable ",
                       "out-of-fold tau-hat; there is no agreement to measure."),
                cleaned$kept_methods),
        x = if (length(cleaned$dropped_methods))
          sprintf("Dropped: %s.",
                  paste(cleaned$dropped_methods, collapse = ", ")) else NULL,
        i = "Compare at least two methods that both fit the data; cf_compare() records per-method failures in metadata$failed_methods."),
      class = "cateomics_stability_single_method"
    )
  }
  if (length(cleaned$dropped_methods) > 0L) {
    rlang::inform(
      sprintf("cf_model_stability: dropped %d failed method(s): %s",
              length(cleaned$dropped_methods),
              paste(cleaned$dropped_methods, collapse = ", ")),
      class = "cateomics_stability_dropped_methods"
    )
  }

  mat <- cleaned$mat
  kept <- cleaned$kept_methods

  # Pairwise correlation matrices.
  cor_matrices <- lapply(cor_methods, function(cm) {
    stats::cor(mat, method = cm, use = "pairwise.complete.obs")
  })
  names(cor_matrices) <- cor_methods

  # Per-method mean off-diagonal agreement under the primary correlation.
  prim <- cor_matrices[[primary_cor]]
  prim_off <- prim
  diag(prim_off) <- NA_real_
  # K >= 2 here: a single column is refused at the entry, so the off-diagonal
  # of each row is non-empty.
  per_method_tau_correlation <- rowMeans(prim_off, na.rm = TRUE)

  # The aggregates of that K-vector are computed where they are shown: there
  # is one place the numbers come from.
  out <- list(
    cor_matrices               = cor_matrices,
    per_method_tau_correlation = per_method_tau_correlation,
    meta                       = list(
      methods_used    = kept,
      methods_dropped = cleaned$dropped_methods,
      primary_cor     = primary_cor,
      n_rows          = cleaned$n_rows_kept,
      n_rows_dropped  = cleaned$n_rows_dropped,
      n_folds         = compare_result$metadata$n_folds
    )
  )
  class(out) <- c("cateomics_model_stability", "list")
  out
}


#' @export
print.cateomics_model_stability <- function(x, ...) {
  v <- x$per_method_tau_correlation
  cat("<cateomics_model_stability>: K = ", length(x$meta$methods_used),
      ", n = ", x$meta$n_rows,
      ", primary_cor = ", x$meta$primary_cor, "\n", sep = "")
  # The reported quantity is a correlation of the out-of-fold estimates, not
  # an agreement statistic in the sense of ICC or kappa, so it is called what
  # it is. With two methods there is one number, so it is printed once, and it
  # says which correlation it is: on its own line the number is otherwise
  # unlabelled.
  if (length(v) <= 2L) {
    cat("  pairwise tau correlation (", x$meta$primary_cor, "): ",
        format(v[1L], digits = 3), "\n", sep = "")
  } else {
    # Every candidate and its number: the mean over the other methods says
    # where each one sits, and a median alone hides which that is. The spread
    # of the pairs follows, because one bad pair does not move a row mean much.
    ord <- order(v, decreasing = TRUE)
    cat("  per-method mean tau correlation: ",
        paste(sprintf("%s %s", names(v)[ord], format(v[ord], digits = 3)),
              collapse = ", "), "\n", sep = "")
    off <- .cf_offdiag_values(x$cor_matrices[[x$meta$primary_cor]])
    cat("  pairwise range (", x$meta$primary_cor, "): ",
        format(min(off, na.rm = TRUE), digits = 3), " to ",
        format(max(off, na.rm = TRUE), digits = 3), "\n", sep = "")
  }
  if (length(x$meta$methods_dropped) > 0L) {
    cat("  methods_dropped: ", paste(x$meta$methods_dropped, collapse = ", "),
        "\n", sep = "")
  }
  cat("  cor_methods: ", paste(names(x$cor_matrices), collapse = ", "), "\n",
      sep = "")
  invisible(x)
}


#' @export
summary.cateomics_model_stability <- function(object, ...) {
  v <- object$per_method_tau_correlation
  # One row per correlation that was computed: its pairs, their range and their
  # median. Otherwise the two matrices that are not the primary one are
  # computed on every call and shown nowhere, and whether the choice of
  # correlation changes the reading cannot be seen.
  pairs <- lapply(object$cor_matrices, function(m) {
    off <- .cf_offdiag_values(m)
    c(min = min(off, na.rm = TRUE), max = max(off, na.rm = TRUE),
      median = stats::median(off, na.rm = TRUE))
  })
  res <- list(
    methods         = object$meta$methods_used,
    methods_dropped = object$meta$methods_dropped,
    n_methods       = length(object$meta$methods_used),
    n_rows          = object$meta$n_rows,
    n_rows_dropped  = object$meta$n_rows_dropped,
    per_method      = v,
    pairwise        = pairs,
    primary_cor     = object$meta$primary_cor
  )
  class(res) <- "summary.cateomics_model_stability"
  res
}


#' The distinct pairs of a symmetric correlation matrix
#'
#' @param m A square, symmetric matrix.
#' @return The upper triangle as a vector: one entry per pair of methods.
#' @keywords internal
#' @noRd
.cf_offdiag_values <- function(m) {
  if (is.null(m) || !is.matrix(m) || nrow(m) < 2L) return(numeric(0))
  m[upper.tri(m)]
}


#' @export
print.summary.cateomics_model_stability <- function(x, ...) {
  cat("cf_model_stability summary\n")
  # `[[` and a length check: an absent count would make the `if` fail with
  # "argument is of length zero" and take the whole summary down.
  dropped <- x[["n_rows_dropped"]]
  cat("  n_methods: ", x[["n_methods"]], ", n_rows: ", x[["n_rows"]],
      if (length(dropped) == 1L && is.finite(dropped) && dropped > 0L)
        sprintf(" (%d row(s) dropped)", as.integer(dropped)) else "",
      "\n", sep = "")
  cat("  primary_cor: ", x$primary_cor, "\n", sep = "")
  cat("  per-method mean tau correlation:\n")
  for (nm in names(x$per_method)) {
    cat("    ", nm, ": ", format(x$per_method[nm], digits = 3), "\n",
        sep = "")
  }
  if (length(x$pairwise) > 0L && length(x$per_method) > 1L) {
    cat("  pairwise tau correlation, by correlation:\n")
    w <- max(nchar(names(x$pairwise)))
    for (nm in names(x$pairwise)) {
      q <- x$pairwise[[nm]]
      cat(sprintf("    %-*s %s to %s (median %s)\n", w, nm,
                  format(q[["min"]], digits = 4),
                  format(q[["max"]], digits = 4),
                  format(q[["median"]], digits = 4)))
    }
  }
  invisible(x)
}
