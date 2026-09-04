#' Aggregate variable importance across multiple CATE methods.
#'
#' @description
#' Cross-method aggregation of \code{cf_interpret()} surrogate-SHAP
#' outputs from \code{\link{cf_compare}}'s cached models via
#' \code{"rra"} (Robust Rank Aggregation, Kolde 2012).
#'
#' The inputs are the surrogate-SHAP outputs from \code{\link{cf_interpret}}
#' (Ridge or XGB). Robust Rank Aggregation is the single cross-method VI
#' consensus.
#'
#' Consensus is Robust Rank Aggregation (Kolde 2012) via
#' \code{RobustRankAggreg::aggregateRanks} over the per-method
#' importance-ranked feature lists -- a rank-based consensus score
#' (RRA p-value), not a calibrated effect estimate.
#'
#' When \code{cmp} was built with \code{cf_compare(base_learners = ...)},
#' aggregation is over the (method x base-learner) configs (keyed by the
#' compound \code{method_base} labels), not the bare methods; pass
#' \code{methods_subset} with those labels to restrict the set.
#'
#' @param cmp \code{cfomics_compare} object.
#' @param top_k Integer (default \code{100L}), top-k cutoff per method
#'   (sorted by SHAP importance).
#' @param methods_subset Optional character vector of config labels to
#'   include (bare method names, or \code{method_base} in grid mode).
#' @param surrogate Character; which cf_interpret() surrogate to read
#'   from. Defaults to \code{"ridge"}; \code{"xgb"} reads the TreeSHAP
#'   output instead.
#' @param ... Forwarded to \code{cf_interpret()}.
#' @return Object of class \code{c("cfomics_aggregate_vi", "list")}.
#' @seealso \code{\link{cf_interpret}}, \code{\link{cf_compare}}.
#' @examples
#' # Two R-native methods so aggregation is meaningful; small MCMC keeps bcf fast.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L)
#' cf_aggregate_vi(cmp, top_k = 20, surrogate = "xgb")
#' @export
cf_aggregate_vi <- function(cmp,
                            top_k          = 100L,
                            methods_subset = NULL,
                            surrogate      = c("ridge", "xgb"),
                            ...) {
  if (!inherits(cmp, "cfomics_compare")) {
    rlang::abort("cmp must be a cfomics_compare.",
                  class = "cfomics_invalid_input")
  }
  method    <- "rra"  # RRA is the only cross-method VI consensus.
  surrogate <- match.arg(surrogate)
  top_k     <- as.integer(top_k)[1L]

  active <- names(cmp$models)[!vapply(cmp$models, is.null, logical(1))]
  if (!is.null(methods_subset)) {
    active <- intersect(active, methods_subset)
  }
  if (length(active) == 0L) {
    rlang::abort("No active models in cmp$models for VI aggregation.",
                  class = "cfomics_aggregate_vi_no_models")
  }

  # Per-method surrogate-SHAP tables. Each entry is the features data.frame
  # (feature_name + importance + direction + rank), sorted by
  # importance so the row order is the per-method ranking.
  per_method <- list()
  failed_methods <- character(0)
  for (m in active) {
    itp <- tryCatch(
      suppressWarnings(suppressMessages(
        cf_interpret(cmp$models[[m]], surrogate = surrogate, top_n = top_k,
                      ...)
      )),
      error = function(e) NULL
    )
    if (is.null(itp) || is.null(itp$vi[[surrogate]]) ||
        is.null(itp$vi[[surrogate]]$features) ||
        nrow(itp$vi[[surrogate]]$features) == 0L) {
      failed_methods <- c(failed_methods, m)
      next
    }
    per_method[[m]] <- itp$vi[[surrogate]]$features
  }

  # Warn when only a subset of methods actually contributed to the
  # aggregate; a silent skip leads to opaque downstream behaviour after
  # reload.
  contributing_n <- length(per_method)
  total_n        <- length(active)
  if (contributing_n == 0L) {
    rlang::abort(
      "All surrogate-SHAP calls failed; nothing to aggregate.",
      class = "cfomics_aggregate_vi_all_failed"
    )
  }
  if (contributing_n < total_n) {
    rlang::warn(
      c(sprintf("cf_aggregate_vi: only %d of %d methods produced a surrogate SHAP table (surrogate = '%s').",
                contributing_n, total_n, surrogate),
        i = sprintf("Failed methods: %s.",
                    paste(failed_methods, collapse = ", "))),
      class = "cfomics_aggregate_vi_partial_methods"
    )
  }

  result <- .aggregate_rra(.ranked_lists_from_shap(per_method))

  structure(
    list(
      results  = result,
      method   = method,
      methods  = names(per_method),
      metadata = list(
        method            = method,
        top_k             = top_k,
        methods           = names(per_method),
        methods_failed    = failed_methods,
        n_methods_used    = contributing_n,
        n_methods_total   = total_n,
        surrogate         = surrogate
      )
    ),
    class = c("cfomics_aggregate_vi", "list")
  )
}


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.ranked_lists_from_shap <- function(per_method) {
  lapply(per_method, function(df) {
    # The surrogate features table is already sorted by abs SHAP descending,
    # so the rank column reflects the per-method ranking. Use feature_name
    # in the existing top-down order.
    df$feature_name
  })
}


#' @keywords internal
#' @noRd
.aggregate_rra <- function(ranked_lists) {
  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    rlang::abort(
      c("Package 'RobustRankAggreg' is required for method = 'rra'.",
        i = "Install via install.packages('RobustRankAggreg')."),
      class = "cfomics_optional_dependency_missing"
    )
  }
  res <- RobustRankAggreg::aggregateRanks(ranked_lists)
  data.frame(
    feature_name = res$Name,
    score        = res$Score,
    rank         = seq_len(nrow(res)),
    stringsAsFactors = FALSE
  )
}


## cf_aggregate_vi returns a `cfomics_aggregate_vi` object. Pathway-level
## interpretation of the ranked features is left to standard Bioconductor
## tools (e.g. fgsea) via the omics-workflow vignette, not a bundled wrapper.

#' @export
print.cfomics_aggregate_vi <- function(x, top_n = 10L, ...) {
  cat("<cfomics_aggregate_vi>: method = ", x$method,
      ", methods aggregated = ", paste(x$methods, collapse = ", "),
      ", features = ", nrow(x$results), "\n", sep = "")
  if (length(x$metadata$methods_failed) > 0L) {
    cat("  methods failed (skipped): ",
        paste(x$metadata$methods_failed, collapse = ", "), "\n", sep = "")
  }
  print(utils::head(x$results, top_n))
  invisible(x)
}
