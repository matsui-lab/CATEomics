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
#' importance-ranked feature lists; \code{score} is its RRA p-value.
#'
#' When \code{cmp} was built with \code{cf_compare(base_learners = ...)},
#' aggregation is over the (method x base-learner) configs (keyed by the
#' compound \code{method_base} labels), not the bare methods; pass
#' \code{methods_subset} with those labels to restrict the set.
#'
#' @param cmp \code{cateomics_compare} object.
#' @param top_k Integer (default \code{100L}), top-k cutoff per method
#'   (sorted by SHAP importance).
#' @param methods_subset Optional character vector of config labels to
#'   include (bare method names, or \code{method_base} in grid mode).
#' @param surrogate Character; which cf_interpret() surrogate to read
#'   from. Defaults to \code{"ridge"}; \code{"xgb"} reads the TreeSHAP
#'   output instead.
#' @param random_state Integer seed (default \code{0L}). It is the seed of
#'   each \code{\link{cf_interpret}} fit this aggregation runs, and the value
#'   used is recorded in \code{metadata}.
#' @param ... Not used; naming an unknown argument is an error.
#' @return Object of class \code{c("cateomics_aggregate_vi", "list")} with
#'   \code{results} (the aggregated table, carrying \code{score} -- the RRA
#'   p-value, so smaller means stronger agreement -- and the consensus rank)
#'   and \code{metadata} (\code{method}, the configs that contributed and
#'   the ones that failed, \code{top_k}, \code{surrogate} and
#'   \code{random_state}). A config whose interpretation warned leaves the
#'   warnings in \code{metadata$method_warnings}, and one that printed
#'   something without raising a caveat leaves it in
#'   \code{metadata$method_messages}: one named character vector per config,
#'   the names being the condition classes and the values the messages. Both
#'   are also reported once as the call returns, and both are \code{NULL}
#'   when there is nothing to report.
#' @seealso \code{\link{cf_interpret}}, \code{\link{cf_compare}}.
#' @examples
#' # Two R-native methods so aggregation is meaningful; small MCMC keeps bcf fast.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)))
#' cf_aggregate_vi(cmp, top_k = 20, surrogate = "xgb")
#' @references
#' Kolde, Laur, Adler & Vilo (2012) Robust rank aggregation for gene list
#' integration and meta-analysis, \emph{Bioinformatics} 28(4):573-580.
#'
#' @export
cf_aggregate_vi <- function(cmp,
                            top_k          = 100L,
                            methods_subset = NULL,
                            surrogate      = c("ridge", "xgb"),
                            random_state   = 0L,
                            ...) {
  if (!inherits(cmp, "cateomics_compare")) {
    rlang::abort("cmp must be a cateomics_compare.",
                  class = "cateomics_invalid_input")
  }
  surrogate <- match.arg(surrogate)
  top_k     <- as.integer(top_k)[1L]

  active <- names(cmp$models)[!vapply(cmp$models, is.null, logical(1))]
  if (!is.null(methods_subset)) {
    active <- intersect(active, methods_subset)
  }
  if (length(active) == 0L) {
    rlang::abort("No active models in cmp$models for VI aggregation.",
                  class = "cateomics_aggregate_vi_no_models")
  }

  # Per-method surrogate-SHAP tables. Each entry is the features data.frame
  # (feature_name + importance + direction + rank), sorted by
  # importance so the row order is the per-method ranking.
  per_method <- list()
  failed_methods <- character(0)
  method_warnings <- list()
  method_messages <- list()
  # A bad dot argument is a caller error, not a per-method failure: catching it
  # below would turn one typo into "every method failed". Surface it once,
  # before the loop, on the same guard cf_interpret uses.
  .cf_reject_dots(..., .fn = "cf_aggregate_vi()")
  for (m in active) {
    # Kept apart: a warning is a caveat about the surrogate, a message is
    # something a backend printed. Reporting the second as the first would say
    # a caveat was raised when none was.
    seen_w <- list(); seen_m <- list()
    itp <- withCallingHandlers(
      tryCatch(
        # The per-method surrogate is refitted here, and it moves with the
        # seed, so the seed this aggregation used is passed down and recorded
        # rather than left at cf_interpret's default.
        cf_interpret(cmp$models[[m]], surrogate = surrogate, top_n = top_k,
                      random_state = random_state, ...),
        error = function(e) NULL
      ),
      # Collected rather than dropped: a surrogate that fitted on part of the
      # rows is still used, and the caller has no other way to learn that.
      warning = function(cnd) {
        seen_w[[.cf_condition_key(cnd)]] <<- stats::setNames(
          conditionMessage(cnd), .cf_condition_class(cnd))
        invokeRestart("muffleWarning")
      },
      message = function(cnd) {
        seen_m[[.cf_condition_key(cnd)]] <<- stats::setNames(
          conditionMessage(cnd), .cf_condition_class(cnd))
        invokeRestart("muffleMessage")
      }
    )
    # `unlist` would paste the de-duplication key onto the class name.
    if (length(seen_w)) method_warnings[[m]] <- do.call(c, unname(seen_w))
    if (length(seen_m)) method_messages[[m]] <- do.call(c, unname(seen_m))
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
      class = "cateomics_aggregate_vi_all_failed"
    )
  }
  if (contributing_n < total_n) {
    rlang::warn(
      c(sprintf("cf_aggregate_vi: only %d of %d methods produced a surrogate SHAP table (surrogate = '%s').",
                contributing_n, total_n, surrogate),
        i = sprintf("Failed methods: %s.",
                    paste(failed_methods, collapse = ", "))),
      class = "cateomics_aggregate_vi_partial_methods"
    )
  }

  # Collected above rather than signalled per method, and said here: a caveat
  # the caller cannot otherwise learn about is not a caveat that has been made.
  # The same shape cf_compare() uses for its candidates.
  if (length(method_warnings) > 0L) {
    rlang::warn(
      sprintf("%d/%d method(s) raised a warning while their surrogate was fitted: %s. See metadata$method_warnings.",
              length(method_warnings), total_n,
              paste(names(method_warnings), collapse = ", ")),
      class = "cateomics_aggregate_vi_method_warnings"
    )
  }
  if (length(method_messages) > 0L) {
    rlang::inform(
      sprintf("%d/%d method(s) printed a message: %s. See metadata$method_messages.",
              length(method_messages), total_n,
              paste(names(method_messages), collapse = ", ")),
      class = "cateomics_aggregate_vi_method_messages"
    )
  }

  result <- .aggregate_rra(.ranked_lists_from_shap(per_method))

  structure(
    list(
      results  = result,
      # What was aggregated, and how, is recorded once, with the other
      # settings. The counts come from the two lists: `length(methods)` is the
      # number that contributed and adding `length(methods_failed)` gives the
      # number that were asked, so neither is stored a second time.
      metadata = list(
        method            = "rra",
        top_k             = top_k,
        methods           = names(per_method),
        methods_failed    = failed_methods,
        # What each config warned about while its surrogate was fitted; the
        # fit went ahead, so this is the only place it is said.
        method_warnings   = if (length(method_warnings) > 0L) method_warnings
                            else NULL,
        method_messages   = if (length(method_messages) > 0L) method_messages
                            else NULL,
        surrogate         = surrogate,
        random_state      = random_state
      )
    ),
    class = c("cateomics_aggregate_vi", "list")
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
      c("Package 'RobustRankAggreg' is required for cf_aggregate_vi().",
        i = "Install via install.packages('RobustRankAggreg')."),
      class = "cateomics_optional_dependency_missing"
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



#' @export
print.cateomics_aggregate_vi <- function(x, top_n = 10L, ...) {
  cat("<cateomics_aggregate_vi>: method = ", x$metadata$method,
      ", methods aggregated = ", paste(x$metadata$methods, collapse = ", "),
      ", features = ", nrow(x$results), "\n", sep = "")
  # The column is an RRA p-value, so the ordering is not the one a column
  # called "score" is usually read with.
  cat("  score = RRA p-value; smaller = stronger agreement across methods\n")
  if (length(x$metadata$methods_failed) > 0L) {
    cat("  methods failed (skipped): ",
        paste(x$metadata$methods_failed, collapse = ", "), "\n", sep = "")
  }
  print(utils::head(x$results, top_n))
  invisible(x)
}
