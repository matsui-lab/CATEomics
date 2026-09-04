#' Select the best CATE method from a cf_compare scoreboard.
#'
#' @description
#' Picks the rank-1 method from a \code{\link{cf_compare}} scoreboard under
#' a pre-registered metric. The selection itself is intentionally thin: it
#' is simply \code{argmin} of the loss column (mirroring EconML's
#' \code{RScorer.best_model}). When the runner-up's between-fold CI
#' overlaps with the winner's, a \code{tie_detected} flag is set and the
#' user is steered toward \code{\link{cf_ensemble}}.
#'
#' The selection is deliberately kept simple: it is argmin of the loss
#' column, without a formal post-selection-inference correction or honest
#' selection CI. The
#' function does \strong{not} refit; it returns the cached full-data
#' model that cf_compare already stored.
#'
#' @param compare_result \code{cfomics_compare} object from
#'   \code{\link{cf_compare}}.
#' @param metric Character; one of the metric columns recorded in
#'   \code{compare_result$scoreboard}.
#' @param top_n Integer (default \code{1L}). When greater than 1, returns
#'   \code{models = list(...)} with the top-N methods. \code{top_n = 1L}
#'   returns a single \code{model = ...} slot.
#' @param tie One of \code{"report"} (default) or \code{"ignore"}. When
#'   \code{"report"}, the runner-up's between-fold CI is compared with the
#'   winner's; overlap sets \code{tie_detected = TRUE} and an informational
#'   message points to \code{cf_ensemble}.
#'
#' @return An object of class \code{c("cfomics_select", "list")} with:
#'   \itemize{
#'     \item \code{selected}: winning method name (or top-N vector);
#'     \item \code{model} (top_n = 1) or \code{models} (top_n > 1): the
#'           cached cf_compare full-data fit(s);
#'     \item \code{selection_stats}: list of
#'           \code{selected_method}, \code{metric}, \code{value},
#'           \code{naive_ci}, \code{runner_up},
#'           \code{runner_up_value},
#'           \code{runner_up_diff_ci},
#'           \code{tie_detected}, \code{tie_method},
#'           \code{K_candidates};
#'     \item \code{top_n}, \code{metric} (for echo).
#'   }
#'
#' @section Tie handling:
#' The two cfomics selection knobs are intentionally parallel:
#' \code{cf_select} returns a \emph{single method}; \code{cf_ensemble}
#' returns a \emph{combined tau hat}. When the runner-up's CI overlaps
#' with the winner's, the spec recommends users consider
#' \code{cf_ensemble} instead of forcing a single pick. The
#' \code{tie_detected} flag exposes that situation.
#'
#' @seealso \code{\link{cf_compare}}, \code{\link{cf_ensemble}},
#'   \code{\link{cf_evaluate}}.
#' @examples
#' # Two R-native methods so selection is meaningful; small MCMC keeps bcf fast.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L)
#' sel <- cf_select(cmp)
#' sel$selection_stats$selected_method
#' @export
cf_select <- function(compare_result,
                      metric = "r_loss",
                      top_n  = 1L,
                      tie    = c("report", "ignore")) {
  if (!inherits(compare_result, "cfomics_compare")) {
    rlang::abort(
      c("cf_select requires output from cf_compare().",
        i = "Single cf_fit output does not contain multiple methods for selection."),
      class = "cfomics_select_input_error"
    )
  }
  tie <- match.arg(tie)

  sb       <- compare_result$scoreboard
  metadata <- compare_result$metadata

  # The scoreboard identity column is `label` in grid mode (method_base) and
  # `method` in the legacy path; key selection / models[[.]] by whichever holds.
  ids <- if ("label" %in% colnames(sb)) sb$label else sb$method

  if (!(metric %in% colnames(sb))) {
    rlang::abort(
      sprintf("metric '%s' not found in compare_result. Available: %s.",
              metric,
              paste(setdiff(colnames(sb), c("label", "method", "base")),
                    collapse = ", ")),
      class = "cfomics_select_metric_missing"
    )
  }

  loss_vec <- sb[[metric]]
  if (all(is.na(loss_vec))) {
    rlang::abort("All methods failed; nothing to select.",
                 class = "cfomics_select_all_failed")
  }

  # honour the metric's optimisation direction. Losses / errors are
  # minimised (lower = better); score-type metrics (r_squared, RATE, AUTOC,
  # QINI, AUC, correlation) are maximised. Sorting ascending for a score metric
  # would pick the WORST method.
  higher_better <- .cf_metric_higher_better(metric)
  ord <- order(loss_vec, decreasing = higher_better, na.last = TRUE)
  ranked_methods <- ids[ord]

  ci_lower <- sb[[paste0(metric, "_ci_lower")]]
  ci_upper <- sb[[paste0(metric, "_ci_upper")]]

  winner_idx <- ord[1L]
  winner     <- ids[winner_idx]
  win_value  <- loss_vec[winner_idx]
  win_naive  <- if (!is.null(ci_lower) && !is.null(ci_upper)) {
    c(ci_lower[winner_idx], ci_upper[winner_idx])
  } else {
    c(NA_real_, NA_real_)
  }

  # Runner-up + tie detection (between-fold CI overlap).
  runner_up         <- NA_character_
  runner_up_value   <- NA_real_
  runner_up_diff_ci <- c(NA_real_, NA_real_)
  tie_detected      <- FALSE

  if (length(ord) >= 2L) {
    ru_idx <- ord[2L]
    runner_up         <- ids[ru_idx]
    runner_up_value   <- loss_vec[ru_idx]
    ru_naive <- if (!is.null(ci_lower) && !is.null(ci_upper)) {
      c(ci_lower[ru_idx], ci_upper[ru_idx])
    } else {
      c(NA_real_, NA_real_)
    }
    if (identical(tie, "report") &&
        all(is.finite(c(win_naive, ru_naive)))) {
      diff_est <- runner_up_value - win_value
      half_w <- (win_naive[2L] - win_naive[1L]) / 2
      half_r <- (ru_naive[2L]  - ru_naive[1L])  / 2
      diff_half <- sqrt(half_w^2 + half_r^2)
      runner_up_diff_ci <- c(diff_est - diff_half, diff_est + diff_half)
      tie_detected <- (runner_up_diff_ci[1L] <= 0) &&
                      (runner_up_diff_ci[2L] >= 0)
    }
  }

  if (isTRUE(tie_detected)) {
    rlang::inform(
      sprintf(paste0("cf_select: '%s' and runner-up '%s' overlap on '%s'. ",
                     "Consider cf_ensemble() to combine them instead of forcing ",
                     "a single pick."),
              winner, runner_up, metric),
      class = "cfomics_select_tie_detected"
    )
  }

  selection_stats <- list(
    selected_method    = winner,
    metric             = metric,
    value              = unname(win_value),
    naive_ci           = unname(win_naive),
    runner_up          = runner_up,
    runner_up_value    = unname(runner_up_value),
    runner_up_diff_ci  = unname(runner_up_diff_ci),
    tie_detected       = isTRUE(tie_detected),
    tie_method         = tie,
    K_candidates       = length(metadata$methods)
  )

  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1L) {
    rlang::abort("top_n must be a positive integer.",
                 class = "cfomics_select_top_n_invalid")
  }

  if (top_n == 1L) {
    out <- list(
      selected        = winner,
      model           = compare_result$models[[winner]],
      # `spec` = list(method, model_outcome, model_final) so the winner is
      # reproducible without parsing the compound label (base_learners grid).
      spec            = compare_result$models[[winner]]$meta$spec,
      selection_stats = selection_stats,
      top_n           = 1L,
      metric          = metric
    )
  } else {
    top_methods <- ranked_methods[seq_len(min(top_n, length(ranked_methods)))]
    out <- list(
      selected        = top_methods,
      models          = compare_result$models[top_methods],
      specs           = lapply(top_methods,
                               function(k) compare_result$models[[k]]$meta$spec),
      selection_stats = selection_stats,
      top_n           = top_n,
      metric          = metric
    )
  }
  class(out) <- c("cfomics_select", "list")
  out
}


#' @export
print.cfomics_select <- function(x, ...) {
  cat("<cfomics_select>: metric = ", x$metric,
      ", winner = ", paste(x$selected, collapse = ", "), "\n", sep = "")
  s <- x$selection_stats
  cat("  value = ", format(s$value, digits = 4), "\n", sep = "")
  if (all(is.finite(s$naive_ci))) {
    cat("  between-fold 95% CI: [",
        format(s$naive_ci[1L], digits = 4), ", ",
        format(s$naive_ci[2L], digits = 4), "]\n", sep = "")
  }
  if (!is.na(s$runner_up)) {
    cat("  runner_up = ", s$runner_up,
        " (value = ", format(s$runner_up_value, digits = 4), ")\n", sep = "")
    if (all(is.finite(s$runner_up_diff_ci))) {
      cat("  runner_up - winner 95% CI: [",
          format(s$runner_up_diff_ci[1L], digits = 4), ", ",
          format(s$runner_up_diff_ci[2L], digits = 4), "]\n", sep = "")
    }
  }
  cat("  tie_detected = ", s$tie_detected,
      " (tie = '", s$tie_method, "')\n", sep = "")
  if (isTRUE(s$tie_detected)) {
    cat("  hint: cf_ensemble() combines methods rather than picking one.\n")
  }
  cat("  K_candidates = ", s$K_candidates, "\n", sep = "")
  invisible(x)
}

#' Whether a selection metric is maximised (higher = better).
#'
#' Losses / errors (r_loss, dr_loss, factual_rmse, ...) are minimised; score
#' metrics (r_squared, RATE, AUTOC, QINI, AUC, correlation) are maximised.
#' @keywords internal
#' @noRd
.cf_metric_higher_better <- function(metric) {
  grepl("score|r_squared|\\br2\\b|\\bauc\\b|\\brate\\b|autoc|qini|correlation",
        metric, ignore.case = TRUE)
}
