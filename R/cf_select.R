#' Select the best CATE method from a cf_compare scoreboard.
#'
#' @description
#' Orders the requested \code{metric} column of a \code{\link{cf_compare}}
#' scoreboard itself -- \code{argmin} for a loss, \code{argmax} for a score
#' -- and returns the best config. The scoreboard's own \code{rank} column
#' follows \code{metric[1]}, so a different \code{metric} can pick a
#' different winner. When the fold-paired difference between the
#' runner-up and the winner covers zero, \code{tie_detected} records that the
#' two are not separated on this metric and this split.
#'
#' The function does \strong{not} refit; it returns the cached full-data
#' model that cf_compare already stored.
#'
#' @param compare_result \code{cateomics_compare} object from
#'   \code{\link{cf_compare}}.
#' @param metric Character; the scoreboard column to select on. The default
#'   \code{NULL} takes the first metric the comparison ranked by
#'   (\code{compare_result$metadata$metric[1]}), which is \code{"r_loss"}
#'   unless \code{cf_compare(metric = )} said otherwise -- so a comparison run
#'   with one non-default metric can be handed straight to this function. Give
#'   it explicitly to select on another of the columns the comparison recorded:
#'   \code{"r_loss"}, \code{"dr_loss"}, \code{"factual_loss"} (minimised) or
#'   \code{"r_squared_score"}, \code{"dr_squared_score"},
#'   \code{"factual_squared_score"} (maximised). A column the comparison did
#'   not record is an error.
#' @param top_n Integer (default \code{1L}). When greater than 1, returns
#'   \code{models = list(...)} with the top-N methods. \code{top_n = 1L}
#'   returns a single \code{model = ...} slot.
#' @param tie One of \code{"report"} (default) or \code{"ignore"}. When
#'   \code{"report"}, the per-fold losses of the runner-up and the winner are
#'   differenced fold by fold and summarised by a t interval on
#'   \code{K - 1} degrees of freedom; an interval covering zero sets
#'   \code{tie_detected = TRUE} and says so. The difference is on the loss
#'   scale. \code{selection_stats$tie_basis} records which comparison was
#'   used; it is \code{NA} when the per-fold losses are unavailable, and the
#'   reason in words when \code{tie = "ignore"} left the question unasked.
#'   \code{tie_detected} is \code{NA} in both of those cases.
#'
#' @return An object of class \code{c("cateomics_select", "list")} with:
#'   \itemize{
#'     \item \code{selected}: winning method name (or top-N vector);
#'     \item \code{model} (top_n = 1) or \code{models} (top_n > 1): the
#'           cached cf_compare full-data fit(s);
#'     \item \code{selection_stats}: list of
#'           \code{selected_method}, \code{metric}, \code{value},
#'           \code{between_fold_ci}, \code{runner_up},
#'           \code{runner_up_value},
#'           \code{runner_up_diff_ci},
#'           \code{tie_detected}, \code{tie_method}, \code{tie_basis},
#'           \code{K_candidates};
#'     \item \code{metric}: the column the selection was made on. The config
#'           behind a winner is that fit's \code{meta$spec}: with
#'           \code{top_n = 1} that is \code{model$meta$spec}, and with
#'           \code{top_n > 1} one per entry of \code{models}.
#'   }
#'
#' @section Tie handling:
#' \code{tie_detected} marks the case where the fold-paired difference
#' between the winner and the runner-up covers zero: the two are not separated
#' by this metric on this split. It is \code{TRUE} then, \code{FALSE} when
#' they are separated or when there is no runner-up (\code{runner_up} is then
#' \code{NA}), and \code{NA} whenever the question was not asked --- which
#' includes \code{tie = "ignore"}, whatever the candidates look like, and a
#' comparison that kept no per-fold losses. \code{tie_basis} says which of
#' those it was.
#'
#' @seealso \code{\link{cf_compare}}, \code{\link{cf_ensemble}},
#'   \code{\link{cf_evaluate}}.
#' @examples
#' # Two R-native methods so selection is meaningful; small MCMC keeps bcf fast.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)))
#' sel <- cf_select(cmp)
#' sel$selection_stats$selected_method
#' @references
#' Nie & Wager (2021) Quasi-Oracle Estimation of Heterogeneous Treatment
#' Effects, \emph{Biometrika} 108(2):299-319 (R-loss).
#'
#' Kennedy (2023) Towards optimal doubly robust estimation of heterogeneous
#' causal effects, \emph{Electron. J. Stat.} 17(2):3008-3049 (DR-loss).
#'
#' Mahajan, Mitliagkas, Neal & Syrgkanis (2024) Empirical analysis of model
#' selection for heterogeneous causal effect estimation, \emph{ICLR}
#' (proxy losses for model selection).
#'
#' @export
cf_select <- function(compare_result,
                      metric = NULL,
                      top_n  = 1L,
                      tie    = c("report", "ignore")) {
  if (!inherits(compare_result, "cateomics_compare")) {
    rlang::abort(
      c("cf_select requires output from cf_compare().",
        i = "Single cf_fit output does not contain multiple methods for selection."),
      class = "cateomics_select_input_error"
    )
  }
  tie <- match.arg(tie)

  sb       <- compare_result$scoreboard
  metadata <- compare_result$metadata

  # The scoreboard identity column is `label` in grid mode (method_base) and
  # `method` in the legacy path; key selection / models[[.]] by whichever holds.
  ids <- if ("label" %in% colnames(sb)) sb$label else sb$method

  # Follow the comparison: it ranked its candidates by the first metric it was
  # asked for, and that is the ranking `print()` and `summary()` report. With a
  # fixed default here, a comparison run with one non-default metric could not
  # be selected from at all -- the column the default names is not in its
  # scoreboard.
  if (is.null(metric)) {
    metric <- metadata$metric[[1L]] %||% "r_loss"
  }
  if (length(metric) != 1L || !is.character(metric)) {
    rlang::abort("`metric` must be a single column name.",
                 class = "cateomics_select_metric_missing")
  }

  if (!(metric %in% colnames(sb))) {
    rlang::abort(
      sprintf("metric '%s' not found in compare_result. Available: %s.",
              metric,
              paste(intersect(names(.cf_metric_directions), colnames(sb)),
                    collapse = ", ")),
      class = "cateomics_select_metric_missing"
    )
  }

  loss_vec <- sb[[metric]]
  if (all(is.na(loss_vec))) {
    rlang::abort("All methods failed; nothing to select.",
                 class = "cateomics_select_all_failed")
  }

  # Losses are minimised, squared-score columns are maximised; the direction
  # comes from the table so an unknown column stops here instead of being
  # sorted the wrong way round.
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

  # Runner-up + tie detection (fold-paired difference of the per-fold
  # losses).
  runner_up         <- NA_character_
  runner_up_value   <- NA_real_
  runner_up_diff_ci <- c(NA_real_, NA_real_)
  # FALSE when there is no runner-up to tie with, NA when a runner-up exists
  # but the paired difference could not be formed.
  tie_detected      <- FALSE
  tie_basis         <- NA_character_

  if (identical(tie, "ignore")) {
    # Not looked at, which is a different thing from looked at and not a tie:
    # `runner_up` is still filled in below, so FALSE here would read as "there
    # is a runner-up and it is separated". Same shape as `evalue_ci`.
    tie_detected <- NA
    tie_basis    <- "not assessed (tie = 'ignore')"
  }
  if (length(ord) >= 2L) {
    ru_idx <- ord[2L]
    runner_up         <- ids[ru_idx]
    runner_up_value   <- loss_vec[ru_idx]
    if (identical(tie, "report")) {
      tie_detected <- NA
      paired <- .cf_select_paired_diff(metadata, metric, winner, runner_up)
      if (!is.null(paired)) {
        runner_up_diff_ci <- paired$ci
        tie_detected      <- (paired$ci[1L] <= 0) && (paired$ci[2L] >= 0)
        tie_basis         <- "fold-paired difference"
      } else {
        rlang::inform(
          c("Tie not assessed: no per-fold losses for this metric.",
            i = "cf_compare() stores them for the metrics it was asked to compute."),
          class = "cateomics_select_tie_not_assessed")
      }
    }
  }

  if (isTRUE(tie_detected)) {
    rlang::inform(
      sprintf(paste0("cf_select: the fold-paired difference between '%s' and ",
                     "runner-up '%s' covers zero on '%s', so this metric and ",
                     "this split do not separate them."),
              winner, runner_up, metric),
      class = "cateomics_select_tie_detected"
    )
  }

  selection_stats <- list(
    selected_method    = winner,
    metric             = metric,
    value              = unname(win_value),
    between_fold_ci    = unname(win_naive),
    runner_up          = runner_up,
    runner_up_value    = unname(runner_up_value),
    runner_up_diff_ci  = unname(runner_up_diff_ci),
    tie_detected       = tie_detected,
    tie_method         = tie,
    tie_basis          = tie_basis,
    K_candidates       = length(metadata$methods)
  )

  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1L) {
    rlang::abort("top_n must be a positive integer.",
                 class = "cateomics_select_top_n_invalid")
  }

  if (top_n == 1L) {
    out <- list(
      selected        = winner,
      model           = compare_result$models[[winner]],
      selection_stats = selection_stats,
      metric          = metric
    )
  } else {
    top_methods <- ranked_methods[seq_len(min(top_n, length(ranked_methods)))]
    out <- list(
      selected        = top_methods,
      models          = compare_result$models[top_methods],
      selection_stats = selection_stats,
      metric          = metric
    )
  }
  class(out) <- c("cateomics_select", "list")
  out
}


#' @export
print.cateomics_select <- function(x, ...) {
  cat("<cateomics_select>: metric = ", x$metric,
      ", winner = ", paste(x$selected, collapse = ", "), "\n", sep = "")
  s <- x$selection_stats
  # A print says what it cannot find rather than stopping part way.
  if (!is.list(s)) {
    cat("  no selection statistics in this object\n")
    return(invisible(x))
  }
  cat("  value = ", format(s$value, digits = 4), "\n", sep = "")
  if (all(is.finite(s$between_fold_ci))) {
    cat("  between-fold 95% CI: [",
        format(s$between_fold_ci[1L], digits = 4), ", ",
        format(s$between_fold_ci[2L], digits = 4), "]\n", sep = "")
  }
  if (length(s$runner_up) == 1L && !is.na(s$runner_up)) {
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
  cat("  K_candidates = ", s$K_candidates, "\n", sep = "")
  invisible(x)
}

# Every scoreboard column cf_compare can be asked to select on, and the
# direction that makes a config better on it. A pattern over the name would
# accept a column that is not there and guess the direction of one that is.
.cf_metric_directions <- c(
  r_loss                = FALSE,
  dr_loss               = FALSE,
  factual_loss          = FALSE,
  r_squared_score       = TRUE,
  dr_squared_score      = TRUE,
  factual_squared_score = TRUE
)

#' Whether a selection metric is maximised (higher = better).
#'
#' The three loss columns are minimised, the three squared-score columns are
#' maximised. Any other name aborts.
#' @keywords internal
#' @noRd
.cf_metric_higher_better <- function(metric) {
  metric <- as.character(metric)[1L]
  if (!metric %in% names(.cf_metric_directions)) {
    rlang::abort(
      c(sprintf("Unknown metric '%s'; its optimisation direction is not defined.",
                metric),
        i = paste0("Known metrics: ",
                   paste(names(.cf_metric_directions), collapse = ", "), ".")),
      class = "cateomics_select_unknown_metric"
    )
  }
  unname(.cf_metric_directions[[metric]])
}


#' Paired per-fold difference between two configs on one metric
#'
#' The two configs were scored on the same folds with the same nuisance, so
#' their losses move together across folds. Differencing within a fold removes
#' that shared movement; adding the two intervals as if independent does not.
#' Returns \code{NULL} when the per-fold losses are unavailable.
#'
#' @param metadata \code{cf_compare} metadata.
#' @param metric Scoreboard column used for selection.
#' @param winner,runner_up Config labels.
#' @return \code{list(estimate, se, ci, k)} on the loss scale: positive means
#'   the runner-up carries the larger loss.
#' @keywords internal
#' @noRd
.cf_select_paired_diff <- function(metadata, metric, winner, runner_up) {
  loss_key <- sub("_squared_score$", "_loss", metric)
  mat  <- metadata$per_fold_loss[[loss_key]]
  if (is.null(mat) || !all(c(winner, runner_up) %in% rownames(mat))) return(NULL)

  # Always difference the loss, including for a score metric: a per-fold
  # score is a ratio, and the mean of per-fold ratios is not the pooled ratio
  # the scoreboard reports, so the interval would not be on either scale.
  vals <- mat[c(winner, runner_up), , drop = FALSE]
  d <- vals[2L, ] - vals[1L, ]
  d <- d[is.finite(d)]
  k <- length(d)
  if (k < 2L) return(NULL)

  est <- mean(d)
  se  <- stats::sd(d) / sqrt(k)
  if (!is.finite(se)) return(NULL)
  tcrit <- stats::qt(0.975, df = k - 1L)
  list(estimate = est, se = se, k = k,
       ci = c(est - tcrit * se, est + tcrit * se))
}
