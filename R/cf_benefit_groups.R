#' Partition samples into benefit groups by predicted CATE
#'
#' @description
#' Assigns each sample to a tau-hat-based benefit group using one of
#' five splits: \code{"median"} (2 groups, default), \code{"tertile"} (3),
#' \code{"quartile"} (4), \code{"extreme_quartile"} (Q1 vs Q4, 2), or
#' \code{"custom"} (user-supplied threshold vector). Operates on the OOF
#' tau-hat from \code{\link{cf_compare}} / \code{cf_fit(partition =
#' "oof")}; falls back to in-sample tau-hat with a partition warning.
#' High-benefit groups collect the individuals with the largest predicted
#' treatment effect (\eqn{\hat\tau}), not the best observed outcome. That
#' rests on CATEomics' convention that a larger \code{Y} is the better outcome
#' (see \code{Y} in \code{\link{cf_fit}}): for an outcome whose smaller
#' values are better, fit the negated outcome, or the group labels read the
#' wrong way round.
#'
#' @section How much the other methods support this split:
#' Given a comparison, the split that was reported can be checked against the
#' other candidates: \code{consensus_source} makes the function draw the same
#' kind of split on every other method's out-of-fold tau-hat and report, for
#' each sample, the share of those methods that put it in the same group. It
#' answers "how stable is this grouping across methods", not "which grouping
#' is right".
#'
#' The share is counted per method as same group or not, so a sample one group
#' away counts like a sample at the other end. With
#' \code{method = "extreme_quartile"} the middle half is in no group: a method
#' that leaves a sample ungrouped does not support it, and a sample the
#' reported split left ungrouped gets \code{NA}. With
#' \code{method = "custom"} the thresholds are absolute, so the same cut
#' points placed on another method's tau-hat can give groups of a different
#' size; the quantile splits are redrawn on each method's own quantiles.
#'
#' The method the groups came from is not one of the voters -- it would agree
#' with itself on every sample -- and it is recognised by its tau-hat matching
#' the one the split used, not by its name (a comparison may label the same
#' method several times, once per base learner).
#'
#' @param model \code{cf_model} (CATE-capable).
#' @param method One of \code{"median"} (default), \code{"tertile"},
#'   \code{"quartile"}, \code{"extreme_quartile"}, \code{"custom"}.
#' @param thresholds Numeric vector of cut points for
#'   \code{method = "custom"}; it is sorted into increasing order before
#'   use. Passing it with any other \code{method} is an error, because the
#'   split would be drawn from quantiles and the thresholds ignored.
#' @param consensus_source A \code{cateomics_compare} object whose
#'   \code{metadata$oof_tau_matrix} the support is counted over. \code{NULL}
#'   (default) skips the calculation and leaves \code{$consensus} empty.
#' @param consensus_methods Character vector of method names to count over
#'   (a subset of \code{colnames(consensus_source$metadata$oof_tau_matrix)});
#'   \code{NULL} (default) uses every method the comparison kept.
#' @param consensus_min_agreement Numeric in [0, 1], or \code{NULL} (default)
#'   to keep every sample. When given, samples whose support falls below it
#'   are written back as \code{NA} in \code{groups}; it requires
#'   \code{consensus_source}.
#' @param ... Not used; naming an unknown argument is an error.
#' @return Object of class \code{c("cateomics_groups", "list")} with
#'   \code{groups} (integer vector), \code{labels} (character),
#'   \code{thresholds}, \code{n_per_group}, \code{tau}
#'   (the effects the split used), \code{consensus} (\code{NULL} unless
#'   \code{consensus_source} was given) and \code{metadata}
#'   (\code{method}, \code{source_method}, and \code{tau_source}, which is
#'   \code{"oof"} or \code{"in-sample"}).
#'
#'   \code{consensus} holds \code{group_support} (one share per sample),
#'   \code{methods_used} (the methods that voted), \code{K_used} (the
#'   denominator) and \code{filter} (the threshold that was applied, if any,
#'   and how many samples it left ungrouped). A method whose out-of-fold
#'   tau-hat is unusable is reported in a message rather than listed, so that
#'   not voting and not working are not the same field.
#' @examples
#' # Split patients by benefit using out-of-fold effects from a causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' grp <- cf_benefit_groups(fit, method = "median")
#' grp$n_per_group
#' @export
cf_benefit_groups <- function(model,
                              method = c("median", "tertile", "quartile",
                                         "extreme_quartile", "custom"),
                              thresholds = NULL,
                              consensus_source        = NULL,
                              consensus_methods       = NULL,
                              consensus_min_agreement = NULL,
                              ...) {
  .cf_reject_dots(..., .fn = "cf_benefit_groups()")
  if (!inherits(model, "cf_model")) {
    rlang::abort("cf_benefit_groups requires a cf_model.",
                 class = "cateomics_invalid_input")
  }
  method <- match.arg(method)
  if (!is.null(thresholds) && !identical(method, "custom")) {
    rlang::abort(
      c(sprintf("`thresholds` was given but method = '%s' cuts on quantiles.",
                method),
        i = "Use method = 'custom' to cut at the values you pass."),
      class = "cateomics_invalid_input"
    )
  }
  if (!is.null(consensus_min_agreement) && is.null(consensus_source)) {
    rlang::abort(
      c("consensus_min_agreement needs something to count agreement over.",
        i = "Pass consensus_source = a cateomics_compare object."),
      class = "cateomics_consensus_source_missing"
    )
  }
  if (!is.null(consensus_min_agreement)) {
    consensus_min_agreement <- as.numeric(consensus_min_agreement)[1L]
    if (!is.finite(consensus_min_agreement) ||
        consensus_min_agreement < 0 || consensus_min_agreement > 1) {
      rlang::abort("consensus_min_agreement must be in [0, 1].",
                   class = "cateomics_consensus_invalid_threshold")
    }
  }
  tau_source <- "oof"
  tau <- model$res$tau_hat_oof
  if (is.null(tau)) {
    # The fallback changes what the groups mean, so it is said on every call;
    # a caller who does not want to hear it can wrap the call.
    rlang::warn(
      c("Model lacks OOF tau-hat; using in-sample tau-hat.",
        i = "Re-fit via cf_compare() for honest groupings."),
      class = "cateomics_partition_none_warning"
    )
    tau <- model$res$tau_hat
    tau_source <- "in-sample"
  }
  if (is.null(tau) || all(is.na(tau))) {
    rlang::abort(
      c("The model has no usable tau-hat to form benefit groups.",
        i = "Fit with cf_fit(..., partition = 'oof'), or use cf_compare()."),
      class = "cateomics_evaluate_no_tau")
  }

  cut_out <- .cf_benefit_cut(tau, method, thresholds)
  groups  <- cut_out$groups
  labels  <- cut_out$labels

  consensus <- NULL
  if (!is.null(consensus_source)) {
    consensus <- .cf_benefit_group_support(
      consensus_source   = consensus_source,
      consensus_methods  = consensus_methods,
      reported_groups    = groups,
      reported_tau       = tau,
      method             = method,
      thresholds         = thresholds,
      reported_ids       = if (isTRUE(model$meta$provenance$sample_names_meaningful))
        rownames(model$data$X_preprocessed %||% model$data$X) else NULL)
    if (!is.null(consensus_min_agreement)) {
      drop <- is.finite(consensus$group_support) &
        consensus$group_support < consensus_min_agreement
      groups[drop] <- NA_integer_
      consensus$filter <- list(min_agreement = consensus_min_agreement,
                               n_dropped = sum(drop))
    }
  }

  structure(
    list(
      groups      = groups,
      labels      = labels,
      thresholds  = cut_out$thresholds,
      n_per_group = tabulate(groups, nbins = length(labels)),
      tau         = tau,
      consensus   = consensus,
      metadata    = list(
        method        = method,
        source_method = model$method,
        tau_source    = tau_source
      )
    ),
    class = c("cateomics_groups", "list")
  )
}


# -------------------------------------------------------------------------
# One cut, used for the reported split and for every method that votes on it
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_benefit_cut <- function(tau, method, thresholds = NULL) {
  if (identical(method, "extreme_quartile")) {
    q <- stats::quantile(tau, probs = c(0.25, 0.75), na.rm = TRUE)
    # Two groups: bottom Q1 vs top Q4; the middle half is in no group.
    groups <- rep(NA_integer_, length(tau))
    groups[tau <= q[1L]] <- 1L
    groups[tau >= q[2L]] <- 2L
    return(list(groups = groups, labels = c("Q1_low", "Q4_high"),
                thresholds = unname(q)))
  }
  if (identical(method, "median")) {
    breaks <- c(-Inf, stats::median(tau, na.rm = TRUE), Inf)
    labels <- c("low_benefit", "high_benefit")
  } else if (identical(method, "tertile")) {
    breaks <- c(-Inf, stats::quantile(tau, probs = c(1/3, 2/3), na.rm = TRUE), Inf)
    labels <- c("low", "mid", "high")
  } else if (identical(method, "quartile")) {
    breaks <- c(-Inf, stats::quantile(tau, probs = c(0.25, 0.5, 0.75),
                                      na.rm = TRUE), Inf)
    labels <- c("Q1", "Q2", "Q3", "Q4")
  } else {                                        # custom
    if (is.null(thresholds) || !is.numeric(thresholds)) {
      rlang::abort("method = 'custom' requires a numeric `thresholds` vector.",
                   class = "cateomics_invalid_input")
    }
    breaks <- c(-Inf, sort(thresholds), Inf)
    labels <- paste0("g", seq_len(length(breaks) - 1L))
  }
  list(groups = as.integer(cut(tau, breaks = breaks, labels = FALSE,
                               include.lowest = TRUE)),
       labels = labels,
       thresholds = unname(breaks[is.finite(breaks)]))
}


# -------------------------------------------------------------------------
# How many of the other methods put each sample in the same group
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_benefit_group_support <- function(consensus_source, consensus_methods,
                                      reported_groups, reported_tau,
                                      method, thresholds,
                                      reported_ids = NULL) {
  if (!inherits(consensus_source, "cateomics_compare")) {
    rlang::abort(
      "consensus_source must be a cateomics_compare object.",
      class = "cateomics_consensus_source_missing"
    )
  }
  oof <- consensus_source$metadata$oof_tau_matrix
  if (is.null(oof) || !is.matrix(oof)) {
    rlang::abort("consensus_source has no oof_tau_matrix.",
                 class = "cateomics_consensus_no_oof")
  }
  n_target <- length(reported_tau)
  if (nrow(oof) != n_target) {
    rlang::abort(
      sprintf(paste0("consensus_source row count (%d) does not match model ",
                     "tau length (%d)."), nrow(oof), n_target),
      class = "cateomics_consensus_size_mismatch"
    )
  }
  # Each voter's estimate is compared with the reported one sample by sample,
  # so the two have to be the same samples in the same order. Equal lengths
  # alone would let a permutation through in silence, which is why the
  # write-back refuses to align by position either. The comparison matrix
  # carries no row names, but its rows are in the order of its candidates'
  # designs, so the names come from one of those.
  oof_ids <- .cf_consensus_row_ids(consensus_source, nrow(oof))
  if (!is.null(oof_ids) && !is.null(reported_ids)) {
    if (!setequal(oof_ids, reported_ids)) {
      rlang::abort(
        c("consensus_source and the model were not fitted on the same samples.",
          i = "Pass the comparison the groups came from."),
        class = "cateomics_consensus_sample_mismatch"
      )
    }
    oof <- oof[match(reported_ids, oof_ids), , drop = FALSE]
  }
  all_methods <- colnames(oof)
  if (is.null(all_methods)) all_methods <- as.character(seq_len(ncol(oof)))
  if (is.null(consensus_methods)) {
    methods_asked <- all_methods
  } else {
    missing <- setdiff(consensus_methods, all_methods)
    if (length(missing) > 0L) {
      rlang::abort(
        sprintf("consensus_methods not in cf_compare: %s",
                paste(missing, collapse = ", ")),
        class = "cateomics_consensus_unknown_methods")
    }
    methods_asked <- consensus_methods
  }

  cleaned <- .cateomics_clean_oof_tau_matrix(oof[, methods_asked, drop = FALSE])
  kept <- cleaned$kept_methods
  unusable <- setdiff(methods_asked, kept)
  if (length(unusable) > 0L) {
    # Said rather than stored: the only other reason a method is not a voter
    # is that it is the method the groups came from, and one field cannot mean
    # both without the reader having to guess which.
    rlang::inform(
      sprintf("Not counting %s: no usable out-of-fold tau-hat.",
              paste(unusable, collapse = ", ")),
      class = "cateomics_consensus_methods_unusable"
    )
  }
  # The method the groups came from is recognised by its estimates, not by its
  # name: a grid comparison labels the same method once per base learner.
  same_as_reported <- vapply(kept, function(k)
    isTRUE(all.equal(unname(oof[, k]), unname(reported_tau))), logical(1))
  voters <- kept[!same_as_reported]

  if (length(voters) == 0L) {
    rlang::warn(
      c("No other method is left to count agreement over.",
        i = "The comparison holds only the method the groups came from."),
      class = "cateomics_consensus_no_clean_methods")
    return(list(group_support = rep(NA_real_, n_target),
                methods_used = character(0),
                K_used = 0L,
                filter = list(min_agreement = NA_real_, n_dropped = 0L)))
  }

  # Each voter draws the same kind of split on its own estimates; a sample it
  # leaves ungrouped (the middle of extreme_quartile) does not support the
  # reported group, so the denominator stays at the number of voters.
  agree <- vapply(voters, function(k) {
    g <- .cf_benefit_cut(oof[, k], method, thresholds)$groups
    same <- !is.na(g) & !is.na(reported_groups) & g == reported_groups
    as.numeric(same)
  }, numeric(n_target))
  if (!is.matrix(agree)) agree <- matrix(agree, nrow = n_target)
  support <- rowMeans(agree)
  support[is.na(reported_groups)] <- NA_real_

  list(group_support   = support,
       methods_used    = voters,
       K_used          = length(voters),
       filter          = list(min_agreement = NA_real_, n_dropped = 0L))
}


#' Sample names for the rows of a comparison's out-of-fold matrix
#'
#' The matrix itself carries none, but its rows are in the order of the
#' candidates' preprocessed designs, which do when the container named its
#' samples. `cf_fit` records whether those names name anything.
#'
#' @param consensus_source A `cateomics_compare` object.
#' @param n_rows The number of rows the names have to cover.
#' @return A character vector of length `n_rows`, or NULL when nothing names
#'   the samples.
#' @keywords internal
#' @noRd
.cf_consensus_row_ids <- function(consensus_source, n_rows) {
  ids <- rownames(consensus_source$metadata$oof_tau_matrix)
  if (is.null(ids)) {
    for (m in consensus_source$models) {
      if (is.null(m)) next
      nm <- rownames(m$data$X_preprocessed %||% m$data$X)
      if (isTRUE(m$meta$provenance$sample_names_meaningful) &&
          length(nm) == n_rows) { ids <- nm; break }
    }
  }
  if (is.null(ids) || length(ids) != n_rows || anyDuplicated(ids)) return(NULL)
  ids
}


#' @export
print.cateomics_groups <- function(x, ...) {
  cat("<cateomics_groups>: method = ", x$metadata$method,
      ", labels = ", paste(x$labels, collapse = "/"), "\n", sep = "")
  cat("  group sizes: ", paste(x$labels, x$n_per_group, sep = "=",
                               collapse = ", "), "\n", sep = "")
  # Report the mean predicted benefit (tau hat) per group so the
  # low/high labels are anchored to an effect scale, not just counts.
  if (!is.null(x$tau) && !is.null(x$groups)) {
    grp_mean <- vapply(seq_along(x$labels), function(g) {
      idx <- which(x$groups == g)
      if (length(idx) == 0L) NA_real_ else mean(x$tau[idx], na.rm = TRUE)
    }, numeric(1))
    cat("  mean tau-hat: ",
        paste(x$labels, format(grp_mean, digits = 3), sep = "=",
              collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$consensus) && x$consensus$K_used > 0L) {
    # With no voter left there is nothing to count, and a heading on its own
    # promises a number that never comes; the warning raised at that point
    # carries the reason.
    cn <- x$consensus
    K <- cn$K_used
    # The support takes K + 1 values, so the counts are the whole story; a
    # median or a mean of them would hide the shape. Full agreement leads,
    # because that is the question being asked.
    counts <- table(factor(round(cn$group_support * K),
                           levels = rev(seq.int(0L, K))))
    cat(sprintf("  consensus (%d other method%s): %s%s\n",
                K, if (K == 1L) "" else "s",
                paste(sprintf("%d/%d agree %d", as.integer(names(counts)), K,
                              as.integer(counts)), collapse = ", "),
                if (anyNA(cn$group_support))
                  sprintf(", ungrouped %d", sum(is.na(cn$group_support))) else ""))
    cat("    (of the other methods, how many put each sample in the same group)\n")
    f <- cn$filter
    if (!is.null(f) && is.finite(f$min_agreement)) {
      sizes <- tabulate(x$groups, nbins = length(x$labels))
      cat(sprintf("    filtered at %s: %d sample(s) left ungrouped%s\n",
                  format(f$min_agreement), f$n_dropped,
                  if (f$n_dropped > 0L)
                    sprintf(" (group sizes become %s)",
                            paste(sizes, collapse = ", ")) else ""))
    }
  }
  invisible(x)
}
