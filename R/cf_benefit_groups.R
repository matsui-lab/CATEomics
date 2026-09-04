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
#' treatment effect (\eqn{\hat\tau}), not the best observed outcome.
#'
#' @param model \code{cf_model} (CATE-capable).
#' @param method One of \code{"median"} (default), \code{"tertile"},
#'   \code{"quartile"}, \code{"extreme_quartile"}, \code{"custom"}.
#' @param thresholds Numeric vector for \code{method = "custom"} (sorted,
#'   in increasing order).
#' @param honest Logical (default \code{TRUE}). When \code{TRUE} use OOF
#'   tau-hat; \code{FALSE} accept in-sample tau-hat without the warning.
#' @param consensus Logical (default \code{FALSE}). When \code{TRUE}, also
#'   compute per-observation cross-method consensus statistics from
#'   \code{consensus_source}. Only valid for select-route (single cf_model)
#'   inputs; passing a \code{cfomics_ensemble} object hard-errors because
#'   the methods have already been combined.
#' @param consensus_source A \code{cfomics_compare} object whose
#'   \code{metadata$oof_tau_matrix} drives the consensus statistics. Ignored
#'   when \code{consensus = FALSE}.
#' @param consensus_methods Character vector of method names to include in
#'   the consensus (must be a subset of
#'   \code{colnames(consensus_source$metadata$oof_tau_matrix)}); \code{NULL}
#'   (default) uses all surviving methods. The selected subset is recorded
#'   in the output metadata.
#' @param consensus_q Numeric in (0, 0.5]; quantile threshold for the
#'   top-q high-benefit definition used by the per-method sign agreement
#'   counts (default 0.25, matching the benefit-group Jaccard convention
#'   used by \code{cf_model_stability}).
#' @param stability_filter Logical (default \code{FALSE}). When \code{TRUE},
#'   individuals flagged as unstable (by \code{consensus$unstable_flag})
#'   are written back as \code{NA} in \code{groups$groups}. Downstream
#'   differential-expression recipes drop \code{NA} rows, so this is the
#'   recommended hook for "stable-only" grouping on the select route.
#'   Requires \code{consensus = TRUE}.
#' @param stability_min_agreement Numeric in [0, 1] (default \code{0.5}).
#'   When \code{stability_filter = TRUE}, individuals whose
#'   \code{sign_agreement} is strictly below this threshold are written
#'   back as \code{NA}.
#' @param ... Ignored.
#' @return Object of class \code{c("cfomics_groups", "list")} with
#'   \code{groups} (integer vector), \code{labels} (character),
#'   \code{thresholds}, \code{metadata}.
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
                                method     = c("median", "tertile",
                                                "quartile",
                                                "extreme_quartile",
                                                "custom"),
                                thresholds = NULL,
                                honest     = TRUE,
                                consensus              = FALSE,
                                consensus_source       = NULL,
                                consensus_methods      = NULL,
                                consensus_q            = 0.25,
                                stability_filter       = FALSE,
                                stability_min_agreement = 0.5,
                                ...) {
  if (isTRUE(stability_filter) && !isTRUE(consensus)) {
    rlang::abort(
      c("stability_filter = TRUE requires consensus = TRUE.",
        i = "Pass consensus_source = a cfomics_compare object."),
      class = "cfomics_stability_requires_consensus"
    )
  }
  stability_min_agreement <- as.numeric(stability_min_agreement)[1L]
  if (!is.finite(stability_min_agreement) ||
      stability_min_agreement < 0 || stability_min_agreement > 1) {
    rlang::abort(
      "stability_min_agreement must be in [0, 1].",
      class = "cfomics_stability_invalid_threshold"
    )
  }
  if (!inherits(model, "cf_model")) {
    rlang::abort("cf_benefit_groups requires a cf_model.",
                  class = "cfomics_invalid_input")
  }
  # Consensus is a select-route concept. An ensemble has already combined
  # the methods, so
  # the "individuals where methods disagree" diagnostic is structurally
  # ill-posed against it.
  if (isTRUE(consensus) && inherits(model, "cfomics_ensemble")) {
    rlang::abort(
      c("consensus = TRUE is not applicable to a cfomics_ensemble result.",
        i = "Pass a single-method cf_model from cf_compare()$models[[winner]] (select route).",
        i = "For ensemble-level interpretation use cf_aggregate_vi() on the member methods instead."),
      class = "cfomics_consensus_not_applicable_to_ensemble"
    )
  }
  method <- match.arg(method)

  tau <- model$res$tau_hat_oof
  if (is.null(tau)) {
    if (isTRUE(honest)) {
      rlang::warn(
        c("Model lacks OOF tau-hat; using in-sample tau-hat.",
          i = "Re-fit via cf_compare() for honest groupings."),
        class = "cfomics_partition_none_warning",
        .frequency = "once",
        .frequency_id = "cf_benefit_groups_partition_none"
      )
    }
    tau <- model$res$tau_hat
  }
  if (is.null(tau) || all(is.na(tau))) {
    rlang::abort(
      c("The model has no usable tau-hat to form benefit groups.",
        i = "Fit with cf_fit(..., partition = 'oof'), or use cf_compare()."),
      class = "cfomics_evaluate_no_tau")
  }

  if (method == "median") {
    breaks <- c(-Inf, stats::median(tau, na.rm = TRUE), Inf)
    labels <- c("low_benefit", "high_benefit")
  } else if (method == "tertile") {
    q <- stats::quantile(tau, probs = c(1/3, 2/3), na.rm = TRUE)
    breaks <- c(-Inf, q, Inf)
    labels <- c("low", "mid", "high")
  } else if (method == "quartile") {
    q <- stats::quantile(tau, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    breaks <- c(-Inf, q, Inf)
    labels <- c("Q1", "Q2", "Q3", "Q4")
  } else if (method == "extreme_quartile") {
    q <- stats::quantile(tau, probs = c(0.25, 0.75), na.rm = TRUE)
    # 2 groups: bottom Q1 vs top Q4, mid samples = NA (excluded downstream).
    groups <- rep(NA_integer_, length(tau))
    groups[tau <= q[1]] <- 1L
    groups[tau >= q[2]] <- 2L
    base_meta_eq <- list(
      method   = method,
      honest   = isTRUE(honest),
      source_method = model$method,
      formula  = model$meta$formula
    )
    consensus_stats_eq <- if (isTRUE(consensus)) {
      .cfomics_benefit_consensus(consensus_source,
                                   consensus_methods, consensus_q,
                                   n_target = length(tau))
    } else {
      NULL
    }
    stab_meta_eq <- list(stability_filter        = isTRUE(stability_filter),
                          stability_min_agreement = NA_real_,
                          n_unstable_dropped       = 0L)
    if (isTRUE(stability_filter) && !is.null(consensus_stats_eq)) {
      flag <- consensus_stats_eq$sign_agreement < stability_min_agreement &
              is.finite(consensus_stats_eq$sign_agreement)
      groups[flag] <- NA_integer_
      stab_meta_eq$stability_min_agreement <- stability_min_agreement
      stab_meta_eq$n_unstable_dropped       <- sum(flag, na.rm = TRUE)
    }
    return(structure(
      list(
        groups       = groups,
        labels       = c("Q1_low", "Q4_high"),
        thresholds   = unname(q),
        n_per_group  = c(sum(groups == 1L, na.rm = TRUE),
                          sum(groups == 2L, na.rm = TRUE)),
        method       = method,
        tau          = tau,
        consensus    = consensus_stats_eq,
        metadata     = c(base_meta_eq,
                         list(consensus       = isTRUE(consensus),
                              consensus_q     = if (isTRUE(consensus))
                                                  consensus_q else NA_real_,
                              consensus_methods = if (isTRUE(consensus))
                                consensus_stats_eq$methods_used else NULL),
                         stab_meta_eq)
      ),
      class = c("cfomics_groups", "list")
    ))
  } else if (method == "custom") {
    if (is.null(thresholds) || !is.numeric(thresholds)) {
      rlang::abort("method = 'custom' requires a numeric `thresholds` vector.",
                    class = "cfomics_invalid_input")
    }
    breaks <- c(-Inf, sort(thresholds), Inf)
    labels <- paste0("g", seq_along(breaks)[-length(breaks)])
  }

  groups <- as.integer(cut(tau, breaks = breaks, labels = FALSE,
                            include.lowest = TRUE))
  n_per <- tabulate(groups, nbins = length(labels))

  base_meta <- list(
    method   = method,
    honest   = isTRUE(honest),
    source_method = model$method,
    formula  = model$meta$formula
  )

  consensus_stats <- if (isTRUE(consensus)) {
    .cfomics_benefit_consensus(consensus_source,
                                 consensus_methods, consensus_q,
                                 n_target = length(tau))
  } else {
    NULL
  }

  stability_meta <- list(stability_filter = isTRUE(stability_filter),
                         stability_min_agreement = NA_real_,
                         n_unstable_dropped       = 0L)
  if (isTRUE(stability_filter) && !is.null(consensus_stats)) {
    flag <- consensus_stats$sign_agreement < stability_min_agreement &
            is.finite(consensus_stats$sign_agreement)
    groups[flag] <- NA_integer_
    n_per <- tabulate(groups, nbins = length(labels))
    stability_meta$stability_min_agreement <- stability_min_agreement
    stability_meta$n_unstable_dropped       <- sum(flag, na.rm = TRUE)
  }

  out_list <- list(
    groups       = groups,
    labels       = labels,
    thresholds   = breaks[is.finite(breaks)],
    n_per_group  = n_per,
    method       = method,
    tau          = tau,
    consensus    = consensus_stats,
    metadata     = c(base_meta,
                     list(consensus       = isTRUE(consensus),
                          consensus_q     = if (isTRUE(consensus))
                                              consensus_q else NA_real_,
                          consensus_methods = if (isTRUE(consensus))
                            consensus_stats$methods_used else NULL),
                     stability_meta)
  )

  structure(out_list, class = c("cfomics_groups", "list"))
}


# -------------------------------------------------------------------------
# consensus computation
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_benefit_consensus <- function(consensus_source,
                                          consensus_methods,
                                          consensus_q,
                                          n_target) {
  if (is.null(consensus_source) ||
      !inherits(consensus_source, "cfomics_compare")) {
    rlang::abort(
      "consensus = TRUE requires consensus_source = a cfomics_compare object.",
      class = "cfomics_consensus_source_missing"
    )
  }
  consensus_q <- as.numeric(consensus_q)[1L]
  if (!is.finite(consensus_q) || consensus_q <= 0 || consensus_q > 0.5) {
    rlang::abort(
      "consensus_q must be in (0, 0.5].",
      class = "cfomics_consensus_invalid_q"
    )
  }

  oof <- consensus_source$metadata$oof_tau_matrix
  if (is.null(oof) || !is.matrix(oof)) {
    rlang::abort(
      "consensus_source has no oof_tau_matrix.",
      class = "cfomics_consensus_no_oof"
    )
  }
  if (nrow(oof) != n_target) {
    rlang::abort(
      sprintf(paste0("consensus_source row count (%d) does not match model ",
                     "tau length (%d)."),
              nrow(oof), n_target),
      class = "cfomics_consensus_size_mismatch"
    )
  }

  all_methods <- colnames(oof)
  if (is.null(all_methods)) {
    all_methods <- as.character(seq_len(ncol(oof)))
  }
  if (is.null(consensus_methods)) {
    methods_used <- all_methods
  } else {
    missing <- setdiff(consensus_methods, all_methods)
    if (length(missing) > 0L) {
      rlang::abort(
        sprintf("consensus_methods not in cf_compare: %s",
                paste(missing, collapse = ", ")),
        class = "cfomics_consensus_unknown_methods"
      )
    }
    methods_used <- consensus_methods
  }

  # Restrict to surviving methods after OOF-matrix cleaning.
  cleaned <- .cfomics_clean_oof_tau_matrix(oof[, methods_used, drop = FALSE])
  if (length(cleaned$kept_methods) == 0L) {
    rlang::warn(
      paste0("No methods survived cleaning for the consensus calculation; ",
             "returning NA stats."),
      class = "cfomics_consensus_no_clean_methods"
    )
    return(list(
      consensus_fraction = rep(NA_real_, n_target),
      sign_agreement     = rep(NA_real_, n_target),
      model_sd           = rep(NA_real_, n_target),
      unstable_flag      = rep(NA, n_target),
      methods_used       = character(0),
      methods_dropped    = methods_used,
      consensus_q        = consensus_q
    ))
  }

  K_used <- length(cleaned$kept_methods)
  tau_mat <- oof[, cleaned$kept_methods, drop = FALSE]

  # Per-method high-benefit set (top-q by tau hat).
  top_n <- max(1L, ceiling(nrow(tau_mat) * consensus_q))
  resp_mat <- matrix(FALSE, nrow = nrow(tau_mat), ncol = K_used,
                     dimnames = list(NULL, cleaned$kept_methods))
  for (k in seq_len(K_used)) {
    col <- tau_mat[, k]
    ord <- order(col, decreasing = TRUE)
    resp_mat[ord[seq_len(top_n)], k] <- TRUE
  }

  consensus_fraction <- rowSums(resp_mat, na.rm = TRUE) / K_used

  # Sign agreement: fraction of methods sharing the sign of the per-row
  # median tau hat (positive vs negative effect).
  med <- apply(tau_mat, 1, stats::median, na.rm = TRUE)
  ref_sign <- sign(med)
  per_sign <- sign(tau_mat)
  match_mat <- per_sign == ref_sign
  match_mat[!is.finite(per_sign)] <- NA
  sign_agreement <- rowMeans(match_mat, na.rm = TRUE)
  # Replace NaN (all-NA rows) with NA_real_.
  sign_agreement[is.nan(sign_agreement)] <- NA_real_

  model_sd <- apply(tau_mat, 1, stats::sd, na.rm = TRUE)

  unstable_flag <- sign_agreement < 0.5 & is.finite(sign_agreement)

  list(
    consensus_fraction = consensus_fraction,
    sign_agreement     = sign_agreement,
    model_sd           = model_sd,
    unstable_flag      = unstable_flag,
    methods_used       = cleaned$kept_methods,
    methods_dropped    = setdiff(methods_used, cleaned$kept_methods),
    consensus_q        = consensus_q
  )
}

#' @export
print.cfomics_groups <- function(x, ...) {
  cat("<cfomics_groups>: method = ", x$method,
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
  if (!is.null(x$consensus)) {
    cn <- x$consensus
    cat("  consensus (top-", cn$consensus_q, "): K_used = ",
        length(cn$methods_used), "\n", sep = "")
    cat("    sign agreement summary: median = ",
        format(stats::median(cn$sign_agreement, na.rm = TRUE), digits = 3),
        ", min = ",
        format(min(cn$sign_agreement, na.rm = TRUE), digits = 3),
        "\n", sep = "")
    cat("    unstable individuals: ", sum(cn$unstable_flag, na.rm = TRUE),
        " of ", length(cn$unstable_flag), "\n", sep = "")
  }
  invisible(x)
}
