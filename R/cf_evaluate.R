#' Evaluate a CATE model with honest metrics
#'
#' @description
#' Compute the default three evaluation metrics for CATE-capable models:
#' \itemize{
#'   \item \strong{RATE} (Yadlowsky, Fleming, Shah, Brunskill & Wager 2025,
#'         JASA). The family is the integral of a weighted Targeting Operator
#'         Characteristic, \eqn{\int \alpha(q)\,TOC(q)\,dq}, with
#'         \eqn{TOC(q) = \bar Y_{DR}(\mathrm{top}\ q) - \bar Y_{DR}}.
#'         \code{rate$toc} is \eqn{TOC(q)} at \code{q = target_quantile},
#'         \code{rate$curve} is \eqn{TOC(q)} over the
#'         \code{rate_quantiles} grid, and \code{rate$autoc} is the grid
#'         mean of that curve (AUTOC, the \eqn{\alpha = 1} weighting).
#'   \item \strong{BLP} (Chernozhukov, Demirer, Duflo & Fernández-Val 2025,
#'         \emph{Econometrica} 93(4); Generic ML). The slope of
#'         \eqn{Y_{DR} \sim \hat\tau_{centered}}, tested against
#'         H0 = 1 (calibration) and H0 = 0 (the heterogeneity loading).
#'   \item \strong{GATES} group-ATE on tau-hat quantile bins plus a
#'         top-minus-bottom group test (non-parametric; degrades gracefully
#'         when tau-hat is near-constant).
#' }
#'
#' Operates on the OOF tau-hat saved by \code{cf_compare} /
#' \code{cf_fit(partition = 'oof')} (\code{model$res$tau_hat_oof}). When
#' the model lacks OOF tau-hat, falls back to in-sample tau-hat
#' (\code{model$res$tau_hat}) and emits \code{cateomics_partition_none_warning}.
#'
#' @details
#' The nuisance for DR pseudo-outcomes / GATES / BLP is the independent shared
#' evaluation nuisance the model was scored with (\code{model$res$nuisance}:
#' cf_compare's shared \code{nuisance_oof}, or the cf_fit-side one). The method
#' is never scored with its own internal nuisance, and the nuisance is not
#' refitted here: it is configured in \code{cf_fit()} / \code{cf_compare()},
#' so the comparison and the evaluation read the same one.
#'
#' TOC, AUTOC and the GATES best-minus-worst contrast all score the units with
#' the \emph{largest} \code{tau_hat}, on CATEomics' convention that a larger
#' \code{Y} is the better outcome (see \code{Y} in \code{\link{cf_fit}}).
#' For an outcome whose smaller values are better, fit the negated outcome
#' rather than reading these metrics the other way round: that keeps
#' \code{\link{cf_benefit_groups}}, \code{\link{cf_interpret}},
#' \code{\link{cf_rule}} and the plots pointing the same way.
#'
#' Standard errors for TOC, AUTOC, GATES and the BLP slope come from a
#' half-sample bootstrap: \code{n_bootstrap} replicates of \code{floor(n / 2)}
#' rows drawn without replacement, with the ranking, the quantile bins and the
#' quantile threshold rebuilt inside each replicate. The reported intervals are
#' \eqn{\hat\theta \pm 1.96\,SE} and the p-values are the normal-approximation
#' two-sided tests on the same SE; \code{metadata$ci_method} names the
#' standard error behind all of them. They cover the variability of these
#' quantities given
#' \code{tau_hat} and the nuisance, not the variability of having estimated
#' those.
#'
#' @param model \code{cf_model} from \code{cf_fit(partition = 'oof')}
#'   or an element of \code{cf_compare()$models}.
#' @param metrics Character vector subset of
#'   \code{c("rate", "blp", "gates")} (default all 3).
#' @param n_groups Integer (default 5), the number of quantile bins the GATES
#'   group effects are taken over. It moves \code{gates} only.
#' @param target_quantile Numeric in (0,1], the treated fraction the
#'   reported \code{rate$toc} is taken at (default 0.2).
#' @param rate_quantiles Numeric grid in (0,1], the treated fractions the
#'   TOC curve is drawn on and AUTOC averages over (default
#'   \code{seq(0.1, 1.0, by = 0.1)}).
#' @param n_bootstrap Integer, half-sample bootstrap replicates for the
#'   metric intervals (default 200).
#' @param random_state Integer seed (default \code{0}) or NULL.
#' @return Object of class \code{c("cateomics_evaluate", "list")} holding
#'   the requested metrics: \code{rate} (\code{toc} at the target quantile
#'   with its SE and interval, \code{autoc}, and \code{curve}, the TOC over
#'   the quantile grid), \code{blp} (the slope with both null hypotheses, the
#'   intercept and the per-sample \code{points} the regression was fitted
#'   on), \code{gates} (the group effects and the top-minus-bottom test, plus
#'   \code{degenerate} and, when that is \code{TRUE}, a \code{message} saying
#'   why the groups could not be formed) and
#'   \code{metadata} (the settings, and \code{ci_method}, the one standard
#'   error construction behind every interval and p-value in the object).
#'
#'   \code{blp$points} has one row per sample, \code{predicted} (the
#'   tau-hat that was evaluated) and \code{observed} (its DR pseudo-outcome);
#'   \code{blp$intercept} is the fitted value at the mean tau-hat, which with
#'   the centred regressor is the DR ATE.
#'
#'   \code{blp$heterogeneity_loading} tests the slope against zero. It is
#'   named after the coefficient rather than after the property because the
#'   slope is zero in two cases: when there is no heterogeneity, and when
#'   \code{tau_hat} is pure noise (Chernozhukov et al. call it the
#'   heterogeneity loading parameter).
#'
#'   Four of the reported quantities carry two p-values from the same estimate
#'   and the same standard error: \code{rate$p_value},
#'   \code{rate$autoc$p_value}, \code{gates$top_bottom_test$p_value} and
#'   \code{blp$heterogeneity_loading$p_value} are two-sided, and each has a
#'   \code{p_value_one_sided} beside it for the alternative that the quantity
#'   is greater than zero. \code{blp$calibration} (H0: slope = 1) has the
#'   two-sided p-value only, because a slope above 1 and a slope below 1 are
#'   both a failure of calibration. The reported intervals
#'   (\code{ci_lower} / \code{ci_upper}) are two-sided throughout, so a
#'   one-sided p-value between 0.025 and 0.05 sits beside an interval that
#'   contains zero. Both are reported; which one answers the question asked is
#'   the reader's to decide.
#' @seealso \code{\link{cf_compare}}, \code{\link{cf_select}},
#'   \code{\link{cf_rule}}.
#' @examples
#' # Evaluate a fast R-native causal forest with an independent nuisance.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' # the fit carries the evaluation nuisance it was scored with
#' ev <- cf_evaluate(fit)
#' ev
#' @references
#' Yadlowsky, Fleming, Shah, Brunskill & Wager (2025) Evaluating treatment
#' prioritization rules via rank-weighted average treatment effects,
#' \emph{JASA} 120(549):38-51 (RATE / TOC / AUTOC).
#'
#' Chernozhukov, Demirer, Duflo & Fernández-Val (2025) Fisher-Schultz Lecture:
#' Generic machine learning inference on heterogeneous treatment effects in
#' randomized experiments, \emph{Econometrica} 93(4):1121-1164 (GATES / BLP).
#'
#' @export
cf_evaluate <- function(model,
                        metrics         = c("rate", "blp", "gates"),
                        n_groups        = 5L,
                        target_quantile = 0.2,
                        rate_quantiles  = seq(0.1, 1.0, by = 0.1),
                        n_bootstrap     = 200L,
                        random_state    = 0L) {
  if (!inherits(model, "cf_model")) {
    rlang::abort(
      "cf_evaluate requires a cf_model from cf_fit() or cf_compare()$models.",
      class = "cateomics_invalid_input"
    )
  }
  metrics <- .cf_match_several(metrics, c("rate", "blp", "gates"),
                               "metrics", "cf_evaluate()")
  # Both quantile arguments carry a 1/q factor and index a top-q slice, so a
  # value outside (0, 1] has no reading: 0 divides by zero and anything above
  # 1 asks for more units than there are. Unchecked, the two failures look
  # nothing alike -- a silent NaN in `rate$toc` and a raw
  # "missing value where TRUE/FALSE needed" from the weighting helper.
  .cf_check_quantile_arg(target_quantile, "target_quantile")
  .cf_check_quantile_arg(rate_quantiles, "rate_quantiles")
  # Which quantities could not be resampled is a fact about this call, so the
  # record starts empty and is reported once on the way out -- naming them,
  # which a warning raised from inside the resampler could not do.
  .cf_se_unavailable$what <- character(0)
  on.exit(.cf_report_se_unavailable(), add = TRUE)

  # tau-hat source: OOF if available (res$tau_hat_oof), else in-sample.
  tau_oof <- model$res$tau_hat_oof
  if (is.null(tau_oof)) {
    rlang::warn(
      c("Model was fit without partition='oof'; using in-sample tau-hat.",
        i = "Evaluation metrics may be optimistically biased.",
        i = "Re-fit via cf_compare() or cf_fit(partition='oof') for honest evaluation."),
      class = "cateomics_partition_none_warning"
    )
    tau_oof <- model$res$tau_hat
    if (is.null(tau_oof) || all(is.na(tau_oof))) {
      rlang::abort(
        "Model lacks usable tau-hat (tau_hat_oof and tau_hat both NA).",
        class = "cateomics_evaluate_no_tau"
      )
    }
  }
  # Source data (Y, T): the model-visible training data. The design matrix is
  # not needed here -- every metric reads the nuisance the model carries.
  if (!is.null(model$data$Y) && !is.null(model$data$T)) {
    Y <- as.numeric(model$data$Y)
    T <- as.integer(model$data$T)
  } else {
    rlang::abort(
      c("cf_evaluate cannot recover (Y, T) from this model.",
        i = "This object carries no training data; refit with cf_fit() or cf_compare()."),
      class = "cateomics_evaluate_missing_data"
    )
  }

  # Nuisance for the DR pseudo-outcome / GATES / BLP regression. This is the
  # independent shared evaluation nuisance the model was scored with
  # (cf_compare's `nuisance_oof`, or the cf_fit-side one), never the method's
  # own internal nuisance -- self-scoring would be optimistic. Reusing exactly
  # that nuisance is what keeps cf_evaluate consistent with the comparison
  # scoreboard, so the learner, its library and a known propensity are set
  # where it is fitted rather than here.
  nuis <- model$res$nuisance
  if (is.null(nuis)) {
    rlang::abort(
      c("This model carries no evaluation nuisance.",
        i = paste("cf_fit() and cf_compare() always fit one; refit the model",
                  "to evaluate it."),
        i = paste("Set `eval_nuisance`, `eval_library`, `eval_params` or",
                  "`eval_propensity_score` there, so the comparison and the",
                  "evaluation read the same nuisance.")),
      class = "cateomics_evaluate_missing_nuisance"
    )
  }
  Y_DR <- .compute_dr_pseudo_outcome(Y, T, nuis$mu0_oof, nuis$mu1_oof,
                                       nuis$pi_oof)

  out <- list()

  if ("rate" %in% metrics) {
    out$rate <- .compute_rate_metric(tau_oof, Y_DR, target_quantile,
                                     n_bootstrap, random_state)
    out$rate$autoc <- .compute_autoc(tau_oof, Y_DR, rate_quantiles,
                                     n_bootstrap, random_state)
    # The curve belongs to the quantity it is the curve of, the way `points`
    # sits inside `blp`.
    out$rate$curve <- .compute_rate_curve(tau_oof, Y_DR, rate_quantiles,
                                          n_bootstrap, random_state)
  }
  if ("blp" %in% metrics) {
    # Calibration slope and BLP-HTE are the SAME regression beta,
    # tested against H0 = 1 (calibration) and H0 = 0 (heterogeneity). Report
    # one BLP metric carrying beta + both tests + the binned calibration data.
    out$blp <- .compute_blp_metric(tau_oof, Y_DR, n_bootstrap, random_state)
  }
  if ("gates" %in% metrics) {
    out$gates <- .compute_gates_metric(tau_oof, Y_DR, n_groups,
                                       n_bootstrap, random_state)
  }

  out$metadata <- list(
    method          = model$method,
    partition       = model$meta$estimator$partition %||% "none",
    n_groups        = as.integer(n_groups),
    target_quantile = target_quantile,
    rate_quantiles  = rate_quantiles,
    # The nuisance always comes from the model, so record where, not which
    # learner: the learner is recorded by the fit that produced it.
    eval_nuisance   = "carried by the model",
    n_bootstrap     = as.integer(n_bootstrap),
    # One construction behind every interval and p-value in the object.
    ci_method       = .cf_ci_method_label(c(
      out$rate$se, out$rate$autoc$se, out$blp$slope$se,
      out$gates$top_bottom_test$se))
  )
  class(out) <- c("cateomics_evaluate", "list")
  out
}

# -------------------------------------------------------------------------
# Metric implementations
# -------------------------------------------------------------------------

#' Half-sample bootstrap for the evaluation metrics
#'
#' RATE, AUTOC, GATES and the BLP slope all build something from the data
#' before averaging -- a ranking, quantile bins, a quantile threshold. Every
#' replicate redoes that construction on its own rows.
#'
#' Replicates draw \code{floor(n / 2)} rows without replacement. The
#' finite-population factor of a half-sample cancels the \code{sqrt(2)}
#' inflation from the smaller sample, so the plain standard deviation across
#' replicates estimates the standard error of the full-sample statistic.
#'
#' \code{stat_fn(idx)} recomputes the statistic on the rows \code{idx} and
#' returns a named numeric vector, the same names every call; a replicate that
#' cannot be computed returns \code{NA} and is dropped.
#'
#' @param n Number of observations.
#' @param stat_fn Statistic of a row subset, evaluated per replicate.
#' @param n_bootstrap Number of replicates.
#' @param random_state Integer seed.
#' @return Named numeric vector of standard errors, \code{NA} when
#'   \code{floor(n / 2)} is under 10 or fewer than 20 replicates were
#'   usable.
#' @keywords internal
#' @noRd
.cf_halfsample_se <- function(n, stat_fn, n_bootstrap = 200L,
                              random_state = NULL,
                                what = "this statistic") {
  m <- n %/% 2L
  template <- stat_fn(seq_len(n))
  if (m < 10L) {
    .cf_se_unavailable$what <- unique(c(.cf_se_unavailable$what, what))
    return(stats::setNames(rep(NA_real_, length(template)),
                           names(template)))
  }

  old <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  # Pin the algorithm, not just the seed: BiocParallel switches the stream to
  # L'Ecuyer-CMRG, which would otherwise give different replicates.
  if (!is.null(random_state))
    set.seed(as.integer(random_state), kind = "Mersenne-Twister")

  reps <- vapply(seq_len(as.integer(n_bootstrap)), function(b) {
    idx <- sample.int(n, m, replace = FALSE)
    out <- tryCatch(stat_fn(idx), error = function(e) NULL)
    if (is.null(out) || length(out) != length(template)) {
      return(stats::setNames(rep(NA_real_, length(template)), names(template)))
    }
    out
  }, numeric(length(template)))
  reps <- matrix(reps, nrow = length(template),
                 dimnames = list(names(template), NULL))

  se <- apply(reps, 1L, function(r) {
    r <- r[is.finite(r)]
    if (length(r) < 20L) NA_real_ else stats::sd(r)
  })
  # Recorded rather than signalled here: this helper runs once per metric, so
  # warning from inside it says the same thing several times for one call and
  # none of the lines says which quantity it is about. `cf_evaluate()` reports
  # the collected names once.
  if (anyNA(se)) .cf_se_unavailable$what <- unique(c(.cf_se_unavailable$what, what))
  se
}

# Which quantities had no usable bootstrap replicates in the call now running.
# Reset by `cf_evaluate()` on entry, so nothing carries over between calls.
.cf_se_unavailable <- new.env(parent = emptyenv())
.cf_se_unavailable$what <- character(0)

#' Report, once for the call, which quantities could not be resampled.
#'
#' @keywords internal
#' @noRd
.cf_report_se_unavailable <- function() {
  what <- .cf_se_unavailable$what
  if (!length(what)) return(invisible(NULL))
  rlang::warn(
    c(sprintf("No bootstrap interval for %d quantit%s: %s.",
              length(what), if (length(what) == 1L) "y" else "ies",
              paste(what, collapse = ", ")),
      i = "Fewer than 20 usable replicates: too few rows to half-sample, or the statistic is undefined on a subsample (a constant tau-hat, an empty bin)."),
    class = "cateomics_bootstrap_se_unavailable")
  invisible(NULL)
}


#' Check a treated-fraction argument
#'
#' @param q The value supplied.
#' @param nm Its argument name, for the message.
#' @return `q`, invisibly.
#' @keywords internal
#' @noRd
.cf_check_quantile_arg <- function(q, nm) {
  ok <- is.numeric(q) && length(q) >= 1L && all(is.finite(q)) &&
    all(q > 0) && all(q <= 1)
  if (!ok) {
    rlang::abort(
      c(sprintf("`%s` must be numeric values in (0, 1].", nm),
        x = sprintf("Supplied: %s.",
                    paste(utils::head(format(q), 5L), collapse = ", ")),
        i = "A treated fraction of 0 has no top-q slice, and one above 1 asks for more units than there are."),
      class = "cateomics_invalid_quantile"
    )
  }
  invisible(q)
}


#' Weight for the targeted top q fraction
#'
#' Units are ranked by `tau`, largest first, and the weight of the tied group
#' at the cut-point is split so exactly `q * n` units are selected. A constant
#' tau-hat therefore gives every unit weight `q`, and the contrast is 0.
#'
#' @param tau Numeric vector.
#' @param q Fraction in (0, 1].
#' @return Numeric weight vector summing to `q * length(tau)`.
#' @keywords internal
#' @noRd
.cf_top_q_weight <- function(tau, q) {
  n <- length(tau)
  k <- q * n
  # The cut-point is the k-th largest tau-hat; units strictly above it take
  # weight 1 and the tied block at the cut-point shares what is left.
  cut <- sort(tau, decreasing = TRUE)[max(1L, ceiling(k))]
  above <- tau > cut
  tied  <- tau == cut
  w <- as.numeric(above)
  n_tied <- sum(tied)
  if (n_tied > 0L) {
    w[tied] <- max(0, k - sum(above)) / n_tied
  }
  w
}

#' One label for every interval in an evaluation.
#'
#' The quantities share one construction, so the object records it once. The
#' label is still read off the standard errors that were produced, so a run
#' whose replicates were unusable does not announce a bootstrap.
#' @keywords internal
#' @noRd
.cf_ci_method_label <- function(ses) {
  if (!length(ses)) return(NA_character_)
  ok <- is.finite(ses)
  if (all(ok)) "half-sample bootstrap SE"
  else if (!any(ok)) "not estimable (no usable half-sample replicates)"
  else "half-sample bootstrap SE (not estimable for some quantities)"
}

#' Normal-approximation p-values for an estimate against a null.
#'
#' One construction for every tested quantity in this file, so the two-sided
#' and the one-sided p-value always read the same estimate and the same
#' standard error. The one-sided alternative is `estimate > h0`: the side the
#' quantity is reported for, since a larger value means the ranking puts more
#' of the treatment effect at the top.
#'
#' Both are `NA` when the standard error is, which is what
#' `.cf_halfsample_se()` returns for a run with no usable replicates.
#'
#' @param estimate The reported point estimate.
#' @param se Its standard error.
#' @param h0 The null value (default 0).
#' @return List with `z`, `two_sided` and `one_sided`.
#' @keywords internal
#' @noRd
.cf_normal_p <- function(estimate, se, h0 = 0) {
  z <- (estimate - h0) / se
  list(z         = unname(z),
       two_sided = unname(2 * stats::pnorm(-abs(z))),
       one_sided = unname(stats::pnorm(-z)))
}

#' Targeting Operator Characteristic at one or more treated fractions.
#'
#' TOC(q) = mean(Y_DR | the top q by predicted effect) - mean(Y_DR), computed
#' through the rank-based weight of `.cf_top_q_weight()`, which sums to
#' `q * n`, hence the 1/q.
#'
#' This is the quantity the curve shows and the one AUTOC integrates, so the
#' point estimate, the curve and the area are on one scale. The RATE family is
#' the integral of a weighted TOC (Yadlowsky et al.): weight 1 gives AUTOC,
#' weight q gives Qini. A value at a single q is a TOC, not a RATE, which is
#' why nothing here is named `rate`.
#'
#' @keywords internal
#' @noRd
.toc_point <- function(tau, Y_DR, quantiles) {
  vapply(quantiles, function(q) {
    w <- .cf_top_q_weight(tau, q)
    (1 / q) * mean(w * Y_DR) - mean(Y_DR)
  }, numeric(1))
}

#' @keywords internal
#' @noRd
.compute_rate_metric <- function(tau, Y_DR, target_quantile,
                                 n_bootstrap = 200L, random_state = NULL) {
  q       <- target_quantile
  toc_est <- unname(.toc_point(tau, Y_DR, q))
  se <- unname(.cf_halfsample_se(
    length(tau),
    stat_fn = function(i) stats::setNames(
      .toc_point(tau[i], Y_DR[i], q), "toc"),
    n_bootstrap = n_bootstrap, random_state = random_state,
    what = "the TOC at the target quantile"))
  pv <- .cf_normal_p(toc_est, se)
  list(
    toc             = toc_est,
    se              = se,
    ci_lower        = toc_est - 1.96 * se,
    ci_upper        = toc_est + 1.96 * se,
    p_value           = pv$two_sided,
    p_value_one_sided = pv$one_sided
  )
}

#' The TOC curve over the grid, whose area is AUTOC.
#' @keywords internal
#' @noRd
.compute_rate_curve <- function(tau, Y_DR, rate_quantiles,
                                n_bootstrap = 200L, random_state = NULL) {
  est <- .toc_point(tau, Y_DR, rate_quantiles)
  se  <- .cf_halfsample_se(
    length(tau),
    stat_fn = function(i) stats::setNames(
      .toc_point(tau[i], Y_DR[i], rate_quantiles),
      paste0("q", seq_along(rate_quantiles))),
    n_bootstrap = n_bootstrap, random_state = random_state,
    what = "the TOC curve")
  data.frame(quantile = rate_quantiles, toc = unname(est),
             ci_lower = unname(est - 1.96 * se),
             ci_upper = unname(est + 1.96 * se))
}

#' Area under the TOC curve (AUTOC).
#'
#' AUTOC = mean over the q-grid of TOC(q), where
#' TOC(q) = (1/q) mean(w_q(tau) * Y_DR) - mean(Y_DR) and `w_q` is the
#' rank-based weight of `.cf_top_q_weight()` (it sums to `q * n`). The grid
#' (q >= smallest rate_quantile) bounds the 1/q weight.
#' @keywords internal
#' @noRd
.autoc_point <- function(tau, Y_DR, rate_quantiles) {
  mean(.toc_point(tau, Y_DR, rate_quantiles))
}

#' @keywords internal
#' @noRd
.compute_autoc <- function(tau, Y_DR, rate_quantiles,
                           n_bootstrap = 200L, random_state = NULL) {
  autoc <- .autoc_point(tau, Y_DR, rate_quantiles)
  se <- unname(.cf_halfsample_se(
    length(tau),
    stat_fn = function(i) stats::setNames(
      .autoc_point(tau[i], Y_DR[i], rate_quantiles), "autoc"),
    n_bootstrap = n_bootstrap, random_state = random_state,
    what = "the AUTOC"))
  pv <- .cf_normal_p(autoc, se)
  list(
    value     = autoc,
    se        = se,
    ci_lower  = autoc - 1.96 * se,
    ci_upper  = autoc + 1.96 * se,
    p_value           = pv$two_sided,
    p_value_one_sided = pv$one_sided
  )
}

#' Safe quantile bins on tau.
#'
#' Near-constant tau makes \code{quantile()} breaks non-unique, which crashes
#' \code{cut()} ("breaks are not unique"). This deduplicates the breaks and
#' reports the (possibly reduced) number of usable groups; when fewer than two
#' distinct bins survive it returns \code{degenerate = TRUE} so callers can
#' skip heterogeneity binning gracefully instead of erroring.
#' @keywords internal
#' @noRd
.cf_safe_quantile_bins <- function(tau, n_groups) {
  br <- unique(stats::quantile(tau, probs = seq(0, 1, length.out = n_groups + 1L),
                               na.rm = TRUE))
  if (length(br) < 3L) {
    return(list(bin = NULL, n_groups = 0L, degenerate = TRUE))
  }
  bin <- cut(tau, breaks = br, include.lowest = TRUE, labels = FALSE)
  list(bin = bin, n_groups = length(br) - 1L, degenerate = FALSE)
}

#' Best Linear Predictor of the treatment effect.
#'
#' The DR-based BLP regresses the DR pseudo-outcome on the centred model
#' tau-hat, \eqn{Y_{DR} \sim \hat\tau_c}, giving a single slope beta tested
#' against two nulls, H0 = 1 and H0 = 0.
#' Chernozhukov et al. 2025 (full citation in \code{?cf_evaluate}). Also returns the binned observed (DR pseudo-outcome) vs
#' predicted tau data for a calibration plot.
#' @keywords internal
#' @noRd
.blp_point <- function(tau, Y_DR) {
  tau_c <- tau - mean(tau)
  stats::setNames(unname(stats::coef(stats::lm(Y_DR ~ tau_c))[2]), "slope")
}

#' @keywords internal
#' @noRd
.compute_blp_metric <- function(tau, Y_DR,
                                n_bootstrap = 200L, random_state = NULL) {
  tau_centered <- tau - mean(tau)
  # Regress the DR pseudo-outcome (E[Y_DR | X] = tau(X)) on the centred model
  # tau-hat; the slope is tested against 1 and against 0. Both coefficients
  # come from the one fit, so the line the calibration view draws is the
  # regression that is reported rather than one rebuilt from the slope alone.
  cf   <- unname(stats::coef(stats::lm(Y_DR ~ tau_centered)))
  beta <- cf[2L]
  boot <- .cf_halfsample_se(
    length(tau),
    stat_fn = function(i) .blp_point(tau[i], Y_DR[i]),
    n_bootstrap = n_bootstrap, random_state = random_state,
    what = "the BLP slope")
  se <- unname(boot[["slope"]])

  # Both tests read the same slope and the same standard error, so the
  # p-values and the interval stay on one construction. The calibration null
  # has no one-sided form that answers the question: a slope above 1 and a
  # slope below 1 are both a failure of calibration, so only the two-sided
  # p-value is reported for it.
  p_cal <- .cf_normal_p(beta, se, h0 = 1)
  p_het <- .cf_normal_p(beta, se, h0 = 0)

  list(
    slope = list(
      estimate  = beta,
      se        = se,
      ci_lower  = beta - 1.96 * se,
      ci_upper  = beta + 1.96 * se
    ),
    # Fitted value at the mean of tau-hat. The regressor is centred, so this
    # is the regression intercept and equals the DR ATE.
    intercept = cf[1L],
    # The pairs the regression was fitted on, one row per sample, so the
    # calibration view can draw them. They are the values already in the
    # object (the out-of-fold tau-hat and the DR pseudo-outcome), carried
    # here because the plot helpers are given the evaluation and nothing else.
    points = data.frame(predicted = as.numeric(tau),
                        observed   = as.numeric(Y_DR)),
    calibration   = list(h0 = 1, z = p_cal$z, p_value = p_cal$two_sided),
    # beta = 0 arises either when there is no heterogeneity or when tau-hat is
    # pure noise, so the field is named after the coefficient the source calls
    # the heterogeneity loading parameter rather than after the property.
    heterogeneity_loading = list(h0 = 0, z = p_het$z, p_value = p_het$two_sided,
                                 p_value_one_sided = p_het$one_sided)
  )
}

#' GATES quantile group ATEs + Best-Worst test.
#' @keywords internal
#' @noRd
.gates_point <- function(tau, Y_DR, n_groups) {
  bins <- .cf_safe_quantile_bins(tau, n_groups)
  est  <- rep(NA_real_, n_groups)
  if (!bins$degenerate) {
    for (k in seq_len(min(bins$n_groups, n_groups))) {
      idx <- which(bins$bin == k)
      if (length(idx)) est[k] <- mean(Y_DR[idx])
    }
  }
  # Best minus worst: the top bin minus the bottom one, on the convention that
  # a larger tau-hat is the better effect (see `Y` in `cf_fit`).
  top <- min(bins$n_groups, n_groups)
  diff <- if (bins$degenerate) NA_real_ else est[top] - est[1L]
  stats::setNames(c(est, diff), c(paste0("g", seq_len(n_groups)), "diff"))
}

#' @keywords internal
#' @noRd
.compute_gates_metric <- function(tau, Y_DR, n_groups,
                                  n_bootstrap = 200L, random_state = NULL) {
  # near-constant tau -> non-unique quantile breaks. Degrade
  # gracefully instead of crashing cut().
  bins <- .cf_safe_quantile_bins(tau, n_groups)
  if (bins$degenerate) {
    return(list(
      group_ate = NULL,
      top_bottom_test = NULL,
      degenerate = TRUE,
      message = "tau-hat is near-constant; GATES heterogeneity is not estimable."
    ))
  }
  bin <- bins$bin
  eff_groups <- bins$n_groups
  boot <- .cf_halfsample_se(
    length(tau),
    stat_fn = function(i) .gates_point(tau[i], Y_DR[i], n_groups),
    n_bootstrap = n_bootstrap, random_state = random_state,
    what = "the GATES group effects")

  rows <- do.call(rbind, lapply(seq_len(eff_groups), function(k) {
    idx <- which(bin == k)
    if (length(idx) == 0L) {
      return(data.frame(group = k, estimate = NA_real_,
                         ci_lower = NA_real_, ci_upper = NA_real_))
    }
    est <- mean(Y_DR[idx])
    se  <- unname(boot[[paste0("g", k)]])
    data.frame(group = k, estimate = est,
                ci_lower = est - 1.96 * se,
                ci_upper = est + 1.96 * se)
  }))

  diff_est <- mean(Y_DR[bin == eff_groups]) - mean(Y_DR[bin == 1L])
  diff_se  <- unname(boot[["diff"]])
  diff_p   <- .cf_normal_p(diff_est, diff_se)

  list(
    group_ate = rows,
    top_bottom_test = list(
      diff      = diff_est,
      se        = diff_se,
      ci_lower  = diff_est - 1.96 * diff_se,
      ci_upper  = diff_est + 1.96 * diff_se,
      p_value           = diff_p$two_sided,
      p_value_one_sided = diff_p$one_sided
    )
  )
}

#' @export
print.cateomics_evaluate <- function(x, ...) {
  # sprintf() with a NULL argument returns character(0), so cat() of it prints
  # nothing and the whole line disappears. An object saved before the one-sided
  # p-values existed would therefore lose four lines silently; NA says the field
  # is absent instead.
  # `sprintf` fed a NULL returns character(0), and `cat()` of that prints
  # nothing: the whole line disappears without a word. Every numeric argument
  # below goes through this, not only the p-values, because an object saved by
  # an older version can be missing any of them.
  num <- function(v) if (is.null(v) || !length(v)) NA_real_ else v
  int <- function(v) if (is.null(v) || !length(v) || !is.finite(v)) NA_integer_ else as.integer(v)
  cat("<cateomics_evaluate>: method = ", x$metadata$method,
      ", partition = ", x$metadata$partition, "\n", sep = "")
  if (!is.null(x$rate)) {
    cat(sprintf("  TOC @ q=%.2f : %.4g +/- %.3g [%.4g, %.4g]  p = %.3g (one-sided %.3g)\n",
                num(x$metadata$target_quantile), num(x$rate$toc),
                num(x$rate$se), num(x$rate$ci_lower), num(x$rate$ci_upper),
                num(x$rate$p_value), num(x$rate$p_value_one_sided)))
    if (!is.null(x$rate$autoc)) {
      a <- x$rate$autoc
      cat(sprintf("  AUTOC       : %.4g +/- %.3g [%.4g, %.4g]  p = %.3g (one-sided %.3g)\n",
                  num(a$value), num(a$se), num(a$ci_lower), num(a$ci_upper),
                  num(a$p_value), num(a$p_value_one_sided)))
    }
  }
  if (!is.null(x$blp)) {
    s <- x$blp$slope
    cat(sprintf("  BLP slope: %.4g +/- %.3g [%.4g, %.4g]\n",
                num(s$estimate), num(s$se), num(s$ci_lower),
                num(s$ci_upper)))
    # Same style as cf_rule's contrasts: the colon stays on the label.
    cat(sprintf("    %-37s p = %.3g\n", "calibration (H0 slope = 1):",
                num(x$blp$calibration$p_value)))
    cat(sprintf("    %-37s p = %.3g (one-sided %.3g)\n",
                "heterogeneity loading (H0 slope = 0):",
                num(x$blp$heterogeneity_loading$p_value),
                num(x$blp$heterogeneity_loading$p_value_one_sided)))
  }
  if (!is.null(x$gates)) {
    if (isTRUE(x$gates$degenerate)) {
      cat("  GATES: not estimable (tau-hat near-constant).\n")
    } else {
      bw <- x$gates$top_bottom_test
      k <- nrow(x$gates$group_ate)
      lab <- sprintf("Q%d - Q1", k)
      cat(sprintf("  GATES %s: %.4g [%.4g, %.4g] (p = %.3g, one-sided %.3g)\n",
                  lab, num(bw$diff), num(bw$ci_lower), num(bw$ci_upper),
                  num(bw$p_value), num(bw$p_value_one_sided)))
    }
  }
  # One label, recorded once; it already reflects whether the replicates were
  # usable, so a run that could not form them does not announce a bootstrap.
  basis <- x$metadata$ci_method
  if (!is.null(basis) && !is.na(basis)) {
    if (grepl("^half-sample bootstrap", basis))
      cat(sprintf("  intervals: %s, R = %d\n", basis,
                  int(x$metadata$n_bootstrap)))
    else
      cat(sprintf("  intervals: %s\n", basis))
  }
  invisible(x)
}
