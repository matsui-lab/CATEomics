#' Evaluate a CATE model with honest metrics
#'
#' @description
#' Compute the default three evaluation metrics for CATE-capable models:
#' \itemize{
#'   \item \strong{RATE} (Yadlowsky, Fleming, Shah, Brunskill & Wager 2025,
#'         JASA). Targeting contrast for the top-\code{target_quantile}
#'         fraction by predicted treatment effect; also returns a per-quantile
#'         curve \code{rate_curve} and the area under it, \code{rate$autoc}
#'         (AUTOC = mean over the q-grid of TOC(q), with an influence-function
#'         SE).
#'   \item \strong{BLP} (Chernozhukov, Demirer, Duflo & Fernández-Val 2025,
#'         \emph{Econometrica} 93(4); Generic ML). The single slope of
#'         \eqn{Y \sim T + (T-\bar{T}) \hat\tau_{centered}}, tested against
#'         H0 = 1 (calibration) and H0 = 0 (heterogeneity).
#'   \item \strong{GATES} group-ATE on tau-hat quantile bins plus a Best-Worst
#'         heterogeneity test (non-parametric; degrades gracefully when tau-hat
#'         is near-constant).
#' }
#'
#' Operates on the OOF tau-hat saved by \code{cf_compare} /
#' \code{cf_fit(partition = 'oof')} (\code{model$res$tau_hat_oof}). When
#' the model lacks OOF tau-hat, falls back to in-sample tau-hat
#' (\code{model$res$tau_hat}) and emits \code{cfomics_partition_none_warning}.
#'
#' @details
#' nuisance for DR pseudo-outcomes / GATES / BLP is the independent
#' shared evaluation nuisance (\code{model$res$nuisance}: cf_compare's shared
#' \code{nuisance_oof} or a cf_fit-side independent eval nuisance), recomputed
#' via \code{eval_nuisance = "rf"} OOF when absent. The method is never scored
#' with its own internal nuisance.
#'
#' @param model \code{cf_model} from \code{cf_fit(partition = 'oof')}
#'   or an element of \code{cf_compare()$models}.
#' @param metrics Character vector subset of
#'   \code{c("rate", "blp", "gates")} (default all 3).
#' @param n_groups Integer, GATES quantile bins (default 5).
#' @param target_quantile Numeric in (0,1], top-q fraction for RATE
#'   point estimate (default 0.2).
#' @param rate_quantiles Numeric grid in (0,1], RATE curve quantiles
#'   (default \code{seq(0.1, 1.0, by = 0.1)}).
#' @param direction Character, \code{"positive"} (treat top-q\% by
#'   tau-hat) or \code{"negative"} (treat bottom-q\%).
#' @param eval_nuisance Base learner for shared nuisance when
#'   recomputation is needed (default \code{"rf"}).
#' @param eval_propensity_score Optional known propensity score (scalar in
#'   (0,1) or length-n vector) for the independent shared evaluation nuisance;
#'   used in a randomised-trial context. \code{NULL} (default)
#'   estimates the propensity internally.
#' @param random_state Integer seed or NULL.
#' @param device One of \code{"auto"} (default), \code{"cpu"}, \code{"gpu"};
#'   broadcast to internal cf_fit calls for any re-fitting required by metrics
#'   (additive; the metric path currently performs no re-fitting, so this is
#'   pass-through).
#' @param ... Ignored.
#' @return Object of class \code{c("cfomics_evaluate", "list")}.
#' @seealso \code{\link{cf_compare}}, \code{\link{cf_select}},
#'   \code{\link{cf_rule}}.
#' @examples
#' # Evaluate a fast R-native causal forest with an independent nuisance.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' # known randomised propensity; independent nuisance for the score
#' ev <- cf_evaluate(fit, eval_nuisance = "rf", eval_propensity_score = 0.5)
#' ev
#' @export
cf_evaluate <- function(model,
                        metrics         = c("rate", "blp", "gates"),
                        n_groups        = 5L,
                        target_quantile = 0.2,
                        rate_quantiles  = seq(0.1, 1.0, by = 0.1),
                        direction       = c("positive", "negative"),
                        eval_nuisance   = "rf",
                        eval_propensity_score = NULL,
                        random_state    = NULL,
                        device          = c("auto", "cpu", "gpu"),
                        ...) {
  device <- match.arg(device)
  if (!inherits(model, "cf_model")) {
    rlang::abort(
      "cf_evaluate requires a cf_model from cf_fit() or cf_compare()$models.",
      class = "cfomics_invalid_input"
    )
  }
  metrics <- match.arg(metrics, c("rate", "blp", "gates"),
                       several.ok = TRUE)
  direction <- match.arg(direction)

  # tau-hat source: OOF if available (res$tau_hat_oof), else in-sample.
  tau_oof <- model$res$tau_hat_oof
  if (is.null(tau_oof)) {
    rlang::warn(
      c("Model was fit without partition='oof'; using in-sample tau-hat.",
        i = "Evaluation metrics may be optimistically biased.",
        i = "Re-fit via cf_compare() or cf_fit(partition='oof') for honest evaluation."),
      class = "cfomics_partition_none_warning",
      .frequency = "once",
      .frequency_id = "cfomics_partition_none_warning"
    )
    tau_oof <- model$res$tau_hat
    if (is.null(tau_oof) || all(is.na(tau_oof))) {
      rlang::abort(
        "Model lacks usable tau-hat (tau_hat_oof and tau_hat both NA).",
        class = "cfomics_evaluate_no_tau"
      )
    }
  }
  if (direction == "negative") {
    tau_oof <- -tau_oof
  }

  # Source data (Y, T, X): the model-visible training data.
  Y <- NULL; T <- NULL; X <- NULL
  if (TRUE) {
    if (!is.null(model$data$Y) && !is.null(model$data$T)) {
      Y <- as.numeric(model$data$Y)
      T <- as.integer(model$data$T)
      ## pre-expand factor / character columns to one-hot numeric before
      ## `as.matrix`,
      ## consistent with cf_compare entry and `.compute_common_nuisance_oof`.
      ## Otherwise `cf_evaluate(model)` on data containing factor columns
      ## (e.g. sex / batch) would fail inside `.compute_common_nuisance_oof`
      ## fallback regression with `NA / character mismatch`.
      raw_X <- model$data$X_preprocessed %||% model$data$X
      X <- suppressWarnings(.expand_factor_columns(raw_X)$X)
    } else {
      rlang::abort(
        c("cf_evaluate cannot recover (Y, T, X) from this model.",
          i = "Models created by cf_compare()$models or cf_fit(..., return_data = TRUE) are recommended."),
        class = "cfomics_evaluate_missing_data"
      )
    }
  }

  # Nuisance for DR pseudo-outcome / GATES / BLP regression. Use the
  # independent shared evaluation nuisance (res$nuisance holds cf_compare's
  # shared nuisance_oof, or a cf_fit-side independent eval nuisance), not the
  # method's own internal nuisance — self-scoring would be optimistic. Falls
  # back to an independent OOF nuisance when absent.
  nuis <- model$res$nuisance
  if (is.null(nuis)) {
    # binary Y = risk difference. When re-deriving the nuisance
    # (the model did not carry one) detect a 0/1 outcome so mu0/mu1 are fit as
    # probability classifiers, consistent with cf_fit / cf_compare. The DR
    # pseudo-outcome / RATE / GATES / BLP formulas are outcome-type agnostic.
    yb <- Y[!is.na(Y)]
    binary_outcome <- length(unique(yb)) == 2L && all(yb %in% c(0, 1))
    folds <- .create_stratified_folds(T, n_folds = 5L,
                                       random_state = random_state)
    nuis  <- .compute_common_nuisance_oof(X, T, Y, folds,
                                           eval_nuisance = eval_nuisance,
                                           eval_propensity_score = eval_propensity_score,
                                           random_state = random_state,
                                           binary_outcome = binary_outcome)
  }
  Y_DR <- .compute_dr_pseudo_outcome(Y, T, nuis$mu0_oof, nuis$mu1_oof,
                                       nuis$pi_oof)

  out <- list()

  if ("rate" %in% metrics) {
    out$rate <- .compute_rate_metric(tau_oof, Y_DR, target_quantile)
    out$rate$autoc <- .compute_autoc(tau_oof, Y_DR, rate_quantiles)
    out$rate_curve <- .compute_rate_curve(tau_oof, Y_DR, rate_quantiles)
  }
  if ("blp" %in% metrics) {
    # Calibration slope and BLP-HTE are the SAME regression beta,
    # tested against H0 = 1 (calibration) and H0 = 0 (heterogeneity). Report
    # one BLP metric carrying beta + both tests + the binned calibration data.
    out$blp <- .compute_blp_metric(tau_oof, n_groups, Y_DR)
  }
  if ("gates" %in% metrics) {
    out$gates <- .compute_gates_metric(tau_oof, Y_DR, n_groups)
  }

  out$metadata <- list(
    method          = model$method,
    partition       = model$meta$partition %||% "none",
    n_groups        = as.integer(n_groups),
    target_quantile = target_quantile,
    rate_quantiles  = rate_quantiles,
    direction       = direction,
    eval_nuisance   = eval_nuisance
  )
  class(out) <- c("cfomics_evaluate", "list")
  out
}

# -------------------------------------------------------------------------
# Metric implementations
# -------------------------------------------------------------------------

#' Yadlowsky 2024 RATE estimator at a single quantile q.
#' @keywords internal
#' @noRd
.compute_rate_metric <- function(tau, Y_DR, target_quantile) {
  q <- target_quantile
  Fhat <- stats::ecdf(tau)
  weight <- as.numeric(Fhat(tau) >= 1 - q)
  rate_est <- mean(weight * Y_DR) - q * mean(Y_DR)
  # Naive SE via empirical variance of (weight * Y_DR - q * Y_DR).
  per_obs <- weight * Y_DR - q * Y_DR
  se      <- stats::sd(per_obs) / sqrt(length(per_obs))
  list(
    value           = rate_est,
    ci_lower        = rate_est - 1.96 * se,
    ci_upper        = rate_est + 1.96 * se,
    target_quantile = q,
    ci_method       = "asymptotic Wald"
  )
}

#' Per-quantile RATE curve.
#' @keywords internal
#' @noRd
.compute_rate_curve <- function(tau, Y_DR, rate_quantiles) {
  rows <- lapply(rate_quantiles, function(q) {
    r <- .compute_rate_metric(tau, Y_DR, q)
    data.frame(quantile = q, rate = r$value,
                ci_lower = r$ci_lower, ci_upper = r$ci_upper)
  })
  do.call(rbind, rows)
}

#' Area under the TOC curve (AUTOC) with an influence-function SE.
#'
#' AUTOC = mean over the q-grid of TOC(q), where
#' TOC(q) = (1/q) E[Y_DR 1{F(tau) >= 1-q}] - E[Y_DR]. The per-observation
#' influence is psi_i = mean_q[(1/q) 1{F(tau_i) >= 1-q} Y_DR_i - Y_DR_i]; AUTOC
#' = mean(psi), SE = sd(psi)/sqrt(n). The grid (q >= smallest rate_quantile)
#' bounds the 1/q weight.
#' @keywords internal
#' @noRd
.compute_autoc <- function(tau, Y_DR, rate_quantiles) {
  p <- stats::ecdf(tau)(tau)
  n <- length(Y_DR)
  per_q <- vapply(rate_quantiles, function(q) {
    w <- as.numeric(p >= 1 - q)
    (1 / q) * (w * Y_DR) - Y_DR          # per-obs TOC(q) contribution
  }, numeric(n))
  psi   <- rowMeans(per_q)                # per-obs AUTOC influence
  autoc <- mean(psi)
  se    <- stats::sd(psi) / sqrt(n)
  list(
    value     = autoc,
    se        = se,
    ci_lower  = autoc - 1.96 * se,
    ci_upper  = autoc + 1.96 * se,
    ci_method = "influence-function"
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
#' against two nulls: H0 = 1 (calibration; is the model correctly scaled?) and
#' H0 = 0 (heterogeneity; does the model capture detectable HTE?).
#' Chernozhukov 2018. Also returns the binned observed (DR pseudo-outcome) vs
#' predicted tau data for a calibration plot.
#' @keywords internal
#' @noRd
.compute_blp_metric <- function(tau, n_groups, Y_DR) {
  tau_centered <- tau - mean(tau)
  # DR-based BLP. A raw lm(Y ~ T + (T - T_bar) * tau_c)
  # uses a constant treatment mean and is only valid under a constant
  # propensity (RCT); under confounding it is biased and inconsistent with the
  # DR-based GATES / RATE. Regress the DR pseudo-outcome (E[Y_DR | X] = tau(X))
  # on the centred model tau-hat; the slope is the best-linear calibration
  # coefficient, tested vs 1 (calibration) and vs 0 (heterogeneity).
  reg <- stats::lm(Y_DR ~ tau_centered)
  beta <- unname(stats::coef(reg)[2])
  se   <- summary(reg)$coefficients[2, 2]
  z_cal <- (beta - 1) / se; p_cal <- 2 * stats::pnorm(-abs(z_cal))
  z_het <- beta / se;       p_het <- 2 * stats::pnorm(-abs(z_het))

  bins <- .cf_safe_quantile_bins(tau, n_groups)
  plot_data <- if (bins$degenerate) NULL else do.call(rbind, lapply(
    seq_len(bins$n_groups), function(k) {
      idx <- which(bins$bin == k)
      if (length(idx) == 0L) return(NULL)
      obs    <- mean(Y_DR[idx])
      se_obs <- stats::sd(Y_DR[idx]) / sqrt(length(idx))
      data.frame(bin = k, predicted_mean = mean(tau[idx]), observed_mean = obs,
                 ci_lower = obs - 1.96 * se_obs, ci_upper = obs + 1.96 * se_obs)
    }))

  list(
    slope = list(
      estimate  = beta,
      se        = se,
      ci_lower  = beta - 1.96 * se,
      ci_upper  = beta + 1.96 * se,
      ci_method = "asymptotic Wald"
    ),
    calibration   = list(h0 = 1, z = unname(z_cal), p_value = unname(p_cal)),
    heterogeneity = list(h0 = 0, z = unname(z_het), p_value = unname(p_het)),
    plot_data = plot_data
  )
}

#' GATES quantile group ATEs + Best-Worst test.
#' @keywords internal
#' @noRd
.compute_gates_metric <- function(tau, Y_DR, n_groups) {
  # near-constant tau -> non-unique quantile breaks. Degrade
  # gracefully instead of crashing cut().
  bins <- .cf_safe_quantile_bins(tau, n_groups)
  if (bins$degenerate) {
    return(list(
      group_ate = NULL,
      best_worst_test = NULL,
      degenerate = TRUE,
      message = "tau-hat is near-constant; GATES heterogeneity is not estimable."
    ))
  }
  bin <- bins$bin
  eff_groups <- bins$n_groups
  rows <- do.call(rbind, lapply(seq_len(eff_groups), function(k) {
    idx <- which(bin == k)
    if (length(idx) == 0L) {
      return(data.frame(group = k, estimate = NA_real_,
                         ci_lower = NA_real_, ci_upper = NA_real_))
    }
    est <- mean(Y_DR[idx])
    se  <- stats::sd(Y_DR[idx]) / sqrt(length(idx))
    data.frame(group = k, estimate = est,
                ci_lower = est - 1.96 * se,
                ci_upper = est + 1.96 * se)
  }))

  idx_top    <- which(bin == eff_groups)
  idx_bottom <- which(bin == 1L)
  diff_est <- mean(Y_DR[idx_top]) - mean(Y_DR[idx_bottom])
  diff_se  <- sqrt(stats::var(Y_DR[idx_top]) / length(idx_top) +
                     stats::var(Y_DR[idx_bottom]) / length(idx_bottom))
  diff_p   <- 2 * stats::pnorm(-abs(diff_est / diff_se))

  list(
    group_ate = rows,
    best_worst_test = list(
      diff      = diff_est,
      ci_lower  = diff_est - 1.96 * diff_se,
      ci_upper  = diff_est + 1.96 * diff_se,
      p_value   = diff_p
    )
  )
}

#' @export
print.cfomics_evaluate <- function(x, ...) {
  cat("<cfomics_evaluate>: method = ", x$metadata$method,
      ", partition = ", x$metadata$partition, "\n", sep = "")
  if (!is.null(x$rate)) {
    cat(sprintf("  RATE @ q=%.2f: %.4g [%.4g, %.4g] (SE not shown)\n",
                x$rate$target_quantile, x$rate$value,
                x$rate$ci_lower, x$rate$ci_upper))
  }
  if (!is.null(x$blp)) {
    s <- x$blp$slope
    cat(sprintf("  BLP slope: %.4g +/- %.3g [%.4g, %.4g]\n",
                s$estimate, s$se, s$ci_lower, s$ci_upper))
    cat(sprintf("    calibration  (H0 slope = 1): p = %.3g\n",
                x$blp$calibration$p_value))
    cat(sprintf("    heterogeneity(H0 slope = 0): p = %.3g\n",
                x$blp$heterogeneity$p_value))
  }
  if (!is.null(x$gates)) {
    if (isTRUE(x$gates$degenerate)) {
      cat("  GATES: not estimable (tau-hat near-constant).\n")
    } else {
      bw <- x$gates$best_worst_test
      cat(sprintf("  GATES Q%d - Q1: %.4g [%.4g, %.4g] (p = %.3g)\n",
                  nrow(x$gates$group_ate), bw$diff, bw$ci_lower, bw$ci_upper,
                  bw$p_value))
    }
  }
  invisible(x)
}
