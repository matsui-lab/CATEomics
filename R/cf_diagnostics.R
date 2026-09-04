#' Post-fit causal diagnostics for a fitted CATE model
#'
#' @description
#' \code{cf_diagnostics()} is the single, on-demand diagnostic entry point
#' for cfomics. It is a \emph{result-interpretation aid} run \emph{after}
#' fitting, not a pre-flight gate. It reports four complementary views of
#' the assumptions behind a fitted model:
#'
#' \enumerate{
#'   \item \strong{Propensity-score overlap} -- the treated/control
#'     propensity-score distributions and their common-support region
#'     (positivity / overlap).
#'   \item \strong{Covariate balance} -- standardized mean differences
#'     (SMD) between treatment groups on the preprocessed covariates,
#'     ranked, top-\code{top_k} (a Love-plot view). In \eqn{n \ll p} omics
#'     the propensity model cannot absorb every feature, so residual
#'     imbalance is the honest thing to surface.
#'   \item \strong{Requested-variable overlap} -- for a user-supplied set
#'     of \code{variables}, the treated/control distributions on the scale
#'     the model actually saw (preprocessed).
#'   \item \strong{E-value} -- VanderWeele & Ding (2017) sensitivity of the
#'     ATE to unmeasured confounding.
#' }
#'
#' It provides a single entry point for propensity overlap, covariate
#' balance, requested-variable overlap, and the E-value in one call.
#'
#' @param model A fitted \code{cf_model} / \code{cfomics_result}.
#' @param variables Optional character vector of covariate names for the
#'   requested-variable overlap view. Must be columns of the model's
#'   preprocessed \code{X}. \code{NULL} (default) skips that view.
#' @param top_k Integer; how many covariates to keep in the SMD Love view
#'   (ranked by \eqn{|SMD|}, default \code{20}).
#' @param smd_threshold Numeric; \eqn{|SMD|} imbalance reference line
#'   (default \code{0.1}).
#' @param evalue_alpha Numeric; significance level for the CI-based E-value
#'   (default \code{0.05}); retained for interface stability.
#' @param ... Ignored.
#'
#' @return An object of class \code{c("cfomics_diagnostics", "list")} with
#'   elements \code{overlap}, \code{smd}, \code{variable_overlap} (or
#'   \code{NULL}), \code{evalue}, and \code{metadata}.
#'
#' @references
#' VanderWeele TJ, Ding P (2017). Sensitivity Analysis in Observational
#' Research: Introducing the E-Value. \emph{Ann Intern Med} 167(4):268-274.
#'
#' @examples
#' # Post-fit assumption checks on a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf", n_features = 40)
#' cf_diagnostics(fit, variables = "age")
#'
#' @export
cf_diagnostics <- function(model,
                           variables     = NULL,
                           top_k         = 20L,
                           smd_threshold = 0.1,
                           evalue_alpha  = 0.05,
                           ...) {
  if (!inherits(model, "cf_model")) {
    rlang::abort(
      "cf_diagnostics requires a fitted cf_model (from cf_fit / cf_compare).",
      class = "cfomics_invalid_input"
    )
  }
  top_k         <- as.integer(top_k)[1L]
  smd_threshold <- as.numeric(smd_threshold)[1L]

  X <- model$data$X
  T_vec <- model$data$T
  if (is.null(X) || is.null(T_vec)) {
    rlang::abort(
      c("cf_diagnostics requires model$data$X and model$data$T.",
        i = "Fit via cf_fit() / cf_compare() so covariates are cached."),
      class = "cfomics_diagnostics_no_data"
    )
  }
  X <- if (is.matrix(X)) X else as.matrix(X)
  T_vec <- as.numeric(T_vec)

  overlap  <- .cf_diag_overlap(model, X, T_vec)
  smd      <- .cf_diag_smd(X, T_vec, top_k, smd_threshold)
  var_ovl  <- if (!is.null(variables))
    .cf_diag_variable_overlap(X, T_vec, variables) else NULL
  evalue   <- .cf_diag_evalue(model, evalue_alpha)

  structure(
    list(
      overlap          = overlap,
      smd              = smd,
      variable_overlap = var_ovl,
      evalue           = evalue,
      metadata         = list(
        method        = model$method %||% NA_character_,
        n             = length(T_vec),
        p             = ncol(X),
        top_k         = top_k,
        smd_threshold = smd_threshold,
        variables     = variables
      )
    ),
    class = c("cfomics_diagnostics", "list")
  )
}


# -------------------------------------------------------------------------
# 1. Propensity-score overlap
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_diag_overlap <- function(model, X, T_vec) {
  # Prefer the cross-fit propensity the model already estimated; fall back
  # to a logistic PS on the preprocessed X only when unavailable.
  ps <- model$res$nuisance$pi_oof
  ps_source <- "model_nuisance"
  if (is.null(ps) || length(ps) != length(T_vec) || all(is.na(ps))) {
    # No cached propensity: fit a quick logistic PS on the preprocessed X.
    # glm.fit has no predict() method, so form the linear predictor directly.
    ps <- tryCatch({
      fit <- stats::glm.fit(cbind(1, X), T_vec, family = stats::binomial())
      as.numeric(stats::binomial()$linkinv(cbind(1, X) %*% fit$coefficients))
    }, error = function(e) rep(mean(T_vec), length(T_vec)))
    ps_source <- "logistic_fallback"
  }
  ps <- as.numeric(ps)
  ps_t1 <- ps[T_vec == 1]; ps_t0 <- ps[T_vec == 0]
  overlap_region <- c(lower = max(min(ps_t0), min(ps_t1)),
                      upper = min(max(ps_t0), max(ps_t1)))
  in_support <- mean(ps >= overlap_region["lower"] &
                     ps <= overlap_region["upper"])
  extreme_idx <- ps < 0.01 | ps > 0.99
  list(
    propensity_scores = ps,
    treatment         = T_vec,
    ps_source         = ps_source,
    overlap_region    = overlap_region,
    support_fraction  = in_support,
    n_extreme         = sum(extreme_idx),
    extreme_pct       = mean(extreme_idx),
    no_overlap        = unname(overlap_region["upper"] <
                                 overlap_region["lower"])
  )
}


# -------------------------------------------------------------------------
# 2. Covariate balance (SMD Love view, top-K)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_diag_smd <- function(X, T_vec, top_k, smd_threshold) {
  feats <- colnames(X)
  if (is.null(feats)) feats <- paste0("X", seq_len(ncol(X)))
  i1 <- T_vec == 1; i0 <- T_vec == 0
  smd <- vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    m1 <- mean(x[i1], na.rm = TRUE); m0 <- mean(x[i0], na.rm = TRUE)
    v1 <- stats::var(x[i1], na.rm = TRUE); v0 <- stats::var(x[i0], na.rm = TRUE)
    sp <- sqrt((v1 + v0) / 2)
    if (is.finite(sp) && sp > 0) (m1 - m0) / sp else 0
  }, numeric(1))
  df <- data.frame(
    variable = feats,
    smd      = smd,
    abs_smd  = abs(smd),
    balanced = abs(smd) < smd_threshold,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$abs_smd, decreasing = TRUE), , drop = FALSE]
  n_total <- nrow(df)
  if (n_total > top_k) df <- df[seq_len(top_k), , drop = FALSE]
  rownames(df) <- NULL
  attr(df, "n_total")     <- n_total
  attr(df, "n_imbalanced") <- sum(!df$balanced)
  attr(df, "threshold")   <- smd_threshold
  df
}


# -------------------------------------------------------------------------
# 3. Requested-variable overlap (preprocessed scale)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_diag_variable_overlap <- function(X, T_vec, variables) {
  feats <- colnames(X)
  missing <- setdiff(variables, feats)
  if (length(missing) > 0L) {
    rlang::abort(
      sprintf("cf_diagnostics: requested variables not in preprocessed X: %s",
              paste(missing, collapse = ", ")),
      class = "cfomics_diagnostics_unknown_variables"
    )
  }
  do.call(rbind, lapply(variables, function(v) {
    x <- X[, v]
    data.frame(
      variable = v,
      value    = x,
      group    = ifelse(T_vec == 1, "Treated", "Control"),
      stringsAsFactors = FALSE
    )
  }))
}


# -------------------------------------------------------------------------
# 4. E-value (VanderWeele-Ding 2017)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_diag_evalue <- function(model, alpha = 0.05) {
  ate      <- model$res$ate
  ci_lower <- model$res$summary$ate_ci_lower
  ci_upper <- model$res$summary$ate_ci_upper
  if (is.null(ate) || !is.finite(ate)) {
    return(list(available = FALSE))
  }

  # The stored ATE CI is a 95% interval; rescale it to the requested (1 - alpha)
  # level so the CI-based E-value reflects alpha. Each half-width is scaled
  # independently by the z-ratio, which preserves any asymmetry in the stored
  # bounds (bcf / cevae / bootstrap percentile CIs are not symmetric about the
  # point estimate) and reproduces the stored bounds exactly at alpha = 0.05.
  # A degenerate or out-of-range alpha, or an unrecoverable interval, keeps the
  # stored bounds.
  if (is.finite(alpha) && alpha > 0 && alpha < 1 &&
      !is.null(ci_lower) && !is.null(ci_upper) &&
      is.finite(ci_lower) && is.finite(ci_upper) && ci_upper > ci_lower) {
    ratio    <- stats::qnorm(1 - alpha / 2) / stats::qnorm(0.975)
    ci_lower <- ate - (ate - ci_lower) * ratio
    ci_upper <- ate + (ci_upper - ate) * ratio
  }

  # standardise the ATE by the SD of the *observed* outcome, then
  # convert the standardized mean difference d to an approximate risk ratio
  # via the EValue-package canonical mapping RR ~= exp(0.91 * d) (the old
  # exp(d) overstated the RR). See VanderWeele & Ding (2017), and the
  # EValue::evalues.MD() documentation.
  y <- model$data$Y
  sd_y <- if (!is.null(y)) stats::sd(as.numeric(y), na.rm = TRUE) else NA_real_
  if (is.null(y) || !is.finite(sd_y) || sd_y == 0) {
    # Fall back to the counterfactual SD only when raw Y is unavailable.
    sd_y <- stats::sd(c(model$res$mu0_hat, model$res$mu1_hat), na.rm = TRUE)
    if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  }

  rr_point <- exp(0.91 * abs(ate) / sd_y)
  rr_ci <- if (!is.null(ci_lower) && is.finite(ci_lower) && ci_lower > 0) {
    exp(0.91 * ci_lower / sd_y)
  } else if (!is.null(ci_upper) && is.finite(ci_upper) && ci_upper < 0) {
    exp(0.91 * abs(ci_upper) / sd_y)
  } else {
    1  # CI crosses the null -> CI E-value is 1.
  }

  evalue_fn <- function(rr) if (rr <= 1) 1 else rr + sqrt(rr * (rr - 1))
  ev_point <- evalue_fn(rr_point)
  ev_ci    <- evalue_fn(rr_ci)

  interp <- if (ev_point > 2) {
    "Robust: an unmeasured confounder would need strong associations with both treatment and outcome to explain away this effect."
  } else if (ev_point > 1.5) {
    "Moderate: a confounder with moderate associations could explain this effect."
  } else {
    "Weak: even a weak unmeasured confounder could explain away this effect."
  }

  list(
    available      = TRUE,
    evalue_point   = ev_point,
    evalue_ci      = ev_ci,
    rr_point       = rr_point,
    rr_ci          = rr_ci,
    ate            = ate,
    # unname so the vector carries clean "lower"/"upper" names rather than
    # compound "lower.estimate" (grf's ATE / CI bounds arrive named "estimate").
    ci             = c(lower = unname(ci_lower %||% NA_real_),
                       upper = unname(ci_upper %||% NA_real_)),
    sd_y           = sd_y,
    interpretation = interp
  )
}


# -------------------------------------------------------------------------
# print (object-lean, one line per element)
# -------------------------------------------------------------------------

#' @export
print.cfomics_diagnostics <- function(x, ...) {
  m <- x$metadata
  cat("<cfomics_diagnostics>: method = ", m$method %||% "<unknown>",
      ", n = ", m$n, ", p = ", m$p, "\n", sep = "")

  ov <- x$overlap
  cat(sprintf("  overlap : common support = [%.3f, %.3f], %.0f%% in support, %d extreme PS (%s)\n",
              ov$overlap_region["lower"], ov$overlap_region["upper"],
              100 * ov$support_fraction, ov$n_extreme, ov$ps_source))
  if (isTRUE(ov$no_overlap))
    cat("            [WARNING] treated/control PS ranges do not overlap.\n")

  sm <- x$smd
  cat(sprintf("  balance : %d/%d covariates |SMD| >= %.2f (top-%d shown)\n",
              attr(sm, "n_imbalanced"), attr(sm, "n_total"),
              attr(sm, "threshold"), nrow(sm)))

  if (!is.null(x$variable_overlap)) {
    cat(sprintf("  variables: overlap view for %d requested covariate(s)\n",
                length(m$variables)))
  }

  ev <- x$evalue
  if (isTRUE(ev$available)) {
    cat(sprintf("  E-value : point = %.2f, CI = %.2f (ATE = %.3f, sd_Y = %.3f)\n",
                ev$evalue_point, ev$evalue_ci, ev$ate, ev$sd_y))
  } else {
    cat("  E-value : not available (no finite ATE).\n")
  }
  invisible(x)
}


# -------------------------------------------------------------------------
# plot (one view per call)
# -------------------------------------------------------------------------

#' Plot one view of a cf_diagnostics object.
#'
#' @param x A \code{cfomics_diagnostics} object (from \code{\link{cf_diagnostics}}).
#' @param type One of \code{"overlap"} (PS density, default),
#'   \code{"smd"} (Love dot plot), \code{"variable"} (requested-variable
#'   distributions), \code{"evalue"} (E-value bias contour).
#' @param ... Ignored.
#' @return A \code{ggplot} object (invisibly), or \code{NULL} if ggplot2 is
#'   unavailable / the requested view has no data.
#' @examples
#' # Overlap (positivity) view from a fast R-native causal forest fit.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf", n_features = 40)
#' dg <- cf_diagnostics(fit)
#' cf_plot_diagnostics(dg, type = "overlap")
#' @export
cf_plot_diagnostics <- function(x,
                                type = c("overlap", "smd",
                                         "variable", "evalue"),
                                ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 is required for cf_diagnostics plots.")
    return(invisible(NULL))
  }
  p <- switch(type,
    overlap  = .cf_diag_plot_overlap(x$overlap),
    smd      = .cf_diag_plot_smd(x$smd),
    variable = .cf_diag_plot_variable(x$variable_overlap),
    evalue   = .cf_diag_plot_evalue(x$evalue)
  )
  if (is.null(p)) return(invisible(NULL))
  p
}

#' @keywords internal
#' @noRd
.cf_diag_plot_overlap <- function(ov) {
  df <- data.frame(ps = ov$propensity_scores,
                   group = ifelse(ov$treatment == 1, "Treated", "Control"))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["ps"]], fill = .data[["group"]])) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = ov$overlap_region,
                        linetype = "dashed", alpha = 0.7) +
    ggplot2::labs(x = "Propensity score", y = "Density",
                  title = "Treatment overlap",
                  subtitle = sprintf("common support [%.3f, %.3f], %.0f%% in support",
                                     ov$overlap_region["lower"],
                                     ov$overlap_region["upper"],
                                     100 * ov$support_fraction),
                  fill = NULL) +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.cf_diag_plot_smd <- function(sm) {
  thr <- attr(sm, "threshold")
  # Love-style dot plot: |SMD| descending, threshold reference line.
  df <- sm
  df$variable <- factor(df$variable, levels = rev(df$variable))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["abs_smd"]],
                                   y = .data[["variable"]])) +
    ggplot2::geom_vline(xintercept = thr, linetype = "dashed",
                        color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(color = .data[["balanced"]]), size = 2.5) +
    ggplot2::labs(x = "|Standardized mean difference|", y = NULL,
                  title = "Covariate balance",
                  subtitle = sprintf("top %d of %d covariates (reference |SMD| = %.2f)",
                                     nrow(df), attr(sm, "n_total"), thr),
                  color = "Balanced") +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.cf_diag_plot_variable <- function(vo) {
  if (is.null(vo)) {
    message("cf_diagnostics: no requested-variable overlap to plot ",
            "(pass `variables =`).")
    return(NULL)
  }
  ggplot2::ggplot(vo, ggplot2::aes(x = .data[["value"]],
                                   fill = .data[["group"]])) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::labs(x = "Value (preprocessed)", y = "Density",
                  title = "Requested-variable overlap", fill = NULL) +
    ggplot2::theme_minimal()
}

#' @keywords internal
#' @noRd
.cf_diag_plot_evalue <- function(ev) {
  if (!isTRUE(ev$available)) {
    message("cf_diagnostics: E-value not available; nothing to plot.")
    return(NULL)
  }
  # canonical VanderWeele-Ding bias contour. For observed risk
  # ratio R, the confounder-treatment RR needed at a given confounder-
  # outcome RR x (x > R) is y = R (x - 1) / (x - R); the E-value point
  # (E, E) lies exactly on this curve.
  R <- ev$rr_point
  E <- ev$evalue_point
  x_seq <- seq(R + 1e-6, max(E * 2, R + 3), length.out = 200)
  curve <- data.frame(x = x_seq, y = R * (x_seq - 1) / (x_seq - R))
  curve <- curve[is.finite(curve$y) & curve$y >= 1, , drop = FALSE]
  ggplot2::ggplot(curve, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(data = data.frame(x = E, y = E),
                        ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                        color = "firebrick", size = 3) +
    ggplot2::annotate("text", x = E, y = E, label = sprintf("E = %.2f", E),
                      vjust = -1, color = "firebrick") +
    ggplot2::labs(x = "Confounder-outcome risk ratio",
                  y = "Confounder-treatment risk ratio",
                  title = "E-value sensitivity contour",
                  subtitle = sprintf("ATE = %.3f, E-value = %.2f (CI %.2f)",
                                     ev$ate, ev$evalue_point, ev$evalue_ci)) +
    ggplot2::theme_minimal()
}
