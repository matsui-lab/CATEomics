#' Post-fit causal diagnostics for a fitted CATE model
#'
#' @description
#' \code{cf_diagnostics()} runs on a fitted model and reports four views of
#' the assumptions behind it:
#'
#' \enumerate{
#'   \item \strong{Propensity-score overlap} -- the treated/control
#'     propensity-score distributions, their common-support region
#'     (positivity / overlap), and how many observations sit at the
#'     \code{clip} bounds. The scores are the ones the fitted model used and
#'     \code{$overlap$ps_type} records where they come from; an ensemble has
#'     no propensity of its own, so it reports the shared evaluation
#'     nuisance (\code{ps_type = "evaluation nuisance"}, floor 1e-4). Methods that
#'     estimate no propensity (slearner, tlearner, tarnet, dragonnet, cevae,
#'     ganite), and a fit whose propensity cannot be read back, get
#'     \code{$overlap = NULL} and a message.
#'   \item \strong{Covariate balance} -- standardized mean differences
#'     (SMD) between treatment groups on the preprocessed covariates,
#'     ranked, top-\code{top_k} (a Love-plot view). How to read it is under
#'     Details.
#'   \item \strong{Requested-variable overlap} -- for a user-supplied set
#'     of \code{variables}, the treated/control distributions on the scale
#'     the model actually saw (preprocessed).
#'   \item \strong{E-value} -- VanderWeele & Ding (2017) sensitivity of the
#'     ATE to unmeasured confounding. The outcome here is continuous, so the
#'     average effect is first standardised into \eqn{d = |ATE| / s} and
#'     turned into an approximate risk ratio with \eqn{RR = \exp(0.91 d)};
#'     \eqn{s} is the pooled within-arm standard deviation of the observed
#'     outcome, which is the spread that conversion is defined against. The
#'     quantity actually used is reported as \code{sd_basis}, because a
#'     treatment arm with fewer than two usable rows leaves only the
#'     whole-sample standard deviation. The interval contributes a second
#'     E-value, taken at the confidence limit nearer the null; an interval
#'     that holds the null gives \code{1}, and a fit that carries no interval
#'     for the average effect gives \code{NA}. \code{ci_basis} says which of
#'     the three it is.
#' }
#'
#' @details
#' An SMD is the between-arm difference in means over the pooled within-arm
#' standard deviation: how far apart the two arms sit on that covariate. Where
#' it is large, one arm thins out at one end of that covariate's range, and
#' \eqn{\mu_0} / \eqn{\mu_1} at those values come from the fitted form rather
#' than from rows of that arm nearby. The view covers the covariates in the
#' model's preprocessed \code{X} -- the ones passed as \code{X_baseline} /
#' \code{omics}; a variable that was not passed cannot appear here. The
#' reported quantiles are taken over every covariate, not only the rows shown.
#'
#' @param model A fitted \code{cf_model}.
#' @param variables Optional character vector of covariate names for the
#'   requested-variable overlap view. Must be columns of the model's
#'   preprocessed \code{X}. \code{NULL} (default) skips that view.
#' @param top_k Integer; how many covariates to keep in the SMD Love view
#'   (ranked by \eqn{|SMD|}, default \code{20}).
#' @param smd_threshold Numeric; the \eqn{|SMD|} reference line drawn on the
#'   balance plot (default \code{0.1}). It is a reference, not a verdict: the
#'   object reports the distribution of \eqn{|SMD|} and labels no covariate
#'   balanced or unbalanced.
#' @param evalue_alpha Numeric in (0, 1); significance level for the CI-based
#'   E-value (default \code{0.05}). The stored 95\% interval is rescaled to
#'   this level, each half-width separately, which keeps an asymmetric
#'   interval asymmetric and reproduces the stored bounds at 0.05.
#' @param ... Not used; naming an unknown argument is an error.
#'
#' @return An object of class \code{c("cateomics_diagnostics", "list")} with
#'   elements \code{overlap} (\code{NULL} when the method estimates no
#'   propensity), \code{smd}, \code{variable_overlap} (or \code{NULL}),
#'   \code{evalue}, and \code{metadata}.
#'
#'   \code{overlap$overlap_region} is the intersection of the two arms'
#'   propensity ranges: its \code{lower} is the larger of the two minima and
#'   its \code{upper} the smaller of the two maxima, so outside it one arm
#'   has no observations at that propensity. \code{support_fraction} is the
#'   share of rows inside it, and is \code{0} exactly when the two ranges do
#'   not meet.
#'
#'   \code{smd} is a data frame of \code{variable}, \code{smd} and
#'   \code{abs_smd} for the top-\code{top_k} covariates, carrying the
#'   \code{median}, \code{p90}, \code{max} and \code{n_total} of
#'   \eqn{|SMD|} over all of them as attributes.
#'
#'   \code{evalue} carries \code{available}, and when that is \code{TRUE}:
#'   \code{evalue_point} and \code{evalue_ci} (the E-value of the point
#'   estimate and of the confidence limit nearer the null; the second is
#'   \code{1} when the interval holds the null and \code{NA} when the fit
#'   carries no interval for the average effect, which is the case when none
#'   was asked for), \code{rr_point} and \code{rr_ci} (the risk ratios they
#'   are computed from), \code{ate} and \code{ci} (the effect and the
#'   interval that were converted, the interval rescaled to
#'   \code{evalue_alpha}), \code{sd_y} and \code{sd_basis} (the spread the
#'   effect was divided by, and which spread that is), and \code{ci_basis}
#'   (what the interval contributed, or why it contributed nothing).
#'
#' @references
#' VanderWeele TJ, Ding P (2017) Sensitivity Analysis in Observational
#' Research: Introducing the E-Value, \emph{Ann Intern Med} 167(4):268-274
#' (the E-value).
#'
#' Chinn S (2000) A simple method for converting an odds ratio to effect size
#' for use in meta-analysis, \emph{Stat Med} 19(22):3127-3131, and
#' VanderWeele TJ (2017) On a square-root transformation of the odds ratio for
#' a common outcome, \emph{Epidemiology} 28(6):e58-e60 (together they give the
#' \code{RR ~ exp(0.91 d)} conversion used for a continuous outcome).
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
      c("cf_diagnostics requires a fitted cf_model (from cf_fit / cf_compare).",
        i = "For a comparison, pass one candidate: cmp$models[[...]]."),
      class = "cateomics_invalid_input"
    )
  }
  .cf_reject_dots(..., .fn = "cf_diagnostics()")
  top_k         <- as.integer(top_k)[1L]
  smd_threshold <- as.numeric(smd_threshold)[1L]
  evalue_alpha  <- as.numeric(evalue_alpha)[1L]
  # A level outside (0, 1) is not a level, and silently keeping the stored 95%
  # interval would report an interval the caller did not ask for.
  if (!is.finite(evalue_alpha) || evalue_alpha <= 0 || evalue_alpha >= 1) {
    rlang::abort("evalue_alpha must be in (0, 1).",
                 class = "cateomics_invalid_input")
  }

  X <- model$data$X_preprocessed %||% model$data$X
  T_vec <- model$data$T
  if (is.null(X) || is.null(T_vec)) {
    rlang::abort(
      c("cf_diagnostics requires model$data$X_preprocessed and model$data$T.",
        i = "Fit via cf_fit() / cf_compare() so covariates are cached."),
      class = "cateomics_diagnostics_no_data"
    )
  }
  X <- if (is.matrix(X)) X else as.matrix(X)
  T_vec <- as.numeric(T_vec)

  # A method without a propensity loses the overlap view, not the whole
  # report: balance, requested variables and the E-value need no propensity.
  overlap <- tryCatch(
    .cf_diag_overlap(model, T_vec),
    cateomics_overlap_no_propensity = function(cnd) {
      rlang::inform(conditionMessage(cnd),
                    class = "cateomics_overlap_view_skipped")
      NULL
    })
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
    class = c("cateomics_diagnostics", "list")
  )
}


# -------------------------------------------------------------------------
# 1. Propensity-score overlap
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cf_diag_overlap <- function(model, T_vec) {
  ps      <- model$res$ps
  ps_type <- model$res$ps_type %||% NA_character_
  clip    <- model$meta$spec$clip %||% c(0.025, 0.975)

  # An ensemble combines per-method tau and holds no backend, so it has no
  # propensity of its own. The shared evaluation nuisance is the propensity its
  # own losses already use; report that, labelled, with the floor that applies
  # to it (`clip` is an estimator-side bound and never touches it).
  if (identical(model$method, "ensemble") &&
      (is.null(ps) || all(is.na(ps)))) {
    ps      <- model$res$nuisance$pi_oof
    ps_type <- "evaluation nuisance"
    clip    <- c(1e-4, 1 - 1e-4)
  }

  if (is.null(ps) || length(ps) != length(T_vec) || anyNA(ps)) {
    m <- model$method %||% "unknown"
    ps_methods <- c("grf", "bcf", "xlearner", "drlearner", "ensemble")
    reason <- if (m %in% ps_methods)
      sprintf("the propensity could not be read back from the '%s' fit", m)
    else
      sprintf("method '%s' estimates no propensity", m)
    rlang::abort(
      c(sprintf("No propensity-score overlap view: %s.", reason),
        i = "Methods that estimate one: grf, bcf, xlearner, drlearner; an ensemble reports the shared evaluation nuisance.",
        i = "The other views are reported as usual; `$overlap` is NULL."),
      class = "cateomics_overlap_no_propensity"
    )
  }

  ps <- as.numeric(ps)
  ps_t1 <- ps[T_vec == 1]; ps_t0 <- ps[T_vec == 0]
  overlap_region <- c(lower = max(min(ps_t0), min(ps_t1)),
                      upper = min(max(ps_t0), max(ps_t1)))
  in_support <- mean(ps >= overlap_region["lower"] &
                     ps <= overlap_region["upper"])

  # Count the scores sitting at the `clip` bounds.
  at_bound <- ps <= clip[1L] + 1e-10 | ps >= clip[2L] - 1e-10

  # No stored flag for "the two ranges do not meet": that is exactly
  # `support_fraction == 0`, and print says so from there.
  list(
    propensity_scores = ps,
    treatment         = T_vec,
    ps_type           = ps_type,
    clip              = clip,
    overlap_region    = overlap_region,
    support_fraction  = in_support,
    n_clipped         = sum(at_bound),
    clipped_pct       = mean(at_bound)
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
    stringsAsFactors = FALSE
  )
  df <- df[order(df$abs_smd, decreasing = TRUE), , drop = FALSE]
  n_total <- nrow(df)
  # The distribution is reported over every covariate, before the table is cut
  # to top_k. No sample is called balanced or not: at n << p an |SMD| of 0.1
  # is what randomisation alone produces for most covariates, and the
  # threshold is kept as a reference line on the plot rather than a verdict.
  quants <- stats::quantile(df$abs_smd, probs = c(0.5, 0.9), na.rm = TRUE)
  if (n_total > top_k) df <- df[seq_len(top_k), , drop = FALSE]
  rownames(df) <- NULL
  attr(df, "n_total")   <- n_total
  attr(df, "median")    <- unname(quants[1L])
  attr(df, "p90")       <- unname(quants[2L])
  attr(df, "max")       <- max(abs(smd), na.rm = TRUE)
  attr(df, "threshold") <- smd_threshold
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
      class = "cateomics_diagnostics_unknown_variables"
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

  # Standardise the ATE, then convert the standardised mean difference d to an
  # approximate risk ratio via RR ~= exp(0.91 * d): ln(OR) = 1.81 d (Chinn
  # 2000) and RR ~= sqrt(OR) for a common outcome (VanderWeele 2017), so
  # 1.81 / 2 = 0.91.
  #
  # The denominator is the pooled within-arm SD, which is the quantity the d of
  # that conversion is defined against. The SD of the observed outcome is not:
  # variance splits exactly into a within-arm part and a between-arm part, and
  # the second is made by the treatment effect itself, so using it puts the
  # effect into its own denominator and shrinks d by an amount that grows with
  # the effect. Each fallback below is named in `sd_basis`, so a caller can
  # tell which quantity the reported E-value was built on.
  y <- if (is.null(model$data$Y)) NULL else as.numeric(model$data$Y)
  w <- if (is.null(model$data$T)) NULL else as.integer(model$data$T)
  sd_y <- NA_real_
  sd_basis <- NA_character_
  if (!is.null(y) && !is.null(w) && length(w) == length(y)) {
    ok  <- is.finite(y) & !is.na(w)
    i1  <- ok & w == 1L
    i0  <- ok & w == 0L
    if (sum(i1) >= 2L && sum(i0) >= 2L) {
      n1 <- sum(i1); n0 <- sum(i0)
      pooled <- ((n1 - 1) * stats::var(y[i1]) + (n0 - 1) * stats::var(y[i0])) /
        (n1 + n0 - 2)
      if (is.finite(pooled) && pooled > 0) {
        sd_y     <- sqrt(pooled)
        sd_basis <- "pooled within-arm SD"
      }
    }
  }
  if (!is.finite(sd_y) || sd_y == 0) {
    # One arm too small, or no treatment vector: the whole-sample SD is the
    # closest thing left, and it is the same number when the arms do not differ.
    sd_y     <- if (!is.null(y)) stats::sd(y, na.rm = TRUE) else NA_real_
    sd_basis <- "SD of the observed outcome"
  }
  if (!is.finite(sd_y) || sd_y == 0) {
    # Fall back to the counterfactual SD only when raw Y is unavailable.
    sd_y     <- stats::sd(c(model$res$mu0_hat, model$res$mu1_hat), na.rm = TRUE)
    sd_basis <- "SD of the modelled potential outcomes"
    if (!is.finite(sd_y) || sd_y == 0) {
      sd_y     <- 1
      sd_basis <- "no usable spread; the effect is reported unstandardised"
    }
  }

  rr_point <- exp(0.91 * abs(ate) / sd_y)
  # The interval contributes the limit nearer the null. An interval that holds
  # the null gives 1; no interval at all gives NA, because "the evidence does
  # not survive the interval" and "there is no interval" are different things
  # and must not be the same number. `ci_basis` says which case this is.
  lo_ok <- !is.null(ci_lower) && is.finite(ci_lower)
  up_ok <- !is.null(ci_upper) && is.finite(ci_upper)
  if (lo_ok && ci_lower > 0) {
    rr_ci    <- exp(0.91 * ci_lower / sd_y)
    ci_basis <- "confidence limit nearer the null"
  } else if (up_ok && ci_upper < 0) {
    rr_ci    <- exp(0.91 * abs(ci_upper) / sd_y)
    ci_basis <- "confidence limit nearer the null"
  } else if (lo_ok && up_ok) {
    rr_ci    <- 1
    ci_basis <- "the interval crosses the null"
  } else {
    rr_ci    <- NA_real_
    ci_basis <- if (isFALSE(model$meta$estimator$bootstrap)) {
      "no interval on the average effect (none was requested: bootstrap_ci = FALSE)"
    } else {
      "no interval on the average effect"
    }
  }

  evalue_fn <- function(rr) {
    if (!is.finite(rr)) NA_real_ else if (rr <= 1) 1 else rr + sqrt(rr * (rr - 1))
  }
  ev_point <- evalue_fn(rr_point)
  ev_ci    <- evalue_fn(rr_ci)

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
    sd_basis       = sd_basis,
    ci_basis       = ci_basis
  )
}


# -------------------------------------------------------------------------
# print (object-lean, one line per element)
# -------------------------------------------------------------------------

#' @export
print.cateomics_diagnostics <- function(x, ...) {
  m <- x$metadata
  cat("<cateomics_diagnostics>: method = ", m$method %||% "<unknown>",
      ", n = ", m$n, ", p = ", m$p, "\n", sep = "")

  ov <- x$overlap
  if (is.null(ov)) {
    cat("  overlap : not available (this method estimates no propensity)\n")
  } else if (length(unique(ov$propensity_scores)) < 2L) {
    # One value for everyone: the support is that value and every row carries
    # it. Printing "100% in support" alone reads as a finding.
    cat(sprintf("  overlap : the propensity takes one value (%.3f); all %d rows carry it (PS: %s)\n",
                ov$propensity_scores[1L], length(ov$propensity_scores),
                ov$ps_type))
  } else {
    cat(sprintf("  overlap : common support = [%.3f, %.3f], %.0f%% in support, %d of %d at the [%.4g, %.4g] bounds (PS: %s)\n",
                ov$overlap_region["lower"], ov$overlap_region["upper"],
                100 * ov$support_fraction, ov$n_clipped, length(ov$propensity_scores),
                ov$clip[1L], ov$clip[2L], ov$ps_type))
    if (isTRUE(ov$support_fraction == 0))
      cat("            [WARNING] treated/control PS ranges do not overlap.\n")
  }

  sm <- x$smd
  if (is.null(sm) || !nrow(sm)) {
    cat("  balance : no table in this object\n")
  } else {
  # The quantiles are over every covariate; the table itself is cut to `top_k`.
  # Saying so when the whole table fits would describe a cut that never
  # happened.
  cut_note <- if (nrow(sm) < (attr(sm, "n_total") %||% nrow(sm)))
    sprintf(" (top-%d shown)", nrow(sm)) else ""
  cat(sprintf("  balance : |SMD| median %.3f, 90th pct %.3f, max %.3f over %d covariates%s\n",
              attr(sm, "median") %||% NA_real_, attr(sm, "p90") %||% NA_real_,
              attr(sm, "max") %||% NA_real_,
              attr(sm, "n_total") %||% NA_integer_, cut_note))
  }

  if (!is.null(x$variable_overlap)) {
    cat(sprintf("  variables: overlap view for %d requested covariate(s)\n",
                length(m$variables)))
  }

  ev <- x$evalue
  if (isTRUE(ev$available)) {
    # Each piece on its own: a field an older object does not carry would
    # otherwise take the whole line, or the whole print, with it.
    nm1 <- function(v) if (length(v) == 1L && is.finite(v)) v else NA_real_
    ci_txt <- if (isTRUE(is.finite(ev$evalue_ci))) sprintf("%.2f", ev$evalue_ci)
              else "not available"
    ate <- nm1(ev$ate); sdy <- nm1(ev$sd_y)
    cat(sprintf("  E-value : point = %.2f, CI = %s [%s] (ATE = %.3f, d = %.3f, %s = %.3f)\n",
                nm1(ev$evalue_point), ci_txt,
                ev$ci_basis %||% "unrecorded", ate,
                abs(ate) / sdy, ev$sd_basis %||% "unrecorded", sdy))
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
#' @param x A \code{cateomics_diagnostics} object (from \code{\link{cf_diagnostics}}).
#' @param type One of \code{"overlap"} (PS density, default),
#'   \code{"smd"} (Love dot plot) or \code{"variable"} (requested-variable
#'   distributions). The E-value is reported as a number by \code{print()}
#'   and in \code{$evalue}; it has no view of its own.
#' @param ... Not used; naming an unknown argument is an error.
#' @return A \code{ggplot} object, or \code{NULL} (invisibly) when the
#'   requested view is not in the object. A view whose panel would be empty is
#'   refused rather than drawn: \code{"overlap"} on a known propensity, which
#'   takes one value for everyone, aborts with \code{cateomics_plot_no_data}
#'   and points at \code{print()} for the support summary.
#'
#'   The object returned is an ordinary \pkg{ggplot2} plot: adding layers,
#'   themes or scales to it restyles the panel, and because the arm is mapped
#'   to the fill, a fill scale added afterwards replaces the one here.
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
                                         "variable"),
                                ...) {
  if (!inherits(x, "cateomics_diagnostics")) {
    rlang::abort("x must be a cateomics_diagnostics from cf_diagnostics().",
                 class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_diagnostics()")
  .cf_require_ggplot("cf_plot_diagnostics")
  type <- match.arg(type)
  p <- switch(type,
    overlap  = .cf_diag_plot_overlap(x$overlap),
    smd      = .cf_diag_plot_smd(x$smd),
    variable = .cf_diag_plot_variable(x$variable_overlap)
  )
  if (is.null(p)) return(invisible(NULL))
  p
}

#' @keywords internal
#' @noRd
.cf_diag_plot_overlap <- function(ov) {
  if (is.null(ov)) return(NULL)
  # A known propensity is one value for everyone, so there are no two
  # distributions to lay over each other and a density has nothing to
  # estimate: the panel would come out empty. The support summary is in
  # print() and in the object, so say that rather than draw nothing.
  if (length(unique(ov$propensity_scores)) < 2L) {
    rlang::abort(
      c("The propensity takes one value, so there is no overlap to draw.",
        i = "This is what a known propensity score looks like; print() reports the support.",
        i = "Estimate the propensity (leave `propensity_score` unset) for a distribution to compare."),
      class = "cateomics_plot_no_data")
  }
  df <- data.frame(ps = ov$propensity_scores,
                   group = ifelse(ov$treatment == 1, "Treated", "Control"))
  # The densities are kernel estimates, so each arm's curve runs past its own
  # data at both ends and both colours show outside the common support even
  # where one arm has no rows. The rug is the scores themselves.
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["ps"]], fill = .data[["group"]])) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_rug(ggplot2::aes(colour = .data[["group"]]), sides = "b",
                      length = ggplot2::unit(0.04, "npc"), alpha = 0.8) +
    ggplot2::geom_vline(xintercept = ov$overlap_region,
                        linetype = "dashed", alpha = 0.7) +
    .CF_ARM_FILL + .CF_ARM_COLOUR +
    ggplot2::labs(x = "Propensity score", y = "Density",
                  title = "Treatment overlap",
                  subtitle = sprintf("common support [%.3f, %.3f], %.0f%% in support; ticks = the scores themselves",
                                     ov$overlap_region["lower"],
                                     ov$overlap_region["upper"],
                                     100 * ov$support_fraction),
                  fill = NULL, colour = NULL) +
    .cf_theme()
}

#' @keywords internal
#' @noRd
.cf_diag_plot_smd <- function(sm) {
  thr <- attr(sm, "threshold")
  # Love-style dot plot: |SMD| descending, with the threshold drawn as a
  # reference line only. Colouring the points by it would have the figure pass
  # a verdict, which is what the threshold is deliberately not used for.
  df <- sm
  df$variable <- factor(df$variable, levels = rev(df$variable))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["abs_smd"]],
                                   y = .data[["variable"]])) +
    ggplot2::geom_vline(xintercept = thr, linetype = "dashed",
                        color = "grey40") +
    ggplot2::geom_point(colour = .CF_OKABE_ITO[1], size = 2.5) +
    ggplot2::labs(x = "|Standardized mean difference|", y = NULL,
                  title = "Covariate balance",
                  subtitle = sprintf("top %d of %d covariates (reference |SMD| = %.2f)",
                                     nrow(df), attr(sm, "n_total"), thr)) +
    .cf_theme()
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
    .CF_ARM_FILL +
    ggplot2::labs(x = "Value (preprocessed)", y = "Density",
                  title = "Requested-variable overlap", fill = NULL) +
    .cf_theme()
}

