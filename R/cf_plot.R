#' Plot helpers for CATEomics objects
#'
#' ggplot2-based plot helpers for the analysis objects
#' (\code{cf_model}, \code{cateomics_compare}, \code{cateomics_evaluate},
#' \code{cateomics_diagnostics}, \code{cateomics_rule},
#' \code{cateomics_rule_grid}, \code{cateomics_interpret},
#' \code{cateomics_model_stability}). Each helper extracts a small data.frame
#' and builds the plot in a frame that holds only that, so the returned
#' ggplot does not capture the source object.
#'
#' @return This is an overview page for the \code{cf_plot_*} family; it does not
#'   define a function and has no return value. Each individual \code{cf_plot_*}
#'   helper returns a \code{ggplot} object.
#' @name cf_plot
NULL

# Okabe-Ito colour-blind-safe qualitative palette.
.CF_OKABE_ITO <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                   "#E69F00", "#56B4E9", "#F0E442", "#999999")

# The views that colour by feature value carry shapviz's upright colour bar,
# which the family's legend position would lay on its side and leave no room
# for the panel. Those keep the bar where shapviz puts it.
.CF_SHAP_COLOUR_BAR <- ggplot2::theme(legend.position = "right")

#' The two treatment arms, wherever they are drawn together.
#'
#' Okabe-Ito green and purple: both are in the colour-blind-safe palette the
#' rest of the family uses, and where the two densities overlap the blend is
#' neutral instead of taking one arm's colour or turning muddy.
#' @keywords internal
#' @noRd
.CF_ARM_COLOURS <- c(Control = "#009E73", Treated = "#CC79A7")

#' @keywords internal
#' @noRd
.CF_ARM_FILL <- ggplot2::scale_fill_manual(values = .CF_ARM_COLOURS)

#' @keywords internal
#' @noRd
.CF_ARM_COLOUR <- ggplot2::scale_colour_manual(values = .CF_ARM_COLOURS)

#' Shared minimal theme.
#' @keywords internal
#' @noRd
.cf_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

#' @keywords internal
#' @noRd
.cf_require_ggplot <- function(fn) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    rlang::abort(
      sprintf("%s requires ggplot2 (Suggests). install.packages('ggplot2').", fn),
      class = "cateomics_optional_dependency_missing"
    )
  }
}


# =========================================================================
# cf_plot_compare: ranked bar of the selection metric.
# =========================================================================

#' Ranked bar chart of a cf_compare scoreboard metric.
#'
#' @param cmp A \code{cateomics_compare} object.
#' @param metric Which scoreboard metric to rank by (default
#'   \code{"r_loss"}): \code{"r_loss"}, \code{"dr_loss"},
#'   \code{"factual_loss"}, \code{"r_squared_score"},
#'   \code{"dr_squared_score"} or \code{"factual_squared_score"}. One this
#'   comparison did not score falls back to one it did; any other name is an
#'   error. The descriptive columns of the scoreboard are not metrics.
#' @param ... Not used; naming an unknown argument is an error.
#' @return A ggplot object (methods ranked best-to-worst by the metric).
#'
#'   The object returned is an ordinary \pkg{ggplot2} plot: adding layers,
#'   themes or scales to it restyles the panel. The single colours are set as
#'   constants rather than mapped, so a colour scale added afterwards does not
#'   reach them.
#' @examples
#' # A one-method comparison keeps the example fast (no Python, no MCMC).
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = "grf", n_features = 40)
#' cf_plot_compare(cmp, metric = "r_loss")
#' @export
cf_plot_compare <- function(cmp, metric = "r_loss", ...) {
  if (!inherits(cmp, "cateomics_compare")) {
    rlang::abort("cmp must be a cateomics_compare.", class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_compare()")
  .cf_require_ggplot("cf_plot_compare")
  sb <- cmp$scoreboard
  # A name that is not one of the six metrics is an error; a metric this
  # comparison did not score falls back to one it did, so the default works
  # whatever cf_compare ran.
  higher_better <- .cf_metric_higher_better(metric)
  if (!metric %in% colnames(sb)) {
    metric <- intersect(names(.cf_metric_directions), colnames(sb))[1L]
    if (is.na(metric)) {
      rlang::abort(
        c("No usable metric column in scoreboard.",
          i = "Pass a metric column of cmp$scoreboard as `metric` (e.g. 'r_loss', 'dr_loss', 'factual_loss')."),
        class = "cateomics_plot_no_metric")
    }
    higher_better <- .cf_metric_higher_better(metric)
  }
  # Bar label = the scoreboard identity: `label` (method_base) in grid mode,
  # else the bare `method`.
  id_vec <- if ("label" %in% colnames(sb)) sb$label else sb$method
  df <- data.frame(method = as.character(id_vec),
                   value  = as.numeric(sb[[metric]]),
                   stringsAsFactors = FALSE)
  df <- df[is.finite(df$value), , drop = FALSE]
  .cf_build_scoreboard(df, metric, higher_better)
}

#' @keywords internal
#' @noRd
.cf_build_scoreboard <- function(df, metric, higher_better) {
  df <- df[order(df$value, decreasing = higher_better), , drop = FALSE]
  # Best method at the top of the horizontal bar chart.
  df$method <- factor(df$method, levels = rev(df$method))
  better <- if (higher_better) "higher = better" else "lower = better"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["value"]],
                                   y = .data[["method"]])) +
    ggplot2::geom_col(fill = .CF_OKABE_ITO[1], width = 0.7) +
    ggplot2::labs(x = metric, y = NULL,
                  title = sprintf("Method ranking by %s", metric),
                  subtitle = better) +
    .cf_theme()
}


# =========================================================================
# cf_plot_evaluate: unified evaluate views by type=.
# rate / gates / blp all read the cateomics_evaluate object.
# =========================================================================

#' Plot a cf_evaluate result: TOC curve, GATES bars, or BLP calibration.
#'
#' @param eval A \code{cateomics_evaluate} object (from \code{cf_evaluate}).
#' @param type Which view: \code{"rate"} (the TOC curve over q, default),
#'   \code{"gates"} (group ATE bars with an overall-ATE line), or
#'   \code{"blp"} (best-linear-predictor calibration: the DR pseudo-outcome
#'   against the predicted effect, one point per sample, with the fitted line
#'   and a slope-1 reference).
#' @param ... Not used; naming an unknown argument is an error.
#' @return A ggplot object.
#'
#'   The object returned is an ordinary \pkg{ggplot2} plot: adding layers,
#'   themes or scales to it restyles the panel. The single colours are set as
#'   constants rather than mapped, so a colour scale added afterwards does not
#'   reach them.
#' @examples
#' # Evaluate a fast R-native causal forest, then plot the TOC curve.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' ev <- cf_evaluate(fit)
#' cf_plot_evaluate(ev, type = "rate")
#' @export
cf_plot_evaluate <- function(eval, type = c("rate", "gates", "blp"), ...) {
  if (!inherits(eval, "cateomics_evaluate")) {
    rlang::abort("eval must be a cateomics_evaluate.", class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_evaluate()")
  .cf_require_ggplot("cf_plot_evaluate")
  type <- match.arg(type)
  switch(type,
         rate  = .cf_plot_eval_rate(eval, ...),
         gates = .cf_plot_eval_gates(eval, ...),
         blp   = .cf_plot_eval_blp(eval, ...))
}

#' BLP calibration view (internal; cf_plot_evaluate(type = "blp")).
#' @keywords internal
#' @noRd
.cf_plot_eval_blp <- function(eval, ...) {
  if (!inherits(eval, "cateomics_evaluate")) {
    rlang::abort("eval must be a cateomics_evaluate.", class = "cateomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  pd <- eval$blp$points
  if (is.null(pd)) {
    rlang::abort(
      c("The BLP calibration has no data to plot.",
        i = "Run cf_evaluate() with 'blp' in `metrics` (the default)."),
      class = "cateomics_plot_no_data")
  }
  df <- data.frame(
    predicted = as.numeric(pd$predicted),
    observed  = as.numeric(pd$observed),
    stringsAsFactors = FALSE
  )
  sl <- eval$blp$slope
  # The regressor is centred, so the stored intercept is the fitted value at
  # the mean tau-hat; shifting it back puts the fitted line on the axes drawn.
  intercept <- eval$blp$intercept - sl$estimate * mean(df$predicted, na.rm = TRUE)
  .cf_build_blp(df, sl$estimate, intercept, sl$ci_lower, sl$ci_upper)
}

#' @keywords internal
#' @noRd
.cf_build_blp <- function(df, slope, intercept, ci_lower = NA_real_,
                          ci_upper = NA_real_) {
  # One point per sample against the pseudo-outcome the slope was fitted on.
  # Binning the points by tau-hat quantile would draw the group effects
  # `cf_plot_evaluate(type = "gates")` already draws, from the same bins. The
  # two axes share a unit, but the pseudo-outcome spreads several times wider
  # than the estimates it is regressed on, so a fixed aspect ratio would make
  # the panel too tall to read; the reference line is named in the subtitle
  # instead of being left to look like a diagonal.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["predicted"]],
                                        y = .data[["observed"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50")
  # The lines whose slope lies inside the reported interval, pivoting through
  # the fitted centroid. It is not a confidence band for the fitted line: the
  # intercept's own uncertainty is not stored, so this is the part of the
  # band the slope accounts for, and the subtitle says so.
  sub <- sprintf("solid = BLP fit (slope = %.2f; 1 = calibrated), dashed = slope 1",
                 slope)
  if (is.finite(ci_lower) && is.finite(ci_upper)) {
    xbar <- mean(df$predicted, na.rm = TRUE)
    ybar <- slope * xbar + intercept
    xs <- seq(min(df$predicted, na.rm = TRUE), max(df$predicted, na.rm = TRUE),
              length.out = 100L)
    band <- data.frame(
      x    = xs,
      ymin = ybar + pmin(ci_lower * (xs - xbar), ci_upper * (xs - xbar)),
      ymax = ybar + pmax(ci_lower * (xs - xbar), ci_upper * (xs - xbar)))
    p <- p + ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["ymin"]],
                   ymax = .data[["ymax"]]),
      inherit.aes = FALSE, fill = .CF_OKABE_ITO[2], alpha = 0.18)
    # Wrapped, not shortened: at the default figure width of seven inches a
    # single line loses its tail, and the clause that goes is the one naming
    # the shading this view added.
    sub <- sprintf(paste0("solid = BLP fit (slope = %.2f [%.2f, %.2f]; ",
                          "1 = calibrated), dashed = slope 1\n",
                          "shaded = slopes inside the interval"),
                   slope, ci_lower, ci_upper)
  }
  p +
    ggplot2::geom_point(colour = .CF_OKABE_ITO[1], size = 1.6, alpha = 0.7) +
    ggplot2::geom_abline(slope = slope, intercept = intercept,
                         colour = .CF_OKABE_ITO[2], linewidth = 0.8) +
    ggplot2::labs(x = "Predicted CATE (out-of-fold)",
                  y = "DR pseudo-outcome",
                  title = "Calibration of predicted treatment effect",
                  subtitle = sub) +
    .cf_theme()
}


# =========================================================================
# .cf_plot_eval_rate: the TOC(q) curve (via cf_plot_evaluate).
# =========================================================================

#' TOC(q) view (internal; cf_plot_evaluate(type = "rate")).
#' @keywords internal
#' @noRd
.cf_plot_eval_rate <- function(eval, ...) {
  if (!inherits(eval, "cateomics_evaluate")) {
    rlang::abort("eval must be a cateomics_evaluate.", class = "cateomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  if (is.null(eval$rate$curve)) {
    rlang::abort(
      c("The evaluation has no TOC curve to plot.",
        i = "Run cf_evaluate(..., metrics = c('rate', ...)) first; 'rate' is included by default."),
      class = "cateomics_plot_no_data")
  }
  rc <- eval$rate$curve
  df <- data.frame(quantile = as.numeric(rc$quantile),
                   toc      = as.numeric(rc$toc),
                   lower    = as.numeric(rc$ci_lower),
                   upper    = as.numeric(rc$ci_upper),
                   stringsAsFactors = FALSE)
  .cf_build_rate(df)
}

#' @keywords internal
#' @noRd
.cf_build_rate <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["quantile"]],
                                   y = .data[["toc"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      alpha = 0.2, fill = .CF_OKABE_ITO[1]) +
    ggplot2::geom_line(colour = .CF_OKABE_ITO[1], linewidth = 0.8) +
    ggplot2::labs(x = "Treated fraction q (top by predicted effect)",
                  y = "TOC(q) = gain over the overall ATE",
                  title = "Targeting operator characteristic (TOC)") +
    .cf_theme()
}


# =========================================================================
# cf_plot_policy_curve: PAPE baseline + treated-fraction axis.
# =========================================================================

#' Plot the policy-value curve from cf_rule(grid=) with a same-budget random
#' baseline (PAPE).
#'
#' @param rule_grid A \code{cateomics_rule_grid} object.
#' @param ... Not used; naming an unknown argument is an error.
#' @return A ggplot object. The x-axis is the treated fraction; a dashed line
#'   marks the same-budget random-assignment policy value (the PAPE baseline).
#'   Adding layers, themes or scales to it restyles the panel. The two lines
#'   and the band are coloured by constants rather than by a mapping, so a
#'   colour scale added afterwards does not reach them.
#' @examples
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' rg <- cf_rule(fit, type = "threshold", grid = TRUE)
#' cf_plot_policy_curve(rg)
#' @export
cf_plot_policy_curve <- function(rule_grid, ...) {
  if (!inherits(rule_grid, "cateomics_rule_grid")) {
    rlang::abort("rule_grid must be a cateomics_rule_grid.",
                 class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_policy_curve()")
  .cf_require_ggplot("cf_plot_policy_curve")
  rd <- rule_grid$results
  n_total <- rule_grid$n %||% max(rd$n_treated, na.rm = TRUE)
  df <- data.frame(
    treated_frac = as.numeric(rd$n_treated) / n_total,
    policy_value = as.numeric(rd$policy_value),
    lower        = as.numeric(rd$ci_lower),
    upper        = as.numeric(rd$ci_upper),
    stringsAsFactors = FALSE
  )
  # Same-budget random baseline (PAPE): treating a random fraction f yields
  # V_none + f * (V_all - V_none).
  df$random <- rule_grid$treat_none %||% NA_real_
  if (is.finite(rule_grid$treat_all %||% NA_real_) &&
      is.finite(rule_grid$treat_none %||% NA_real_)) {
    df$random <- rule_grid$treat_none +
      df$treated_frac * (rule_grid$treat_all - rule_grid$treat_none)
  }
  .cf_build_policy_curve(df)
}

#' @keywords internal
#' @noRd
.cf_build_policy_curve <- function(df) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["treated_frac"]],
                                        y = .data[["policy_value"]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      alpha = 0.2, fill = .CF_OKABE_ITO[1]) +
    ggplot2::geom_line(colour = .CF_OKABE_ITO[1], linewidth = 0.8)
  if (any(is.finite(df$random))) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(y = .data[["random"]]), linetype = "dashed",
      colour = .CF_OKABE_ITO[2], na.rm = TRUE)
  }
  p +
    ggplot2::labs(x = "Treated fraction",
                  y = "Policy value (DR)",
                  title = "Policy value vs treated fraction",
                  subtitle = "dashed = same-budget random assignment (PAPE baseline)") +
    .cf_theme()
}


# =========================================================================
# .cf_plot_eval_gates: group ATE bars + overall ATE line.
# =========================================================================

#' GATES group-ATE view (internal; cf_plot_evaluate(type = "gates")).
#' @keywords internal
#' @noRd
.cf_plot_eval_gates <- function(eval, ...) {
  if (!inherits(eval, "cateomics_evaluate")) {
    rlang::abort("eval must be a cateomics_evaluate.", class = "cateomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  if (is.null(eval$gates) || isTRUE(eval$gates$degenerate) ||
      is.null(eval$gates$group_ate)) {
    rlang::abort(
      c("The GATES plot has no group ATEs.",
        i = "Run cf_evaluate() with 'gates' in `metrics` (the default); near-constant estimated effects leave no distinct benefit groups."),
      class = "cateomics_plot_no_data")
  }
  ga <- eval$gates$group_ate
  df <- data.frame(group = as.integer(ga$group),
                   estimate = as.numeric(ga$estimate),
                   lower = as.numeric(ga$ci_lower),
                   upper = as.numeric(ga$ci_upper),
                   stringsAsFactors = FALSE)
  # Overall-ATE reference line: with equal-count quantile bins the
  # unweighted mean of the per-group DR-ATEs equals the overall DR-ATE, and
  # ties that collapse a break make it an approximation.
  overall <- mean(df$estimate, na.rm = TRUE)
  bw_p <- eval$gates$top_bottom_test$p_value %||% NA_real_
  .cf_build_gates(df, overall, bw_p)
}

#' @keywords internal
#' @noRd
.cf_build_gates <- function(df, overall, bw_p) {
  df$group <- factor(df$group)
  sub <- if (is.finite(bw_p))
    sprintf("dashed = overall ATE; top - bottom: p = %.3g", bw_p)
  else "dashed = overall ATE"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["group"]],
                                   y = .data[["estimate"]])) +
    ggplot2::geom_col(fill = .CF_OKABE_ITO[1], width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      width = 0.2, colour = "grey30") +
    ggplot2::geom_hline(yintercept = overall, linetype = "dashed",
                        colour = .CF_OKABE_ITO[2]) +
    ggplot2::labs(x = "Predicted-effect quantile group (low to high)",
                  y = "Group ATE (DR)",
                  title = "Group ATEs across predicted-effect quantiles",
                  subtitle = sub) +
    .cf_theme()
}


# =========================================================================
# cf_plot_interpret: unified SHAP views by type=, via shapviz.
# =========================================================================

#' Plot a cf_interpret result: SHAP importance bar, beeswarm, or dependence.
#'
#' Coerces the \code{cf_interpret} result to a \pkg{shapviz} object
#' internally and delegates to \pkg{shapviz}. \code{type = "bar"} /
#' \code{"beeswarm"} call \code{shapviz::sv_importance} (the mean|SHAP| ranking
#' or the per-observation distribution); \code{type = "dependence"} calls
#' \code{shapviz::sv_dependence} for a single \code{feature}. Requires the
#' \pkg{shapviz} package.
#'
#' @param itp A \code{cateomics_interpret} object (from \code{cf_interpret}).
#' @param type One of \code{"bar"} (mean|SHAP| importance, default),
#'   \code{"beeswarm"} (per-observation SHAP distribution), or
#'   \code{"dependence"} (SHAP vs feature value for one \code{feature}).
#' @param feature Character; the feature for \code{type = "dependence"}
#'   (required for that view, ignored otherwise).
#' @param top_n Integer; maximum features for the bar / beeswarm view
#'   (default \code{20}).
#' @param surrogate Which surrogate to plot ("ridge" / "xgb"); \code{NULL}
#'   (default) picks the first available. For the linear surrogate
#'   (\code{"ridge"}) the dependence view is a straight line by construction:
#'   its SHAP value for a feature is that feature's coefficient times its
#'   standardised value, so the relation it draws has no curvature to show.
#'   The points are coloured by a second feature that \pkg{shapviz} picks
#'   with its own interaction heuristic, not one that was asked for, and with
#'   no interaction in the data the candidates score almost alike.
#' @param ... Passed to the underlying \pkg{shapviz} function.
#' @return A ggplot object, titled and themed like the other
#'   \code{cf_plot_*} views; adding \code{ggplot2::labs(title = ...)} or a
#'   \pkg{ggplot2} theme to it replaces either. Adding layers or scales
#'   restyles the panel, but the single colours are set as constants rather
#'   than mapped, so a colour scale added afterwards does not reach them.
#' @examples
#' # Surrogate SHAP importances from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' itp <- cf_interpret(fit, surrogate = "xgb")
#' if (requireNamespace("shapviz", quietly = TRUE))
#'   cf_plot_interpret(itp, type = "bar")
#' @export
cf_plot_interpret <- function(itp, type = c("bar", "beeswarm", "dependence"),
                              feature = NULL, top_n = 20L, surrogate = NULL,
                              ...) {
  if (!inherits(itp, "cateomics_interpret")) {
    rlang::abort("itp must be a cateomics_interpret (from cf_interpret).",
                 class = "cateomics_invalid_input")
  }
  type <- match.arg(type)
  sv <- as_shapviz(itp, surrogate = surrogate)
  if (identical(type, "dependence")) {
    if (is.null(feature) || !is.character(feature) || length(feature) != 1L) {
      rlang::abort(
        c("type = 'dependence' requires a single `feature` name.",
          i = "Pass exactly one covariate name to `feature`, e.g. feature = 'GENE1'."),
        class = "cateomics_invalid_input")
    }
    return(shapviz::sv_dependence(sv, v = feature, ...) +
             ggplot2::labs(title = "SHAP dependence") +
             .cf_theme() + .CF_SHAP_COLOUR_BAR)
  }
  # shapviz titles none of its plots and draws them on its own theme, so the
  # title and the look the rest of the family carries are added here. The
  # SHAP values are the surrogate model's, which the title says so that the
  # panel is not read as the estimator's own. The bar is the one mark whose
  # colour carries nothing, so it takes the family's first colour; the
  # beeswarm and the dependence plot colour by feature value and are left as
  # shapviz draws them.
  args <- list(sv, kind = type, max_display = as.integer(top_n)[1L], ...)
  if (identical(type, "bar") && is.null(args$fill)) {
    args$fill <- .CF_OKABE_ITO[1]
  }
  p <- do.call(shapviz::sv_importance, args) +
    ggplot2::labs(title = "SHAP importance (surrogate model)") +
    .cf_theme()
  if (!identical(type, "bar")) p <- p + .CF_SHAP_COLOUR_BAR
  p
}


# =========================================================================
# cf_plot_stability: correlation heatmap only, with cor_method.
# =========================================================================

#' Correlation heatmap of a cf_model_stability object.
#'
#' @param stability A \code{cateomics_model_stability} object.
#' @param cor_method Which stored correlation matrix to show; defaults to the
#'   object's primary correlation method.
#' @param ... Not used; naming an unknown argument is an error.
#' @return A ggplot object (K x K pairwise correlation of OOF tau-hat).
#'
#'   The object returned is an ordinary \pkg{ggplot2} plot: adding layers,
#'   themes or scales to it restyles the panel, and because the fill is mapped,
#'   a fill scale added afterwards replaces the one here.
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf) so the pairwise
#' # correlation is non-degenerate; known randomised propensity.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   # 1 chain keeps bcf fast
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)),
#'                   propensity_score = 0.5, eval_propensity_score = 0.5,
#'                   random_state = 1)
#' stab <- cf_model_stability(cmp)
#' cf_plot_stability(stab)
#' @export
cf_plot_stability <- function(stability, cor_method = NULL, ...) {
  if (!inherits(stability, "cateomics_model_stability")) {
    rlang::abort("stability must be a cateomics_model_stability object.",
                 class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_stability()")
  .cf_require_ggplot("cf_plot_stability")
  if (is.null(cor_method)) cor_method <- stability$meta$primary_cor
  mat <- stability$cor_matrices[[cor_method]]
  if (is.null(mat) || nrow(mat) == 0L) {
    rlang::abort(
      c(sprintf("No '%s' correlation matrix is stored on this stability object.", cor_method),
        i = sprintf("Available: %s. Pass one as `cor_method`.",
                    paste(names(stability$cor_matrices), collapse = ", "))),
      class = "cateomics_plot_stability_no_data")
  }
  rn <- rownames(mat); K <- nrow(mat)
  df <- data.frame(row = rep(rn, times = K),
                   col = rep(colnames(mat), each = K),
                   val = as.numeric(mat),
                   stringsAsFactors = FALSE)
  .cf_build_stability(df, rn, colnames(mat), cor_method)
}

#' @keywords internal
#' @noRd
.cf_build_stability <- function(df, rn, cn, cor_method) {
  # A discrete y axis puts its first level at the bottom, so the levels are
  # reversed to put the diagonal where a correlation matrix is read from:
  # top left to bottom right.
  df$row <- factor(df$row, levels = rev(rn))
  df$col <- factor(df$col, levels = cn)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["col"]], y = .data[["row"]],
                                   fill = .data[["val"]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = format(.data[["val"]], digits = 2)),
                       size = 3) +
    ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white",
                                  high = "#2166AC", midpoint = 0,
                                  limits = c(-1, 1), name = "correlation") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Pairwise tau correlation across methods",
                  subtitle = sprintf("%s correlation of out-of-fold tau-hat",
                                     tools::toTitleCase(cor_method))) +
    .cf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}


# =========================================================================
# cf_plot_rule (policy tree): renderer for cf_rule(type = "tree").
# =========================================================================

#' Plot the policy tree from cf_rule(type = "tree").
#'
#' @param rule A \code{cateomics_rule} produced by \code{cf_rule(type = "tree")}.
#' @param ... Not used; naming an unknown argument is an error.
#' @return The \code{DiagrammeR} htmlwidget for the tree, with the leaves
#'   labelled by what they do and a caption naming the scale the splits are
#'   on (the preprocessed covariates). When \code{DiagrammeR} is not
#'   installed, the text summary is printed and \code{NULL} is returned
#'   invisibly.
#' @examples
#' # Draw the policy tree from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' rule <- cf_rule(fit, type = "tree")
#' cf_plot_rule(rule)
#' @export
cf_plot_rule <- function(rule, ...) {
  if (!inherits(rule, "cateomics_rule")) {
    rlang::abort("rule must be a cateomics_rule produced by cf_rule().",
                 class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_rule()")
  if (!identical(rule$type, "tree") || is.null(rule$tree)) {
    rlang::abort(
      "cf_plot_rule() requires a cf_rule(type = 'tree') with a fitted tree.",
      class = "cateomics_invalid_input")
  }
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    # policytree's plot method builds a DiagrammeR widget and returns it, so
    # the widget has to be handed back for anything to be rendered outside an
    # auto-printing console. Its leaves are numbered actions and it carries no
    # caption, so the graph alone does not say which action treats, nor that
    # the splits are on the preprocessed scale. Both are added to the DOT the
    # widget holds rather than redrawn.
    w <- plot(rule$tree)
    dot <- w$x$diagram
    actions <- rule$tree$action.names %||% c("control", "treat")
    for (k in seq_along(actions)) {
      dot <- gsub(sprintf("action = %d", k),
                  if (identical(actions[k], "treat")) "treat" else "do not treat",
                  dot, fixed = TRUE)
    }
    dot <- sub("{", paste0("{\n labelloc=\"t\"; fontsize=11; ",
                           "label=\"Policy tree: who to treat",
                           "\\n(splits are on the preprocessed covariate scale)\";"),
               dot, fixed = TRUE)
    w$x$diagram <- dot
    w
  } else {
    rlang::inform(
      c("DiagrammeR not installed; printing the policy_tree text summary.",
        i = "install.packages('DiagrammeR') for the graph."),
      class = "cateomics_policy_tree_text_fallback")
    print(rule$tree)
    invisible(NULL)
  }
}


# =========================================================================
# cf_plot_cate: distribution of the estimated CATE for a fitted model.
# =========================================================================

#' Plot the distribution of estimated CATE (tau-hat) for a fitted model.
#'
#' Histogram of the tau-hat the fit carries. An ensemble has no in-sample
#' estimate, so its out-of-fold one is shown.
#'
#' @param fit A \code{cf_model} object from \code{\link{cf_fit}}.
#' @param bins Integer; number of histogram bins (default \code{30}).
#' @param ... Not used; naming an unknown argument is an error.
#' @return A ggplot object (histogram of predicted CATE), or \code{NULL}
#'   invisibly if no per-sample effects are available.
#'
#'   The object returned is an ordinary \pkg{ggplot2} plot: adding layers,
#'   themes or scales to it restyles the panel. The single colours are set as
#'   constants rather than mapped, so a colour scale added afterwards does not
#'   reach them.
#' @examples
#' # Distribution of estimated CATE from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf", n_features = 40)
#' cf_plot_cate(fit)
#' @export
cf_plot_cate <- function(fit, bins = 30L, ...) {
  if (!inherits(fit, "cf_model")) {
    rlang::abort("fit must be a cf_model from cf_fit().",
                 class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_plot_cate()")
  .cf_require_ggplot("cf_plot_cate")
  tau <- tryCatch(predict(fit, type = "tau"), error = function(e) NULL)
  if (is.null(tau) || all(is.na(tau))) {
    rlang::inform("No per-sample CATE available; nothing to plot.")
    return(invisible(NULL))
  }
  # Build in a frame that holds only the data frame and the bin count: a
  # ggplot keeps a reference to the environment its aesthetics were made in,
  # so building here would carry the whole fit along with the plot.
  .cf_build_cate(data.frame(cate = as.numeric(tau)), as.integer(bins)[1L],
                 ate = fit$res$ate %||% NA_real_)
}

#' @keywords internal
#' @noRd
.cf_build_cate <- function(df, bins, ate = NA_real_) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["cate"]])) +
    ggplot2::geom_histogram(bins = bins,
                            fill = .CF_OKABE_ITO[1], colour = "white")
  # The average effect, so the spread can be read against something: a
  # near-constant tau-hat and a wide one otherwise look alike.
  sub <- NULL
  if (is.finite(ate)) {
    p <- p + ggplot2::geom_vline(xintercept = ate, linetype = "dashed",
                                 colour = "grey40")
    sub <- sprintf("dashed = average effect (%.3g)", ate)
  }
  p +
    ggplot2::labs(x = "Predicted CATE (tau-hat)", y = "Count",
                  title = "Distribution of estimated treatment effects",
                  subtitle = sub) +
    .cf_theme()
}
