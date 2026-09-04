#' Plot helpers for cfomics objects
#'
#' ggplot2-based plot helpers for the analysis objects
#' (\code{cfomics_compare}, \code{cfomics_evaluate}, \code{cfomics_rule_grid},
#' \code{cfomics_interpret}, \code{cfomics_model_stability}). Each helper
#' extracts a small data.frame and builds the plot in a clean environment so
#' the returned ggplot never captures the (potentially large) source object.
#'
#' @return This is an overview page for the \code{cf_plot_*} family; it does not
#'   define a function and has no return value. Each individual \code{cf_plot_*}
#'   helper returns a \code{ggplot} object.
#' @name cf_plot
NULL

# Okabe-Ito colour-blind-safe qualitative palette.
.CF_OKABE_ITO <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                   "#E69F00", "#56B4E9", "#F0E442", "#999999")

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
      class = "cfomics_optional_dependency_missing"
    )
  }
}


# =========================================================================
# cf_plot_compare: ranked bar of the selection metric.
# =========================================================================

#' Ranked bar chart of a cf_compare scoreboard metric.
#'
#' @param cmp A \code{cfomics_compare} object.
#' @param metric Which scoreboard column to rank by (default \code{"r_loss"}).
#' @param ... Ignored.
#' @return A ggplot object (methods ranked best-to-worst by the metric).
#' @examples
#' # A one-method comparison keeps the example fast (no Python, no MCMC).
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = "grf", n_features = 40)
#' cf_plot_compare(cmp, metric = "r_loss")
#' @export
cf_plot_compare <- function(cmp, metric = "r_loss", ...) {
  if (!inherits(cmp, "cfomics_compare")) {
    rlang::abort("cmp must be a cfomics_compare.", class = "cfomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_compare")
  sb <- cmp$scoreboard
  if (!metric %in% colnames(sb)) {
    metric <- intersect(c("r_loss", "dr_loss", "factual_loss"), colnames(sb))[1L]
    if (is.na(metric)) {
      rlang::abort(
        c("No usable metric column in scoreboard.",
          i = "Pass a loss column of cmp$scoreboard as `metric` (e.g. 'r_loss', 'dr_loss', 'factual_loss')."),
        class = "cfomics_plot_no_metric")
    }
  }
  higher_better <- .cf_metric_higher_better(metric)
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
# rate / gates / blp all read the cfomics_evaluate object.
# =========================================================================

#' Plot a cf_evaluate result: RATE curve, GATES bars, or BLP calibration.
#'
#' @param eval A \code{cfomics_evaluate} object (from \code{cf_evaluate}).
#' @param type Which view: \code{"rate"} (TOC(q) targeting-operator curve,
#'   default), \code{"gates"} (group ATE bars with an overall-ATE line), or
#'   \code{"blp"} (best-linear-predictor calibration of predicted vs observed).
#' @param ... Ignored.
#' @return A ggplot object.
#' @examples
#' # Evaluate a fast R-native causal forest, then plot the RATE curve.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' ev <- cf_evaluate(fit, eval_nuisance = "rf", eval_propensity_score = 0.5)
#' cf_plot_evaluate(ev, type = "rate")
#' @export
cf_plot_evaluate <- function(eval, type = c("rate", "gates", "blp"), ...) {
  if (!inherits(eval, "cfomics_evaluate")) {
    rlang::abort("eval must be a cfomics_evaluate.", class = "cfomics_invalid_input")
  }
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
  if (!inherits(eval, "cfomics_evaluate")) {
    rlang::abort("eval must be a cfomics_evaluate.", class = "cfomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  pd <- eval$blp$plot_data
  if (is.null(pd)) {
    rlang::abort(
      c("The BLP calibration has no data to plot.",
        i = "Run cf_evaluate() with 'blp' in `metrics` (the default); near-constant estimated effects collapse the bins and produce none."),
      class = "cfomics_plot_no_data")
  }
  df <- data.frame(
    predicted = as.numeric(pd$predicted_mean),
    observed  = as.numeric(pd$observed_mean),
    lower     = as.numeric(pd$ci_lower),
    upper     = as.numeric(pd$ci_upper),
    stringsAsFactors = FALSE
  )
  slope <- eval$blp$slope$estimate
  .cf_build_blp(df, slope)
}

#' @keywords internal
#' @noRd
.cf_build_blp <- function(df, slope) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["predicted"]],
                                   y = .data[["observed"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      width = 0, colour = "grey60") +
    ggplot2::geom_point(colour = .CF_OKABE_ITO[1], size = 2.5) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Predicted CATE (bin mean)",
                  y = "Observed effect (DR, bin mean)",
                  title = "Calibration of predicted treatment effect",
                  subtitle = sprintf("BLP slope = %.2f (1 = calibrated)", slope)) +
    .cf_theme()
}


# =========================================================================
# .cf_plot_eval_rate: TOC(q) curve (via cf_plot_evaluate).
# =========================================================================

#' RATE / TOC(q) view (internal; cf_plot_evaluate(type = "rate")).
#' @keywords internal
#' @noRd
.cf_plot_eval_rate <- function(eval, ...) {
  if (!inherits(eval, "cfomics_evaluate")) {
    rlang::abort("eval must be a cfomics_evaluate.", class = "cfomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  if (is.null(eval$rate_curve)) {
    rlang::abort(
      c("The evaluation has no RATE curve to plot.",
        i = "Run cf_evaluate(..., metrics = c('rate', ...)) first; 'rate' is included by default."),
      class = "cfomics_plot_no_data")
  }
  rc <- eval$rate_curve
  df <- data.frame(quantile = as.numeric(rc$quantile),
                   rate     = as.numeric(rc$rate),
                   lower    = as.numeric(rc$ci_lower),
                   upper    = as.numeric(rc$ci_upper),
                   stringsAsFactors = FALSE)
  .cf_build_rate(df)
}

#' @keywords internal
#' @noRd
.cf_build_rate <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["quantile"]],
                                   y = .data[["rate"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      alpha = 0.2, fill = .CF_OKABE_ITO[1]) +
    ggplot2::geom_line(colour = .CF_OKABE_ITO[1], linewidth = 0.8) +
    ggplot2::labs(x = "Treated fraction q (top by predicted effect)",
                  y = "TOC(q): gain over the overall ATE",
                  title = "Targeting operator characteristic (RATE)") +
    .cf_theme()
}


# =========================================================================
# cf_plot_policy_curve: PAPE baseline + treated-fraction axis.
# =========================================================================

#' Plot the policy-value curve from cf_rule(grid=) with a same-budget random
#' baseline (PAPE).
#'
#' @param rule_grid A \code{cfomics_rule_grid} object.
#' @param ... Ignored.
#' @return A ggplot object. The x-axis is the treated fraction; a dashed line
#'   marks the same-budget random-assignment policy value (the PAPE baseline).
#' @examples
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' rg <- cf_rule(fit, type = "threshold", grid = TRUE)
#' cf_plot_policy_curve(rg)
#' @export
cf_plot_policy_curve <- function(rule_grid, ...) {
  if (!inherits(rule_grid, "cfomics_rule_grid")) {
    rlang::abort("rule_grid must be a cfomics_rule_grid.",
                 class = "cfomics_invalid_input")
  }
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
  if (!inherits(eval, "cfomics_evaluate")) {
    rlang::abort("eval must be a cfomics_evaluate.", class = "cfomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_evaluate")
  if (is.null(eval$gates) || isTRUE(eval$gates$degenerate) ||
      is.null(eval$gates$group_ate)) {
    rlang::abort(
      c("The GATES plot has no group ATEs.",
        i = "Run cf_evaluate() with 'gates' in `metrics` (the default); near-constant estimated effects leave no distinct benefit groups."),
      class = "cfomics_plot_no_data")
  }
  ga <- eval$gates$group_ate
  df <- data.frame(group = as.integer(ga$group),
                   estimate = as.numeric(ga$estimate),
                   lower = as.numeric(ga$ci_lower),
                   upper = as.numeric(ga$ci_upper),
                   stringsAsFactors = FALSE)
  # Overall-ATE reference line: GATES uses equal-count quantile bins, so the
  # unweighted mean of the per-group DR-ATEs equals the overall DR-ATE.
  overall <- mean(df$estimate, na.rm = TRUE)
  bw_p <- eval$gates$best_worst_test$p_value %||% NA_real_
  .cf_build_gates(df, overall, bw_p)
}

#' @keywords internal
#' @noRd
.cf_build_gates <- function(df, overall, bw_p) {
  df$group <- factor(df$group)
  sub <- if (is.finite(bw_p))
    sprintf("dashed = overall ATE; best-worst p = %.3g", bw_p)
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
#' @param itp A \code{cfomics_interpret} object (from \code{cf_interpret}).
#' @param type One of \code{"bar"} (mean|SHAP| importance, default),
#'   \code{"beeswarm"} (per-observation SHAP distribution), or
#'   \code{"dependence"} (SHAP vs feature value for one \code{feature}).
#' @param feature Character; the feature for \code{type = "dependence"}
#'   (required for that view, ignored otherwise).
#' @param top_n Integer; maximum features for the bar / beeswarm view
#'   (default \code{20}).
#' @param surrogate Which surrogate to plot ("ridge" / "xgb"); \code{NULL}
#'   (default) picks the first available.
#' @param ... Passed to the underlying \pkg{shapviz} function.
#' @return A ggplot object.
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
  if (!inherits(itp, "cfomics_interpret")) {
    rlang::abort("itp must be a cfomics_interpret (from cf_interpret).",
                 class = "cfomics_invalid_input")
  }
  type <- match.arg(type)
  sv <- as_shapviz(itp, surrogate = surrogate)
  if (identical(type, "dependence")) {
    if (is.null(feature) || !is.character(feature) || length(feature) != 1L) {
      rlang::abort(
        c("type = 'dependence' requires a single `feature` name.",
          i = "Pass exactly one covariate name to `feature`, e.g. feature = 'GENE1'."),
        class = "cfomics_invalid_input")
    }
    return(shapviz::sv_dependence(sv, v = feature, ...))
  }
  shapviz::sv_importance(sv, kind = type,
                         max_display = as.integer(top_n)[1L], ...)
}


# =========================================================================
# cf_plot_stability: correlation heatmap only, with cor_method.
# =========================================================================

#' Correlation heatmap of a cf_model_stability object.
#'
#' @param stability A \code{cfomics_model_stability} object.
#' @param cor_method Which stored correlation matrix to show; defaults to the
#'   object's primary correlation method.
#' @param ... Ignored.
#' @return A ggplot object (K x K pairwise correlation of OOF tau-hat).
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf) so the pairwise
#' # correlation is non-degenerate; known randomised propensity.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L,  # 1 chain keeps bcf fast
#'                   propensity_score = 0.5, random_state = 1)
#' stab <- cf_model_stability(cmp)
#' cf_plot_stability(stab)
#' @export
cf_plot_stability <- function(stability, cor_method = NULL, ...) {
  if (!inherits(stability, "cfomics_model_stability")) {
    rlang::abort("cf_plot_stability requires a cfomics_model_stability object.",
                 class = "cfomics_plot_stability_invalid_input")
  }
  .cf_require_ggplot("cf_plot_stability")
  if (is.null(cor_method)) cor_method <- stability$meta$primary_cor
  mat <- stability$cor_matrices[[cor_method]]
  if (is.null(mat) || nrow(mat) == 0L) {
    rlang::abort(
      c(sprintf("No '%s' correlation matrix is stored on this stability object.", cor_method),
        i = sprintf("Available: %s. Pass one as `cor_method`.",
                    paste(names(stability$cor_matrices), collapse = ", "))),
      class = "cfomics_plot_stability_no_data")
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
  df$row <- factor(df$row, levels = rn)
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
                  title = "Pairwise agreement of predicted effects",
                  subtitle = sprintf("%s correlation of OOF tau-hat", cor_method)) +
    .cf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}


# =========================================================================
# cf_plot_rule (policy tree): renderer for cf_rule(type = "tree").
# =========================================================================

#' Plot the policy tree from cf_rule(type = "tree").
#'
#' @param rule A \code{cfomics_rule} produced by \code{cf_rule(type = "tree")}.
#' @param backend One of \code{"auto"}, \code{"rpart.plot"}, \code{"base"}.
#' @param ... Ignored.
#' @return Rendered to the active device (invisible \code{NULL}).
#' @examples
#' # Draw the policy tree from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' rule <- cf_rule(fit, type = "tree")
#' cf_plot_rule(rule)
#' @export
cf_plot_rule <- function(rule, backend = c("auto", "rpart.plot", "base"), ...) {
  backend <- match.arg(backend)
  if (!inherits(rule, "cfomics_rule")) {
    rlang::abort("rule must be a cfomics_rule produced by cf_rule().",
                 class = "cfomics_invalid_input")
  }
  if (!identical(rule$type, "tree") || is.null(rule$tree)) {
    rlang::abort(
      "cf_plot_rule() requires a cf_rule(type = 'tree') with a fitted tree.",
      class = "cfomics_invalid_input")
  }
  if (inherits(rule$tree, "policy_tree")) {
    if (requireNamespace("DiagrammeR", quietly = TRUE)) {
      plot(rule$tree)
    } else {
      rlang::inform(
        c("DiagrammeR not installed; printing the policy_tree text summary.",
          i = "install.packages('DiagrammeR') for the graph."),
        class = "cfomics_policy_tree_text_fallback")
      print(rule$tree)
    }
    return(invisible(NULL))
  }
  use_rpartplot <- (backend %in% c("auto", "rpart.plot")) &&
                   requireNamespace("rpart.plot", quietly = TRUE)
  if (use_rpartplot) {
    rpart.plot::rpart.plot(rule$tree, type = 2L, extra = 101L,
                           fallen.leaves = TRUE, main = "Treatment policy tree")
  } else if (identical(backend, "rpart.plot")) {
    rlang::abort(
      c("rpart.plot is required for backend='rpart.plot'.",
        i = "install.packages('rpart.plot') or use backend='base'."),
      class = "cfomics_optional_dependency_missing")
  } else {
    plot(rule$tree, uniform = TRUE, main = "Treatment policy tree")
    graphics::text(rule$tree, use.n = TRUE, all = TRUE, cex = 0.7)
  }
  invisible(NULL)
}


# =========================================================================
# cf_plot_cate: distribution of the estimated CATE for a fitted model.
# =========================================================================

#' Plot the distribution of estimated CATE (tau-hat) for a fitted model.
#'
#' Histogram of the in-sample predicted conditional average treatment effects
#' from a \code{cf_model} -- a quick look at how much treatment-effect
#' heterogeneity the model found.
#'
#' @param fit A \code{cf_model} object from \code{\link{cf_fit}}.
#' @param bins Integer; number of histogram bins (default \code{30}).
#' @param ... Ignored.
#' @return A ggplot object (histogram of predicted CATE), or \code{NULL}
#'   invisibly if no per-sample effects are available.
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
                 class = "cfomics_invalid_input")
  }
  .cf_require_ggplot("cf_plot_cate")
  tau <- tryCatch(predict(fit, type = "tau"), error = function(e) NULL)
  if (is.null(tau) || all(is.na(tau))) {
    rlang::inform("No per-sample CATE available; nothing to plot.")
    return(invisible(NULL))
  }
  df <- data.frame(cate = as.numeric(tau))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["cate"]])) +
    ggplot2::geom_histogram(bins = as.integer(bins)[1L],
                            fill = .CF_OKABE_ITO[1], colour = "white") +
    ggplot2::labs(x = "Predicted CATE (tau-hat)", y = "Count",
                  title = "Distribution of estimated treatment effects") +
    .cf_theme()
}
