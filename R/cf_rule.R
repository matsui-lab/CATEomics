#' Build a treatment-decision rule from a CATE model
#'
#' @description
#' Builds a binary policy \eqn{d(X) \in \{0, 1\}} from the OOF tau-hat
#' of a CATE-capable \code{cf_model} using one of three rule types:
#' \itemize{
#'   \item \code{"threshold"}: \eqn{d(X) = 1\{\hat\tau(X) > threshold\}}.
#'   \item \code{"topk"}: \eqn{d(X) = 1\{rank(\hat\tau(X)) \in top\_K\}}.
#'   \item \code{"tree"}: \code{policytree::policy_tree} with reward
#'     \code{Gamma = [control = 0, treat = tau-hat]}, depth
#'     \code{tree_depth}. The tree maximises \eqn{\sum d(X)\,\hat\tau(X)}
#'     welfare, reflecting the model's tau-hat directly (like threshold /
#'     topk) rather than re-learning from DR scores.
#' }
#' When \code{compute_policy_value = TRUE} (default) the honest DR policy value
#' (independent shared nuisance) is computed and attached, together with
#' treat-all / treat-none / same-budget-random (PAPE) baselines.
#'
#' @section Policy tree engine:
#' The \code{type = "tree"} engine uses \code{policytree::policy_tree}
#' (welfare-maximising), not a descriptive \code{rpart} surrogate.
#'
#' @param model \code{cf_model} (CATE-capable).
#' @param type \code{"threshold"} / \code{"topk"} / \code{"tree"}.
#' @param threshold Numeric, threshold for \code{type = "threshold"}.
#' @param topk_frac Numeric in (0, 1], fraction to treat for
#'   \code{type = "topk"} (default 0.2).
#' @param grid Optional grid sweep for \code{type = "threshold"} / \code{"topk"}:
#'   \code{NULL} (default) fits a single rule; \code{TRUE} sweeps the default
#'   quantile (threshold) / equispaced (topk) grid; a numeric vector sweeps
#'   those values. A grid sweep returns a \code{cfomics_rule_grid} object
#'   (used by \code{\link{cf_plot_policy_curve}}).
#' @param direction \code{"positive"} (treat positive tau-hat) or
#'   \code{"negative"} (treat negative tau-hat). Default \code{"positive"}.
#' @param tree_depth Integer, tree depth for \code{type = "tree"}
#'   (passed to \code{policytree::policy_tree(depth = ...)}, default 2L
#'   = policy_tree default).
#' @param min_node_size Integer, \code{policy_tree} \code{min.node.size}
#'   (default 1L = policy_tree default).
#' @param compute_policy_value Logical, compute DR policy value.
#' @param ... Ignored.
#' @return Object of class \code{c("cfomics_rule", "list")}.
#' @references
#' Zhou, Athey & Wager (2023) Offline Multi-Action Policy Learning:
#' Generalization and Optimization, \emph{Operations Research} 71(1):148–183
#' (policy-tree learning; \pkg{policytree}: Sverdrup, Kanodia, Zhou, Athey &
#' Wager (2020) \emph{JOSS} 5(50):2232).
#' Imai & Li (2023) Experimental Evaluation of Individualized Treatment Rules,
#' \emph{JASA} 118(541):242–256 (PAPE baseline).
#' Athey & Wager (2021) Policy Learning With Observational Data,
#' \emph{Econometrica} 89(1):133–161 (DR policy value).
#' @examples
#' # Build a treat/don't-treat rule from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' rule <- cf_rule(fit, type = "threshold")
#' rule
#' # a policy tree over the covariates
#' tree_rule <- cf_rule(fit, type = "tree")
#' @export
cf_rule <- function(model,
                    type             = c("threshold", "topk", "tree"),
                    threshold        = 0,
                    topk_frac        = 0.2,
                    grid             = NULL,
                    direction        = c("positive", "negative"),
                    tree_depth       = 2L,
                    min_node_size    = 1L,
                    compute_policy_value = TRUE,
                    ...) {
  if (!inherits(model, "cf_model")) {
    rlang::abort("cf_rule requires a cf_model.",
                  class = "cfomics_invalid_input")
  }
  type      <- match.arg(type)
  direction <- match.arg(direction)

  # A non-NULL `grid` sweeps threshold / topk and returns a rule grid
  # (grid = TRUE for the default grid, or a numeric vector of thresholds /
  # topk fractions).
  if (!is.null(grid)) {
    if (identical(type, "tree")) {
      rlang::abort("`grid` is only for type = 'threshold' / 'topk'.",
                   class = "cfomics_invalid_input")
    }
    return(.cf_rule_grid_impl(model, type = type, grid = grid,
                              compute_policy_value = compute_policy_value,
                              direction = direction, ...))
  }

  tau <- model$res$tau_hat_oof
  if (is.null(tau)) {
    rlang::warn(
      c("Model lacks OOF tau-hat; using in-sample tau-hat.",
        i = "Re-fit via cf_compare() for honest rule learning."),
      class = "cfomics_partition_none_warning",
      .frequency = "once",
      .frequency_id = "cfomics_partition_none_warning"
    )
    tau <- model$res$tau_hat
  }
  if (direction == "negative") tau <- -tau

  rule_fn  <- NULL
  tree_obj <- NULL
  d_vec    <- integer(length(tau))

  if (type == "threshold") {
    d_vec   <- as.integer(tau > threshold)
    rule_fn <- function(x_input) as.integer(x_input > threshold)
  } else if (type == "topk") {
    # NA-robust top-k. Rank only the finite tau-hat; treat NA
    # units as not-treated rather than letting sort()/comparison propagate
    # NA into the whole decision vector.
    tau_ok <- tau[is.finite(tau)]
    if (length(tau_ok) == 0L) {
      rlang::abort("cf_rule(type = 'topk'): no finite tau-hat to rank.",
                   class = "cfomics_rule_no_tau")
    }
    k <- max(1L, floor(length(tau_ok) * topk_frac))
    cutoff <- sort(tau_ok, decreasing = TRUE)[k]
    d_vec  <- as.integer(is.finite(tau) & tau >= cutoff)
    rule_fn <- function(x_input) as.integer(is.finite(x_input) & x_input >= cutoff)
  } else if (type == "tree") {
    if (!requireNamespace("policytree", quietly = TRUE)) {
      rlang::abort(
        c("policytree is required for cf_rule(type = 'tree').",
          i = "install.packages('policytree')"),
        class = "cfomics_optional_dependency_missing"
      )
    }
    X <- model$data$X_preprocessed %||% model$data$X
    if (is.null(X)) {
      rlang::abort(
        c("cf_rule(type = 'tree') needs the training X.",
          i = paste("Use a model from cf_compare()$models (X cached)",
                    "or refit via cf_fit() with return_data = TRUE.")),
        class = "cfomics_evaluate_missing_data"
      )
    }
    X_mat <- as.matrix(X)
    if (!is.numeric(X_mat)) {
      rlang::abort(
        c("policytree requires a numeric covariate matrix.",
          i = "Convert factors / characters to numeric before cf_rule."),
        class = "cfomics_rule_X_not_numeric"
      )
    }
    # the tree rewards the model's tau-hat directly
    # (Gamma = [control = 0, treat = tau]) so it maximises sum(d * tau) welfare,
    # faithfully reflecting tau-hat (like threshold / topk) instead of
    # re-learning from DR scores. `direction` is already applied to `tau`, so
    # no separate action flip is needed. Honest evaluation still uses the
    # independent DR policy value.
    Gamma <- cbind(control = rep(0, length(tau)), treat = tau)
    tree_obj <- policytree::policy_tree(
      X             = X_mat,
      Gamma         = Gamma,
      depth         = as.integer(tree_depth),
      min.node.size = as.integer(min_node_size),
      verbose       = FALSE
    )
    # predict() returns action index (1 = control, 2 = treat) by default.
    d_vec <- as.integer(stats::predict(tree_obj, X_mat)) - 1L
    rule_fn <- function(x_input) {
      as.integer(stats::predict(tree_obj, as.matrix(x_input))) - 1L
    }
  }

  out <- list(
    type             = type,
    threshold        = if (type == "threshold") threshold else NA_real_,
    topk_frac        = if (type == "topk") topk_frac else NA_real_,
    direction        = direction,
    decisions        = d_vec,
    rule_fn          = rule_fn,
    tree             = tree_obj,
    treated_frac     = mean(d_vec),
    policy_value     = NULL
  )

  if (isTRUE(compute_policy_value)) {
    out$policy_value <- .compute_policy_value(model, d_vec)
  }
  class(out) <- c("cfomics_rule", "list")
  out
}

#' Internal: sweep cf_rule over a threshold / topk grid.
#'
#' Called by \code{cf_rule(..., grid = TRUE)} or a numeric \code{grid}. The
#' default grid is quantile-based for \code{type = "threshold"} and equispaced
#' for \code{type = "topk"}. Returns \code{c("cfomics_rule_grid", "list")}.
#' @keywords internal
#' @noRd
.cf_rule_grid_impl <- function(model,
                               type                 = c("threshold", "topk"),
                               grid                 = TRUE,
                               compute_policy_value = TRUE,
                               ...) {
  type <- match.arg(type)
  if (isTRUE(grid)) {
    tau <- model$res$tau_hat_oof %||% model$res$tau_hat
    grid <- if (identical(type, "threshold")) {
      unname(stats::quantile(tau, probs = seq(0, 1, by = 0.05),
                              na.rm = TRUE))
    } else {
      seq(0.05, 1.0, by = 0.05)
    }
  }

  # capture the treat-all / treat-none policy values once so the
  # policy-value curve can draw the same-budget random-assignment baseline
  # (they do not depend on the rule, only on the shared DR pseudo-outcomes).
  treat_all <- NA_real_; treat_none <- NA_real_
  results <- do.call(rbind, lapply(grid, function(g) {
    rule <- if (identical(type, "threshold")) {
      cf_rule(model, type = "threshold", threshold = g,
              compute_policy_value = compute_policy_value, ...)
    } else {
      cf_rule(model, type = "topk", topk_frac = g,
              compute_policy_value = compute_policy_value, ...)
    }
    pv <- rule$policy_value
    if (is.null(pv)) pv <- list(estimate = NA_real_, ci_lower = NA_real_,
                                  ci_upper = NA_real_)
    if (is.na(treat_all) && !is.null(pv$baseline_treat_all))
      treat_all <<- pv$baseline_treat_all$estimate
    if (is.na(treat_none) && !is.null(pv$baseline_treat_none))
      treat_none <<- pv$baseline_treat_none$estimate
    data.frame(
      parameter        = g,
      n_treated        = sum(rule$decisions),
      policy_value     = pv$estimate,
      ci_lower         = pv$ci_lower,
      ci_upper         = pv$ci_upper
    )
  }))

  out <- list(
    type       = type,
    grid       = grid,
    results    = results,
    n          = length(model$res$tau_hat_oof %||% model$res$tau_hat),
    treat_all  = treat_all,
    treat_none = treat_none
  )
  class(out) <- c("cfomics_rule_grid", "list")
  out
}

# -------------------------------------------------------------------------
# Helper: DR policy value (Athey-Wager 2017)
# -------------------------------------------------------------------------

.compute_policy_value <- function(model, d_vec) {
  Y <- model$data$Y
  T <- model$data$T
  nuis <- model$res$nuisance   # independent shared eval nuisance
  if (is.null(Y) || is.null(T) || is.null(nuis)) {
    return(list(
      estimate              = NA_real_,
      ci_lower              = NA_real_,
      ci_upper              = NA_real_,
      baseline_treat_all    = list(estimate = NA_real_),
      baseline_treat_none   = list(estimate = NA_real_),
      baseline_random       = list(estimate = NA_real_),
      improvement_over_all  = list(estimate = NA_real_),
      improvement_over_none = list(estimate = NA_real_),
      improvement_over_random = list(estimate = NA_real_),
      type                  = "dr"
    ))
  }
  Y <- as.numeric(Y)
  T <- as.numeric(T)
  pi_oof  <- pmin(pmax(nuis$pi_oof,  1e-4), 1 - 1e-4)
  mu0_oof <- nuis$mu0_oof
  mu1_oof <- nuis$mu1_oof

  Y_DR_1 <- mu1_oof + (T / pi_oof) * (Y - mu1_oof)
  Y_DR_0 <- mu0_oof + ((1 - T) / (1 - pi_oof)) * (Y - mu0_oof)

  pv_per_obs <- d_vec * Y_DR_1 + (1 - d_vec) * Y_DR_0
  pv_est     <- mean(pv_per_obs)
  pv_se      <- stats::sd(pv_per_obs) / sqrt(length(pv_per_obs))

  treat_all_est <- mean(Y_DR_1)
  treat_all_se  <- stats::sd(Y_DR_1) / sqrt(length(Y_DR_1))
  treat_none_est <- mean(Y_DR_0)
  treat_none_se  <- stats::sd(Y_DR_0) / sqrt(length(Y_DR_0))

  imp_all_est <- pv_est - treat_all_est
  imp_all_se  <- stats::sd(pv_per_obs - Y_DR_1) / sqrt(length(Y_DR_1))
  imp_none_est <- pv_est - treat_none_est
  imp_none_se  <- stats::sd(pv_per_obs - Y_DR_0) / sqrt(length(Y_DR_0))

  # same-budget random baseline (evalITR PAPE): treat a random
  # `treated_frac` fraction. Its DR value is the treated_frac-weighted mix of
  # treat-all / treat-none. improvement_over_random is the headline deployment
  # gain over a rule with the same treatment budget but no targeting.
  treated_frac    <- mean(d_vec)
  random_per_obs  <- treated_frac * Y_DR_1 + (1 - treated_frac) * Y_DR_0
  random_est      <- mean(random_per_obs)
  random_se       <- stats::sd(random_per_obs) / sqrt(length(random_per_obs))
  imp_random_est  <- pv_est - random_est
  imp_random_se   <- stats::sd(pv_per_obs - random_per_obs) / sqrt(length(pv_per_obs))

  list(
    estimate              = pv_est,
    ci_lower              = pv_est - 1.96 * pv_se,
    ci_upper              = pv_est + 1.96 * pv_se,
    baseline_treat_all    = list(
      estimate = treat_all_est,
      ci_lower = treat_all_est - 1.96 * treat_all_se,
      ci_upper = treat_all_est + 1.96 * treat_all_se
    ),
    baseline_treat_none   = list(
      estimate = treat_none_est,
      ci_lower = treat_none_est - 1.96 * treat_none_se,
      ci_upper = treat_none_est + 1.96 * treat_none_se
    ),
    baseline_random       = list(   # PAPE same-budget random
      estimate = random_est,
      ci_lower = random_est - 1.96 * random_se,
      ci_upper = random_est + 1.96 * random_se
    ),
    improvement_over_random = list(  # headline deployment gain
      estimate = imp_random_est,
      ci_lower = imp_random_est - 1.96 * imp_random_se,
      ci_upper = imp_random_est + 1.96 * imp_random_se
    ),
    improvement_over_all  = list(
      estimate = imp_all_est,
      ci_lower = imp_all_est - 1.96 * imp_all_se,
      ci_upper = imp_all_est + 1.96 * imp_all_se
    ),
    improvement_over_none = list(
      estimate = imp_none_est,
      ci_lower = imp_none_est - 1.96 * imp_none_se,
      ci_upper = imp_none_est + 1.96 * imp_none_se
    ),
    type                  = "dr"
  )
}

#' @export
print.cfomics_rule <- function(x, ...) {
  cat("<cfomics_rule>: type = ", x$type, ", direction = ", x$direction,
      "\n", sep = "")
  cat(sprintf("  treated fraction: %.3f\n", x$treated_frac))
  if (!is.null(x$policy_value)) {
    pv <- x$policy_value
    cat(sprintf("  policy value (DR): %.4g [%.4g, %.4g]\n",
                pv$estimate, pv$ci_lower, pv$ci_upper))
    ir <- pv$improvement_over_random
    if (!is.null(ir) && is.finite(ir$estimate)) {
      cat(sprintf("  improvement over same-budget random (PAPE): %.4g [%.4g, %.4g]\n",
                  ir$estimate, ir$ci_lower, ir$ci_upper))
    }
  }
  invisible(x)
}

#' @export
print.cfomics_rule_grid <- function(x, ...) {
  cat("<cfomics_rule_grid>: type = ", x$type,
      ", grid points = ", length(x$grid), "\n", sep = "")
  print(utils::head(x$results, 5))
  if (nrow(x$results) > 5) cat("  ...\n")
  invisible(x)
}
