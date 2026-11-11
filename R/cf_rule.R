#' Build a treatment-decision rule from a CATE model
#'
#' @description
#' Builds a binary policy \eqn{d(X) \in \{0, 1\}} from the OOF tau-hat
#' of a CATE-capable \code{cf_model} using one of three rule types. Every
#' type treats the units with the \emph{largest} \code{tau_hat}, on CATEomics'
#' convention that a larger \code{Y} is the better outcome (see \code{Y} in
#' \code{\link{cf_fit}}); for an outcome whose smaller values are better, fit
#' the negated outcome.
#' \itemize{
#'   \item \code{"threshold"}: \eqn{d(X) = 1\{\hat\tau(X) > threshold\}}.
#'   \item \code{"topk"}: \eqn{d(X) = 1\{rank(\hat\tau(X)) \in top\_K\}}.
#'   \item \code{"tree"}: \code{policytree::policy_tree} with reward
#'     \code{Gamma = [control = 0, treat = tau-hat]}, depth
#'     \code{tree_depth}. The tree maximises \eqn{\sum d(X)\,\hat\tau(X)}
#'     welfare, so it is built from the model's tau-hat, as threshold and
#'     topk are.
#' }
#' When \code{compute_policy_value = TRUE} (default) the DR policy value is
#' computed from the model's shared evaluation nuisance (out-of-fold when the
#' fit used \code{partition = "oof"}) and attached, together with
#' treat-all / treat-none / same-budget-random (PAPE) baselines.
#'
#' @section Policy tree engine:
#' The \code{type = "tree"} engine is \code{policytree::policy_tree}. The tree
#' is fitted once, on all rows, and the decisions come from predicting on
#' those same rows, so the tree's structure is in-sample; the rewards it
#' maximises are the out-of-fold \code{tau_hat} and the policy value is scored
#' with the out-of-fold DR pseudo-outcome. The bootstrap behind the policy
#' value does not refit it: it re-applies the same tree to each replicate.
#'
#' @param model \code{cf_model} (CATE-capable).
#' @param type \code{"threshold"} (default) / \code{"topk"} / \code{"tree"}.
#' @param threshold Numeric (default \code{0}), the effect above which
#'   \code{type = "threshold"} treats. When every estimated effect has the
#'   same sign the default treats everyone (or no one): the rule is then one
#'   of the fixed policies it is compared with, and \code{print()} says so.
#' @param topk_frac Numeric in (0, 1], fraction to treat for
#'   \code{type = "topk"} (default 0.2).
#' @param grid Optional grid sweep for \code{type = "threshold"} / \code{"topk"}:
#'   \code{NULL} (default) fits a single rule; \code{TRUE} sweeps the default
#'   quantile (threshold) / equispaced (topk) grid; a numeric vector sweeps
#'   those values. A grid sweep returns a \code{cateomics_rule_grid} object
#'   (used by \code{\link{cf_plot_policy_curve}}).
#' @param tree_depth Integer, tree depth for \code{type = "tree"}
#'   (passed to \code{policytree::policy_tree(depth = ...)}, default 2L
#'   = policy_tree default).
#' @param min_node_size Integer, \code{policy_tree} \code{min.node.size}
#'   (default 1L = policy_tree default).
#' @param compute_policy_value Logical (default \code{TRUE}); compute the
#'   DR policy value and its baselines.
#' @param n_bootstrap Integer, half-sample bootstrap replicates behind the
#'   policy-value standard errors (default 200).
#' @param random_state Integer seed for that bootstrap (default \code{0}).
#' @param ... Not used; naming an unknown argument is an error.
#' @return Object of class \code{c("cateomics_rule", "list")} with
#'   \code{type}, \code{threshold} (the cut-point given for
#'   \code{type = "threshold"}, \code{NA} otherwise), \code{topk_frac},
#'   \code{cutoff} (the effect at which \code{type = "topk"} ended up cutting
#'   this sample, \code{NA} otherwise), \code{decisions} (0/1 per sample),
#'   \code{tree} (the fitted \pkg{policytree} object for
#'   \code{type = "tree"}; fitted on all rows and predicted on those same
#'   rows, so its structure is in-sample), \code{treated_frac} and, unless
#'   \code{compute_policy_value = FALSE}, \code{policy_value} with its
#'   treat-all / treat-none / same-budget-random baselines. With
#'   \code{grid}, an object of class \code{c("cateomics_rule_grid", "list")}
#'   instead.
#'
#' @section What \code{policy_value} reports:
#'   Four levels -- the rule's own value and the treat-all, treat-none and
#'   same-budget-random baselines -- and the three differences of the rule from
#'   those baselines. Each carries \code{estimate}, \code{se},
#'   \code{ci_lower}, \code{ci_upper} and \code{ci_method}; the interval is
#'   \eqn{\hat\theta \pm 1.96\,SE}.
#'
#'   \code{ci_method} names the standard error that component actually used.
#'   Every component whose value depends on the rule, or on the budget the rule
#'   spends, is resampled with the half-sample bootstrap: the replicates draw
#'   rows and apply the rule as it is stated, which for \code{"topk"} means
#'   reading its cut-point in the replicate and for \code{"threshold"} and
#'   \code{"tree"} means applying the same fixed rule. Treat-all and
#'   treat-none are plain means of fixed per-row values, so they keep the
#'   influence-function standard error. A component with fewer than 20 usable
#'   replicates keeps that closed form too, and \code{ci_method} says so.
#'
#'   The three differences also carry \code{p_value}, the two-sided
#'   normal-approximation p-value against \eqn{H_0 = 0} on the reported
#'   standard error: that the rule buys nothing over the baseline.
#'   \code{improvement_over_random} and \code{improvement_over_none} add
#'   \code{p_value_one_sided} against the alternative that the rule does
#'   better. \code{improvement_over_all} has the two-sided p-value only,
#'   because both of its directions are substantive. The four levels carry no
#'   p-value: \eqn{H_0 = 0} would read "the mean outcome under this policy is
#'   zero", which is a point on the outcome scale rather than a statement about
#'   the rule.
#' @references
#' Zhou, Athey & Wager (2023) Offline Multi-Action Policy Learning:
#' Generalization and Optimization, \emph{Operations Research} 71(1):148-183
#' (policy-tree learning; \pkg{policytree}: Sverdrup, Kanodia, Zhou, Athey &
#' Wager (2020) \emph{JOSS} 5(50):2232).
#'
#' Imai & Li (2023) Experimental Evaluation of Individualized Treatment Rules,
#' \emph{JASA} 118(541):242-256 (PAPE baseline).
#'
#' Athey & Wager (2021) Policy Learning With Observational Data,
#' \emph{Econometrica} 89(1):133-161 (DR policy value).
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
                    tree_depth       = 2L,
                    min_node_size    = 1L,
                    compute_policy_value = TRUE,
                    n_bootstrap      = 200L,
                    random_state     = 0L,
                    ...) {
  if (!inherits(model, "cf_model")) {
    rlang::abort("cf_rule requires a cf_model.",
                  class = "cateomics_invalid_input")
  }
  .cf_reject_dots(..., .fn = "cf_rule()")
  type      <- match.arg(type)
  # The fraction is read as "the top q of the estimates", so it has to be a
  # quantile. Checked with the same helper `cf_evaluate()` uses for its own
  # quantile arguments.
  if (identical(type, "topk")) .cf_check_quantile_arg(topk_frac, "topk_frac")

  # A non-NULL `grid` sweeps threshold / topk and returns a rule grid
  # (grid = TRUE for the default grid, or a numeric vector of thresholds /
  # topk fractions).
  if (!is.null(grid)) {
    if (identical(type, "tree")) {
      rlang::abort("`grid` is only for type = 'threshold' / 'topk'.",
                   class = "cateomics_invalid_input")
    }
    return(.cf_rule_grid_impl(model, type = type, grid = grid,
                              compute_policy_value = compute_policy_value,
                              n_bootstrap = n_bootstrap,
                              random_state = random_state, ...))
  }

  tau <- model$res$tau_hat_oof
  if (is.null(tau)) {
    rlang::warn(
      c("Model lacks OOF tau-hat; using in-sample tau-hat.",
        i = "Re-fit via cf_compare() for honest rule learning."),
      class = "cateomics_partition_none_warning"
    )
    tau <- model$res$tau_hat
  }

  cutoff   <- NA_real_
  tree_obj <- NULL
  d_vec    <- integer(length(tau))

  if (type == "threshold") {
    d_vec   <- as.integer(tau > threshold)
  } else if (type == "topk") {
    # NA-robust top-k. Rank only the finite tau-hat; treat NA
    # units as not-treated rather than letting sort()/comparison propagate
    # NA into the whole decision vector.
    tau_ok <- tau[is.finite(tau)]
    if (length(tau_ok) == 0L) {
      rlang::abort("cf_rule(type = 'topk'): no finite tau-hat to rank.",
                   class = "cateomics_rule_no_tau")
    }
    k <- max(1L, floor(length(tau_ok) * topk_frac))
    cutoff <- sort(tau_ok, decreasing = TRUE)[k]
    d_vec  <- as.integer(is.finite(tau) & tau >= cutoff)
  } else if (type == "tree") {
    if (!requireNamespace("policytree", quietly = TRUE)) {
      rlang::abort(
        c("policytree is required for cf_rule(type = 'tree').",
          i = "install.packages('policytree')"),
        class = "cateomics_optional_dependency_missing"
      )
    }
    X <- model$data$X_preprocessed %||% model$data$X
    if (is.null(X)) {
      rlang::abort(
        c("cf_rule(type = 'tree') needs the training X.",
          i = paste("This object carries no X; refit with cf_fit() or",
                    "cf_compare().")),
        class = "cateomics_evaluate_missing_data"
      )
    }
    X_mat <- as.matrix(X)
    if (!is.numeric(X_mat)) {
      rlang::abort(
        c("policytree requires a numeric covariate matrix.",
          i = "Convert factors / characters to numeric before cf_rule."),
        class = "cateomics_rule_X_not_numeric"
      )
    }
    # Gamma = [control = 0, treat = tau], so the tree maximises
    # sum(d * tau).
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
  }

  out <- list(
    type             = type,
    threshold        = if (type == "threshold") threshold else NA_real_,
    topk_frac        = if (type == "topk") topk_frac else NA_real_,
    decisions        = d_vec,
    cutoff           = cutoff,
    tree             = tree_obj,
    treated_frac     = mean(d_vec),
    policy_value     = NULL
  )

  if (isTRUE(compute_policy_value)) {
    # The replicates resample rows; the rule stays the one that was reported.
    # "topk" is stated as a fraction of whoever is at hand, so applying it to a
    # subsample means reading its cut-point there. "threshold" and "tree" are
    # stated as a fixed rule, so a row's decision does not depend on who else
    # is in the replicate and the stored decisions are read directly.
    decide_fn <- if (type == "topk") {
      force(topk_frac)
      function(i) {
        tau_sub <- tau[i]
        ok <- tau_sub[is.finite(tau_sub)]
        if (length(ok) == 0L) return(rep(0L, length(tau_sub)))
        kk <- max(1L, floor(length(ok) * topk_frac))
        as.integer(is.finite(tau_sub) & tau_sub >= sort(ok, decreasing = TRUE)[kk])
      }
    } else {
      function(i) d_vec[i]
    }
    out$policy_value <- .compute_policy_value(
      model, d_vec, decide_fn = decide_fn,
      n_bootstrap = n_bootstrap, random_state = random_state)
  }
  class(out) <- c("cateomics_rule", "list")
  out
}

#' Internal: sweep cf_rule over a threshold / topk grid.
#'
#' Called by \code{cf_rule(..., grid = TRUE)} or a numeric \code{grid}. The
#' default grid is quantile-based for \code{type = "threshold"} and equispaced
#' for \code{type = "topk"}. Returns \code{c("cateomics_rule_grid", "list")}.
#' @keywords internal
#' @noRd
.cf_rule_grid_impl <- function(model,
                               type                 = c("threshold", "topk"),
                               grid                 = TRUE,
                               compute_policy_value = TRUE,
                               n_bootstrap          = 200L,
                               random_state         = NULL,
                               ...) {
  type <- match.arg(type)
  if (isTRUE(grid)) {
    # The default threshold grid is a set of tau-hat quantiles, taken on the
    # scale each grid point is compared against in `cf_rule()`.
    tau <- model$res$tau_hat_oof %||% model$res$tau_hat
    grid <- if (identical(type, "threshold")) {
      unname(stats::quantile(tau, probs = seq(0, 1, by = 0.05),
                              na.rm = TRUE))
    } else {
      seq(0.05, 1.0, by = 0.05)
    }
  }
  # Checked before the sweep starts: a bad grid point would otherwise abort
  # part way through, after the earlier points had already been fitted.
  if (identical(type, "topk")) .cf_check_quantile_arg(grid, "grid")

  # capture the treat-all / treat-none policy values once so the
  # policy-value curve can draw the same-budget random-assignment baseline
  # (they do not depend on the rule, only on the shared DR pseudo-outcomes).
  treat_all <- NA_real_; treat_none <- NA_real_
  # Every grid point builds a rule from the same model, so a caveat about that
  # model is the same caveat each time: the sweep says it once.
  results <- .cf_dedup_conditions(do.call(rbind, lapply(grid, function(g) {
    rule <- if (identical(type, "threshold")) {
      cf_rule(model, type = "threshold", threshold = g,
              compute_policy_value = compute_policy_value,
              n_bootstrap = n_bootstrap, random_state = random_state, ...)
    } else {
      cf_rule(model, type = "topk", topk_frac = g,
              compute_policy_value = compute_policy_value,
              n_bootstrap = n_bootstrap, random_state = random_state, ...)
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
  })))

  out <- list(
    type       = type,
    grid       = grid,
    results    = results,
    n          = length(model$res$tau_hat_oof %||% model$res$tau_hat),
    treat_all  = treat_all,
    treat_none = treat_none
  )
  class(out) <- c("cateomics_rule_grid", "list")
  out
}

# -------------------------------------------------------------------------
# Helper: DR policy value (Athey & Wager 2021)
# -------------------------------------------------------------------------

.compute_policy_value <- function(model, d_vec, decide_fn = NULL,
                                  n_bootstrap = 200L,
                                  random_state = NULL) {
  Y <- model$data$Y
  T <- model$data$T
  nuis <- model$res$nuisance   # independent shared eval nuisance
  if (is.null(Y) || is.null(T) || is.null(nuis)) {
    na <- NA_real_
    return(list(
      estimate              = na,
      se                    = na,
      ci_lower              = na,
      ci_upper              = na,
      baseline_treat_all    = .cf_pv_component(na, na, NA_character_, "none"),
      baseline_treat_none   = .cf_pv_component(na, na, NA_character_, "none"),
      baseline_random       = .cf_pv_component(na, na, NA_character_, "none"),
      improvement_over_random = .cf_pv_component(na, na, NA_character_, "both"),
      improvement_over_all  = .cf_pv_component(na, na, NA_character_, "two_sided"),
      improvement_over_none = .cf_pv_component(na, na, NA_character_, "both"),
      type                  = "dr",
      ci_method             = NA_character_
    ))
  }
  Y <- as.numeric(Y)
  T <- as.numeric(T)
  # The evaluation nuisance carries an absolute floor only, not the
  # estimator-side `clip`; this reapplies that floor.
  pi_oof  <- .cf_clip_ps(nuis$pi_oof)
  mu0_oof <- nuis$mu0_oof
  mu1_oof <- nuis$mu1_oof

  Y_DR_1 <- mu1_oof + (T / pi_oof) * (Y - mu1_oof)
  Y_DR_0 <- mu0_oof + ((1 - T) / (1 - pi_oof)) * (Y - mu0_oof)

  pv_per_obs <- d_vec * Y_DR_1 + (1 - d_vec) * Y_DR_0
  pv_est     <- mean(pv_per_obs)

  treat_all_est  <- mean(Y_DR_1)
  treat_none_est <- mean(Y_DR_0)
  imp_all_est    <- pv_est - treat_all_est
  imp_none_est   <- pv_est - treat_none_est

  # same-budget random baseline (evalITR PAPE): treat a random
  # `treated_frac` fraction. Its DR value is the treated_frac-weighted mix of
  # treat-all / treat-none. improvement_over_random is the difference from a
  # rule with the same treatment budget but no targeting.
  treated_frac    <- mean(d_vec)
  random_per_obs  <- treated_frac * Y_DR_1 + (1 - treated_frac) * Y_DR_0
  random_est      <- mean(random_per_obs)
  imp_random_est  <- pv_est - random_est

  # Standard errors. The closed form below is exact for a quantity that is a
  # plain mean of per-row values fixed by the rule, which treat-all and
  # treat-none always are and which the rule's own value is whenever the rule
  # is stated as a fixed threshold or a fitted tree. It is not exact for the
  # two components defined through the budget the rule spends, because that
  # budget is realised rather than given. The replicates are therefore run for
  # everything except the two fixed policies, and each component keeps the
  # closed form only when the replicates could not serve it.
  n_obs <- length(pv_per_obs)
  se <- c(pv         = stats::sd(pv_per_obs) / sqrt(n_obs),
          all        = stats::sd(Y_DR_1) / sqrt(n_obs),
          none       = stats::sd(Y_DR_0) / sqrt(n_obs),
          random     = stats::sd(random_per_obs) / sqrt(n_obs),
          imp_all    = stats::sd(pv_per_obs - Y_DR_1) / sqrt(n_obs),
          imp_none   = stats::sd(pv_per_obs - Y_DR_0) / sqrt(n_obs),
          imp_random = stats::sd(pv_per_obs - random_per_obs) / sqrt(n_obs))
  boot_used <- character(0)
  if (is.function(decide_fn)) {
    # A rule that treats nobody or everybody makes the same-budget random
    # policy one of the two fixed policies, so both components defined through
    # that budget are identically zero or identically a fixed baseline and the
    # closed form is exact for them. Leave them out of the replicates rather
    # than have every replicate come back unusable.
    budget_informative <- treated_frac > 0 && treated_frac < 1
    boot <- .cf_halfsample_se(
      length(Y),
      stat_fn = function(i) {
        d_i   <- decide_fn(i)
        pv_i  <- d_i * Y_DR_1[i] + (1 - d_i) * Y_DR_0[i]
        frac  <- mean(d_i)
        out   <- c(pv       = mean(pv_i),
                   imp_all  = mean(pv_i) - mean(Y_DR_1[i]),
                   imp_none = mean(pv_i) - mean(Y_DR_0[i]))
        if (!budget_informative) return(out)
        # A single replicate can still land on a degenerate budget; drop that
        # replicate for the two budget components. The per-component NA path in
        # .cf_halfsample_se falls back to the closed form if too many do.
        degenerate <- !is.finite(frac) || frac <= 0 || frac >= 1
        rnd_i <- frac * Y_DR_1[i] + (1 - frac) * Y_DR_0[i]
        c(out,
          random     = if (degenerate) NA_real_ else mean(rnd_i),
          imp_random = if (degenerate) NA_real_ else mean(pv_i) - mean(rnd_i))
      },
      n_bootstrap = n_bootstrap, random_state = random_state)
    for (nm in names(boot)) if (is.finite(boot[[nm]])) {
      se[[nm]] <- boot[[nm]]
      boot_used <- c(boot_used, nm)
    }
  }
  # The label names the standard error that was actually used, per component.
  # Every component whose value depends on the rule or on the budget the rule
  # spends is resampled; treat-all and treat-none are plain means of fixed
  # per-row values, so they keep the closed form. Below 20 usable replicates a
  # component keeps the closed form too, and the label says so.
  .se_label <- function(nm) if (nm %in% boot_used)
    "half-sample bootstrap SE" else "influence function SE"
  pv_se <- unname(se[["pv"]]); treat_all_se <- unname(se[["all"]])
  treat_none_se <- unname(se[["none"]]); random_se <- unname(se[["random"]])
  imp_all_se <- unname(se[["imp_all"]]); imp_none_se <- unname(se[["imp_none"]])
  imp_random_se <- unname(se[["imp_random"]])

  list(
    estimate              = pv_est,
    se                    = pv_se,
    ci_lower              = pv_est - 1.96 * pv_se,
    ci_upper              = pv_est + 1.96 * pv_se,
    baseline_treat_all    = .cf_pv_component(
      treat_all_est,  treat_all_se,  .se_label("all"),        "none"),
    baseline_treat_none   = .cf_pv_component(
      treat_none_est, treat_none_se, .se_label("none"),       "none"),
    baseline_random       = .cf_pv_component(   # PAPE same-budget random
      random_est,     random_se,     .se_label("random"),     "none"),
    improvement_over_random = .cf_pv_component(
      imp_random_est, imp_random_se, .se_label("imp_random"), "both"),
    improvement_over_all  = .cf_pv_component(
      imp_all_est,    imp_all_se,    .se_label("imp_all"),    "two_sided"),
    improvement_over_none = .cf_pv_component(
      imp_none_est,   imp_none_se,   .se_label("imp_none"),   "both"),
    type                  = "dr",
    ci_method             = .se_label("pv")
  )
}

#' Internal: one component of the policy value
#'
#' Builds the \code{estimate} / \code{se} / interval / p-value / label list that
#' every baseline and improvement component of \code{policy_value} carries, so
#' that the fields are the same whether or not the quantities could be computed.
#'
#' A p-value is attached to the three improvement components only. For the
#' levels the null \eqn{H_0 = 0} would read "the mean outcome under this policy
#' is zero", which is a point on the outcome scale and not a hypothesis about
#' the rule. \code{improvement_over_all} gets the two-sided p-value alone
#' because both directions are substantive: above zero says that treating a
#' subset beats treating everybody, below zero that it loses benefit. The other
#' two improvements also get the one-sided p-value against the alternative that
#' the rule does better than the baseline.
#'
#' @param est,se The estimate and the standard error reported for it.
#' @param label The \code{ci_method} string naming that standard error.
#' @param p One of \code{"none"} / \code{"two_sided"} / \code{"both"}.
#' @return A named list; the interval is the normal approximation
#'   \eqn{\hat\theta \pm 1.96\,SE} and the p-values use the same SE.
#' @keywords internal
#' @noRd
.cf_pv_component <- function(est, se, label, p = "none") {
  usable <- is.finite(est) && is.finite(se) && se > 0
  pp <- if (usable) .cf_normal_p(est, se, 0) else
    list(two_sided = NA_real_, one_sided = NA_real_)
  out <- list(estimate = est, se = se,
              ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
  if (p != "none") out$p_value <- pp$two_sided
  if (p == "both") out$p_value_one_sided <- pp$one_sided
  out$ci_method <- label
  out
}

#' @export
print.cateomics_rule <- function(x, ...) {
  cat("<cateomics_rule>: type = ", x$type, "\n", sep = "")
  # A rule that treats everyone (or no one) is one of the fixed policies it is
  # being compared with, so two of the three contrasts are zero by definition
  # rather than by measurement. Saying which case it is keeps a structural
  # zero from reading as "the rule gained nothing".
  frac <- x$treated_frac
  all_treated  <- isTRUE(is.finite(frac) && frac >= 1)
  none_treated <- isTRUE(is.finite(frac) && frac <= 0)
  no_split <- identical(x$type, "tree") &&
    !is.null(x$tree) && isTRUE(length(x$tree$nodes) <= 1L)
  # A stump assigns one action to everyone, and which action that is depends
  # on the rewards: saying only that everyone was treated would leave the
  # reader of a control-assigning stump without the one fact the note exists
  # for, that the tree never split.
  note <- if (no_split && all_treated)
    "  (the tree did not split: this rule is treat-all)"
  else if (no_split && none_treated)
    "  (the tree did not split: this rule is treat-none)"
  else if (all_treated) "  (every sample treated: this rule is treat-all)"
  else if (none_treated) "  (no sample treated: this rule is treat-none)"
  else ""
  cat(sprintf("  treated fraction: %.3f%s\n", frac, note))
  if (!is.null(x$policy_value)) {
    pv <- x$policy_value
    cat(sprintf("  policy value (DR): %.4g [%.4g, %.4g]\n",
                pv$estimate, pv$ci_lower, pv$ci_upper))
    # All three contrasts, not one: an improvement over treat-all can be the
    # one that matters and was not shown.
    lines <- list(
      c("improvement over same-budget random (PAPE)", "improvement_over_random",
        if (all_treated || none_treated) "zero by construction" else ""),
      c("improvement over treat-all", "improvement_over_all",
        if (all_treated) "zero by construction" else ""),
      c("improvement over treat-none", "improvement_over_none",
        if (none_treated) "zero by construction" else ""))
    # The colon stays attached to its label and the padding goes after it, so
    # the numbers line up without the colon floating away from the words.
    width <- max(vapply(lines, function(l) nchar(l[1L]), 0L)) + 1L
    for (l in lines) {
      it <- pv[[l[2L]]]
      # Length as well as value: a component whose estimate is absent gives
      # `logical(0)` to the test, which stops the print instead of skipping
      # the line.
      if (is.null(it) || length(it$estimate) != 1L ||
          !is.finite(it$estimate)) next
      lab <- paste0(l[1L], ":")
      if (nzchar(l[3L])) {
        cat(sprintf("  %-*s %.4g  (%s)\n", width, lab, it$estimate, l[3L]))
      } else {
        cat(sprintf("  %-*s %.4g [%.4g, %.4g]\n", width, lab,
                    it$estimate, it$ci_lower, it$ci_upper))
      }
    }
  }
  invisible(x)
}

#' @export
print.cateomics_rule_grid <- function(x, ...) {
  cat("<cateomics_rule_grid>: type = ", x$type,
      ", grid points = ", length(x$grid), "\n", sep = "")
  if (is.null(x$results) || !nrow(x$results)) {
    cat("  no sweep results in this object\n")
    return(invisible(x))
  }
  print(utils::head(x$results, 5))
  if (nrow(x$results) > 5) cat("  ...\n")
  invisible(x)
}
