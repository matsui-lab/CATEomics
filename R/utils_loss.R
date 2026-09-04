#' Resolve eval_nuisance into a list(outcome=, propensity=) spec.
#'
#' eval_nuisance can be a single string ("rf", "xgboost", etc.) or a named
#' list with elements \code{outcome} and \code{propensity}. A string is
#' broadcast: list(outcome = s, propensity = s).
#' When a list is provided, missing elements fall back to the partner. The
#' special value \code{propensity = "known"} signals that an externally
#' supplied propensity score should be used in place of a fitted PS model
#' (RCT context).
#'
#' Allowed learner names match \code{.fit_base_learner()} (linear, logistic,
#' lasso, ridge, elastic_net, rf, xgboost), plus the sentinel "known" for
#' the propensity slot only. Validation is strict.
#'
#' @keywords internal
#' @noRd
.cfomics_resolve_eval_nuisance <- function(eval_nuisance) {
  allowed <- c("linear", "logistic", "lasso", "ridge", "elastic_net",
               "rf", "xgboost")
  if (is.character(eval_nuisance) && length(eval_nuisance) == 1L &&
      !is.null(eval_nuisance) && !is.na(eval_nuisance)) {
    if (!(eval_nuisance %in% allowed)) {
      rlang::abort(
        sprintf("eval_nuisance '%s' is not in the allowed pool: %s",
                eval_nuisance, paste(allowed, collapse = ", ")),
        class = "cfomics_eval_nuisance_unknown"
      )
    }
    return(list(outcome = eval_nuisance, propensity = eval_nuisance))
  }
  if (is.list(eval_nuisance)) {
    out <- eval_nuisance$outcome
    pro <- eval_nuisance$propensity
    if (is.null(out) && is.null(pro)) {
      rlang::abort(
        paste0("eval_nuisance list must have at least one of `outcome` ",
               "or `propensity`."),
        class = "cfomics_eval_nuisance_empty_list"
      )
    }
    # Broadcast the sibling, but "known" is only valid on the propensity
    # slot. When the user supplied `propensity = "known"` without an outcome
    # learner, fall back to "rf" (the eval default) instead of copying
    # "known" into the outcome slot.
    if (is.null(out)) out <- if (identical(pro, "known")) "rf" else pro
    if (is.null(pro)) pro <- out
    if (!(out %in% allowed)) {
      rlang::abort(
        sprintf("eval_nuisance$outcome '%s' is not allowed: %s",
                out, paste(allowed, collapse = ", ")),
        class = "cfomics_eval_nuisance_unknown"
      )
    }
    if (!(pro %in% c(allowed, "known"))) {
      rlang::abort(
        sprintf(paste0("eval_nuisance$propensity '%s' is not allowed: ",
                       "%s or 'known'."),
                pro, paste(allowed, collapse = ", ")),
        class = "cfomics_eval_nuisance_unknown"
      )
    }
    return(list(outcome = out, propensity = pro))
  }
  rlang::abort(
    "eval_nuisance must be a single learner name or a named list.",
    class = "cfomics_eval_nuisance_invalid_type"
  )
}


#' Broadcast a known propensity score to length n with finite clipping.
#' @keywords internal
#' @noRd
.cfomics_broadcast_known_ps <- function(known_ps, n) {
  if (is.null(known_ps)) {
    rlang::abort(
      c(paste0("eval_nuisance$propensity = 'known' requires a numeric ",
               "`eval_propensity_score` (scalar or length-n vector)."),
        i = "Supply eval_propensity_score = 0.5 for a randomised trial, or a length-n vector of known scores."),
      class = "cfomics_eval_known_ps_missing"
    )
  }
  known_ps <- as.numeric(known_ps)
  if (length(known_ps) == 1L) {
    out <- rep(known_ps, n)
  } else if (length(known_ps) == n) {
    out <- known_ps
  } else {
    rlang::abort(
      sprintf(paste0("`eval_propensity_score` length (%d) must be 1 ",
                     "(scalar) or n (= %d)."),
              length(known_ps), n),
      class = "cfomics_eval_known_ps_length_mismatch"
    )
  }
  if (any(!is.finite(out) | out <= 0 | out >= 1)) {
    rlang::abort(
      "`eval_propensity_score` must be in (0, 1).",
      class = "cfomics_eval_known_ps_out_of_range"
    )
  }
  pmin(pmax(out, 1e-4), 1 - 1e-4)
}


#' Common-nuisance K-fold OOF estimates (mu0, mu1, pi).
#'
#' Fits shared nuisances (mu_0, mu_1, e) via K-fold OOF. Used both by
#' cf_compare (selection loss) and cf_evaluate (DR pseudo-outcome).
#'
#' \code{eval_nuisance} may be a list \code{list(outcome=, propensity=)} so
#' the propensity and outcome models can be specified independently. When
#' \code{propensity = "known"} the \code{eval_propensity_score} argument is
#' used instead of fitting a PS model (RCT context).
#'
#' @param X Numeric covariate matrix.
#' @param T Binary treatment vector.
#' @param Y Continuous outcome vector.
#' @param folds List of length K of test-fold indices (1-based).
#' @param eval_nuisance Base learner name OR
#'   \code{list(outcome=, propensity=)}. \code{propensity = "known"} uses
#'   \code{eval_propensity_score}.
#' @param random_state Integer seed or NULL.
#' @param eval_propensity_score Scalar in (0,1) or length-n vector; required
#'   when \code{eval_nuisance$propensity = "known"}, otherwise ignored.
#' @return List with elements \code{mu0_oof}, \code{mu1_oof}, \code{pi_oof}
#'   (numeric vectors of length \code{n}).
#' @keywords internal
#' @noRd
.compute_common_nuisance_oof <- function(X, T, Y, folds,
                                          eval_nuisance = "rf",
                                          random_state = NULL,
                                          eval_propensity_score = NULL,
                                          binary_outcome = FALSE) {
  spec <- .cfomics_resolve_eval_nuisance(eval_nuisance)
  # Binary outcome -> mu0 / mu1 OOF nuisances are probability classifiers.
  otype <- if (isTRUE(binary_outcome)) "classification" else "regression"
  n <- length(Y)
  ## factor / character / logical cols in X (e.g., colData(SE)$sex factor
  ## M/F) must be expanded to one-hot numeric before `.coerce_to_matrix`
  ## inside `.fit_base_learner`. Skip the warning emitted by
  ## `.expand_factor_columns` (single-level / high-cardinality) since this
  ## helper is called from the cf_compare oof path, not directly from user
  ## code.
  if (is.data.frame(X)) {
    X_mat <- suppressWarnings(.expand_factor_columns(X)$X)
  } else {
    X_mat <- as.matrix(X)
  }
  if (!is.numeric(X_mat)) {
    storage.mode(X_mat) <- "double"
  }
  T_vec <- as.integer(T)
  Y_vec <- as.numeric(Y)

  # Save/restore the caller's `.Random.seed` (consistent with the sister
  # helper `.compute_full_data_nuisance`) to prevent RNG-state leakage. The
  # caller's RNG state is preserved even when invoked via cf_compare, keeping
  # reproducibility intact.
  if (!is.null(random_state)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  mu0_oof <- numeric(n)
  mu1_oof <- numeric(n)
  pi_oof  <- numeric(n)
  m_oof   <- numeric(n)   # marginal outcome regression E[Y|X], for the R-loss

  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(n), test_idx)

    X_train <- X_mat[train_idx, , drop = FALSE]
    Y_train <- Y_vec[train_idx]
    T_train <- T_vec[train_idx]
    X_test  <- X_mat[test_idx,  , drop = FALSE]

    # propensity: "known" uses the user-supplied eval_propensity_score.
    if (identical(spec$propensity, "known")) {
      known_ps_full <- .cfomics_broadcast_known_ps(eval_propensity_score, n)
      pi_oof[test_idx] <- known_ps_full[test_idx]
    } else {
      ps_model <- .fit_base_learner(X_train, T_train, type = "classification",
                                     learner = spec$propensity,
                                     random_state = random_state)
      pi_test <- .predict_base_learner(ps_model, newdata = X_test)
      pi_oof[test_idx] <- pmin(pmax(pi_test, 1e-4), 1 - 1e-4)
    }

    # outcome regressions per arm
    idx_t0 <- which(T_train == 0L)
    idx_t1 <- which(T_train == 1L)
    if (length(idx_t0) > 1L) {
      m0 <- .fit_base_learner(X_train[idx_t0, , drop = FALSE],
                              Y_train[idx_t0], type = otype,
                              learner = spec$outcome,
                              random_state = random_state)
      mu0_oof[test_idx] <- .predict_base_learner(m0, newdata = X_test)
    }
    if (length(idx_t1) > 1L) {
      m1 <- .fit_base_learner(X_train[idx_t1, , drop = FALSE],
                              Y_train[idx_t1], type = otype,
                              learner = spec$outcome,
                              random_state = random_state)
      mu1_oof[test_idx] <- .predict_base_learner(m1, newdata = X_test)
    }

    # marginal outcome regression E[Y|X] on the pooled fold (both arms), used
    # as the R-loss residualiser. Estimated directly rather than reconstructed
    # from the per-arm models so it does not inherit the propensity estimate.
    m_fit <- .fit_base_learner(X_train, Y_train, type = otype,
                               learner = spec$outcome,
                               random_state = random_state)
    m_oof[test_idx] <- .predict_base_learner(m_fit, newdata = X_test)
  }

  list(mu0_oof = mu0_oof, mu1_oof = mu1_oof, pi_oof = pi_oof, m_oof = m_oof)
}


#' Compute full-data (in-sample) evaluation nuisance.
#'
#' When partition = "none", `eval_nuisance` is fit once on the full data and
#' the in-sample predictions populate
#' `fit$nuisance = list(mu0_oof, mu1_oof, pi_oof)`. cf_evaluate / cf_interpret
#' use this in-sample nuisance (with a warning).
#'
#' Field names reuse the K-fold `*_oof` convention (for a unified downstream
#' interface), although the values are in-sample.
#'
#' @param X Numeric matrix (n × p).
#' @param T Binary treatment vector (0/1).
#' @param Y Numeric outcome vector.
#' @param eval_nuisance Base-learner name (default `"rf"`).
#' @param random_state Optional integer seed; the caller's RNG state is
#'   saved/restored inside `.fit_base_learner`.
#'
#' @return List(mu0_oof, mu1_oof, pi_oof) — full-data in-sample predictions.
#'
#' @keywords internal
#' @noRd
.compute_full_data_nuisance <- function(X, T, Y, eval_nuisance = "rf",
                                        random_state = NULL,
                                        eval_propensity_score = NULL,
                                        binary_outcome = FALSE) {
  spec <- .cfomics_resolve_eval_nuisance(eval_nuisance)
  # Binary outcome -> the mu0 / mu1 nuisance regressions become probability
  # classifiers (P(Y = 1 | X, T = t)), matching the outcome scale.
  otype <- if (isTRUE(binary_outcome)) "classification" else "regression"
  n <- length(Y)
  ## Pre-expand factor / character columns to one-hot before bare
  ## `as.matrix()`. This is the partition='none' fallback path called by
  ## cf_fit -> cf_evaluate when no K-fold nuisance is available. Matches the
  ## `.compute_common_nuisance_oof` pre-expansion, so DR pseudo-outcome /
  ## GATES / BLP loss helpers see a numeric matrix.
  X_mat <- suppressWarnings(.expand_factor_columns(X)$X)
  T_vec <- as.integer(T)
  Y_vec <- as.numeric(Y)

  # Save/restore the caller's `.Random.seed` to prevent RNG-state leakage
  # when this helper is invoked via cf_fit's partition='none' path.
  if (!is.null(random_state)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  # propensity: "known" uses the user-supplied eval_propensity_score.
  if (identical(spec$propensity, "known")) {
    pi_in <- .cfomics_broadcast_known_ps(eval_propensity_score, n)
  } else {
    ps_model <- .fit_base_learner(X_mat, T_vec, type = "classification",
                                   learner = spec$propensity,
                                   random_state = random_state)
    pi_in <- .predict_base_learner(ps_model, newdata = X_mat)
    pi_in <- pmin(pmax(pi_in, 1e-4), 1 - 1e-4)
  }

  mu0_in <- numeric(n)
  mu1_in <- numeric(n)
  idx_t0 <- which(T_vec == 0L)
  idx_t1 <- which(T_vec == 1L)
  if (length(idx_t0) > 1L) {
    m0 <- .fit_base_learner(X_mat[idx_t0, , drop = FALSE], Y_vec[idx_t0],
                            type = otype, learner = spec$outcome,
                            random_state = random_state)
    mu0_in <- .predict_base_learner(m0, newdata = X_mat)
  }
  if (length(idx_t1) > 1L) {
    m1 <- .fit_base_learner(X_mat[idx_t1, , drop = FALSE], Y_vec[idx_t1],
                            type = otype, learner = spec$outcome,
                            random_state = random_state)
    mu1_in <- .predict_base_learner(m1, newdata = X_mat)
  }

  list(mu0_oof = mu0_in, mu1_oof = mu1_in, pi_oof = pi_in)
}


#' DR pseudo-outcome from OOF nuisance + observed Y/T.
#' @keywords internal
#' @noRd
.compute_dr_pseudo_outcome <- function(Y, T, mu0_oof, mu1_oof, pi_oof) {
  T_vec <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  mu1_oof + (T_vec / pi_oof) * (Y_vec - mu1_oof) -
    (mu0_oof + ((1 - T_vec) / (1 - pi_oof)) * (Y_vec - mu0_oof))
}

#' R-loss (Nie-Wager 2021) per observation.
#'
#' \code{r_loss_i = ((Y_i - m_oof(X_i)) - tau_hat(X_i) * (T_i - pi_oof(X_i)))^2}
#' where \code{m_oof(X_i)} is the out-of-fold marginal outcome regression
#' \eqn{E[Y|X]}.
#' @keywords internal
#' @noRd
.compute_r_loss_per_obs <- function(Y, T, m_oof, pi_oof, tau_hat) {
  Yres <- as.numeric(Y) - m_oof
  Tres <- as.numeric(T) - pi_oof
  (Yres - tau_hat * Tres)^2
}

#' DR-loss per observation: (tau_hat - Y_DR)^2 = the DR pseudo-outcome
#' regression loss (Kennedy 2023, \emph{EJS} 17(2):3008-3049); used for CATE
#' model selection as in Saito & Yasui 2020 (Counterfactual Cross-Validation).
#' @keywords internal
#' @noRd
.compute_dr_loss_per_obs <- function(Y, T, mu0_oof, mu1_oof, pi_oof, tau_hat) {
  Y_DR <- .compute_dr_pseudo_outcome(Y, T, mu0_oof, mu1_oof, pi_oof)
  (tau_hat - Y_DR)^2
}

#' Factual loss per observation: (Y - mu_T(X))^2 using mu_t OOF.
#' @keywords internal
#' @noRd
.compute_factual_loss_per_obs <- function(Y, T, mu0_oof, mu1_oof) {
  T_vec <- as.integer(T)
  mu_T  <- ifelse(T_vec == 1L, mu1_oof, mu0_oof)
  (as.numeric(Y) - mu_T)^2
}

#' Aggregate per-obs loss into per-fold means.
#' @keywords internal
#' @noRd
.aggregate_loss_by_fold <- function(per_obs, fold_id) {
  out <- vapply(sort(unique(fold_id)), function(k) {
    mean(per_obs[fold_id == k], na.rm = TRUE)
  }, numeric(1))
  out
}

#' Between-fold variance t-CI for a loss vector of length K.
#' @keywords internal
#' @noRd
.between_fold_ci <- function(per_fold_loss, alpha = 0.05) {
  K <- length(per_fold_loss)
  value <- mean(per_fold_loss)
  if (K < 2L) {
    return(list(value = value, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  se  <- stats::sd(per_fold_loss) / sqrt(K)
  tcrit <- stats::qt(1 - alpha / 2, df = K - 1)
  list(
    value    = value,
    ci_lower = value - tcrit * se,
    ci_upper = value + tcrit * se
  )
}

#' R-squared analog: 1 - loss / base_loss.
#' \code{base_loss} is the loss under a constant-CATE model (theta = ATE_OOF).
#' @keywords internal
#' @noRd
.compute_r_squared_score <- function(per_obs, base_per_obs) {
  # Numerator and denominator are averaged over the SAME finite observations so
  # the ratio is well-defined even when per_obs carries method-specific NAs
  # (e.g. a method that only surfaces potential outcomes on some folds). This is
  # a no-op when per_obs is fully finite (the usual case).
  ok <- is.finite(per_obs) & is.finite(base_per_obs)
  if (!any(ok)) return(NA_real_)
  loss      <- mean(per_obs[ok])
  base_loss <- mean(base_per_obs[ok])
  if (!is.finite(base_loss) || base_loss <= 0) return(NA_real_)
  1 - loss / base_loss
}


# -------------------------------------------------------------------------
# Per-observation loss recomputation helpers
# -------------------------------------------------------------------------
#
# cf_compare stores per-method scalar losses in `scoreboard` but discards
# per-obs residuals / DR scores after aggregation. cf_ensemble's
# r_stacking / dr_stacking require per-obs material (since the loss of the
# weighted-average tau is NOT the weighted average of per-method losses).
# Cf_select with non-scoreboard metrics also benefits from on-the-fly
# recomputation.
#
# The required ingredients are all already cached in cf_compare's output:
#
#   * `cmp$metadata$nuisance_oof` = list(mu0_oof, mu1_oof, pi_oof)
#   * `cmp$models[[m]]$data$Y` / `$T` (any surviving model is enough)
#   * `cmp$metadata$oof_tau_matrix` = n x K tau hats
#
# All recomputation is pure arithmetic (no nuisance refit) — this is the
# entire point of these helpers.


#' Extract the cf_compare cache needed by per-obs loss recomputation.
#'
#' Pulls Y / T from the first surviving model, nuisance_oof from metadata,
#' and exposes the oof_tau_matrix. Guards against the edge case where every
#' method failed (no model retains $data$Y / $data$T).
#'
#' @param compare_result A `cfomics_compare` object.
#' @return A list with `Y`, `T`, `mu0_oof`, `mu1_oof`, `pi_oof`,
#'   `oof_tau_matrix`, and `methods` (column names of oof_tau_matrix).
#'   Returns `NULL` with a class-tagged error when the cache is unusable.
#' @keywords internal
#' @noRd
.cfomics_extract_compare_cache <- function(compare_result) {
  if (!inherits(compare_result, "cfomics_compare")) {
    rlang::abort(
      "compare_result must be a cfomics_compare object.",
      class = "cfomics_invalid_compare_result"
    )
  }

  meta <- compare_result$metadata
  nuis <- meta$nuisance_oof
  if (is.null(nuis) ||
      is.null(nuis$mu0_oof) || is.null(nuis$mu1_oof) || is.null(nuis$pi_oof)) {
    rlang::abort(
      paste0("cf_compare metadata is missing nuisance_oof; ",
             "ensure partition = 'oof' was used."),
      class = "cfomics_missing_nuisance_oof"
    )
  }

  oof_tau_matrix <- meta$oof_tau_matrix
  if (is.null(oof_tau_matrix) || !is.matrix(oof_tau_matrix)) {
    rlang::abort(
      "cf_compare metadata is missing oof_tau_matrix.",
      class = "cfomics_missing_oof_tau_matrix"
    )
  }

  # Find any surviving model that still carries the cached training pair.
  models <- compare_result$models
  YT <- NULL
  if (is.list(models) && length(models) > 0L) {
    for (m in names(models)) {
      mdl <- models[[m]]
      if (!is.null(mdl) && is.list(mdl$data) &&
          !is.null(mdl$data$Y) && !is.null(mdl$data$T)) {
        YT <- list(Y = mdl$data$Y, T = mdl$data$T)
        break
      }
    }
  }

  if (is.null(YT)) {
    rlang::abort(
      paste0("cf_compare cache has no surviving model with Y/T; ",
             "all methods may have failed."),
      class = "cfomics_compare_cache_empty"
    )
  }

  list(
    Y              = as.numeric(YT$Y),
    T              = as.integer(YT$T),
    mu0_oof        = as.numeric(nuis$mu0_oof),
    mu1_oof        = as.numeric(nuis$mu1_oof),
    pi_oof         = as.numeric(nuis$pi_oof),
    # m_oof must be the direct OOF regression E[Y|X] (Work B). Never fall back
    # to mu1*pi + mu0*(1-pi): that propensity-weighted reconstruction is the
    # residualiser R-loss canonicalization deliberately removed. The shared eval
    # nuisance always supplies m_oof, so a missing value is a contract breach.
    m_oof          = if (!is.null(nuis$m_oof)) as.numeric(nuis$m_oof)
                     else rlang::abort(
                       "cf_compare cache nuisance is missing m_oof (E[Y|X]).",
                       class = "cfomics_missing_m_oof"),
    oof_tau_matrix = oof_tau_matrix,
    methods        = colnames(oof_tau_matrix)
  )
}


#' Compute method-independent stacking materials.
#'
#' Returns the per-observation Y / T residuals (R-stacking core) and the
#' DR pseudo-outcome (DR-stacking core). These vectors do NOT depend on
#' which CATE method is being evaluated, so they only need to be computed
#' once per cf_compare result and can be reused across stacking solvers.
#'
#' For r_stacking, the objective is:
#'
#'   min_w  sum_i ((Yres_i) - (T_res_i) * sum_k w_k tau_k(x_i))^2
#'
#' which factors as a least-squares problem in the rows of the
#' oof_tau_matrix scaled by Tres, against the target Yres.
#'
#' For dr_stacking, the objective is:
#'
#'   min_w  sum_i (sum_k w_k tau_k(x_i) - Y_DR_i)^2
#'
#' i.e. ordinary least squares with target Y_DR.
#'
#' @param compare_result A `cfomics_compare` object.
#' @return A list with `Yres`, `Tres`, `Y_DR`, plus the underlying
#'   `Y`, `T`, `mu0_oof`, `mu1_oof`, `pi_oof`, and `oof_tau_matrix` for
#'   convenience.
#' @keywords internal
#' @noRd
.cfomics_compute_stacking_materials <- function(compare_result) {
  cache <- .cfomics_extract_compare_cache(compare_result)

  Yres <- cache$Y - cache$m_oof
  Tres <- as.numeric(cache$T) - cache$pi_oof
  Y_DR <- .compute_dr_pseudo_outcome(cache$Y, cache$T,
                                     cache$mu0_oof, cache$mu1_oof,
                                     cache$pi_oof)

  c(cache, list(
    Yres   = Yres,
    Tres   = Tres,
    Y_DR   = Y_DR
  ))
}


#' Compute per-observation loss for a method from cf_compare cache.
#'
#' Thin wrapper that selects `tau_hat` for `method` from
#' `compare_result$metadata$oof_tau_matrix`, extracts the cached
#' Y / T / nuisance via `.cfomics_extract_compare_cache`, and dispatches to
#' the underlying per-obs primitive (`.compute_r_loss_per_obs`,
#' `.compute_dr_loss_per_obs`, `.compute_factual_loss_per_obs`).
#'
#' Used by `cf_select` (when the requested metric is not in `scoreboard`)
#' and by `cf_ensemble` (softmax / stacking / cba for legacy compatibility
#' paths). Pure arithmetic, no nuisance refit.
#'
#' @param compare_result A `cfomics_compare` object.
#' @param method Character; one of `colnames(oof_tau_matrix)`.
#' @param loss_type Character; one of `"r_loss"`, `"dr_loss"`,
#'   `"factual_loss"`.
#' @return Numeric vector of length `n` with the per-observation loss.
#' @keywords internal
#' @noRd
.cfomics_compute_per_obs_loss_from_compare <- function(compare_result, method,
                                                      loss_type = c("r_loss",
                                                                    "dr_loss",
                                                                    "factual_loss")) {
  loss_type <- match.arg(loss_type)
  cache     <- .cfomics_extract_compare_cache(compare_result)

  if (!is.character(method) || length(method) != 1L) {
    rlang::abort("`method` must be a single character string.",
                 class = "cfomics_invalid_method")
  }
  if (!method %in% cache$methods) {
    rlang::abort(
      sprintf("method '%s' not found in oof_tau_matrix columns: %s",
              method, paste(cache$methods, collapse = ", ")),
      class = "cfomics_method_not_in_compare"
    )
  }

  tau_hat <- cache$oof_tau_matrix[, method]

  switch(loss_type,
    r_loss       = .compute_r_loss_per_obs(cache$Y, cache$T,
                                           cache$m_oof, cache$pi_oof, tau_hat),
    dr_loss      = .compute_dr_loss_per_obs(cache$Y, cache$T,
                                            cache$mu0_oof, cache$mu1_oof,
                                            cache$pi_oof, tau_hat),
    factual_loss = .compute_factual_loss_per_obs(cache$Y, cache$T,
                                                 cache$mu0_oof, cache$mu1_oof)
  )
}


#' Compute per-observation loss for an *ensemble* tau hat from cf_compare.
#'
#' Convenience wrapper for cf_ensemble: given a length-n combined tau hat
#' vector (e.g. weighted average over methods), compute the per-obs loss
#' using the cached nuisance + Y/T. Mirrors
#' `.cfomics_compute_per_obs_loss_from_compare` but takes the tau vector
#' directly so the caller can plug in an ensemble tau without first having
#' to register a synthetic method column.
#'
#' @param compare_result A `cfomics_compare` object.
#' @param tau_hat Numeric vector of length `n` (ensemble tau hat).
#' @param loss_type Character; `"r_loss"` / `"dr_loss"` / `"factual_loss"`.
#' @return Numeric vector of per-obs losses (length `n`).
#' @keywords internal
#' @noRd
.cfomics_compute_per_obs_loss_for_tau <- function(compare_result, tau_hat,
                                                  loss_type = c("r_loss",
                                                                "dr_loss",
                                                                "factual_loss")) {
  loss_type <- match.arg(loss_type)
  cache     <- .cfomics_extract_compare_cache(compare_result)

  tau_hat <- as.numeric(tau_hat)
  n <- length(cache$Y)
  if (length(tau_hat) != n) {
    rlang::abort(
      sprintf("tau_hat length (%d) does not match cf_compare n (%d).",
              length(tau_hat), n),
      class = "cfomics_tau_length_mismatch"
    )
  }

  switch(loss_type,
    r_loss       = .compute_r_loss_per_obs(cache$Y, cache$T,
                                           cache$m_oof, cache$pi_oof, tau_hat),
    dr_loss      = .compute_dr_loss_per_obs(cache$Y, cache$T,
                                            cache$mu0_oof, cache$mu1_oof,
                                            cache$pi_oof, tau_hat),
    factual_loss = .compute_factual_loss_per_obs(cache$Y, cache$T,
                                                 cache$mu0_oof, cache$mu1_oof)
  )
}
