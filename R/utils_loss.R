#' Resolve eval_nuisance into a list(outcome=, propensity=) spec.
#'
#' eval_nuisance can be a single string ("super_learner", "rf", etc.) or a named
#' list with elements \code{outcome} and \code{propensity}. A string is
#' broadcast: list(outcome = s, propensity = s).
#' When a list is provided, missing elements fall back to the partner. The
#' special value \code{propensity = "known"} signals that an externally
#' supplied propensity score should be used in place of a fitted PS model
#' (RCT context).
#'
#' Allowed learner names match \code{.fit_base_learner()} (linear, logistic,
#' lasso, ridge, elastic_net, rf, xgboost), plus \code{"super_learner"} (the
#' default, a weighted combination of \code{eval_library}) and the sentinel
#' "known" for the propensity slot only. Validation is strict.
#'
#' @keywords internal
#' @noRd
.cateomics_resolve_eval_nuisance <- function(eval_nuisance) {
  allowed <- .CF_EVAL_LEARNERS
  if (is.character(eval_nuisance) && length(eval_nuisance) == 1L &&
      !is.null(eval_nuisance) && !is.na(eval_nuisance)) {
    if (!(eval_nuisance %in% allowed)) {
      rlang::abort(
        sprintf("eval_nuisance '%s' is not in the allowed pool: %s",
                eval_nuisance, paste(allowed, collapse = ", ")),
        class = "cateomics_eval_nuisance_unknown"
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
        class = "cateomics_eval_nuisance_empty_list"
      )
    }
    # Broadcast the sibling, but "known" is only valid on the propensity
    # slot, so it is never copied into the outcome slot. The public entries
    # refuse `"known"` as a learner name, so this arm is reached only by the
    # spec the entries build from `eval_propensity_score`.
    if (is.null(out)) out <- if (identical(pro, "known")) .CF_SL_TOKEN else pro
    if (is.null(pro)) pro <- out
    if (!(out %in% allowed)) {
      rlang::abort(
        sprintf("eval_nuisance$outcome '%s' is not allowed: %s",
                out, paste(allowed, collapse = ", ")),
        class = "cateomics_eval_nuisance_unknown"
      )
    }
    if (!(pro %in% c(allowed, "known"))) {
      rlang::abort(
        sprintf(paste0("eval_nuisance$propensity '%s' is not allowed: ",
                       "%s or 'known'."),
                pro, paste(allowed, collapse = ", ")),
        class = "cateomics_eval_nuisance_unknown"
      )
    }
    return(list(outcome = out, propensity = pro))
  }
  rlang::abort(
    "eval_nuisance must be a single learner name or a named list.",
    class = "cateomics_eval_nuisance_invalid_type"
  )
}

#' Validate the public `eval_nuisance` argument.
#'
#' Two shapes are public, matching what `eval_library` and `eval_params`
#' already accept: a single learner name from `.CF_EVAL_LEARNERS`, applying to
#' both slots, or `list(outcome = , propensity = )` naming a learner per slot,
#' with a missing slot falling back to its partner. A known propensity is
#' expressed by passing a value to `eval_propensity_score`, so the `"known"`
#' sentinel stays internal: a learner name carries no score to use.
#'
#' Called from the `cf_fit` and `cf_compare` entries so both refuse the same
#' values with the same condition.
#'
#' @keywords internal
#' @noRd
.cf_check_eval_nuisance <- function(eval_nuisance) {
  bad <- function(what, got) {
    rlang::abort(
      c(sprintf("%s must be one of: %s.", what,
                paste(shQuote(.CF_EVAL_LEARNERS), collapse = ", ")),
        x = sprintf("Received: %s.", paste(deparse(got), collapse = " ")),
        i = paste("Pass one name for both slots, or",
                  "list(outcome = , propensity = ) for one each."),
        i = "A known propensity goes to `eval_propensity_score` as a value."),
      class = "cateomics_invalid_eval_nuisance"
    )
  }
  if (is.list(eval_nuisance)) {
    slots <- c("outcome", "propensity")
    nms <- names(eval_nuisance)
    if (is.null(nms) || !length(eval_nuisance) ||
        !all(nms %in% slots) || anyDuplicated(nms)) {
      rlang::abort(
        c("`eval_nuisance` list must be named `outcome` and / or `propensity`.",
          x = sprintf("Received names: %s.",
                      paste(deparse(nms), collapse = " "))),
        class = "cateomics_invalid_eval_nuisance"
      )
    }
    for (sl in intersect(slots, nms)) {
      v <- eval_nuisance[[sl]]
      if (!is.character(v) || length(v) != 1L || !(v %in% .CF_EVAL_LEARNERS)) {
        bad(sprintf("`eval_nuisance$%s`", sl), v)
      }
    }
    return(invisible(eval_nuisance))
  }
  if (!is.character(eval_nuisance) || length(eval_nuisance) != 1L ||
      !(eval_nuisance %in% .CF_EVAL_LEARNERS)) {
    bad("`eval_nuisance`", eval_nuisance)
  }
  invisible(eval_nuisance)
}

#' Build the internal evaluation-nuisance spec.
#'
#' A numeric `eval_propensity_score` switches the propensity slot to the
#' supplied value instead of fitting a model; `NULL` (the default) leaves the
#' spec as the caller gave it, so the estimated nuisance is unchanged. The
#' outcome slot always keeps the learner the caller named, whether that came
#' from a single name or from the `outcome` slot of a list.
#'
#' Naming a propensity learner and supplying a known score at the same time is
#' a contradiction, the same one `propensity_score` and `model_propensity`
#' refuse on the estimator side, so it aborts rather than one silently winning.
#'
#' @keywords internal
#' @noRd
.cf_eval_nuisance_spec <- function(eval_nuisance, eval_propensity_score) {
  if (is.null(eval_propensity_score)) return(eval_nuisance)
  if (is.list(eval_nuisance)) {
    if (!is.null(eval_nuisance$propensity)) {
      rlang::abort(
        c(paste("`eval_nuisance$propensity` and `eval_propensity_score`",
                "cannot both be given."),
          i = paste("Drop the propensity slot to use the known score, or",
                    "drop `eval_propensity_score` to fit the learner.")),
        class = "cateomics_eval_propensity_args_conflict"
      )
    }
    out <- eval_nuisance$outcome %||% .CF_SL_TOKEN
    return(list(outcome = out, propensity = "known"))
  }
  list(outcome = eval_nuisance, propensity = "known")
}

#' Report a known propensity supplied to one side only.
#'
#' `propensity_score` is the estimator-side propensity and
#' `eval_propensity_score` the one the proxy losses and the evaluation
#' metrics read. They are separate settings, so supplying one leaves the other
#' estimated; this says so once rather than letting the asymmetry be silent.
#'
#' @keywords internal
#' @noRd
.cf_inform_known_ps_asymmetry <- function(propensity_score,
                                          eval_propensity_score) {
  has_fit  <- !is.null(propensity_score)
  has_eval <- !is.null(eval_propensity_score)
  if (has_fit == has_eval) return(invisible(NULL))
  if (has_fit) {
    rlang::inform(
      c("`propensity_score` is set but `eval_propensity_score` is not.",
        i = paste("The evaluation nuisance still estimates a propensity.",
                  "Pass `eval_propensity_score` as well to use the known",
                  "score in the losses and evaluation metrics.")),
      class = "cateomics_known_ps_one_side"
    )
  } else {
    rlang::inform(
      c("`eval_propensity_score` is set but `propensity_score` is not.",
        i = paste("The estimator still estimates its own propensity.",
                  "Pass `propensity_score` as well to use the known score",
                  "in the estimator.")),
      class = "cateomics_known_ps_one_side"
    )
  }
  invisible(NULL)
}



#' Broadcast a known propensity score to length n with finite clipping.
#' @keywords internal
#' @noRd
.cateomics_broadcast_known_ps <- function(known_ps, n) {
  if (is.null(known_ps)) {
    rlang::abort(
      c(paste0("eval_nuisance$propensity = 'known' requires a numeric ",
               "`eval_propensity_score` (scalar or length-n vector)."),
        i = "Supply eval_propensity_score = 0.5 for a randomised trial, or a length-n vector of known scores."),
      class = "cateomics_eval_known_ps_missing"
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
      class = "cateomics_eval_known_ps_length_mismatch"
    )
  }
  if (any(!is.finite(out) | out <= 0 | out >= 1)) {
    rlang::abort(
      "`eval_propensity_score` must be in (0, 1).",
      class = "cateomics_eval_known_ps_out_of_range"
    )
  }
  .cf_clip_ps(out)
}


# The propensity clip bounds. `cf_fit` / `cf_compare` leave the argument NULL
# so the dispatcher can tell "the caller set it" from "the caller did not",
# which is what lets an unreachable `clip` be reported; the value itself lives
# here and in the formals of the four fitters that bound a propensity.
.CF_PS_CLIP <- c(0.025, 0.975)

#' Validate the public `clip` argument.
#'
#' Called from the `cf_fit` and `cf_compare` entries so both refuse the same
#' values: `clip` bounds the estimator-side propensity, and cf_compare forwards
#' it to every config.
#'
#' @keywords internal
#' @noRd
.cf_check_clip <- function(clip) {
  if (is.null(clip)) return(invisible(NULL))
  if (!is.numeric(clip) || length(clip) != 2L ||
      clip[1] >= clip[2] || clip[1] <= 0 || clip[2] >= 1) {
    rlang::abort(
      sprintf("`clip` must be length-2 numeric with 0 < clip[1] < clip[2] < 1 (got: %s).",
              paste(deparse(clip), collapse = " ")),
      class = "cateomics_invalid_clip"
    )
  }
  invisible(clip)
}

#' Keep a propensity away from 0 and 1
#'
#' The evaluation propensity carries this absolute floor only; the
#' estimator-side `clip` is not applied to it.
#'
#' @keywords internal
#' @noRd
.cf_clip_ps <- function(p, clip = c(1e-4, 1 - 1e-4)) {
  lo <- max(as.numeric(clip)[1L], 1e-4)
  hi <- min(as.numeric(clip)[2L], 1 - 1e-4)
  pmin(pmax(p, lo), hi)
}


#' Common-nuisance K-fold OOF estimates (mu0, mu1, pi).
#'
#' Fits shared nuisances (mu_0, mu_1, e) via K-fold OOF. Called by
#' cf_compare for the selection loss and by the fold-wise path in
#' `utils_kfold.R`; the downstream evaluators read what these produce rather
#' than refitting.
#'
#' \code{eval_nuisance} may be a list \code{list(outcome=, propensity=)} so
#' the propensity and outcome models can be specified independently. When
#' \code{propensity = "known"} the \code{eval_propensity_score} argument is
#' used instead of fitting a PS model (RCT context).
#'
#' @param X Numeric covariate matrix.
#' @param T Binary treatment vector.
#' @param Y Outcome vector; a class label when `binary_outcome` is TRUE.
#' @param folds List of length K of test-fold indices (1-based).
#' @param eval_nuisance Base learner name OR
#'   \code{list(outcome=, propensity=)}. \code{propensity = "known"} uses
#'   \code{eval_propensity_score}.
#' @param random_state Integer seed or NULL.
#' @param eval_propensity_score Scalar in (0,1) or length-n vector; required
#'   when \code{eval_nuisance$propensity = "known"}, otherwise ignored.
#' @param propensity_X Character vector naming the covariates the propensity
#'   slot is fitted on, or NULL for all of them. The outcome slots always see
#'   the full design matrix.
#' @return List with elements \code{mu0_oof}, \code{mu1_oof}, \code{pi_oof}
#'   and \code{m_oof} (numeric vectors of length \code{n}; \code{m_oof} is
#'   the pooled OOF regression of Y on X that the R-loss residualises with),
#'   plus \code{sl_weights} when the Super Learner was used.
#' @keywords internal
#' @noRd
.compute_common_nuisance_oof <- function(X, T, Y, folds,
                                          eval_nuisance = "super_learner",
                                          random_state = NULL,
                                          eval_propensity_score = NULL,
                                          binary_outcome = FALSE,
                                          eval_params = list(),
                                          eval_library = NULL,
                                          propensity_X = NULL) {
  spec <- .cateomics_resolve_eval_nuisance(eval_nuisance)
  out_is_sl <- .cf_sl_is_token(spec$outcome)
  ps_is_sl  <- .cf_sl_is_token(spec$propensity)
  # The outcome slot is fitted on the outcome scale and the propensity slot is
  # always a classifier, so each is validated against the arguments its own
  # fitter reads.
  .otype_for_params <- if (isTRUE(binary_outcome)) "classification" else "regression"
  .sl_out_lib <- if (out_is_sl) .cf_sl_resolve_library(eval_library, "outcome") else NULL
  .sl_ps_lib  <- if (ps_is_sl)  .cf_sl_resolve_library(eval_library, "propensity") else NULL
  .eval_out_p <- if (out_is_sl) {
    .cf_eval_params_members(eval_params, "outcome", .sl_out_lib,
                            type = .otype_for_params)
  } else {
    .cf_eval_params_slot(eval_params, "outcome", spec$outcome,
                         type = .otype_for_params)
  }
  .eval_ps_p  <- if (ps_is_sl) {
    .cf_eval_params_members(eval_params, "propensity", .sl_ps_lib,
                            type = "classification")
  } else {
    .cf_eval_params_slot(eval_params, "propensity", spec$propensity,
                         type = "classification")
  }
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

  # A missing seed reads as 0 here, the same as in the dispatch.
  random_state <- as.integer(random_state %||% 0L)
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

  # `propensity_X` restricts the propensity slot to a subset of the covariates,
  # as it does on the estimator side; the outcome slots always see all of X.
  # NULL (the default) returns X unchanged, so the estimates do not move.
  X_ps <- .resolve_propensity_X(X_mat, propensity_X, colnames(X_mat))$X_ps

  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(n), test_idx)

    X_train <- X_mat[train_idx, , drop = FALSE]
    Y_train <- Y_vec[train_idx]
    T_train <- T_vec[train_idx]
    X_test  <- X_mat[test_idx,  , drop = FALSE]
    X_ps_train <- X_ps[train_idx, , drop = FALSE]
    X_ps_test  <- X_ps[test_idx,  , drop = FALSE]

    # propensity: "known" uses the user-supplied eval_propensity_score.
    if (identical(spec$propensity, "known")) {
      known_ps_full <- .cateomics_broadcast_known_ps(eval_propensity_score, n)
      pi_oof[test_idx] <- known_ps_full[test_idx]
    } else if (!ps_is_sl) {
      ps_model <- .fit_base_learner(X_ps_train, T_train,
                                     type = "classification",
                                     learner = spec$propensity,
                                     random_state = random_state,
                                     .params = .eval_ps_p)
      pi_test <- .predict_base_learner(ps_model, newdata = X_ps_test)
      pi_oof[test_idx] <- .cf_clip_ps(pi_test)
    }

    if (!out_is_sl) {
      # outcome regressions per arm
      idx_t0 <- which(T_train == 0L)
      idx_t1 <- which(T_train == 1L)
      if (length(idx_t0) > 1L) {
        m0 <- .fit_base_learner(X_train[idx_t0, , drop = FALSE],
                                Y_train[idx_t0], type = otype,
                                learner = spec$outcome,
                                random_state = random_state,
                                .params = .eval_out_p)
        mu0_oof[test_idx] <- .predict_base_learner(m0, newdata = X_test)
      }
      if (length(idx_t1) > 1L) {
        m1 <- .fit_base_learner(X_train[idx_t1, , drop = FALSE],
                                Y_train[idx_t1], type = otype,
                                learner = spec$outcome,
                                random_state = random_state,
                                .params = .eval_out_p)
        mu1_oof[test_idx] <- .predict_base_learner(m1, newdata = X_test)
      }

      # marginal outcome regression E[Y|X] on the pooled fold (both arms), used
      # as the R-loss residualiser. Estimated directly rather than reconstructed
      # from the per-arm models so it does not inherit the propensity estimate.
      m_fit <- .fit_base_learner(X_train, Y_train, type = otype,
                                 learner = spec$outcome,
                                 random_state = random_state,
                                 .params = .eval_out_p)
      m_oof[test_idx] <- .predict_base_learner(m_fit, newdata = X_test)
    }
  }

  # The Super Learner needs the whole out-of-fold matrix before it can weight
  # its members, so those quantities are assembled after the loop rather than
  # inside it. Same folds, same held-out rows.
  sl_weights <- list()
  if (ps_is_sl) {
    r <- .cf_sl_oof_quantity(X_ps, T_vec, folds, .sl_ps_lib,
                             type = "classification",
                             random_state = random_state,
                             member_params = .eval_ps_p,
                             where = "evaluation propensity")
    pi_oof <- .cf_clip_ps(r$pred)
    sl_weights$pi <- r$weights
  }
  if (out_is_sl) {
    rows_t0 <- which(T_vec == 0L)
    rows_t1 <- which(T_vec == 1L)
    r0 <- .cf_sl_oof_quantity(X_mat, Y_vec, folds, .sl_out_lib, type = otype,
                              random_state = random_state,
                              member_params = .eval_out_p,
                              train_rows = rows_t0, weight_rows = rows_t0,
                              where = "evaluation outcome (control arm)")
    r1 <- .cf_sl_oof_quantity(X_mat, Y_vec, folds, .sl_out_lib, type = otype,
                              random_state = random_state,
                              member_params = .eval_out_p,
                              train_rows = rows_t1, weight_rows = rows_t1,
                              where = "evaluation outcome (treated arm)")
    rm_ <- .cf_sl_oof_quantity(X_mat, Y_vec, folds, .sl_out_lib, type = otype,
                               random_state = random_state,
                               member_params = .eval_out_p,
                               where = "evaluation outcome (pooled)")
    mu0_oof <- r0$pred
    mu1_oof <- r1$pred
    m_oof   <- rm_$pred
    sl_weights$mu0 <- r0$weights
    sl_weights$mu1 <- r1$weights
    sl_weights$m   <- rm_$weights
  }

  list(mu0_oof = mu0_oof, mu1_oof = mu1_oof, pi_oof = pi_oof, m_oof = m_oof,
       sl_weights = if (length(sl_weights)) sl_weights else NULL)
}


#' Compute full-data (in-sample) evaluation nuisance.
#'
#' When partition = "none" the in-sample predictions populate
#' `fit$res$nuisance = list(mu0_oof, mu1_oof, pi_oof, m_oof)`, which `cf_evaluate`
#' and `cf_rule` then read (with a warning). A single learner is fit once on
#' the full data; the Super Learner solves its weights on an internal split
#' and refits every member on all rows.
#'
#' Field names reuse the K-fold `*_oof` convention (for a unified downstream
#' interface), although the values are in-sample.
#'
#' @param X Numeric matrix (n × p).
#' @param T Binary treatment vector (0/1).
#' @param Y Numeric outcome vector.
#' @param eval_nuisance Base-learner name, or `"super_learner"` (the default).
#' @param random_state Optional integer seed; the caller's RNG state is
#'   saved/restored inside `.fit_base_learner`.
#' @param eval_library Super Learner library for this fit.
#' @param propensity_X Character vector naming the covariates the propensity
#'   slot is fitted on, or NULL for all of them. The outcome slots always see
#'   the full design matrix.
#' @param n_folds Fold count for the internal split the Super Learner weights
#'   are solved on. This path has no outer folds of its own.
#'
#' @return List(mu0_oof, mu1_oof, pi_oof, m_oof), plus `sl_weights` when the
#'   Super Learner was used — full-data in-sample predictions. Same four
#'   elements as the OOF helper, so a nuisance has one shape whichever
#'   partition produced it.
#'
#' @keywords internal
#' @noRd
.compute_full_data_nuisance <- function(X, T, Y, eval_nuisance = "super_learner",
                                        random_state = NULL,
                                        eval_propensity_score = NULL,
                                        binary_outcome = FALSE,
                                        eval_params = list(),
                                        eval_library = NULL,
                                        n_folds = 5L,
                                        propensity_X = NULL) {
  spec <- .cateomics_resolve_eval_nuisance(eval_nuisance)
  out_is_sl <- .cf_sl_is_token(spec$outcome)
  ps_is_sl  <- .cf_sl_is_token(spec$propensity)
  # The outcome slot is fitted on the outcome scale and the propensity slot is
  # always a classifier, so each is validated against the arguments its own
  # fitter reads.
  .otype_for_params <- if (isTRUE(binary_outcome)) "classification" else "regression"
  .sl_out_lib <- if (out_is_sl) .cf_sl_resolve_library(eval_library, "outcome") else NULL
  .sl_ps_lib  <- if (ps_is_sl)  .cf_sl_resolve_library(eval_library, "propensity") else NULL
  .eval_out_p <- if (out_is_sl) {
    .cf_eval_params_members(eval_params, "outcome", .sl_out_lib,
                            type = .otype_for_params)
  } else {
    .cf_eval_params_slot(eval_params, "outcome", spec$outcome,
                         type = .otype_for_params)
  }
  .eval_ps_p  <- if (ps_is_sl) {
    .cf_eval_params_members(eval_params, "propensity", .sl_ps_lib,
                            type = "classification")
  } else {
    .cf_eval_params_slot(eval_params, "propensity", spec$propensity,
                         type = "classification")
  }
  # Binary outcome -> the mu0 / mu1 nuisance regressions become probability
  # classifiers (P(Y = 1 | X, T = t)), matching the outcome scale.
  otype <- if (isTRUE(binary_outcome)) "classification" else "regression"
  n <- length(Y)
  ## Pre-expand factor / character columns to one-hot before bare
  ## `as.matrix()`. This is the `partition = "none"` path, which `cf_fit`
  ## calls instead of the K-fold helper.
  ## `.compute_common_nuisance_oof` pre-expansion, so DR pseudo-outcome /
  ## GATES / BLP loss helpers see a numeric matrix.
  X_mat <- suppressWarnings(.expand_factor_columns(X)$X)
  T_vec <- as.integer(T)
  Y_vec <- as.numeric(Y)
  # `propensity_X` restricts the propensity slot to a subset of the
  # covariates, as it does on the estimator side; the outcome slots always
  # see all of X. NULL (the default) returns X unchanged.
  X_ps  <- .resolve_propensity_X(X_mat, propensity_X, colnames(X_mat))$X_ps

  # A missing seed is read as 0, as the dispatch does: see the sister helper.
  random_state <- as.integer(random_state %||% 0L)
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

  sl_weights <- list()

  # propensity: "known" uses the user-supplied eval_propensity_score.
  if (identical(spec$propensity, "known")) {
    pi_in <- .cateomics_broadcast_known_ps(eval_propensity_score, n)
  } else if (ps_is_sl) {
    r <- .cf_sl_full_quantity(X_ps, T_vec, .sl_ps_lib, type = "classification",
                              strata = T_vec, n_folds = n_folds,
                              random_state = random_state,
                              member_params = .eval_ps_p,
                              where = "evaluation propensity")
    pi_in <- .cf_clip_ps(r$pred)
    sl_weights$pi <- r$weights
  } else {
    ps_model <- .fit_base_learner(X_ps, T_vec, type = "classification",
                                   learner = spec$propensity,
                                   random_state = random_state,
                                   .params = .eval_ps_p)
    pi_in <- .predict_base_learner(ps_model, newdata = X_ps)
    pi_in <- .cf_clip_ps(pi_in)
  }

  mu0_in <- numeric(n)
  mu1_in <- numeric(n)
  # The pooled regression E[Y|X] that the R-loss residualises with. The OOF
  # helper returns it, so the in-sample path returns it too and a nuisance has
  # the same four elements whichever partition produced it.
  m_in   <- numeric(n)
  idx_t0 <- which(T_vec == 0L)
  idx_t1 <- which(T_vec == 1L)
  # With no folds the group's outcome model is fitted on all of its rows, so
  # the only requirement is that there is something to fit. Below two rows the
  # fit was skipped and the quantity stayed at its zero initialisation, which
  # then entered the doubly-robust pseudo-outcome as if it were an estimate.
  short <- c("control", "treated")[c(length(idx_t0) < 2L, length(idx_t1) < 2L)]
  if (length(short)) {
    rlang::abort(
      c(sprintf("The %s group has fewer than 2 rows, so its outcome model cannot be fitted.",
                paste(short, collapse = " and the ")),
        i = "Every doubly-robust quantity needs a conditional mean for both groups.",
        i = "Collect more of that group, or score only the quantities that do not use it."),
      class = "cateomics_group_too_small"
    )
  }
  if (out_is_sl) {
    if (length(idx_t0) > 1L) {
      r0 <- .cf_sl_full_quantity(X_mat, Y_vec, .sl_out_lib, type = otype,
                                 strata = T_vec, n_folds = n_folds,
                                 random_state = random_state,
                                 member_params = .eval_out_p,
                                 train_rows = idx_t0, weight_rows = idx_t0,
                                 where = "evaluation outcome (control arm)")
      mu0_in <- r0$pred
      sl_weights$mu0 <- r0$weights
    }
    if (length(idx_t1) > 1L) {
      r1 <- .cf_sl_full_quantity(X_mat, Y_vec, .sl_out_lib, type = otype,
                                 strata = T_vec, n_folds = n_folds,
                                 random_state = random_state,
                                 member_params = .eval_out_p,
                                 train_rows = idx_t1, weight_rows = idx_t1,
                                 where = "evaluation outcome (treated arm)")
      mu1_in <- r1$pred
      sl_weights$mu1 <- r1$weights
    }
    rm_ <- .cf_sl_full_quantity(X_mat, Y_vec, .sl_out_lib, type = otype,
                                strata = T_vec, n_folds = n_folds,
                                random_state = random_state,
                                member_params = .eval_out_p,
                                where = "evaluation outcome (pooled)")
    m_in <- rm_$pred
    sl_weights$m <- rm_$weights
  } else {
    if (length(idx_t0) > 1L) {
      m0 <- .fit_base_learner(X_mat[idx_t0, , drop = FALSE], Y_vec[idx_t0],
                              type = otype, learner = spec$outcome,
                              random_state = random_state,
                              .params = .eval_out_p)
      mu0_in <- .predict_base_learner(m0, newdata = X_mat)
    }
    if (length(idx_t1) > 1L) {
      m1 <- .fit_base_learner(X_mat[idx_t1, , drop = FALSE], Y_vec[idx_t1],
                              type = otype, learner = spec$outcome,
                              random_state = random_state,
                              .params = .eval_out_p)
      mu1_in <- .predict_base_learner(m1, newdata = X_mat)
    }
    mm <- .fit_base_learner(X_mat, Y_vec, type = otype,
                            learner = spec$outcome,
                            random_state = random_state,
                            .params = .eval_out_p)
    m_in <- .predict_base_learner(mm, newdata = X_mat)
  }

  list(mu0_oof = mu0_in, mu1_oof = mu1_in, pi_oof = pi_in, m_oof = m_in,
       sl_weights = if (length(sl_weights)) sl_weights else NULL)
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
#' regression loss (Kennedy 2023); used for CATE model selection (Mahajan et al.
#' 2024). Full citations are in \code{?cf_compare}.
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
#' \code{base_per_obs} is the loss of the reference model the caller passes:
#' the R-loss-optimal constant for r_loss, the mean of the DR pseudo-outcome
#' for dr_loss, and the per-arm sample means for factual_loss.
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
# r_stacking / dr_stacking need per-obs material (since the loss of the
# weighted-average tau is NOT the weighted average of per-method losses);
# they take it from `.cateomics_compute_stacking_materials` below.
#
# The required ingredients are all already cached in cf_compare's output:
#
#   * `cmp$metadata$nuisance_oof` = list(mu0_oof, mu1_oof, pi_oof, m_oof)
#   * `cmp$models[[m]]$data$Y` / `$T` (any surviving model is enough)
#   * `cmp$metadata$oof_tau_matrix` = n x K tau hats
#
# All recomputation is pure arithmetic (no nuisance refit).


#' Extract the cf_compare cache needed by per-obs loss recomputation.
#'
#' Pulls Y / T from the first surviving model, nuisance_oof from metadata,
#' and exposes the oof_tau_matrix. Guards against the edge case where every
#' method failed (no model retains $data$Y / $data$T).
#'
#' @param compare_result A `cateomics_compare` object.
#' @return A list with `Y`, `T`, `mu0_oof`, `mu1_oof`, `pi_oof`,
#'   `oof_tau_matrix`, and `methods` (column names of oof_tau_matrix).
#'   Aborts when the cache is unusable (`cateomics_missing_nuisance_oof` /
#'   `cateomics_missing_oof_tau_matrix` / `cateomics_compare_cache_empty` /
#'   `cateomics_missing_m_oof`).
#' @keywords internal
#' @noRd
.cateomics_extract_compare_cache <- function(compare_result) {
  if (!inherits(compare_result, "cateomics_compare")) {
    rlang::abort(
      "compare_result must be a cateomics_compare object.",
      class = "cateomics_invalid_compare_result"
    )
  }

  meta <- compare_result$metadata
  nuis <- meta$nuisance_oof
  if (is.null(nuis) ||
      is.null(nuis$mu0_oof) || is.null(nuis$mu1_oof) || is.null(nuis$pi_oof)) {
    rlang::abort(
      paste0("cf_compare metadata is missing nuisance_oof; ",
             "this object was not produced by cf_compare(), or its ",
             "metadata was modified."),
      class = "cateomics_missing_nuisance_oof"
    )
  }

  oof_tau_matrix <- meta$oof_tau_matrix
  if (is.null(oof_tau_matrix) || !is.matrix(oof_tau_matrix)) {
    rlang::abort(
      "cf_compare metadata is missing oof_tau_matrix.",
      class = "cateomics_missing_oof_tau_matrix"
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
      class = "cateomics_compare_cache_empty"
    )
  }

  list(
    Y              = as.numeric(YT$Y),
    T              = as.integer(YT$T),
    mu0_oof        = as.numeric(nuis$mu0_oof),
    mu1_oof        = as.numeric(nuis$mu1_oof),
    pi_oof         = as.numeric(nuis$pi_oof),
    # m_oof is the direct OOF regression E[Y|X]; it is never reconstructed as
    # mu1*pi + mu0*(1-pi). The shared eval nuisance always supplies it, so a
    # missing value is a contract breach.
    m_oof          = if (!is.null(nuis$m_oof)) as.numeric(nuis$m_oof)
                     else rlang::abort(
                       "cf_compare cache nuisance is missing m_oof (E[Y|X]).",
                       class = "cateomics_missing_m_oof"),
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
#' @param compare_result A `cateomics_compare` object.
#' @return A list with `Yres`, `Tres`, `Y_DR`, plus the underlying
#'   `Y`, `T`, `mu0_oof`, `mu1_oof`, `pi_oof`, and `oof_tau_matrix` for
#'   convenience.
#' @keywords internal
#' @noRd
.cateomics_compute_stacking_materials <- function(compare_result) {
  cache <- .cateomics_extract_compare_cache(compare_result)

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
#' Y / T / nuisance via `.cateomics_extract_compare_cache`, and dispatches to
#' the underlying per-obs primitive (`.compute_r_loss_per_obs`,
#' `.compute_dr_loss_per_obs`, `.compute_factual_loss_per_obs`).
#'
#' Used by `cf_ensemble(method = "softmax")` when the requested loss column
#' is absent from the scoreboard. Pure arithmetic, no nuisance refit.
#'
#' @param compare_result A `cateomics_compare` object.
#' @param method Character; one of `colnames(oof_tau_matrix)`.
#' @param loss_type Character; one of `"r_loss"`, `"dr_loss"`,
#'   `"factual_loss"`.
#' @return Numeric vector of length `n` with the per-observation loss.
#' @keywords internal
#' @noRd
.cateomics_compute_per_obs_loss_from_compare <- function(compare_result, method,
                                                      loss_type = c("r_loss",
                                                                    "dr_loss",
                                                                    "factual_loss")) {
  loss_type <- match.arg(loss_type)
  cache     <- .cateomics_extract_compare_cache(compare_result)

  if (!is.character(method) || length(method) != 1L) {
    rlang::abort("`method` must be a single character string.",
                 class = "cateomics_invalid_method")
  }
  if (!method %in% cache$methods) {
    rlang::abort(
      sprintf("method '%s' not found in oof_tau_matrix columns: %s",
              method, paste(cache$methods, collapse = ", ")),
      class = "cateomics_method_not_in_compare"
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
