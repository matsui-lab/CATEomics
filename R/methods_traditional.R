#' Fit Generalized Random Forest (GRF) causal forest model
#'
#' Two-stage orthogonalized causal forest following Athey, Tibshirani & Wager
#' (2019). Nuisance surfaces \eqn{\hat W(X)} and \eqn{\hat Y(X)} are fit by
#' \code{grf::regression_forest()}; the treatment-effect surface
#' \eqn{\hat\tau(X)} is then fit by \code{grf::causal_forest()} on the
#' residualized scores.
#'
#' Returns the full \code{cateomics_result} contract (\code{tau_hat},
#' \code{mu0_hat}, \code{mu1_hat}, \code{ate} and summary). The ATE and its
#' 95\% Wald CI come
#' from \code{grf::average_treatment_effect()} (doubly-robust score), not
#' from \eqn{\bar{\hat\tau}(X)}.
#'
#' Binary (0/1) outcomes are modelled on the native continuous scale, so the
#' estimand is the risk difference; the counterfactual means may fall outside
#' \eqn{[0, 1]}.
#'
#' @section Potential outcomes -- Class C, Robinson residualization:
#' \code{mu0_hat = Y.hat - W.hat * tau_hat}, \code{mu1_hat = Y.hat +
#' (1 - W.hat) * tau_hat}, where \code{Y.hat = E[Y|X]} is a single outcome
#' regression and \code{W.hat} is the propensity estimate (Robinson
#' 1988 decomposition). \code{mu1_hat - mu0_hat == tau_hat} by construction.
#'
#' \strong{NOTE}: \code{mu_t} are reconstructed from the single
#' \code{Y.hat} model + residualisation — \emph{not} separate outcome
#' regressions per arm. The DR-score ATE may differ from
#' \code{mean(ite)}; CATEomics uses \code{grf::average_treatment_effect()}.
#'
#' @details
#' The potential-outcome surfaces are reconstructed from the residualized
#' decomposition:
#' \deqn{\hat y_0 = \hat Y(X) - \hat W(X)\,\hat\tau(X)}
#' \deqn{\hat y_1 = \hat Y(X) + (1 - \hat W(X))\,\hat\tau(X)}
#' which satisfies \eqn{\hat y_1 - \hat y_0 = \hat\tau(X)} exactly by
#' construction.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param random_state Integer random seed. When set, yields bit-identical
#'   results across runs with the same inputs. The caller's global
#'   \code{.Random.seed} is saved on entry and restored on exit so seeding
#'   does not leak into the user's session.
#' @param num.trees Integer, number of trees for \code{causal_forest}
#'   (default \code{2000}).
#' @param tune.parameters Character, which GRF parameters to tune;
#'   default \code{"none"} (the grf default). Pass \code{"all"} for grf
#'   cross-validation tuning. Forwarded to both
#'   \code{regression_forest} (W.hat / Y.hat nuisance) and
#'   \code{causal_forest} (CATE forest).
#' @param estimate.variance Logical, whether \code{predict.causal_forest}
#'   should return variance estimates (default \code{TRUE}).
#' @param params Named list of extra arguments. \code{params$nuisance}
#'   configures the orthogonalization forests; everything else goes to
#'   \code{grf::causal_forest()}. Names are checked on the R side against the
#'   arguments of this fitter and of \code{grf::causal_forest()}, so an
#'   unknown name aborts before the forest is grown.
#' @param propensity_score Known propensity score: a scalar in (0, 1) or a
#'   length-n vector. \code{NULL} (default) estimates it internally.
#' @param propensity_X Optional character vector of covariate names to
#'   restrict the propensity model to. Reached through \code{\link{cf_fit}} or
#'   \code{\link{cf_compare}}, leaving it unset fits the propensity on the
#'   clinical covariates alone; called directly, \code{NULL} means every
#'   column of the design matrix.
#' @param clip Length-2 numeric, the bounds the propensity is held within
#'   (default \code{c(0.025, 0.975)}). Values are replaced by the bound; no
#'   sample is dropped.
#' @return A list with \code{model} (the fitted \code{grf::causal_forest}
#'   object) and \code{res} conforming to the \code{cateomics_result} contract.
#' @references
#' Athey, Tibshirani & Wager (2019) Generalized Random Forests,
#' \emph{Ann. Statist.} 47(2):1148-1178.
#'
#' Robinson (1988) Root-N-Consistent Semiparametric Regression,
#' \emph{Econometrica} 56(4):931-954 (the residualisation the potential
#' outcomes are reconstructed with).
#'
#' @keywords internal
cf_fit_grf <- function(X, T, Y,
                       covariate_names = NULL,
                       random_state = 0L,
                       num.trees = 2000L,
                       tune.parameters = "none",
                       estimate.variance = TRUE,
                       propensity_score = NULL,
                       propensity_X = NULL,
                       clip = .CF_PS_CLIP,    # unified propensity clip bounds
                       params = list()) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    rlang::abort("Package 'grf' is required for method='grf'. Please install it.", class = "cateomics_missing_dependency")
  }

  # A user-supplied scalar / vector propensity_score is forwarded to
  # grf::causal_forest as W.hat; when NULL, W.hat is estimated internally by
  # grf::regression_forest.
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cateomics_invalid_propensity_score"
      )
    }
    if (length(propensity_score) > 1L && length(propensity_score) != length(T)) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match length(T) = %d.",
                length(propensity_score), length(T)),
        class = "cateomics_invalid_propensity_score"
      )
    }
  }

  if (!is.null(random_state)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    # Pin the RNG algorithm (not just the seed) so the forest is reproducible
    # regardless of the caller's RNG kind. cf_compare fits per method inside
    # BiocParallel::bplapply, which switches the stream to L'Ecuyer-CMRG even
    # under SerialParam; without this, cf_fit (Mersenne-Twister) and cf_compare
    # would draw different grf forests from the same seed. The saved
    # `.Random.seed` (which encodes the kind) is restored on exit.
    set.seed(random_state, kind = "Mersenne-Twister")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(W)

  # `params$nuisance` configures the two orthogonalization forests; every other
  # entry is a `grf::causal_forest()` argument. Both were validated upstream.
  .extra_grf     <- params[setdiff(names(params), "nuisance")]

  # Match grf's own orthogonalization settings (its `args.orthog`): a quarter
  # of the trees, a fixed minimum node size, honesty on, no confidence groups,
  # and the main forest's splitting arguments carried over.
  # `seed` is refused in `params` (random_state is the one knob), so it is
  # not listed here: there is nothing to carry over.
  .grf_carried <- c("sample.fraction", "mtry", "alpha", "imbalance.penalty",
                    "honesty.prune.leaves", "sample.weights", "clusters",
                    "equalize.cluster.weights", "num.threads")
  .grf_split_args <- utils::modifyList(
    list(sample.fraction      = 0.5,
         mtry                 = min(ceiling(sqrt(ncol(X_mat)) + 20), ncol(X_mat)),
         alpha                = 0.05,
         imbalance.penalty    = 0,
         honesty.prune.leaves = TRUE),
    .extra_grf[intersect(names(.extra_grf), .grf_carried)])
  .nuisance_args <- utils::modifyList(
    c(.grf_split_args,
      list(num.trees        = max(50, as.integer(num.trees) / 4),
           min.node.size    = 5,
           honesty          = TRUE,
           honesty.fraction = 0.5,
           ci.group.size    = 1,
           tune.parameters  = tune.parameters)),
    params[["nuisance"]] %||% list())

  # grf builds both nuisance forests on the same X, so its mtry always fits.
  # With `propensity_X` the treatment forest sees fewer columns, and an mtry
  # larger than that turns it into near-exhaustive bagging.
  .fit_nuisance_forest <- function(X_in, Y_in) {
    args <- .nuisance_args
    if (!is.null(args$mtry)) args$mtry <- min(args$mtry, ncol(X_in))
    do.call(grf::regression_forest, c(list(X = X_in, Y = Y_in), args))
  }

  # With a user-supplied PS, skip forest.W estimation and use W.hat directly.
  if (!is.null(propensity_score)) {
    W.hat <- if (length(propensity_score) == 1L) {
      rep(as.numeric(propensity_score), length.out = n)
    } else {
      as.numeric(propensity_score)
    }
    # Apply the clip to (lower, upper).
    W.hat <- pmin(pmax(W.hat, clip[1L]), clip[2L])
    # No propensity forest is kept: the propensity is known, so there is
    # nothing to fit. `predict(newdata =)` can replay a constant value but not
    # a per-sample one, and says so when it cannot.
    forest.W <- NULL
  } else {
    # propensity_X can subset the X used to estimate W.hat, while the outcome
    # forest.Y still uses the full X (DML standard practice: the outcome model
    # sees full X, the PS model sees the user-selected X). Small-sample failure
    # under opt-in `tune.parameters = "all"` is left to grf to handle.
    .ps_res <- .resolve_propensity_X(X_mat, propensity_X, covariate_names)
    X_for_ps <- .ps_res$X_ps
    forest.W <- .fit_nuisance_forest(X_for_ps, W)
    # Clip the estimated propensity on the same bounds as a supplied one, so
    # `clip` means the same thing whoever produced the score.
    W.hat <- pmin(pmax(predict(forest.W)$predictions, clip[1L]), clip[2L])
  }

  forest.Y <- .fit_nuisance_forest(X_mat, Y_vec)
  Y.hat <- predict(forest.Y)$predictions

  tau.forest <- do.call(grf::causal_forest, c(list(
    X = X_mat, Y = Y_vec, W = W,
    W.hat = W.hat,
    Y.hat = Y.hat,
    num.trees = as.integer(num.trees),
    tune.parameters = tune.parameters
  ), .extra_grf))

  tau.pred <- predict(tau.forest, estimate.variance = estimate.variance)
  ite <- tau.pred$predictions

  ate <- grf::average_treatment_effect(tau.forest)
  ate_estimate <- ate["estimate"]
  ate_se <- ate["std.err"]

  mu0_hat <- Y.hat - W.hat * ite
  mu1_hat <- Y.hat + (1 - W.hat) * ite
  
  ate_ci_lower_val <- as.numeric(ate_estimate - 1.96 * ate_se)
  ate_ci_upper_val <- as.numeric(ate_estimate + 1.96 * ate_se)
  ate_ci <- c(ate_ci_lower_val, ate_ci_upper_val)

  # Per-row ite_ci from GRF variance estimates (Wald, IJ variance via
  # predict(causal_forest, estimate.variance = TRUE)).
  ite_ci <- NULL
  if (estimate.variance && "variance.estimates" %in% names(tau.pred)) {
    se_ite <- sqrt(pmax(tau.pred$variance.estimates, 0))
    ite_ci <- data.frame(
      lower = as.numeric(ite) - 1.96 * se_ite,
      upper = as.numeric(ite) + 1.96 * se_ite
    )
  }

  # The one builder every entry uses, so a fit, a comparison row and a direct
  # call to this fitter describe the effect the same way.
  summary_stats <- .cf_build_summary(ate_estimate, ate_ci, ite)

  list(
    model    = tau.forest,
    forest.Y = forest.Y,   # persist nuisance forests for newdata mu_t
    forest.W = forest.W,   # persist nuisance forests for newdata mu_t
    ps_clip  = clip,       # so newdata bounds the propensity the same way
    res = list(
      tau_hat = ite,
      ate = ate_estimate,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat,
      # W.hat as an external inspection surface. predict(regression_forest)
      # without newdata returns out-of-bag predictions, so an estimated W.hat
      # is already held out from each row's own tree ensemble.
      ps      = W.hat,
      ps_type = if (!is.null(propensity_score)) "known" else "oob",
      ate_ci  = ate_ci,
      tau_ci  = ite_ci,
      summary = summary_stats,
      metadata = list(ci_method = "asymptotic Wald")
    )
  )
}

#' Predict from GRF model
#'
#' @description
#' \code{type = "all"} returns the flat list
#' \code{$ate / $ate_ci / $tau_hat / $tau_ci / $mu0_hat / $mu1_hat / $ps /
#' $ps_type}. On newdata the per-row \code{tau_ci} comes from the IJ
#' variance estimate (Wald) and \code{ate_ci} is \code{NA}. The
#' counterfactual potential outcomes \eqn{\hat\mu_0}, \eqn{\hat\mu_1} on
#' newdata are reconstructed from the persisted nuisance forests
#' \code{forest.Y} / \code{forest.W} via the Robinson decomposition.
#'
#' @param object \code{cf_model} object with \code{method = "grf"}.
#' @param newdata Optional new covariate matrix. When \code{NULL}, the
#'   stored training-sample predictions are returned.
#' @param type One of \code{"all"}, \code{"ate"}, \code{"ate_ci"},
#'   \code{"tau"}, \code{"tau_ci"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param estimate.variance Logical, forwarded to
#'   \code{predict.causal_forest}. Default \code{TRUE}.
#' @param ... Ignored.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_grf <- function(object, newdata = NULL, type = "all",
                           estimate.variance = TRUE, ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }

  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.", class = "cateomics_unsupported_predict_type")
  }

  tau.forest <- object$model$model
  X_new <- as.matrix(newdata)

  # grf's causal_forest predict aligns columns positionally, so newdata must be
  # replayed in the exact training covariate order. When the newdata carries
  # column names, select + reorder to the training covariate_names (dropping any
  # extra columns); error on missing ones.
  cov_names <- object$meta$covariate_names
  if (!is.null(cov_names) && !is.null(colnames(X_new))) {
    missing_cov <- setdiff(cov_names, colnames(X_new))
    if (length(missing_cov) > 0L) {
      rlang::abort(
        sprintf("grf newdata is missing training covariate(s): %s.",
                paste(missing_cov, collapse = ", ")),
        class = "cateomics_grf_newdata_missing_columns"
      )
    }
    X_new <- X_new[, cov_names, drop = FALSE]
  }

  tau.pred <- predict(tau.forest, X_new, estimate.variance = estimate.variance)
  ite_new <- tau.pred$predictions

  # Wald CI from IJ variance (estimate.variance = TRUE)
  var_new <- if (!is.null(tau.pred$variance.estimates)) tau.pred$variance.estimates else NULL
  ite_ci_new <- NULL
  ate_ci_new <- NULL
  if (!is.null(var_new)) {
    se_new <- sqrt(pmax(var_new, 0))
    ite_ci_new <- data.frame(
      lower = as.numeric(ite_new) - 1.96 * se_new,
      upper = as.numeric(ite_new) + 1.96 * se_new
    )
    # No honest ATE CI on newdata: grf has no outcomes there, and the SE of the
    # mean of correlated per-point IJ variances is not sqrt(mean(var_i)) (that
    # ignores the 1/n averaging). Leave ate_ci_new = NA; the per-point
    # ite_ci_new (Wald from IJ variance) remains valid.
  }
  ate_new <- mean(ite_new)

  # Reconstruct mu_t = Y.hat +/- W.hat-weighted tau on newdata.
  #   Y.hat comes from the persisted outcome forest (forest.Y).
  #   W.hat comes from the persisted propensity forest (forest.W); when a known
  #   scalar PS was used at fit time forest.W is NULL, but the constant W.hat is
  #   stored on the causal forest and can be replayed for new rows. A per-obs
  #   (vector) known PS cannot be recovered for new rows, so mu is left NA.
  forest_Y <- object$model$forest.Y
  forest_W <- object$model$forest.W
  Y_hat_new <- if (!is.null(forest_Y)) predict(forest_Y, X_new)$predictions else NULL

  W_hat_new <- if (!is.null(forest_W)) {
    # When propensity_X is set, forest_W was trained on X[, propensity_X] so
    # newdata must be subset to the same columns before predicting W.hat.
    px <- object$meta$propensity_X
    if (!is.null(px)) {
      px_idx <- match(px, colnames(X_new))
      if (anyNA(px_idx)) {
        rlang::abort(sprintf(
          "newdata is missing propensity_X column(s): %s",
          paste(px[is.na(px_idx)], collapse = ", ")
        ), class = "cateomics_newdata_missing_column")
      }
      predict(forest_W, X_new[, px_idx, drop = FALSE])$predictions
    } else {
      predict(forest_W, X_new)$predictions
    }
  } else {
    # No propensity forest: replay a constant known PS if that is what was used.
    w_train <- object$model$model$W.hat
    if (!is.null(w_train) && length(unique(round(w_train, 12))) == 1L) {
      rep(w_train[1L], nrow(X_new))
    } else {
      # A propensity given per sample is a property of the rows it was given
      # for, and there is no value for a row the fit never saw. The effect is
      # still predicted; the counterfactual means are what need the propensity.
      rlang::warn(
        c("The known propensity was supplied per sample, so it cannot be carried to new rows.",
          i = "`mu0_hat` and `mu1_hat` come back as NA for `newdata`; `tau_hat` is unaffected.",
          i = "Leave `propensity_score` unset if the counterfactual means are needed on new data."),
        class = "cateomics_newdata_propensity_unavailable"
      )
      NULL
    }
  }

  # Bound the newdata propensity on the same interval as the fit did,
  # otherwise predict(fit, newdata = training X) would not reproduce
  # predict(fit) and the held-out mu on the OOF path would be on a different
  # scale from the full-model mu.
  if (!is.null(W_hat_new)) {
    ps_clip <- object$model$ps_clip %||% c(0.025, 0.975)
    W_hat_new <- pmin(pmax(W_hat_new, ps_clip[1L]), ps_clip[2L])
  }

  if (!is.null(Y_hat_new) && !is.null(W_hat_new)) {
    y0_new <- Y_hat_new - W_hat_new * ite_new
    y1_new <- Y_hat_new + (1 - W_hat_new) * ite_new
  } else {
    y0_new <- rep(NA_real_, length(ite_new))
    y1_new <- rep(NA_real_, length(ite_new))
  }

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    ate_ci  = ate_ci_new,
    tau     = as.numeric(ite_new),
    tau_ci  = ite_ci_new,
    mu0_hat = y0_new,
    mu1_hat = y1_new,
    ps      = if (!is.null(W_hat_new)) as.numeric(W_hat_new) else
                rep(NA_real_, length(ite_new)),
    ps_type = if (is.null(W_hat_new)) NA_character_ else
              if (is.null(forest_W)) "known" else "out-of-sample"
  )
}

