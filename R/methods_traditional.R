#' Fit Generalized Random Forest (GRF) causal forest model
#'
#' Two-stage orthogonalized causal forest following Athey, Tibshirani & Wager
#' (2019). Nuisance surfaces \eqn{\hat W(X)} and \eqn{\hat Y(X)} are fit by
#' \code{grf::regression_forest()}; the treatment-effect surface
#' \eqn{\hat\tau(X)} is then fit by \code{grf::causal_forest()} on the
#' residualized scores.
#'
#' Returns the full \code{cfomics_result} contract (\code{tau_hat},
#' \code{mu0_hat}, \code{mu1_hat}, \code{ate} and summary). The ATE and its
#' 95\% Wald CI come
#' from \code{grf::average_treatment_effect()} (doubly-robust score), not
#' from \eqn{\bar{\hat\tau}(X)}.
#'
#' Binary (0/1) outcomes are modelled on the native continuous scale, so the
#' estimand is the risk difference; the counterfactual means may fall outside
#' \eqn{[0, 1]}. For a proper binary likelihood use a meta/DR learner
#' (classifier nuisance) or a neural method.
#'
#' @section Potential outcomes (Class C: Robinson residualization):
#' \code{mu0_hat = Y.hat - W.hat * tau_hat}, \code{mu1_hat = Y.hat +
#' (1 - W.hat) * tau_hat}, where \code{Y.hat = E[Y|X]} is a single outcome
#' regression and \code{W.hat} is the propensity estimate (Robinson
#' 1988 decomposition). \code{mu1_hat - mu0_hat == tau_hat} by construction.
#'
#' \strong{NOTE}: \code{mu_t} are reconstructed from the single
#' \code{Y.hat} model + residualisation — \emph{not} separate outcome
#' regressions per arm. The DR-score ATE may differ from
#' \code{mean(ite)}; cfomics uses \code{grf::average_treatment_effect()}.
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
#'   default \code{"none"}. This
#'   matches the grf upstream default and keeps cfomics' method
#'   comparison framework "tuning-OFF for all methods" (only DML /
#'   DR / catenets / etc. parametric backends select their
#'   `model_*` base learner, which is hyperparameter-fixed lasso by
#'   default). Pass \code{"all"} to opt into grf cross-validation
#'   tuning; note that for small training folds (n <~ 200) auto-tune
#'   may abort with grf's "honesty fraction is too close to 1 or 0"
#'   error — refer to grf docs for guidance. Forwarded to both
#'   \code{regression_forest} (W.hat / Y.hat nuisance) and
#'   \code{causal_forest} (CATE forest).
#' @param estimate.variance Logical, whether \code{predict.causal_forest}
#'   should return variance estimates (default \code{TRUE}).
#' @param ... Additional arguments forwarded to
#'   \code{grf::causal_forest()}.
#' @return A list with \code{model} (the fitted \code{grf::causal_forest}
#'   object) and \code{res} conforming to the \code{cfomics_result} contract.
#' @keywords internal
cf_fit_grf <- function(X, T, Y,
                       covariate_names = NULL,
                       random_state = 0L,
                       num.trees = 2000L,
                       tune.parameters = "none",
                       estimate.variance = TRUE,
                       propensity_score = NULL,
                       propensity_X = NULL,
                       trim = c(0.025, 0.975),    # unified propensity clip bounds
                       ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    rlang::abort("Package 'grf' is required for method='grf'. Please install it.", class = "cfomics_missing_dependency")
  }

  # A user-supplied scalar / vector propensity_score is forwarded to
  # grf::causal_forest as W.hat; when NULL, W.hat is estimated internally by
  # grf::regression_forest.
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cfomics_invalid_propensity_score"
      )
    }
    if (length(propensity_score) > 1L && length(propensity_score) != length(T)) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match length(T) = %d.",
                length(propensity_score), length(T)),
        class = "cfomics_invalid_propensity_score"
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

  # With a user-supplied PS, skip forest.W estimation and use W.hat directly.
  if (!is.null(propensity_score)) {
    W.hat <- if (length(propensity_score) == 1L) {
      rep(as.numeric(propensity_score), length.out = n)
    } else {
      as.numeric(propensity_score)
    }
    # Apply the trim clip to (lower, upper).
    W.hat <- pmin(pmax(W.hat, trim[1L]), trim[2L])
    # forest.W = NULL is stored in res (cannot be reused for newdata; warns).
    forest.W <- NULL
  } else {
    # propensity_X can subset the X used to estimate W.hat, while the outcome
    # forest.Y still uses the full X (DML standard practice: the outcome model
    # sees full X, the PS model sees the user-selected X). Small-sample failure
    # under opt-in `tune.parameters = "all"` is left to grf to handle.
    .ps_res <- .resolve_propensity_X(X_mat, propensity_X, covariate_names)
    X_for_ps <- .ps_res$X_ps
    forest.W <- grf::regression_forest(X_for_ps, W, tune.parameters = tune.parameters)
    W.hat <- predict(forest.W)$predictions
  }

  forest.Y <- grf::regression_forest(X_mat, Y_vec, tune.parameters = tune.parameters)
  Y.hat <- predict(forest.Y)$predictions

  ## Via cf_compare, `...` may carry arguments meant for other methods (e.g.
  ## n_iter for catenets, num_epochs for cevae), so intersect with
  ## grf::causal_forest's formals and filter to the accepted names.
  .extra_args <- list(...)
  .grf_formals <- names(formals(grf::causal_forest))
  .extra_grf <- .extra_args[intersect(names(.extra_args), .grf_formals)]
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
  
  ite_summary <- list(
    mean = mean(ite),
    std = stats::sd(ite),
    quantile_05 = stats::quantile(ite, 0.05),
    quantile_50 = stats::quantile(ite, 0.50),
    quantile_95 = stats::quantile(ite, 0.95)
  )
  
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

  summary_stats <- list(
    ate = ate_estimate,
    ate_ci_lower = ate_ci_lower_val,
    ate_ci_upper = ate_ci_upper_val,
    ite_mean = ite_summary$mean,
    ite_std = ite_summary$std,
    ite_quantiles = list(
      q05 = ite_summary$quantile_05,
      q50 = ite_summary$quantile_50,
      q95 = ite_summary$quantile_95
    )
  )

  list(
    model    = tau.forest,
    forest.Y = forest.Y,   # persist nuisance forests for newdata mu_t
    forest.W = forest.W,   # persist nuisance forests for newdata mu_t
    res = list(
      tau_hat = ite,
      ate = ate_estimate,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat,
      propensity = W.hat,  # W.hat as an external inspection surface
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
#' \code{$ate / $ate_ci / $tau_hat / $tau_ci / $mu0_hat / $mu1_hat}. CI on
#' newdata is supported via the IJ variance estimate (Wald). The
#' counterfactual potential outcomes \eqn{\hat y_0}, \eqn{\hat y_1} on
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
#' @return Depending on `type`: the full list of predictions when
#'   `type = "all"`, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_grf <- function(object, newdata = NULL, type = "all",
                           estimate.variance = TRUE, ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }

  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.", class = "cfomics_unsupported_predict_type")
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
        class = "cfomics_grf_newdata_missing_columns"
      )
    }
    X_new <- X_new[, cov_names, drop = FALSE]
  }

  tau.pred <- predict(tau.forest, X_new, estimate.variance = estimate.variance)
  ite_new <- tau.pred$predictions

  # Wald CI from IJ variance (estimate.variance = TRUE)
  var_new <- if (!is.null(tau.pred$variance.estimates)) tau.pred$variance.estimates else NULL
  ite_ci_new <- NA_real_
  ate_ci_new <- NA_real_
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
        ), class = "cfomics_newdata_missing_column")
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
      NULL
    }
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
    mu1_hat = y1_new
  )
}

