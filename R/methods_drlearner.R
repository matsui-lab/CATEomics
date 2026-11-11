#' Fit a DR-Learner via econml
#'
#' Thin R wrapper around \code{econml.dr.DRLearner} using the shared CATEomics
#' base-learner vocabulary for nuisance estimators (propensity, outcome,
#' final CATE).
#'
#' Returns the full \code{cateomics_result} contract. The ITE is taken from
#' the orthogonalized final stage (\code{dr$effect}); potential outcomes
#' \eqn{\hat y_0, \hat y_1} are obtained separately by averaging the
#' internal \code{models_regression} nuisance predictions over Monte Carlo
#' iterations and cross-validation folds (see
#' \code{inst/python/metalearner_helpers.py::drlearner_potential_outcomes}).
#'
#' @section Execution model:
#' All Python work runs inside the CATEomics basilisk environment via
#' \code{.cf_py_run()}. The fitted \code{DRLearner}
#' cannot be returned as a live object, so it is serialised with
#' \code{joblib.dump()} to a cache file during fit and reloaded at predict
#' time. The \code{model} pillar of the returned \code{cf_model} therefore
#' holds \code{list(model_path, backend)} rather than a live estimator.
#'
#' @section Potential outcomes -- Class B prime, outcome regression as an intermediate nuisance:
#' \code{mu0_hat} / \code{mu1_hat} are the cross-fit nuisance regressions
#' (\code{models_regression}) averaged over the Monte Carlo iterations and the cross-validation folds.
#' \code{tau_hat} is the orthogonalised CATE projected by the final-stage
#' regression on pseudo-outcomes.
#'
#' \strong{NOTE}: \code{mu1_hat - mu0_hat != tau_hat}. The Neyman-orthogonal
#' final stage regresses the pseudo-outcome rather than differencing the
#' plug-in nuisance fits, so the two quantities come from different stages.
#' \code{ate} is reported as \code{mean(tau_hat)} (orthogonal score), not
#' \code{mean(mu1_hat - mu0_hat)}.
#'
#' @section Confidence intervals:
#' When \code{model_final = "linear"} the final stage uses econml's
#' \code{StatsModelsLinearRegression} so that \code{ate_interval()} /
#' \code{effect_interval()} return asymptotic Wald CIs. \code{res$ate_ci}
#' (length-2 numeric) and \code{res$tau_ci} (data.frame with \code{lower} /
#' \code{upper} columns) are populated. For non-linear \code{model_final},
#' set \code{bootstrap_ci = TRUE} to obtain percentile-method CIs from
#' \code{n_bootstrap} refits (computed inside the same basilisk process).
#'
#' @details
#' \eqn{\hat y_1 - \hat y_0} is the plug-in
#' contrast of the nuisance outcome models, while \code{tau_hat} is the
#' Neyman-orthogonal CATE from the final stage; the two do not agree
#' numerically.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param model_propensity Base-learner name for the propensity model.
#'   Default \code{"logistic"}.
#' @param model_outcome Base-learner name for the outcome regression used
#'   internally by econml. Default \code{"ridge"}.
#' @param model_final Base-learner name for the final CATE regression.
#'   Default \code{"ridge"}. It does \emph{not} track \code{model_outcome},
#'   unlike \code{\link{cf_fit_xlearner}}'s Stage-2 model. \code{"linear"}
#'   enables asymptotic Wald CIs via
#'   \code{StatsModelsLinearRegression}.
#' @param cv Integer, number of cross-fitting folds. Default \code{5}.
#' @param bootstrap_ci Logical. When \code{TRUE}, override any native CI
#'   path and compute percentile-method CIs from \code{n_bootstrap} refits.
#'   The resampling is implemented in CATEomics and is the same for every
#'   method that accepts it. Default \code{FALSE}.
#' @param n_bootstrap Integer, number of bootstrap replicates when
#'   \code{bootstrap_ci = TRUE}. Default \code{100}.
#' @param clip Length-2 numeric propensity clipping bounds.
#' @param propensity_score Optional known propensity score (RCT context);
#'   scalar or length-n vector. Mutually exclusive with \code{model_propensity}.
#' @param propensity_X Optional character vector of covariate names the
#'   propensity model is restricted to.
#'   With \code{propensity_X} unset this is the clinical covariates alone;
#'   see \code{\link{cf_fit}}.
#' @param random_state Integer random seed. When set, yields reproducible
#'   results. Forwarded to each sklearn nuisance learner and to numpy's RNG
#'   inside the basilisk worker. That worker can be the calling R process
#'   itself, so the Python-side generators are snapshotted and restored around
#'   the call and the seeding does not leak into the caller's session.
#' @param binary_outcome Logical. When \code{TRUE} the stage-1 outcome model
#'   becomes a probability classifier and DRLearner runs with
#'   \code{discrete_outcome = TRUE}.
#' @param params Named list with optional \code{outcome} / \code{final}
#'   blocks of hyperparameters for the base learners. A block reaches the
#'   scikit-learn estimator built in \code{inst/python/learner_factory.py},
#'   which raises on a hyperparameter it does not take. A flat entry is
#'   refused on the R side, because this method reads the blocks only.
#'   A value the Python side requires to be whole is written as an R
#'   integer (\code{10L}, not \code{10}), which R otherwise makes a double.
#' @param propensity_params Named list of hyperparameters for the propensity model.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @references
#' Kennedy (2023) Towards optimal doubly robust estimation of heterogeneous
#' causal effects, \emph{Electron. J. Stat.} 17(2):3008-3049.
#'
#' Battocchi, Dillon, Hei, Lewis, Oka, Oprescu & Syrgkanis (2019) EconML: A
#' Python Package for ML-Based Heterogeneous Treatment Effects Estimation,
#' \url{https://github.com/py-why/EconML} (implementation).
#'
#' @keywords internal
cf_fit_drlearner <- function(X, T, Y,
                             covariate_names  = NULL,
                             model_propensity = "logistic",
                             model_outcome = "ridge",
                             model_final      = "ridge",
                             cv               = 5L,
                             bootstrap_ci     = FALSE,
                             n_bootstrap      = 100L,
                             clip             = .CF_PS_CLIP,
                             propensity_score = NULL,
                             propensity_X     = NULL,
                             random_state     = NULL,
                             binary_outcome   = FALSE,
                             params            = list(),
                             propensity_params = list()) {
  # ---- R-side validation / resolution (outside the basilisk process) --------
  .validate_learner_spec(model_propensity)
  .validate_learner_spec(model_outcome)
  .validate_learner_spec(model_final)
  # econml's StatsModelsLinearRegression has nothing to configure, so a
  # `final` block here would be accepted and never used.
  if (identical(model_final, "linear") &&
      length(params[["final"]] %||% list())) {
    rlang::abort(
      c("`params$final` has no effect with model_final = \"linear\".",
        i = "That stage is econml's StatsModelsLinearRegression, which takes no hyperparameters.",
        i = "Use another model_final, or drop the block."),
      class = "cateomics_unused_params"
    )
  }
  if (!is.null(propensity_score) && length(propensity_params)) {
    rlang::abort(
      c("`propensity_params` has no effect when `propensity_score` is supplied.",
        i = "A known propensity is used as given; no propensity model is fitted."),
      class = "cateomics_propensity_params_unused"
    )
  }

  X_mat <- as.matrix(X)
  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }
  colnames(X_mat) <- covariate_names
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)

  # propensity_X → 0-based column indices (resolved in the main process where
  # CATEomics is loaded; passed to the worker as a plain integer vector).
  ps_indices_0b <- NULL
  if (!is.null(propensity_X)) {
    ps_indices_0b <- as.integer(
      .resolve_propensity_X(X_mat, propensity_X, covariate_names)$indices_0b
    )
  }

  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  model_path <- .cf_new_model_path("drlearner", "joblib")
  pypath     <- system.file("python", package = "CATEomics")

  # ---- Python fit + serialize (inside the basilisk process) -----------------
  # The worker is a self-contained closure: it references only its arguments,
  # base R, and reticulate (CATEomics is NOT loaded in the worker, so the
  # python helper path is passed in explicitly). `convert = FALSE` is required
  # so numpy arrays keep their Python identity (spike-proven).
  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, covariate_names,
             model_propensity, model_outcome, model_final, cv, clip,
             propensity_score, ps_indices_0b, seed_py, binary_outcome,
             bootstrap_ci, n_bootstrap,
             outcome_params, final_params, propensity_params) {
      sys     <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      clipped <- reticulate::import_from_path("clipped_propensity",
                                              path = pypath, convert = FALSE)
      interact <- reticulate::import_from_path("interaction_augmented",
                                               path = pypath, convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      DRLearner <- reticulate::import("econml.dr", convert = FALSE)$DRLearner

      if (!is.null(seed_py)) np$random$seed(seed_py)

      # A fresh DRLearner build (econml clones the sklearn learners per fit,
      # so a new build per bootstrap replicate is the clean path).
      build_dr <- function() {
        # econml fits this on hstack([XW, T]) and then predicts at T = 0 / 1.
        # Wrapping it expands that design to [XW, T, XW * T], so the two
        # counterfactual surfaces can differ in slope and not only in level.
        ml_g <- interact$make_interaction_augmented(
          if (isTRUE(binary_outcome)) {
            factory$make_classifier(model_outcome, random_state = seed_py,
                                    params = outcome_params)
          } else {
            factory$make_regressor(model_outcome, random_state = seed_py,
                                   params = outcome_params)
          })
        if (is.null(propensity_score)) {
          ml_m <- factory$make_classifier(model_propensity, random_state = seed_py,
                                          params = propensity_params)
          if (is.null(ps_indices_0b)) {
            ml_m <- clipped$ClippedPropensityModel(
              ml_m, lower = as.numeric(clip[1]), upper = as.numeric(clip[2]))
          } else {
            ml_m <- clipped$ColumnSubsetPropensityModel(
              base_estimator = ml_m,
              column_indices = as.integer(ps_indices_0b),
              lower = as.numeric(clip[1]), upper = as.numeric(clip[2]))
          }
        } else {
          ml_m <- clipped$FixedPropensityModel(as.numeric(propensity_score),
                                       lower = as.numeric(clip[1]),
                                       upper = as.numeric(clip[2]))
        }
        ml_final <- if (identical(model_final, "linear")) {
          factory$make_drlearner_final_linear(random_state = seed_py)
        } else {
          factory$make_regressor(model_final, random_state = seed_py,
                                 params = final_params)
        }
        DRLearner(model_propensity = ml_m, model_regression = ml_g,
                  model_final = ml_final, cv = as.integer(cv),
                  discrete_outcome = isTRUE(binary_outcome))
      }

      X_py <- np$asarray(X, dtype = "float64")
      Y_py <- np$asarray(Y, dtype = "float64")
      T_py <- np$asarray(as.numeric(T), dtype = "float64")

      dr <- build_dr()
      dr$fit(Y_py, T_py, X = X_py)

      tau <- as.numeric(reticulate::py_to_r(dr$effect(X_py, T0 = 0L, T1 = 1L)))
      po  <- reticulate::py_to_r(helpers$drlearner_potential_outcomes(dr, X_py))
      mu0 <- as.numeric(po[[1]])
      mu1 <- as.numeric(po[[2]])
      ate <- mean(tau)

      ate_ci <- NULL; tau_lo <- NULL; tau_hi <- NULL; ci_method <- "unavailable"
      if (isTRUE(bootstrap_ci)) {
        n <- length(Y)
        B <- as.integer(n_bootstrap)
        ate_boot <- numeric(B)
        tau_boot <- matrix(NA_real_, nrow = B, ncol = n)
        # Draw the resamples with numpy so a single np seed governs all
        # worker randomness (the worker has no seeded R RNG): one seeded
        # stream for the whole resampling, drawn from in turn.
        if (!is.null(seed_py)) np$random$seed(as.integer(seed_py))
        for (b in seq_len(B)) {
          idx <- as.integer(reticulate::py_to_r(
            np$random$choice(as.integer(n), size = as.integer(n), replace = TRUE))) + 1L
          Xb <- np$asarray(X[idx, , drop = FALSE], dtype = "float64")
          Yb <- np$asarray(Y[idx], dtype = "float64")
          Tb <- np$asarray(as.numeric(T)[idx], dtype = "float64")
          drb <- tryCatch({ d <- build_dr(); d$fit(Yb, Tb, X = Xb); d },
                          error = function(e) NULL)
          if (is.null(drb)) { ate_boot[b] <- NA_real_; next }
          eff <- as.numeric(reticulate::py_to_r(drb$effect(X_py, T0 = 0L, T1 = 1L)))
          tau_boot[b, ] <- eff
          ate_boot[b]   <- mean(eff)
        }
        ate_ci <- as.numeric(stats::quantile(ate_boot, c(0.025, 0.975),
                                              names = FALSE, na.rm = TRUE))
        tau_lo <- apply(tau_boot, 2, stats::quantile, probs = 0.025,
                        names = FALSE, na.rm = TRUE)
        tau_hi <- apply(tau_boot, 2, stats::quantile, probs = 0.975,
                        names = FALSE, na.rm = TRUE)
        ci_method <- "nonparametric bootstrap"
      } else if (identical(model_final, "linear")) {
        aci <- tryCatch({
          ci <- reticulate::py_to_r(dr$ate_interval(X_py, T0 = 0L, T1 = 1L, alpha = 0.05))
          c(as.numeric(ci[[1]]), as.numeric(ci[[2]]))
        }, error = function(e) NULL)
        eci <- tryCatch({
          e <- reticulate::py_to_r(dr$effect_interval(X_py, T0 = 0L, T1 = 1L, alpha = 0.05))
          list(lo = as.numeric(e[[1]]), hi = as.numeric(e[[2]]))
        }, error = function(e) NULL)
        ate_ci <- aci
        if (!is.null(eci)) { tau_lo <- eci$lo; tau_hi <- eci$hi }
        ci_method <- if (is.null(aci)) "unavailable" else "asymptotic Wald"
      }

      joblib$dump(dr, model_path)
      # Averaged over the cross-fit propensity models: econml does not record
      # which rows each fold held out, so a row cannot be paired with the one
      # model that did not train on it.
      ps_hat <- tryCatch(
        as.numeric(reticulate::py_to_r(helpers$drlearner_propensity(dr, X_py))),
        error = function(e) NULL)
      list(tau_hat = tau, ate = ate, mu0_hat = mu0, mu1_hat = mu1,
           ate_ci = ate_ci, tau_ci_lower = tau_lo, tau_ci_upper = tau_hi,
           ci_method = ci_method, ps = ps_hat)
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, Y = Y_vec, T = T_vec, covariate_names = covariate_names,
    model_propensity = model_propensity, model_outcome = model_outcome,
    model_final = model_final, cv = cv, clip = clip,
    propensity_score = propensity_score, ps_indices_0b = ps_indices_0b,
    seed_py = seed_py, binary_outcome = binary_outcome,
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap,
    outcome_params = params[["outcome"]] %||% list(),
    final_params   = params[["final"]] %||% list(),
    propensity_params = propensity_params
  )

  # ---- Assemble the cateomics_result (main process) ---------------------------
  tau_ci <- if (!is.null(fit_out$tau_ci_lower)) {
    data.frame(lower = fit_out$tau_ci_lower, upper = fit_out$tau_ci_upper)
  } else NULL

  # Surface the propensity for cf_diagnostics: the user-supplied PS when
  # known, otherwise the one extracted from econml above. NA only if that
  # extraction failed, which is what `ps_type` then records.
  n_obs <- length(Y_vec)
  ps_known <- !is.null(propensity_score)
  propensity_hat <- if (ps_known) {
    if (length(propensity_score) == 1L) rep(as.numeric(propensity_score), n_obs)
    else as.numeric(propensity_score)
  } else {
    .cf_validate_ps(fit_out$ps, n_obs) %||% rep(NA_real_, n_obs)
  }
  ps_type_fit <- if (ps_known) "known" else
                 if (anyNA(propensity_hat)) NA_character_ else "cross-fit (approx)"

  list(
    model = list(
      model_path = model_path,
      backend    = "econml.dr.DRLearner"
    ),
    res = list(
      tau_hat    = fit_out$tau_hat,
      ate        = fit_out$ate,
      mu0_hat    = fit_out$mu0_hat,
      mu1_hat    = fit_out$mu1_hat,
      ps         = propensity_hat,
      ps_type    = ps_type_fit,
      ate_ci     = fit_out$ate_ci,
      tau_ci     = tau_ci,
      # res$metadata holds small predict-config read by predict_cf_drlearner
      # (model_final selects the newdata CI path); cf_fit() lifts ci_method to
      # res$ci_method + meta$estimator and strips it from here.
      metadata = list(
        model_propensity = model_propensity,
        model_outcome = model_outcome,
        model_final      = model_final,
        cv               = as.integer(cv),
        bootstrap_ci     = isTRUE(bootstrap_ci),
        ci_method        = fit_out$ci_method,
        backend          = "econml.dr.DRLearner",
        note             = "mu1_hat - mu0_hat != tau_hat by design (Neyman-orthogonal CATE final stage vs plug-in nuisance regression)"
      )
    )
  )
}

#' Predict from a DRLearner
#'
#' @param object cf_model with \code{method = "drlearner"}.
#' @param newdata Optional new covariates.
#' @param type One of \code{"all"} (the default), \code{"ate"},
#'   \code{"ate_ci"}, \code{"tau"}, \code{"tau_ci"}, \code{"mu0"},
#'   \code{"mu1"}, \code{"ps"}, \code{"summary"}.
#' @param ... Ignored.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_drlearner <- function(object, newdata = NULL, type = "all", ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }

  if (type == "summary") {
    stop("type 'summary' is not supported with newdata.")
  }

  model_final  <- object$res$metadata$model_final %||% "ridge"
  ci_supported <- identical(model_final, "linear")
  ci_requested <- type %in% c("ate_ci", "tau_ci")
  if (ci_requested && !ci_supported) {
    rlang::abort(
      c(
        "DR-Learner CI on newdata requires `model_final = 'linear'`.",
        i = sprintf("Current model_final is '%s'.", model_final),
        i = "Refit with model_final = 'linear' (econml StatsModelsInference).",
        i = "Other CATE methods with newdata CI: grf (IJ variance), bcf (posterior)."
      ),
      class = "cateomics_drlearner_ci_unsupported_nonlinear_final"
    )
  }

  X_new <- as.matrix(newdata)
  covariate_names <- object$meta$covariate_names
  if (!is.null(covariate_names) && all(covariate_names %in% colnames(X_new))) {
    X_new <- X_new[, covariate_names, drop = FALSE]
  }

  pred <- .cf_py_run(
    function(pypath, model_path, X_new, want_ci) {
      sys    <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np     <- reticulate::import("numpy", convert = FALSE)
      joblib <- reticulate::import("joblib", convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      dr   <- joblib$load(model_path)
      X_py <- np$asarray(X_new, dtype = "float64")

      tau <- as.numeric(reticulate::py_to_r(dr$effect(X_py, T0 = 0L, T1 = 1L)))
      po  <- reticulate::py_to_r(helpers$drlearner_potential_outcomes(dr, X_py))
      mu0 <- as.numeric(po[[1]])
      mu1 <- as.numeric(po[[2]])
      ps  <- tryCatch(
        as.numeric(reticulate::py_to_r(helpers$drlearner_propensity(dr, X_py))),
        error = function(e) NULL)

      tau_lo <- NULL; tau_hi <- NULL
      if (isTRUE(want_ci)) {
        eci <- tryCatch({
          e <- reticulate::py_to_r(dr$effect_interval(X_py, T0 = 0L, T1 = 1L, alpha = 0.05))
          list(lo = as.numeric(e[[1]]), hi = as.numeric(e[[2]]))
        }, error = function(e) NULL)
        if (!is.null(eci)) { tau_lo <- eci$lo; tau_hi <- eci$hi }
      }
      list(tau = tau, mu0 = mu0, mu1 = mu1, tau_lo = tau_lo, tau_hi = tau_hi,
           ps = ps)
    },
    pypath = system.file("python", package = "CATEomics"),
    model_path = object$model$model$model_path,
    X_new = X_new, want_ci = ci_supported
  )

  ate_new    <- mean(pred$tau)
  ate_ci_new <- NULL
  tau_ci_new <- NULL
  if (!is.null(pred$tau_lo)) {
    tau_ci_new <- data.frame(lower = pred$tau_lo, upper = pred$tau_hi)
    ate_ci_new <- c(mean(pred$tau_lo), mean(pred$tau_hi))
  }

  n_new  <- nrow(X_new)
  ps_new <- .cf_validate_ps(pred$ps, n_new)
  # See predict_cf_xlearner: a per-observation known propensity does not carry
  # over to new rows.
  known_fit <- identical(object$res$ps_type %||% NA_character_, "known")
  if (known_fit && length(object$res$ps) > 1L &&
      length(unique(round(as.numeric(object$res$ps), 12L))) > 1L) {
    ps_new <- NULL
  }

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    ate_ci  = ate_ci_new,
    tau     = pred$tau,
    tau_ci  = tau_ci_new,
    mu0_hat = pred$mu0,
    mu1_hat = pred$mu1,
    ps      = ps_new %||% rep(NA_real_, n_new),
    ps_type = if (is.null(ps_new)) NA_character_ else
              if (known_fit) "known" else "out-of-sample"
  )
}
