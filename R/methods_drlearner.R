#' Fit a DR-Learner via econml
#'
#' Thin R wrapper around \code{econml.dr.DRLearner} using the shared cfomics
#' base-learner vocabulary for nuisance estimators (propensity, outcome,
#' final CATE).
#'
#' Returns the full \code{cfomics_result} contract. The ITE is taken from
#' the orthogonalized final stage (\code{dr$effect}); potential outcomes
#' \eqn{\hat y_0, \hat y_1} are obtained separately by averaging the
#' internal \code{models_regression} nuisance predictions over Monte Carlo
#' iterations and cross-validation folds (see
#' \code{inst/python/metalearner_helpers.py::drlearner_potential_outcomes}).
#'
#' @section Execution model:
#' All Python work runs inside the cfomics basilisk environment via
#' \code{.cf_py_run()} (a separate process). The fitted \code{DRLearner}
#' cannot be returned as a live object, so it is serialised with
#' \code{joblib.dump()} to a cache file during fit and reloaded at predict
#' time. The \code{model} pillar of the returned \code{cf_model} therefore
#' holds \code{list(model_path, backend)} rather than a live estimator.
#'
#' @section Potential outcomes (Class B': Outcome regression as intermediate nuisance):
#' \code{mu0_hat} / \code{mu1_hat} are the cross-fit nuisance regressions
#' (\code{models_regression}) averaged over Monte Carlo iterations.
#' \code{tau_hat} is the orthogonalised CATE projected by the final-stage
#' regression on pseudo-outcomes.
#'
#' \strong{NOTE}: \code{mu1_hat - mu0_hat != tau_hat} by design — the
#' Neyman-orthogonal final stage decouples CATE estimation from the
#' plug-in nuisance contrast. \code{ate} is reported as \code{mean(tau_hat)}
#' (orthogonal score), not \code{mean(mu1_hat - mu0_hat)}.
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
#' \strong{Option B semantics.} \eqn{\hat y_1 - \hat y_0} is the plug-in
#' contrast of the nuisance outcome models, while \code{tau_hat} is the
#' Neyman-orthogonal CATE from the final stage; these two quantities are
#' \emph{not} required to agree numerically.
#'
#' \strong{Linear base learners produce a constant CATE.} With
#' \code{model_outcome = model_final = "linear"}, the final-stage CATE
#' projection collapses onto a single-coefficient contrast; use
#' \code{"lasso"} / \code{"rf"} / \code{"elastic_net"} for at least the
#' final stage to recover heterogeneity. See the family-level discussion
#' in \code{\link{metalearners}}.
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
#'   Default \code{"ridge"} (a regularised linear final stage), matching
#'   econml's \code{DRLearner} convention of a linear \code{model_final}
#'   (\code{StatsModelsLinearRegression}). \strong{Asymmetry with
#'   \code{\link{cf_fit_xlearner}}}: this \code{model_final} does \emph{not}
#'   track \code{model_outcome} (it stays linear unless overridden), whereas
#'   xlearner's Stage-2 model defaults to \code{model_outcome}. Both mirror
#'   their respective econml estimator. Pass \code{model_final} explicitly
#'   (e.g. \code{"rf"}) to give the DR final stage a non-linear form.
#'   Default \code{"ridge"}. Use \code{"linear"} for asymptotic Wald CIs
#'   via \code{StatsModelsLinearRegression}.
#' @param cv Integer, number of cross-fitting folds. Default \code{5} (econml's
#'   own default is 2; cfomics uses 5 for lower-variance nuisance in the n < p
#'   regime).
#' @param bootstrap_ci Logical. When \code{TRUE}, override any native CI
#'   path and compute percentile-method CIs from \code{n_bootstrap} refits.
#'   Default \code{FALSE}.
#' @param n_bootstrap Integer, number of bootstrap replicates when
#'   \code{bootstrap_ci = TRUE}. Default \code{100}.
#' @param trim Length-2 numeric propensity clipping bounds.
#' @param propensity_score Optional known propensity score (RCT context);
#'   scalar or length-n vector. Mutually exclusive with \code{model_propensity}.
#' @param propensity_X Optional covariate subset (names or indices) the
#'   propensity model is restricted to.
#' @param random_state Integer random seed. When set, yields reproducible
#'   results. Forwarded to each sklearn nuisance learner and to numpy's RNG
#'   inside the basilisk process. Because that process is ephemeral, seeding
#'   cannot leak into the caller's session.
#' @param binary_outcome Logical. When \code{TRUE} the stage-1 outcome model
#'   becomes a probability classifier and DRLearner runs with
#'   \code{discrete_outcome = TRUE}.
#' @param ... Ignored (reserved for future arguments).
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cfomics_result} contract.
#' @keywords internal
cf_fit_drlearner <- function(X, T, Y,
                             covariate_names  = NULL,
                             model_propensity = "logistic",
                             model_outcome = "ridge",
                             model_final      = "ridge",
                             cv               = 5L,
                             bootstrap_ci     = FALSE,
                             n_bootstrap      = 100L,
                             trim             = c(0.025, 0.975),
                             propensity_score = NULL,
                             propensity_X     = NULL,
                             random_state     = NULL,
                             binary_outcome   = FALSE,
                             ...) {
  # ---- R-side validation / resolution (outside the basilisk process) --------
  .validate_learner_spec(model_propensity)
  .validate_learner_spec(model_outcome)
  .validate_learner_spec(model_final)

  X_mat <- as.matrix(X)
  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }
  colnames(X_mat) <- covariate_names
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)

  # propensity_X → 0-based column indices (resolved in the main process where
  # cfomics is loaded; passed to the worker as a plain integer vector).
  ps_indices_0b <- NULL
  if (!is.null(propensity_X)) {
    ps_indices_0b <- as.integer(
      .resolve_propensity_X(X_mat, propensity_X, covariate_names)$indices_0b
    )
  }

  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  model_path <- .cf_new_model_path("drlearner", "joblib")
  pypath     <- system.file("python", package = "cfomics")

  # ---- Python fit + serialize (inside the basilisk process) -----------------
  # The worker is a self-contained closure: it references only its arguments,
  # base R, and reticulate (cfomics is NOT loaded in the worker, so the
  # python helper path is passed in explicitly). `convert = FALSE` is required
  # so numpy arrays keep their Python identity (spike-proven).
  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, covariate_names,
             model_propensity, model_outcome, model_final, cv, trim,
             propensity_score, ps_indices_0b, seed_py, binary_outcome,
             bootstrap_ci, n_bootstrap) {
      sys     <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      clipped <- reticulate::import_from_path("clipped_propensity",
                                              path = pypath, convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      DRLearner <- reticulate::import("econml.dr", convert = FALSE)$DRLearner

      if (!is.null(seed_py)) np$random$seed(seed_py)

      # A fresh DRLearner build (econml clones the sklearn learners per fit,
      # so a new build per bootstrap replicate is the clean path).
      build_dr <- function() {
        ml_g <- if (isTRUE(binary_outcome)) {
          factory$make_classifier(model_outcome, random_state = seed_py)
        } else {
          factory$make_regressor(model_outcome, random_state = seed_py)
        }
        if (is.null(propensity_score)) {
          ml_m <- factory$make_classifier(model_propensity, random_state = seed_py)
          if (is.null(ps_indices_0b)) {
            ml_m <- clipped$ClippedPropensityModel(
              ml_m, lower = as.numeric(trim[1]), upper = as.numeric(trim[2]))
          } else {
            ml_m <- clipped$ColumnSubsetPropensityModel(
              base_estimator = ml_m,
              column_indices = as.integer(ps_indices_0b),
              lower = as.numeric(trim[1]), upper = as.numeric(trim[2]))
          }
        } else {
          ml_m <- clipped$FixedPropensityModel(as.numeric(propensity_score))
        }
        ml_final <- if (identical(model_final, "linear")) {
          factory$make_drlearner_final_linear(random_state = seed_py)
        } else {
          factory$make_regressor(model_final, random_state = seed_py)
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
        for (b in seq_len(B)) {
          # Draw the resample with numpy so a single np seed governs all
          # worker randomness (the ephemeral process has no seeded R RNG).
          if (!is.null(seed_py)) np$random$seed(as.integer(seed_py) + b)
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
      list(tau_hat = tau, ate = ate, mu0_hat = mu0, mu1_hat = mu1,
           ate_ci = ate_ci, tau_ci_lower = tau_lo, tau_ci_upper = tau_hi,
           ci_method = ci_method)
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, Y = Y_vec, T = T_vec, covariate_names = covariate_names,
    model_propensity = model_propensity, model_outcome = model_outcome,
    model_final = model_final, cv = cv, trim = trim,
    propensity_score = propensity_score, ps_indices_0b = ps_indices_0b,
    seed_py = seed_py, binary_outcome = binary_outcome,
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap
  )

  # ---- Assemble the cfomics_result (main process) ---------------------------
  tau_ci <- if (!is.null(fit_out$tau_ci_lower)) {
    data.frame(lower = fit_out$tau_ci_lower, upper = fit_out$tau_ci_upper)
  } else NULL

  # Surface a `propensity` field for cf_diagnostics. econml does not expose a
  # single cross-fit PS vector, so this surfaces the user-supplied PS when
  # known, otherwise NA.
  n_obs <- length(Y_vec)
  propensity_hat <- if (!is.null(propensity_score)) {
    if (length(propensity_score) == 1L) rep(as.numeric(propensity_score), n_obs)
    else as.numeric(propensity_score)
  } else {
    rep(NA_real_, n_obs)
  }

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
      propensity = propensity_hat,
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
#' @param type One of \code{"ate"}, \code{"tau"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param ... Ignored.
#' @return Depending on `type`: the full list of predictions when
#'   `type = "all"`, otherwise the requested component (a scalar ATE, a numeric
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
      class = "cfomics_drlearner_ci_unsupported_nonlinear_final"
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

      tau_lo <- NULL; tau_hi <- NULL
      if (isTRUE(want_ci)) {
        eci <- tryCatch({
          e <- reticulate::py_to_r(dr$effect_interval(X_py, T0 = 0L, T1 = 1L, alpha = 0.05))
          list(lo = as.numeric(e[[1]]), hi = as.numeric(e[[2]]))
        }, error = function(e) NULL)
        if (!is.null(eci)) { tau_lo <- eci$lo; tau_hi <- eci$hi }
      }
      list(tau = tau, mu0 = mu0, mu1 = mu1, tau_lo = tau_lo, tau_hi = tau_hi)
    },
    pypath = system.file("python", package = "cfomics"),
    model_path = object$model$model$model_path,
    X_new = X_new, want_ci = ci_supported
  )

  ate_new    <- mean(pred$tau)
  ate_ci_new <- NA_real_
  tau_ci_new <- NA_real_
  if (!is.null(pred$tau_lo)) {
    tau_ci_new <- data.frame(lower = pred$tau_lo, upper = pred$tau_hi)
    ate_ci_new <- c(mean(pred$tau_lo), mean(pred$tau_hi))
  }

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    ate_ci  = ate_ci_new,
    tau     = pred$tau,
    tau_ci  = tau_ci_new,
    mu0_hat = pred$mu0,
    mu1_hat = pred$mu1
  )
}
