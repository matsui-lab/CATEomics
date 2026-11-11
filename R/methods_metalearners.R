#' Metalearners (S/T/X-learner)
#'
#' S-learner fits a CATEomics interaction-augmented sklearn model directly;
#' T-learner and X-learner are thin R wrappers around
#' \code{econml.metalearners.TLearner} and \code{XLearner}.
#' Nuisance models are resolved through the shared CATEomics base-learner
#' vocabulary (see \code{inst/python/learner_factory.py}).
#'
#' All three metalearners return \code{tau_hat}, \code{mu0_hat},
#' \code{mu1_hat} and \code{ate}; \code{res$summary} is assembled centrally
#' by the entry that called them.
#'
#' @section Execution model:
#' All Python work runs inside the CATEomics basilisk environment via
#' \code{.cf_py_run()}. The fitted estimator is
#' serialised with \code{joblib.dump()} during fit and reloaded at predict
#' time, so the \code{model} pillar holds \code{list(model_path, backend)}
#' rather than a live estimator. The worker can be the calling R process
#' itself, so \code{.cf_py_run()} snapshots the Python-side generators and
#' restores them afterwards; seeding numpy inside the fit does not leak into
#' the caller's session.
#'
#' @section Base-learner backends:
#' Python-side methods (metalearners, \code{drlearner}) map the
#' shared base-learner names to sklearn via
#' \code{inst/python/learner_factory.py}. Numerical outputs are therefore not
#' bit-identical to the R-side \code{.fit_base_learner()} path for the same
#' token (CV fold assignment, regularization paths, solver conventions).
#'
#' @name metalearners
#' @keywords internal
NULL

# -------------------------------------------------------------------------
# Internal utilities (R-side)
# -------------------------------------------------------------------------

.build_meta_X <- function(X, covariate_names) {
  X_mat <- as.matrix(X)
  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }
  colnames(X_mat) <- covariate_names
  X_mat
}

# Assemble the two-pillar (model + res) return list shared by every
# metalearner. `fit_out` is the R-serialisable list returned by the basilisk
# fit worker; `extra_meta` are method-specific metadata fields.
.metalearner_assemble <- function(fit_out, model_path, backend,
                                   bootstrap_ci, extra_meta = list(),
                                   propensity = NULL,
                                   ps_type = NA_character_) {
  tau_ci <- if (!is.null(fit_out$tau_ci_lower)) {
    data.frame(lower = fit_out$tau_ci_lower, upper = fit_out$tau_ci_upper)
  } else NULL

  res <- list(
    tau_hat = fit_out$tau_hat,
    ate     = fit_out$ate,
    mu0_hat = fit_out$mu0_hat,
    mu1_hat = fit_out$mu1_hat,
    ate_ci  = fit_out$ate_ci,
    tau_ci  = tau_ci,
    metadata = c(
      extra_meta,
      list(
        bootstrap_ci = isTRUE(bootstrap_ci),
        ci_method    = fit_out$ci_method,
        backend      = backend
      )
    )
  )
  if (!is.null(propensity)) {
    res$ps      <- propensity
    res$ps_type <- ps_type
  }

  list(
    model = list(model_path = model_path, backend = backend),
    res   = res
  )
}

# -------------------------------------------------------------------------
# S-Learner
# -------------------------------------------------------------------------

#' Fit an S-learner
#'
#' Fits a single outcome model on the interaction-augmented
#' \code{[X | T | X*T]} matrix, so the fitted surface is not constant in
#' \code{X} even for a linear base learner. Potential outcomes are read back
#' by setting the \code{T} column to 1 / 0.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param model_outcome Base-learner name for the single outcome model.
#'   Default \code{"ridge"}.
#' @param bootstrap_ci Logical. When \code{TRUE}, compute percentile CIs from
#'   \code{n_bootstrap} refits (inside the basilisk process). The resampling
#'   is implemented in CATEomics and is the same for every method that accepts
#'   it. Default \code{FALSE}.
#' @param n_bootstrap Integer bootstrap replicates. Default \code{100}.
#' @param random_state Integer seed forwarded to the sklearn learner and
#'   numpy inside the basilisk process.
#' @param binary_outcome Logical. When \code{TRUE} the outcome model is a
#'   probability classifier so \code{mu1 - mu0} is a risk difference.
#' @param params Named list with an optional \code{outcome} (and, for
#'   x-learner, \code{final}) block of hyperparameters for the base learners.
#'   A block reaches the scikit-learn estimator built in
#'   \code{inst/python/learner_factory.py}, which raises on a hyperparameter
#'   it does not take. A flat entry is refused on the R side, because these
#'   methods read the blocks only and it would have nowhere to go.
#'   A value the Python side requires to be whole is written as an R
#'   integer (\code{10L}, not \code{10}), which R otherwise makes a double.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @references
#' Künzel, Sekhon, Bickel & Yu (2019) Metalearners for estimating heterogeneous
#' treatment effects using machine learning, \emph{PNAS} 116(10):4156-4165.
#'
#' @keywords internal
cf_fit_slearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome   = "ridge",
                            bootstrap_ci    = FALSE,
                            n_bootstrap     = 100L,
                            random_state    = NULL,
                            binary_outcome  = FALSE,
                            params          = list()) {
  .validate_learner_spec(model_outcome)
  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)
  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  outcome_params <- params[["outcome"]] %||% list()
  model_path <- .cf_new_model_path("slearner", "joblib")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, binary_outcome,
             seed_py, bootstrap_ci, n_bootstrap, outcome_params) {
      sys     <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      interact <- reticulate::import_from_path("interaction_augmented",
                                               path = pypath, convert = FALSE)
      if (!is.null(seed_py)) np$random$seed(seed_py)

      # The S-learner design is [X, T, X * T]; the wrapper builds it from
      # [X, T], which is the same expansion the DR outcome nuisance now uses.
      build_fit <- function(Xm, Yv, Tv) {
        est <- interact$make_interaction_augmented(
          factory$make_outcome_model(model_outcome,
                                     binary = isTRUE(binary_outcome),
                                     random_state = seed_py,
                                     params = outcome_params))
        Xp  <- np$asarray(Xm, dtype = "float64")
        Tp  <- np$asarray(matrix(as.numeric(Tv), ncol = 1L), dtype = "float64")
        est$fit(np$hstack(list(Xp, Tp)),
                np$asarray(as.numeric(Yv), dtype = "float64"))
        est
      }
      po_of <- function(est, Xm) {
        Xp <- np$asarray(Xm, dtype = "float64")
        m  <- nrow(Xm)
        ones  <- np$ones(reticulate::r_to_py(list(m, 1L)),  dtype = "float64")
        zeros <- np$zeros(reticulate::r_to_py(list(m, 1L)), dtype = "float64")
        list(mu0 = as.numeric(reticulate::py_to_r(
               est$predict(np$hstack(list(Xp, zeros))))),
             mu1 = as.numeric(reticulate::py_to_r(
               est$predict(np$hstack(list(Xp, ones))))))
      }

      tau_on  <- function(est, Xm) { p <- po_of(est, Xm); p$mu1 - p$mu0 }
      boot_ci <- function() {
        n <- length(Y); B <- as.integer(n_bootstrap)
        ab <- numeric(B); tb <- matrix(NA_real_, B, n)
        # One seeded stream for the whole resampling, drawn from in turn.
        if (!is.null(seed_py)) np$random$seed(as.integer(seed_py))
        for (b in seq_len(B)) {
          idx <- as.integer(reticulate::py_to_r(np$random$choice(
            as.integer(n), size = as.integer(n), replace = TRUE))) + 1L
          eb <- tryCatch(build_fit(X[idx, , drop = FALSE], Y[idx], T[idx]),
                         error = function(e) NULL)
          if (is.null(eb)) { ab[b] <- NA_real_; next }
          tt <- tau_on(eb, X); tb[b, ] <- tt; ab[b] <- mean(tt)
        }
        list(ate_ci = as.numeric(stats::quantile(ab, c(0.025, 0.975),
                                                 names = FALSE, na.rm = TRUE)),
             tau_lo = apply(tb, 2, stats::quantile, probs = 0.025,
                            names = FALSE, na.rm = TRUE),
             tau_hi = apply(tb, 2, stats::quantile, probs = 0.975,
                            names = FALSE, na.rm = TRUE))
      }

      est <- build_fit(X, Y, T)
      po  <- po_of(est, X)
      tau <- po$mu1 - po$mu0
      ate <- mean(tau)

      ate_ci <- NULL; tau_lo <- NULL; tau_hi <- NULL; ci_method <- "unavailable"
      if (isTRUE(bootstrap_ci)) {
        bb <- boot_ci()
        ate_ci <- bb$ate_ci; tau_lo <- bb$tau_lo; tau_hi <- bb$tau_hi
        ci_method <- "nonparametric bootstrap"
      }

      joblib$dump(est, model_path)
      list(tau_hat = tau, ate = ate, mu0_hat = po$mu0, mu1_hat = po$mu1,
           ate_ci = ate_ci, tau_ci_lower = tau_lo, tau_ci_upper = tau_hi,
           ci_method = ci_method)
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, Y = Y_vec, T = T_vec, model_outcome = model_outcome,
    binary_outcome = binary_outcome, seed_py = seed_py,
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap,
    outcome_params = outcome_params
  )

  .metalearner_assemble(
    fit_out, model_path,
    backend = "CATEomics.slearner_python_sklearn_interactions",
    bootstrap_ci = bootstrap_ci,
    extra_meta = list(model_outcome = model_outcome)
  )
}

#' Predict from an S-learner
#' @param object cf_model with \code{method = "slearner"}.
#' @param newdata Optional new covariates.
#' @param type One of \code{"all"} (the default), \code{"ate"},
#'   \code{"ate_ci"}, \code{"tau"}, \code{"tau_ci"}, \code{"mu0"},
#'   \code{"mu1"}, \code{"ps"}, \code{"summary"}.
#' @param ... Ignored.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_slearner <- function(object, newdata = NULL, type = "all", ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }
  if (type %in% c("ate_ci", "tau_ci")) {
    .abort_ite_ci_unsupported("slearner")
  }
  nd <- as.matrix(newdata)
  covariate_names <- object$meta$covariate_names
  if (!is.null(covariate_names) && all(covariate_names %in% colnames(nd))) {
    nd <- nd[, covariate_names, drop = FALSE]
  }

  pred <- .cf_py_run(
    function(pypath, model_path, X) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np     <- reticulate::import("numpy", convert = FALSE)
      joblib <- reticulate::import("joblib", convert = FALSE)
      # The saved estimator is the interaction-augmented wrapper, so it takes
      # [X, T] and builds the products itself.
      est <- joblib$load(model_path)
      Xp  <- np$asarray(X, dtype = "float64")
      m   <- nrow(X)
      ones  <- np$ones(reticulate::r_to_py(list(m, 1L)),  dtype = "float64")
      zeros <- np$zeros(reticulate::r_to_py(list(m, 1L)), dtype = "float64")
      list(mu0 = as.numeric(reticulate::py_to_r(
             est$predict(np$hstack(list(Xp, zeros))))),
           mu1 = as.numeric(reticulate::py_to_r(
             est$predict(np$hstack(list(Xp, ones))))))
    },
    pypath = system.file("python", package = "CATEomics"),
    model_path = object$model$model$model_path, X = nd
  )
  ite_new <- pred$mu1 - pred$mu0
  .assemble_predict_output(
    type = type, ate = mean(ite_new), ate_ci = NULL,
    tau = ite_new, tau_ci = NULL, mu0_hat = pred$mu0, mu1_hat = pred$mu1
  )
}

# -------------------------------------------------------------------------
# T-Learner
# -------------------------------------------------------------------------

#' Fit a T-learner via econml.metalearners
#'
#' Fits two separate outcome regressions and reports
#' \eqn{\hat\mu_1(X) - \hat\mu_0(X)} as the CATE.
#'
#' @param X,T,Y,covariate_names,bootstrap_ci,n_bootstrap,random_state,binary_outcome,params
#'   See \code{\link{cf_fit_slearner}}.
#' @param model_outcome Base-learner name for both arms. Default \code{"ridge"}.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @references
#' Künzel, Sekhon, Bickel & Yu (2019) Metalearners for estimating heterogeneous
#' treatment effects using machine learning, \emph{PNAS} 116(10):4156-4165.
#'
#' Battocchi, Dillon, Hei, Lewis, Oka, Oprescu & Syrgkanis (2019) EconML: A
#' Python Package for ML-Based Heterogeneous Treatment Effects Estimation,
#' \url{https://github.com/py-why/EconML} (implementation).
#'
#' @keywords internal
cf_fit_tlearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome   = "ridge",
                            bootstrap_ci    = FALSE,
                            n_bootstrap     = 100L,
                            random_state    = NULL,
                            binary_outcome  = FALSE,
                            params          = list()) {
  .validate_learner_spec(model_outcome)
  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)
  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  outcome_params <- params[["outcome"]] %||% list()
  model_path <- .cf_new_model_path("tlearner", "joblib")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, binary_outcome,
             seed_py, bootstrap_ci, n_bootstrap, outcome_params) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      metal   <- reticulate::import("econml.metalearners", convert = FALSE)
      if (!is.null(seed_py)) np$random$seed(seed_py)

      build_fit <- function(Xm, Yv, Tv) {
        arm <- factory$make_outcome_model(model_outcome,
                                          binary = isTRUE(binary_outcome),
                                          random_state = seed_py,
                                          params = outcome_params)
        est <- metal$TLearner(models = arm)
        est$fit(np$asarray(as.numeric(Yv), dtype = "float64"),
                np$asarray(as.integer(Tv), dtype = "float64"),
                X = np$asarray(Xm, dtype = "float64"))
        est
      }
      po_of <- function(est, Xm) {
        Xp  <- np$asarray(Xm, dtype = "float64")
        tau <- as.numeric(reticulate::py_to_r(est$effect(Xp, T0 = 0L, T1 = 1L)))
        po  <- reticulate::py_to_r(helpers$tlearner_potential_outcomes(est, Xp))
        list(tau = tau, mu0 = as.numeric(po[[1]]), mu1 = as.numeric(po[[2]]))
      }
      boot_ci <- function() {
        n <- length(Y); B <- as.integer(n_bootstrap)
        ab <- numeric(B); tb <- matrix(NA_real_, B, n)
        # One seeded stream for the whole resampling, drawn from in turn.
        if (!is.null(seed_py)) np$random$seed(as.integer(seed_py))
        for (b in seq_len(B)) {
          idx <- as.integer(reticulate::py_to_r(np$random$choice(
            as.integer(n), size = as.integer(n), replace = TRUE))) + 1L
          eb <- tryCatch(build_fit(X[idx, , drop = FALSE], Y[idx], T[idx]),
                         error = function(e) NULL)
          if (is.null(eb)) { ab[b] <- NA_real_; next }
          tt <- po_of(eb, X)$tau; tb[b, ] <- tt; ab[b] <- mean(tt)
        }
        list(ate_ci = as.numeric(stats::quantile(ab, c(0.025, 0.975),
                                                 names = FALSE, na.rm = TRUE)),
             tau_lo = apply(tb, 2, stats::quantile, probs = 0.025,
                            names = FALSE, na.rm = TRUE),
             tau_hi = apply(tb, 2, stats::quantile, probs = 0.975,
                            names = FALSE, na.rm = TRUE))
      }

      est <- build_fit(X, Y, T)
      p0  <- po_of(est, X)
      ate <- mean(p0$tau)

      ate_ci <- NULL; tau_lo <- NULL; tau_hi <- NULL; ci_method <- "unavailable"
      if (isTRUE(bootstrap_ci)) {
        bb <- boot_ci()
        ate_ci <- bb$ate_ci; tau_lo <- bb$tau_lo; tau_hi <- bb$tau_hi
        ci_method <- "nonparametric bootstrap"
      }

      joblib$dump(est, model_path)
      list(tau_hat = p0$tau, ate = ate, mu0_hat = p0$mu0, mu1_hat = p0$mu1,
           ate_ci = ate_ci, tau_ci_lower = tau_lo, tau_ci_upper = tau_hi,
           ci_method = ci_method)
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, Y = Y_vec, T = T_vec, model_outcome = model_outcome,
    binary_outcome = binary_outcome, seed_py = seed_py,
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap,
    outcome_params = outcome_params
  )

  .metalearner_assemble(
    fit_out, model_path, backend = "econml.metalearners.TLearner",
    bootstrap_ci = bootstrap_ci, extra_meta = list(model_outcome = model_outcome)
  )
}

#' Predict from a T-learner
#' @inheritParams predict_cf_slearner
#' @param object cf_model with \code{method = "tlearner"}.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_tlearner <- function(object, newdata = NULL, type = "all", ...) {
  .predict_metalearner_generic(object, newdata, type,
                               po_fn_name = "tlearner_potential_outcomes",
                               method_name = "tlearner")
}

# -------------------------------------------------------------------------
# X-Learner
# -------------------------------------------------------------------------

#' Fit an X-learner via econml.metalearners
#'
#' Two-stage metalearner (Künzel et al. 2019). \code{mu0_hat} / \code{mu1_hat}
#' are the Stage-1 outcome fits; \code{tau_hat} is the Stage-3
#' propensity-weighted blend of the two imputed-effect fits, which is not the
#' difference of the Stage-1 fits, so \code{mu1 - mu0 != tau_hat}.
#'
#' @param X,T,Y,covariate_names,bootstrap_ci,n_bootstrap,random_state,binary_outcome,params
#'   See \code{\link{cf_fit_slearner}}.
#' @param propensity_params Hyperparameters for the propensity model.
#' @param model_outcome Base-learner for Stage-1 outcome regressions. Default
#'   \code{"ridge"}.
#' @param model_final Base-learner for the Stage-2 CATE model.
#'   Default \code{NULL}, which resolves to \code{model_outcome} — i.e. the
#'   Stage-2 CATE learner \strong{tracks the base learner}.
#'   \strong{Asymmetry with \code{\link{cf_fit_drlearner}}}:
#'   xlearner's final stage follows \code{model_outcome}, whereas drlearner's
#'   \code{model_final} does not and defaults to \code{"ridge"}. So passing
#'   only \code{model_outcome} changes xlearner's final stage but not
#'   drlearner's.

#' @param model_propensity Base-learner for the propensity model used in the
#'   final blend.
#'   Default \code{"logistic"}.
#' @param clip Length-2 propensity clipping bounds.
#' @param propensity_score Optional known PS (scalar or length-n).
#' @param propensity_X Optional covariate subset for the propensity model.
#'   With \code{propensity_X} unset this is the clinical covariates alone;
#'   see \code{\link{cf_fit}}.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @references
#' Künzel, Sekhon, Bickel & Yu (2019) Metalearners for estimating heterogeneous
#' treatment effects using machine learning, \emph{PNAS} 116(10):4156-4165.
#'
#' Battocchi, Dillon, Hei, Lewis, Oka, Oprescu & Syrgkanis (2019) EconML: A
#' Python Package for ML-Based Heterogeneous Treatment Effects Estimation,
#' \url{https://github.com/py-why/EconML} (implementation).
#'
#' @keywords internal
cf_fit_xlearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome    = "ridge",
                            model_final      = NULL,
                            model_propensity = "logistic",
                            bootstrap_ci     = FALSE,
                            n_bootstrap      = 100L,
                            clip             = .CF_PS_CLIP,
                            propensity_score = NULL,
                            propensity_X     = NULL,
                            random_state     = NULL,
                            binary_outcome   = FALSE,
                            params           = list(),
                            propensity_params = list()) {
  model_final <- model_final %||% model_outcome
  .validate_learner_spec(model_outcome)
  .validate_learner_spec(model_final)
  .validate_learner_spec(model_propensity)
  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)
  seed_py <- if (is.null(random_state)) NULL else as.integer(random_state)
  ps_indices_0b <- NULL
  if (!is.null(propensity_X)) {
    ps_indices_0b <- as.integer(
      .resolve_propensity_X(X_mat, propensity_X, colnames(X_mat))$indices_0b)
  }
  outcome_params <- params[["outcome"]] %||% list()
  final_params   <- params[["final"]] %||% list()
  if (!is.null(propensity_score) && length(propensity_params)) {
    rlang::abort(
      c("`propensity_params` has no effect when `propensity_score` is supplied.",
        i = "A known propensity is used as given; no propensity model is fitted."),
      class = "cateomics_propensity_params_unused"
    )
  }
  model_path <- .cf_new_model_path("xlearner", "joblib")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, model_final,
             model_propensity, clip, propensity_score, ps_indices_0b,
             binary_outcome, seed_py, bootstrap_ci, n_bootstrap,
             outcome_params, final_params, propensity_params) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      clipped <- reticulate::import_from_path("clipped_propensity",
                                              path = pypath, convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      metal   <- reticulate::import("econml.metalearners", convert = FALSE)
      if (!is.null(seed_py)) np$random$seed(seed_py)

      build_ps <- function() {
        if (is.null(propensity_score)) {
          pm <- factory$make_classifier(model_propensity, random_state = seed_py,
                                        params = propensity_params)
          if (is.null(ps_indices_0b)) {
            clipped$ClippedPropensityModel(pm, lower = as.numeric(clip[1]),
                                           upper = as.numeric(clip[2]))
          } else {
            clipped$ColumnSubsetPropensityModel(
              base_estimator = pm,
              column_indices = as.integer(ps_indices_0b),
              lower = as.numeric(clip[1]), upper = as.numeric(clip[2]))
          }
        } else {
          clipped$FixedPropensityModel(as.numeric(propensity_score),
                                       lower = as.numeric(clip[1]),
                                       upper = as.numeric(clip[2]))
        }
      }
      build_fit <- function(Xm, Yv, Tv) {
        outcome <- factory$make_outcome_model(model_outcome,
                                              binary = isTRUE(binary_outcome),
                                              random_state = seed_py,
                                              params = outcome_params)
        cate <- factory$make_regressor(model_final, random_state = seed_py,
                                       params = final_params)
        est <- metal$XLearner(models = outcome, cate_models = cate,
                              propensity_model = build_ps())
        est$fit(np$asarray(as.numeric(Yv), dtype = "float64"),
                np$asarray(as.integer(Tv), dtype = "float64"),
                X = np$asarray(Xm, dtype = "float64"))
        est
      }
      po_of <- function(est, Xm) {
        Xp  <- np$asarray(Xm, dtype = "float64")
        tau <- as.numeric(reticulate::py_to_r(est$effect(Xp, T0 = 0L, T1 = 1L)))
        po  <- reticulate::py_to_r(helpers$xlearner_potential_outcomes(est, Xp))
        list(tau = tau, mu0 = as.numeric(po[[1]]), mu1 = as.numeric(po[[2]]))
      }
      boot_ci <- function() {
        n <- length(Y); B <- as.integer(n_bootstrap)
        ab <- numeric(B); tb <- matrix(NA_real_, B, n)
        # One seeded stream for the whole resampling, drawn from in turn.
        if (!is.null(seed_py)) np$random$seed(as.integer(seed_py))
        for (b in seq_len(B)) {
          idx <- as.integer(reticulate::py_to_r(np$random$choice(
            as.integer(n), size = as.integer(n), replace = TRUE))) + 1L
          eb <- tryCatch(build_fit(X[idx, , drop = FALSE], Y[idx], T[idx]),
                         error = function(e) NULL)
          if (is.null(eb)) { ab[b] <- NA_real_; next }
          tt <- po_of(eb, X)$tau; tb[b, ] <- tt; ab[b] <- mean(tt)
        }
        list(ate_ci = as.numeric(stats::quantile(ab, c(0.025, 0.975),
                                                 names = FALSE, na.rm = TRUE)),
             tau_lo = apply(tb, 2, stats::quantile, probs = 0.025,
                            names = FALSE, na.rm = TRUE),
             tau_hi = apply(tb, 2, stats::quantile, probs = 0.975,
                            names = FALSE, na.rm = TRUE))
      }

      est <- build_fit(X, Y, T)
      p0  <- po_of(est, X)
      ate <- mean(p0$tau)

      ate_ci <- NULL; tau_lo <- NULL; tau_hi <- NULL; ci_method <- "unavailable"
      if (isTRUE(bootstrap_ci)) {
        bb <- boot_ci()
        ate_ci <- bb$ate_ci; tau_lo <- bb$tau_lo; tau_hi <- bb$tau_hi
        ci_method <- "nonparametric bootstrap"
      }

      joblib$dump(est, model_path)
      # The X-learner weights its two imputed-effect surfaces by this PS, so it
      # is the propensity the fitted model actually used.
      ps_hat <- tryCatch(
        as.numeric(reticulate::py_to_r(
          helpers$xlearner_propensity(est, np$asarray(X, dtype = "float64")))),
        error = function(e) NULL)
      list(tau_hat = p0$tau, ate = ate, mu0_hat = p0$mu0, mu1_hat = p0$mu1,
           ate_ci = ate_ci, tau_ci_lower = tau_lo, tau_ci_upper = tau_hi,
           ci_method = ci_method, ps = ps_hat)
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, Y = Y_vec, T = T_vec, model_outcome = model_outcome,
    model_final = model_final, model_propensity = model_propensity,
    clip = clip, propensity_score = propensity_score,
    ps_indices_0b = ps_indices_0b, binary_outcome = binary_outcome,
    seed_py = seed_py, bootstrap_ci = isTRUE(bootstrap_ci),
    n_bootstrap = n_bootstrap, outcome_params = outcome_params,
    final_params = final_params, propensity_params = propensity_params
  )

  # Known PS when supplied, otherwise the PS the fitted X-learner estimated
  # (in-sample: the propensity model saw every row).
  n_obs <- length(Y_vec)
  ps_known <- !is.null(propensity_score)
  propensity_hat <- if (ps_known) {
    if (length(propensity_score) == 1L) rep(as.numeric(propensity_score), n_obs)
    else as.numeric(propensity_score)
  } else {
    .cf_validate_ps(fit_out$ps, n_obs) %||% rep(NA_real_, n_obs)
  }
  ps_type_fit <- if (ps_known) "known" else
                 if (anyNA(propensity_hat)) NA_character_ else "in-sample"

  .metalearner_assemble(
    fit_out, model_path, backend = "econml.metalearners.XLearner",
    bootstrap_ci = bootstrap_ci,
    extra_meta = list(
      model_outcome = model_outcome, model_propensity = model_propensity,
      note = "mu1_hat - mu0_hat != tau_hat by design (propensity-weighted Stage 3 blend, not direct Stage 1 difference)"
    ),
    propensity = propensity_hat,
    ps_type    = ps_type_fit
  )
}

#' Predict from an X-learner
#' @inheritParams predict_cf_slearner
#' @param object cf_model with \code{method = "xlearner"}.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_xlearner <- function(object, newdata = NULL, type = "all", ...) {
  .predict_metalearner_generic(object, newdata, type,
                               po_fn_name = "xlearner_potential_outcomes",
                               ps_fn_name = "xlearner_propensity",
                               method_name = "xlearner")
}

# -------------------------------------------------------------------------
# Shared predict dispatcher (T / X, econml effect + potential outcomes)
# -------------------------------------------------------------------------

.predict_metalearner_generic <- function(object, newdata, type, po_fn_name,
                                          ps_fn_name = NULL,
                                          method_name = NULL) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }
  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.",
                 class = "cateomics_unsupported_predict_type")
  }
  if (type %in% c("ate_ci", "tau_ci")) {
    .abort_ite_ci_unsupported(method_name %||% "metalearner")
  }

  X_new <- as.matrix(newdata)
  covariate_names <- object$meta$covariate_names
  if (!is.null(covariate_names) && all(covariate_names %in% colnames(X_new))) {
    X_new <- X_new[, covariate_names, drop = FALSE]
  }

  pred <- .cf_py_run(
    function(pypath, model_path, X, po_fn_name, ps_fn_name) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      helpers <- reticulate::import_from_path("metalearner_helpers",
                                              path = pypath, convert = FALSE)
      est <- joblib$load(model_path)
      Xp  <- np$asarray(X, dtype = "float64")
      tau <- as.numeric(reticulate::py_to_r(est$effect(Xp, T0 = 0L, T1 = 1L)))
      po  <- reticulate::py_to_r(helpers[[po_fn_name]](est, Xp))
      ps  <- if (is.null(ps_fn_name)) NULL else tryCatch(
        as.numeric(reticulate::py_to_r(helpers[[ps_fn_name]](est, Xp))),
        error = function(e) NULL)
      list(tau = tau, mu0 = as.numeric(po[[1]]), mu1 = as.numeric(po[[2]]),
           ps = ps)
    },
    pypath = system.file("python", package = "CATEomics"),
    model_path = object$model$model$model_path, X = X_new,
    po_fn_name = po_fn_name, ps_fn_name = ps_fn_name
  )

  n_new  <- nrow(X_new)
  ps_new <- .cf_validate_ps(pred$ps, n_new)

  # A per-observation known propensity cannot be mapped to new rows, and the
  # fixed-propensity model would hand back the training vector whenever the row
  # counts happen to agree, so drop it instead of pairing it with the wrong
  # rows. A scalar known propensity is the same value everywhere and stays.
  fit_ps      <- object$res$ps
  known_fit   <- identical(object$res$ps_type %||% NA_character_, "known")
  known_vector <- known_fit && length(fit_ps) > 1L &&
                  length(unique(round(as.numeric(fit_ps), 12L))) > 1L
  if (known_vector) ps_new <- NULL

  .assemble_predict_output(
    type = type, ate = mean(pred$tau), ate_ci = NULL,
    tau = pred$tau, tau_ci = NULL, mu0_hat = pred$mu0, mu1_hat = pred$mu1,
    ps      = ps_new %||% rep(NA_real_, n_new),
    ps_type = if (is.null(ps_new)) NA_character_ else
              if (known_fit) "known" else "out-of-sample"
  )
}
