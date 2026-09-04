#' Metalearners (S/T/X-learner) via econml
#'
#' Thin R wrappers around \code{econml.metalearners.SLearner} (a cfomics
#' interaction-augmented sklearn model), \code{TLearner}, and \code{XLearner}.
#' Nuisance models are resolved through the shared cfomics base-learner
#' vocabulary (see \code{inst/python/learner_factory.py}).
#'
#' All three metalearners return the full \code{cfomics_result} contract:
#' \code{tau_hat}, \code{mu0_hat}, \code{mu1_hat} plus \code{ate} and summary.
#'
#' @section Execution model:
#' All Python work runs inside the cfomics basilisk environment via
#' \code{.cf_py_run()} (a separate process). The fitted estimator is
#' serialised with \code{joblib.dump()} during fit and reloaded at predict
#' time, so the \code{model} pillar holds \code{list(model_path, backend)}
#' rather than a live estimator. Because the worker process is ephemeral,
#' seeding numpy inside it cannot leak into the caller's session and no
#' RNG save/restore is needed.
#'
#' @section CATE constancy with linear base learners:
#' When every nuisance model is \code{"linear"}, the resulting CATE surface
#' collapses to a constant. If you need heterogeneous CATE, use
#' \code{"lasso"} / \code{"rf"} / \code{"elastic_net"} for at least the
#' outcome stage, or switch to \code{grf}.
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
                                   propensity = NULL) {
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
  if (!is.null(propensity)) res$propensity <- propensity

  list(
    model = list(model_path = model_path, backend = backend),
    res   = res
  )
}

# -------------------------------------------------------------------------
# S-Learner
# -------------------------------------------------------------------------

#' Fit an S-learner via econml
#'
#' Fits a single outcome model on the interaction-augmented
#' \code{[X | T | X*T]} matrix (a cfomics S-learner variant that guarantees
#' heterogeneous CATE for any base learner). Potential outcomes are read back
#' by setting the \code{T} column to 1 / 0.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param model_outcome Base-learner name for the single outcome model.
#'   Default \code{"ridge"}.
#' @param bootstrap_ci Logical. When \code{TRUE}, compute percentile CIs from
#'   \code{n_bootstrap} refits (inside the basilisk process). Default
#'   \code{FALSE}.
#' @param n_bootstrap Integer bootstrap replicates. Default \code{100}.
#' @param random_state Integer seed forwarded to the sklearn learner and
#'   numpy inside the basilisk process.
#' @param binary_outcome Logical. When \code{TRUE} the outcome model is a
#'   probability classifier so \code{mu1 - mu0} is a risk difference.
#' @param ... Ignored.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cfomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @keywords internal
cf_fit_slearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome   = "ridge",
                            bootstrap_ci    = FALSE,
                            n_bootstrap     = 100L,
                            random_state    = NULL,
                            binary_outcome  = FALSE,
                            ...) {
  .validate_learner_spec(model_outcome)
  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)
  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  model_path <- .cf_new_model_path("slearner", "joblib")
  pypath     <- system.file("python", package = "cfomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, binary_outcome,
             seed_py, bootstrap_ci, n_bootstrap) {
      sys     <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      np      <- reticulate::import("numpy", convert = FALSE)
      joblib  <- reticulate::import("joblib", convert = FALSE)
      factory <- reticulate::import_from_path("learner_factory",
                                              path = pypath, convert = FALSE)
      if (!is.null(seed_py)) np$random$seed(seed_py)

      build_fit <- function(Xm, Yv, Tv) {
        est <- factory$make_outcome_model(model_outcome,
                                          binary = isTRUE(binary_outcome),
                                          random_state = seed_py)
        Xp  <- np$asarray(Xm, dtype = "float64")
        Tp  <- np$asarray(matrix(as.numeric(Tv), ncol = 1L), dtype = "float64")
        Xint <- np$hstack(list(Xp, Tp, np$multiply(Xp, Tp)))
        est$fit(Xint, np$asarray(as.numeric(Yv), dtype = "float64"))
        est
      }
      po_of <- function(est, Xm) {
        Xp <- np$asarray(Xm, dtype = "float64")
        m  <- nrow(Xm)
        ones  <- np$ones(reticulate::r_to_py(list(m, 1L)),  dtype = "float64")
        zeros <- np$zeros(reticulate::r_to_py(list(m, 1L)), dtype = "float64")
        X1 <- np$hstack(list(Xp, ones,  np$multiply(Xp, ones)))
        X0 <- np$hstack(list(Xp, zeros, np$multiply(Xp, zeros)))
        list(mu0 = as.numeric(reticulate::py_to_r(est$predict(X0))),
             mu1 = as.numeric(reticulate::py_to_r(est$predict(X1))))
      }

      tau_on  <- function(est, Xm) { p <- po_of(est, Xm); p$mu1 - p$mu0 }
      boot_ci <- function() {
        n <- length(Y); B <- as.integer(n_bootstrap)
        ab <- numeric(B); tb <- matrix(NA_real_, B, n)
        for (b in seq_len(B)) {
          if (!is.null(seed_py)) np$random$seed(as.integer(seed_py) + b)
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
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap
  )

  .metalearner_assemble(
    fit_out, model_path,
    backend = "cfomics.slearner_python_sklearn_interactions",
    bootstrap_ci = bootstrap_ci,
    extra_meta = list(model_outcome = model_outcome)
  )
}

#' Predict from an S-learner
#' @param object cf_model with \code{method = "slearner"}.
#' @param newdata Optional new covariates.
#' @param type One of \code{"ate"}, \code{"tau"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param ... Ignored.
#' @return Depending on `type`: the full list of predictions when
#'   `type = "all"`, otherwise the requested component (a scalar ATE, a numeric
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
      est <- joblib$load(model_path)
      Xp  <- np$asarray(X, dtype = "float64")
      m   <- nrow(X)
      ones  <- np$ones(reticulate::r_to_py(list(m, 1L)),  dtype = "float64")
      zeros <- np$zeros(reticulate::r_to_py(list(m, 1L)), dtype = "float64")
      X1 <- np$hstack(list(Xp, ones,  np$multiply(Xp, ones)))
      X0 <- np$hstack(list(Xp, zeros, np$multiply(Xp, zeros)))
      list(mu0 = as.numeric(reticulate::py_to_r(est$predict(X0))),
           mu1 = as.numeric(reticulate::py_to_r(est$predict(X1))))
    },
    pypath = system.file("python", package = "cfomics"),
    model_path = object$model$model$model_path, X = nd
  )
  ite_new <- pred$mu1 - pred$mu0
  .assemble_predict_output(
    type = type, ate = mean(ite_new), ate_ci = NA_real_,
    tau = ite_new, tau_ci = NA_real_, mu0_hat = pred$mu0, mu1_hat = pred$mu1
  )
}

# -------------------------------------------------------------------------
# T-Learner
# -------------------------------------------------------------------------

#' Fit a T-learner via econml
#'
#' Fits two separate outcome regressions and reports
#' \eqn{\hat\mu_1(X) - \hat\mu_0(X)} as the CATE.
#'
#' @param X,T,Y,covariate_names,bootstrap_ci,n_bootstrap,random_state,binary_outcome,...
#'   See \code{\link{cf_fit_slearner}}.
#' @param model_outcome Base-learner name for both arms. Default \code{"ridge"}.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cfomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @keywords internal
cf_fit_tlearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome   = "ridge",
                            bootstrap_ci    = FALSE,
                            n_bootstrap     = 100L,
                            random_state    = NULL,
                            binary_outcome  = FALSE,
                            ...) {
  .validate_learner_spec(model_outcome)
  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)
  seed_py    <- if (is.null(random_state)) NULL else as.integer(random_state)
  model_path <- .cf_new_model_path("tlearner", "joblib")
  pypath     <- system.file("python", package = "cfomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, binary_outcome,
             seed_py, bootstrap_ci, n_bootstrap) {
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
                                          random_state = seed_py)
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
        for (b in seq_len(B)) {
          if (!is.null(seed_py)) np$random$seed(as.integer(seed_py) + b)
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
    bootstrap_ci = isTRUE(bootstrap_ci), n_bootstrap = n_bootstrap
  )

  .metalearner_assemble(
    fit_out, model_path, backend = "econml.metalearners.TLearner",
    bootstrap_ci = bootstrap_ci, extra_meta = list(model_outcome = model_outcome)
  )
}

#' Predict from a T-learner
#' @inheritParams predict_cf_slearner
#' @return Depending on `type`: the full list of predictions when
#'   `type = "all"`, otherwise the requested component (a scalar ATE, a numeric
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

#' Fit an X-learner via econml
#'
#' Two-stage metalearner (Künzel et al. 2019). \code{mu0_hat} / \code{mu1_hat}
#' are the Stage-1 outcome fits; \code{tau_hat} is the Stage-3
#' propensity-weighted blend (so \code{mu1 - mu0 != tau_hat} by design).
#'
#' @param X,T,Y,covariate_names,bootstrap_ci,n_bootstrap,random_state,binary_outcome,...
#'   See \code{\link{cf_fit_slearner}}.
#' @param model_outcome Base-learner for Stage-1 outcome regressions. Default
#'   \code{"ridge"}.
#' @param model_final Base-learner for the Stage-2 CATE model.
#'   Default \code{NULL}, which resolves to \code{model_outcome} — i.e. the
#'   Stage-2 CATE learner \strong{tracks the base learner}. This mirrors
#'   econml's \code{XLearner}, whose \code{cate_models} default to the outcome
#'   models when unset. \strong{Asymmetry with \code{\link{cf_fit_drlearner}}}:
#'   xlearner's final stage follows \code{model_outcome}, whereas drlearner's
#'   \code{model_final} defaults to a fixed linear (\code{"ridge"}). So passing
#'   only \code{model_outcome} changes xlearner's final stage but not
#'   drlearner's.
#'   \code{NULL} reuses \code{model_outcome}.
#' @param model_propensity Base-learner for the Stage-2 propensity model.
#'   Default \code{"logistic"}.
#' @param trim Length-2 propensity clipping bounds.
#' @param propensity_score Optional known PS (scalar or length-n).
#' @param propensity_X Optional covariate subset for the propensity model.
#' @return A list with \code{model} (\code{list(model_path, backend)}) and
#'   \code{res} conforming to the \code{cfomics_result} contract.
#' @seealso \code{\link{metalearners}}.
#' @keywords internal
cf_fit_xlearner <- function(X, T, Y,
                            covariate_names = NULL,
                            model_outcome    = "ridge",
                            model_final      = NULL,
                            model_propensity = "logistic",
                            bootstrap_ci     = FALSE,
                            n_bootstrap      = 100L,
                            trim             = c(0.025, 0.975),
                            propensity_score = NULL,
                            propensity_X     = NULL,
                            random_state     = NULL,
                            binary_outcome   = FALSE,
                            ...) {
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
  model_path <- .cf_new_model_path("xlearner", "joblib")
  pypath     <- system.file("python", package = "cfomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, Y, T, model_outcome, model_final,
             model_propensity, trim, propensity_score, ps_indices_0b,
             binary_outcome, seed_py, bootstrap_ci, n_bootstrap) {
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
          pm <- factory$make_classifier(model_propensity, random_state = seed_py)
          if (is.null(ps_indices_0b)) {
            clipped$ClippedPropensityModel(pm, lower = as.numeric(trim[1]),
                                           upper = as.numeric(trim[2]))
          } else {
            clipped$ColumnSubsetPropensityModel(
              base_estimator = pm,
              column_indices = as.integer(ps_indices_0b),
              lower = as.numeric(trim[1]), upper = as.numeric(trim[2]))
          }
        } else {
          clipped$FixedPropensityModel(as.numeric(propensity_score))
        }
      }
      build_fit <- function(Xm, Yv, Tv) {
        outcome <- factory$make_outcome_model(model_outcome,
                                              binary = isTRUE(binary_outcome),
                                              random_state = seed_py)
        cate <- factory$make_regressor(model_final, random_state = seed_py)
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
        for (b in seq_len(B)) {
          if (!is.null(seed_py)) np$random$seed(as.integer(seed_py) + b)
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
    model_final = model_final, model_propensity = model_propensity,
    trim = trim, propensity_score = propensity_score,
    ps_indices_0b = ps_indices_0b, binary_outcome = binary_outcome,
    seed_py = seed_py, bootstrap_ci = isTRUE(bootstrap_ci),
    n_bootstrap = n_bootstrap
  )

  # Surface the known PS when supplied, otherwise NA. econml does not expose a
  # single cross-fit PS vector.
  n_obs <- length(Y_vec)
  propensity_hat <- if (!is.null(propensity_score)) {
    if (length(propensity_score) == 1L) rep(as.numeric(propensity_score), n_obs)
    else as.numeric(propensity_score)
  } else {
    rep(NA_real_, n_obs)
  }

  .metalearner_assemble(
    fit_out, model_path, backend = "econml.metalearners.XLearner",
    bootstrap_ci = bootstrap_ci,
    extra_meta = list(
      model_outcome = model_outcome, model_propensity = model_propensity,
      note = "mu1_hat - mu0_hat != tau_hat by design (propensity-weighted Stage 3 blend, not direct Stage 1 difference)"
    ),
    propensity = propensity_hat
  )
}

#' Predict from an X-learner
#' @inheritParams predict_cf_slearner
#' @return Depending on `type`: the full list of predictions when
#'   `type = "all"`, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_xlearner <- function(object, newdata = NULL, type = "all", ...) {
  .predict_metalearner_generic(object, newdata, type,
                               po_fn_name = "xlearner_potential_outcomes",
                               method_name = "xlearner")
}

# -------------------------------------------------------------------------
# Shared predict dispatcher (T / X, econml effect + potential outcomes)
# -------------------------------------------------------------------------

.predict_metalearner_generic <- function(object, newdata, type, po_fn_name,
                                          method_name = NULL) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }
  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.",
                 class = "cfomics_unsupported_predict_type")
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
    function(pypath, model_path, X, po_fn_name) {
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
      list(tau = tau, mu0 = as.numeric(po[[1]]), mu1 = as.numeric(po[[2]]))
    },
    pypath = system.file("python", package = "cfomics"),
    model_path = object$model$model$model_path, X = X_new,
    po_fn_name = po_fn_name
  )

  .assemble_predict_output(
    type = type, ate = mean(pred$tau), ate_ci = NA_real_,
    tau = pred$tau, tau_ci = NA_real_, mu0_hat = pred$mu0, mu1_hat = pred$mu1
  )
}
