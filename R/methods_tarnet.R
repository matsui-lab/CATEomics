#' Fit a TARNet CATE model via catenets (PyTorch)
#'
#' Thin R wrapper around \pkg{catenets} (van der Schaar lab, v0.2.x)
#' \code{TARNet} (Shalit-Johansson-Sontag, ICML 2017). With
#' \code{penalty_disc > 0} the model becomes CFRNet (counterfactual
#' regression with an integral probability metric penalty on the shared
#' representation). CATEomics ships the wrapper at
#' \code{inst/python/catenets_wrapper.py}; the Python module is loaded
#' lazily via reticulate.
#'
#' @section Potential outcomes -- Class A, direct potential-outcome output:
#' \code{mu0_hat} / \code{mu1_hat} are the two output heads' predictions on
#' the shared representation. By construction
#' \code{mu1_hat - mu0_hat == ite} (plug-in CATE).
#'
#' @section RNG bridge (catenets reproducibility):
#' catenets v0.2.x exposes a \code{seed} constructor argument but only
#' uses it for train/validation splitting. Per-batch shuffling and weight
#' initialisation use the \emph{global} numpy / torch / Python RNGs. The wrapper
#' at \code{inst/python/catenets_wrapper.py} saves and restores those
#' globals around the fit, calling \code{torch.manual_seed},
#' \code{np.random.seed}, and \code{random.seed} for bit-exact
#' reproducibility.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector (continuous; \code{binary_y = FALSE} default).
#' @param covariate_names Character vector of covariate names.
#' @param random_state Integer seed (forwarded to numpy / torch / Python
#'   RNGs inside the wrapper). When \code{NULL}, no seeding is performed.
#'   The seed does not reach the thread count: a fit repeats bit-exactly
#'   while the thread count is the same, and its low-order digits change
#'   when that count changes. The thread count is set outside R, for
#'   instance through \code{OMP_NUM_THREADS}.
#' @param penalty_disc Numeric. CFRNet IPM penalty on the shared
#'   representation. Default \code{0} (plain TARNet); set to a positive
#'   value (e.g. \code{0.1}) to enable CFRNet-style discrepancy
#'   regularisation.
#' @param n_iter Integer number of training iterations. Default
#'   \code{10000} (catenets upstream default, with early stopping at
#'   \code{n_iter_min = 200} + \code{patience = 10}).
#' @param val_split_prop Numeric fraction of the supplied rows catenets holds
#'   out to decide when to stop. Default \code{0.3} (catenets upstream
#'   default), so the weights are fitted on the other 70 percent; under
#'   \code{partition = "oof"} that is 70 percent of each fold's training part.
#'   The held-out rows score the factual outcome only. \code{0} makes the
#'   training rows the held-out set as well, so every row trains; catenets
#'   keeps stopping early, now on the training loss, since that is governed by
#'   a separate argument CATEomics does not set.
#' @param batch_size Integer mini-batch size. Default \code{100}
#'   (catenets upstream default).
#' @param learning_rate Numeric learning rate. Default \code{1e-4}
#'   (catenets upstream default).
#' @param binary_y Logical, whether the outcome head is binary. Default
#'   \code{FALSE} (continuous Y, matches catenets upstream default). Binary
#'   outcomes are normally set automatically via \code{cf_fit(outcome_type =)}.
#' @param params Named list passed straight to the catenets \code{TARNet}
#'   constructor, so its names are catenets' own. CATEomics does not check them
#'   here: an unknown name raises a \code{TypeError} from catenets.
#'   A value the Python side requires to be whole is written as an R
#'   integer (\code{10L}, not \code{10}), which R otherwise makes a double.
#' @param device Character, \code{"auto"} / \code{"cpu"} / \code{"gpu"}.
#' @param binary_outcome Logical; when \code{TRUE} the outcome is treated as
#'   a class label and the effect is a risk difference.
#' @return A list with \code{model} (the path the fitted estimator was
#'   serialised to, and the backend name) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @references
#' Shalit, Johansson & Sontag (2017) Estimating individual treatment effect:
#' generalization bounds and algorithms, \emph{PMLR} 70:3076-3085.
#'
#' Curth & van der Schaar (2021) Nonparametric Estimation of Heterogeneous
#' Treatment Effects: From Theory to Learning Algorithms, \emph{AISTATS}
#' (catenets, implementation).
#'
#' @keywords internal
cf_fit_tarnet <- function(X, T, Y,
                          covariate_names = NULL,
                          random_state    = NULL,
                          penalty_disc    = 0,
                          n_iter          = 10000L,
                          batch_size      = 100L,
                          learning_rate   = 1e-4,
                          val_split_prop  = 0.3,
                          binary_y        = FALSE,
                          device          = "cpu",
                          binary_outcome  = FALSE,
                          params          = list()) {
  binary_y <- isTRUE(binary_y) || isTRUE(binary_outcome)

  X_mat <- .build_meta_X(X, covariate_names)
  Y_vec <- as.numeric(Y)
  T_vec <- as.integer(T)

  # X is standardised centrally in preprocessing; catenets standardises neither
  # X nor Y internally, so CATEomics additionally z-scores Y R-side, fits catenets
  # in the scaled space, and inverse-transforms the returned mu / ite back.
  z <- .catenets_scale_y(X_mat, Y_vec)
  if (isTRUE(binary_y)) {
    # keep the raw 0/1 outcome (sigmoid / BCE head); skip Y z-score.
    z$Y_scaled <- Y_vec
    z$y_zscore <- list(Y_mean = 0, Y_sd = 1)
  }

  model_path <- .cf_new_model_path("tarnet", "pt")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, T, Y, random_state, binary_y,
             penalty_disc, n_iter, batch_size, learning_rate, val_split_prop,
             device, extra) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      cw <- reticulate::import_from_path("catenets_wrapper",
                                         path = pypath, convert = FALSE)
      out <- do.call(cw$fit_tarnet, c(list(
        X = X, T = as.integer(T), Y = as.numeric(Y),
        random_state  = random_state,
        binary_y      = isTRUE(binary_y),
        penalty_disc  = as.numeric(penalty_disc),
        n_iter        = as.integer(n_iter),
        batch_size    = as.integer(batch_size),
        learning_rate = as.numeric(learning_rate),
        val_split_prop = as.numeric(val_split_prop),
        device        = device), extra))
      cw$save_model(out$model, model_path)
      list(
        ite = as.numeric(reticulate::py_to_r(out$ite)),
        mu0 = as.numeric(reticulate::py_to_r(out$mu0_hat)),
        mu1 = as.numeric(reticulate::py_to_r(out$mu1_hat)),
        device_effective = tryCatch(
          as.character(reticulate::py_to_r(out$device_effective)),
          error = function(e) NA_character_))
    },
    pypath = pypath, model_path = model_path,
    X = z$X_scaled, T = T_vec, Y = z$Y_scaled,
    random_state = if (is.null(random_state)) NULL else as.integer(random_state),
    binary_y = isTRUE(binary_y), penalty_disc = as.numeric(penalty_disc),
    n_iter = as.integer(n_iter), batch_size = as.integer(batch_size),
    learning_rate = as.numeric(learning_rate),
    val_split_prop = as.numeric(val_split_prop), device = device,
    extra = params
  )

  inv <- .catenets_inverse_predictions(fit_out$mu0, fit_out$mu1, fit_out$ite,
                                       z$y_zscore)

  list(
    model = list(model_path = model_path,
                 backend = "catenets.models.torch.TARNet"),
    res = list(
      tau_hat = inv$ite,
      ate     = inv$ate,
      mu0_hat = inv$mu0_hat,
      mu1_hat = inv$mu1_hat,
      ate_ci  = NULL,
      tau_ci  = NULL,
      # cf_fit moves x_zscore / y_zscore into meta; cf_predict re-applies the
      # train statistics to newdata (data-leakage defence).
      x_zscore = z$x_zscore,
      y_zscore = z$y_zscore,
      device_effective = fit_out$device_effective,
      metadata = list(
        ci_method     = "unavailable",
        backend       = "catenets.models.torch.TARNet",
        binary_y      = isTRUE(binary_y),
        penalty_disc  = as.numeric(penalty_disc),
        n_iter        = as.integer(n_iter),
        batch_size    = as.integer(batch_size),
        learning_rate = as.numeric(learning_rate),
        val_split_prop = as.numeric(val_split_prop)
      )
    )
  )
}

#' Predict from a TARNet model
#'
#' @param object \code{cf_model} with \code{method = "tarnet"}.
#' @param newdata Optional new covariates. When \code{NULL} the cached
#'   training predictions are returned.
#' @param type One of \code{"all"} (the default), \code{"ate"},
#'   \code{"ate_ci"}, \code{"tau"}, \code{"tau_ci"}, \code{"mu0"},
#'   \code{"mu1"}, \code{"ps"}, \code{"summary"}.
#' @param ... Ignored.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_tarnet <- function(object, newdata = NULL, type = "all", ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }
  if (type %in% c("ate_ci", "tau_ci")) {
    .abort_ite_ci_unsupported("tarnet")
  }
  X_new <- as.matrix(newdata)
  covariate_names <- object$meta$covariate_names
  if (!is.null(covariate_names) && all(covariate_names %in% colnames(X_new))) {
    X_new <- X_new[, covariate_names, drop = FALSE]
  }
  # X z-score is applied upstream by cf_predict.R via meta$x_zscore
  # (double-z-score guard); predict_cf_tarnet receives scaled X and only
  # inverse-transforms Y here.
  pred <- .cf_py_run(
    function(pypath, model_path, X, device) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      cw <- reticulate::import_from_path("catenets_wrapper",
                                         path = pypath, convert = FALSE)
      m <- cw$load_model(model_path, device = device)
      p <- cw$predict_with(m, X, device = device)
      list(ite = as.numeric(reticulate::py_to_r(p$ite)),
           mu0 = as.numeric(reticulate::py_to_r(p$mu0_hat)),
           mu1 = as.numeric(reticulate::py_to_r(p$mu1_hat)))
    },
    pypath = system.file("python", package = "CATEomics"),
    model_path = object$model$model$model_path, X = X_new,
    # Predict where the fit ran. Without this the reload takes whatever GPU
    # is visible, so the same saved model answers differently per host.
    device = .cf_predict_device(object)
  )
  ite_scaled <- pred$ite
  mu0_scaled <- pred$mu0
  mu1_scaled <- pred$mu1
  ystats <- object$meta$y_zscore
  if (!is.null(ystats)) {
    inv <- .catenets_inverse_predictions(mu0_scaled, mu1_scaled, ite_scaled,
                                         ystats)
    y0_new  <- inv$mu0_hat
    y1_new  <- inv$mu1_hat
    ite_new <- inv$ite
  } else {
    y0_new  <- mu0_scaled
    y1_new  <- mu1_scaled
    ite_new <- ite_scaled
  }
  .assemble_predict_output(
    type    = type,
    ate     = mean(ite_new),
    # NULL, not NA: an interval that was not computed is absent, and the
    # stored path and every other method's newdata path say it the same way.
    ate_ci  = NULL,
    tau     = ite_new,
    tau_ci  = NULL,
    mu0_hat = y0_new,
    mu1_hat = y1_new
  )
}
