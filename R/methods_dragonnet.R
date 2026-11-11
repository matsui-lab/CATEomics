#' Fit a DragonNet CATE model via catenets (PyTorch)
#'
#' Thin R wrapper around \pkg{catenets} (van der Schaar lab, v0.2.x)
#' \code{DragonNet} (Shi-Blei-Veitch, NeurIPS 2019). DragonNet shares a
#' representation backbone with TARNet and adds an auxiliary propensity head,
#' whose binary cross-entropy enters the training loss weighted by
#' \code{prop_loss_multiplier}. catenets v0.2.x implements no targeted
#' regularisation.
#'
#' @section Potential outcomes -- Class A, direct potential-outcome output:
#' \code{mu0_hat} / \code{mu1_hat} are the two output heads' predictions on
#' the shared representation. By construction
#' \code{mu1_hat - mu0_hat == ite} (plug-in CATE).
#'
#' @section RNG bridge (catenets reproducibility):
#' See \code{\link{cf_fit_tarnet}}; identical RNG handling applies.
#'
#' @inheritParams cf_fit_tarnet
#' @param params Named list passed straight to the catenets
#'   \code{DragonNet} constructor, so its names are catenets' own. CATEomics
#'   does not check them here: an unknown name raises a \code{TypeError}
#'   from catenets.
#'   A value the Python side requires to be whole is written as an R
#'   integer (\code{10L}, not \code{10}), which R otherwise makes a double.
#' @param device Character, \code{"auto"} / \code{"cpu"} / \code{"gpu"}.
#' @param binary_outcome Logical; when \code{TRUE} the outcome is treated as
#'   a class label and the effect is a risk difference.
#' @return A list with \code{model} (the path the fitted estimator was
#'   serialised to, and the backend name) and
#'   \code{res} conforming to the \code{cateomics_result} contract.
#' @references
#' Shi, Blei & Veitch (2019) Adapting neural networks for the estimation of
#' treatment effects, \emph{NeurIPS} 32.
#'
#' Curth & van der Schaar (2021) Nonparametric Estimation of Heterogeneous
#' Treatment Effects: From Theory to Learning Algorithms, \emph{AISTATS}
#' (catenets, implementation).
#'
#' @keywords internal
cf_fit_dragonnet <- function(X, T, Y,
                             covariate_names = NULL,
                             random_state    = NULL,
                             # Default matches the catenets upstream
                             # BasicDragonNet.
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
  # X nor Y internally, so CATEomics additionally z-scores Y R-side (shared helper
  # with TARNet).
  z <- .catenets_scale_y(X_mat, Y_vec)
  if (isTRUE(binary_y)) {
    z$Y_scaled <- Y_vec
    z$y_zscore <- list(Y_mean = 0, Y_sd = 1)
  }

  model_path <- .cf_new_model_path("dragonnet", "pt")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, T, Y, random_state, binary_y,
             n_iter, batch_size, learning_rate, val_split_prop,
             device, extra) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      cw <- reticulate::import_from_path("catenets_wrapper",
                                         path = pypath, convert = FALSE)
      out <- do.call(cw$fit_dragonnet, c(list(
        X = X, T = as.integer(T), Y = as.numeric(Y),
        random_state  = random_state,
        binary_y      = isTRUE(binary_y),
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
    binary_y = isTRUE(binary_y), n_iter = as.integer(n_iter),
    batch_size = as.integer(batch_size),
    learning_rate = as.numeric(learning_rate),
    val_split_prop = as.numeric(val_split_prop), device = device,
    extra = params
  )

  inv <- .catenets_inverse_predictions(fit_out$mu0, fit_out$mu1, fit_out$ite,
                                       z$y_zscore)

  list(
    model = list(model_path = model_path,
                 backend = "catenets.models.torch.DragonNet"),
    res = list(
      tau_hat = inv$ite,
      ate     = inv$ate,
      mu0_hat = inv$mu0_hat,
      mu1_hat = inv$mu1_hat,
      ate_ci  = NULL,
      tau_ci  = NULL,
      x_zscore = z$x_zscore,
      y_zscore = z$y_zscore,
      device_effective = fit_out$device_effective,
      metadata = list(
        ci_method     = "unavailable",
        backend       = "catenets.models.torch.DragonNet",
        binary_y      = isTRUE(binary_y),
        n_iter        = as.integer(n_iter),
        batch_size    = as.integer(batch_size),
        learning_rate = as.numeric(learning_rate),
        val_split_prop = as.numeric(val_split_prop)
      )
    )
  )
}

#' Predict from a DragonNet model
#'
#' @inheritParams predict_cf_tarnet
#' @param object \code{cf_model} with \code{method = "dragonnet"}.
#' @return Depending on \code{type}: the full list of predictions when
#'   \code{type = "all"}, otherwise the requested component (a scalar ATE, a numeric
#'   vector of per-sample effects or potential outcomes, or a summary list).
#' @keywords internal
predict_cf_dragonnet <- function(object, newdata = NULL, type = "all", ...) {
  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }
  if (type %in% c("ate_ci", "tau_ci")) {
    .abort_ite_ci_unsupported("dragonnet")
  }
  X_new <- as.matrix(newdata)
  covariate_names <- object$meta$covariate_names
  if (!is.null(covariate_names) && all(covariate_names %in% colnames(X_new))) {
    X_new <- X_new[, covariate_names, drop = FALSE]
  }
  # X z-score applied upstream via meta$x_zscore (double guard).
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
