#' @title GANITE Method Implementation
#' @description Internal functions for the GANITE causal inference method.
#'
#' GANITE (Generative Adversarial Nets for Inference of Individualized
#' Treatment Effects) uses GANs to estimate individual treatment effects.
#' The implementation follows the loss formulation of Yoon et al. (2018,
#' Section 4), with the sigmoid output replaced by a linear head and BCE on
#' Y replaced by MSE to support continuous outcomes. The paper's four loss
#' components are all retained:
#' \itemize{
#'   \item Counterfactual block: \eqn{V_{CF}} (BCE on t) + \eqn{\alpha \cdot L_S^G}
#'     (factual MSE).
#'   \item ITE block (k = 2 binary T): \eqn{V_{ITE}} (BCE distinguishing
#'     \eqn{\bar y} from \eqn{\hat y}) + \eqn{\beta \cdot L_S^I} with
#'     \eqn{L_S^I = ((\bar y_1 - \bar y_0) - (\hat y_1 - \hat y_0))^2}.
#' }
#'
#' Reference: Jinsung Yoon, James Jordon, Mihaela van der Schaar,
#' "GANITE: Estimation of Individualized Treatment Effects using Generative
#' Adversarial Nets", ICLR 2018. \url{https://openreview.net/forum?id=ByKWUeWA-}

#' Get GANITE Python module
#'
#' Uses PyPI ganite 0.1.2 (Yoon 2018 official PyTorch implementation) via
#' `ganite_wrapper.py`, which imports `from ganite import Ganite`. The former
#' self-hosted TensorFlow backend was removed so it cannot shadow the PyPI
#' package on `sys.path`.
#'
#' Backend: PyPI ganite (PyTorch), with native Blackwell sm_120 support. A
#' 3-source RNG (numpy + torch + random) save/restore around the fit prevents
#' global state leak; a module-global DEVICE override respects the caller's
#' device choice.
#'
#' @return The Python module object for the ganite backend.
#' @keywords internal
.get_ganite_module <- function() {
  reticulate::import_from_path(
    module = "ganite_wrapper",
    path = system.file("python", package = "cfomics")
  )
}

#' Fit GANITE model
#'
#' Internal fit used by \code{cf_fit(method = "ganite")}. \code{Y} is z-score
#' normalized on the R side (\code{X} is standardised centrally in
#' preprocessing) and the data is handed to the Python GANITE module. The
#' trained Inference network is serialized to disk (\code{torch.save}) under the
#' package model cache; the returned object holds only the model path, and
#' \code{predict(newdata = ...)} reloads it to score fresh samples.
#'
#' @param X Numeric matrix of covariates.
#' @param T Binary treatment vector (0/1).
#' @param Y Numeric outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param random_state Integer random seed for reproducibility. The backend
#'   (PyTorch) seeds the numpy, torch and Python \code{random} RNGs for the fit
#'   inside the ephemeral basilisk worker process, so seeding does not leak into
#'   the user's session.
#' @param h_dim Integer hidden-layer width (default 100).
#' @param alpha Numeric weight of the supervised factual loss
#'   \eqn{L_S^G} in the generator objective (default \code{2.0}, the
#'   IHDP-tuned value from the GANITE paper, Appendix Table 7).
#' @param beta Numeric weight of the supervised PEHE-type loss
#'   \eqn{L_S^I} in the inference objective, which balances GAN
#'   (\eqn{V_{ITE}}) and supervised terms (default \code{5.0}, IHDP-tuned).
#' @param iterations Integer number of training iterations per phase
#'   (default \code{5000}, the PyPI ganite upstream default; the GANITE paper
#'   uses 10000).
#' @param batch_size Integer mini-batch size (default 256).
#' @param verbose Logical, print training progress.
#' @param ... Additional arguments forwarded to Python (currently ignored).
#' @return A \code{cf_model_ganite} object.
#' @keywords internal
cf_fit_ganite <- function(X, T, Y,
                          covariate_names = NULL,
                          random_state = 0L,
                          h_dim = 100L,
                          # iterations = 5000L is the PyPI ganite 0.1.2 upstream
                          # default and within the range recommended in Yoon 2018
                          # Section 4. alpha = 2.0 / beta = 5.0 are the original
                          # values from Yoon et al. ICLR 2018 (used for the IHDP /
                          # Jobs benchmarks in the paper's Table 4), kept as a
                          # policy choice; the PyPI library defaults (alpha=0.1,
                          # beta=0) are the van der Schaar lab's conservative
                          # defaults, which depart from the paper recommendation.
                          alpha = 2.0,
                          beta = 5.0,
                          iterations = 5000L,
                          batch_size = 256L,
                          verbose = FALSE,
                          device = "cpu",
                          binary_outcome = FALSE,
                          ...) {
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "double"
  Y_vec <- as.numeric(Y)
  T_vec <- as.numeric(T)

  # X is standardised centrally in `.cf_preprocess_input` (binary/one-hot
  # skipped), so feed it through unchanged. Only Y is standardised here (the
  # GANITE backend expects a normalised outcome) and inverse-transformed after
  # fit.
  X_norm <- X_mat

  y_mean <- mean(Y_vec)
  y_sd <- stats::sd(Y_vec)
  if (!is.finite(y_sd) || y_sd < 1e-8) y_sd <- 1
  Y_norm <- (Y_vec - y_mean) / y_sd
  # binary outcome -> feed the raw 0/1 outcome (sigmoid head); skip the
  # Y z-score so inverse-transform is identity.
  if (isTRUE(binary_outcome)) {
    y_mean <- 0; y_sd <- 1
    Y_norm <- Y_vec
  }

  model_path <- .cf_new_model_path("ganite", "pt")
  pypath     <- system.file("python", package = "cfomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X_norm, T, Y_norm, covariate_names,
             random_state, h_dim, alpha, beta, iterations, batch_size,
             verbose, device) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      torch <- reticulate::import("torch", convert = FALSE)
      mod   <- reticulate::import_from_path("ganite_wrapper",
                                            path = pypath, convert = FALSE)
      cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
      if (identical(device, "auto")) device <- if (cuda_ok) "cuda" else "cpu"
      else if (identical(device, "gpu")) device <- if (cuda_ok) "cuda" else "cpu"

      ret <- mod$run_ganite_pypi(
        X = X_norm, T = as.numeric(T), Y = as.numeric(Y_norm),
        variable_names = covariate_names, random_state = as.integer(random_state),
        h_dim = as.integer(h_dim), alpha = as.numeric(alpha),
        beta = as.numeric(beta), iterations = as.integer(iterations),
        batch_size = as.integer(batch_size), verbose = verbose, device = device)
      mod$save_ganite(ret$model, model_path)
      res_norm <- ret$res
      list(
        mu0 = as.numeric(reticulate::py_to_r(res_norm$mu0_hat)),
        mu1 = as.numeric(reticulate::py_to_r(res_norm$mu1_hat)),
        device_effective = tryCatch(
          as.character(reticulate::py_to_r(res_norm$device_effective)),
          error = function(e) NA_character_))
    },
    pypath = pypath, model_path = model_path,
    X_norm = X_norm, T = T_vec, Y_norm = Y_norm,
    covariate_names = covariate_names, random_state = random_state,
    h_dim = h_dim, alpha = alpha, beta = beta, iterations = iterations,
    batch_size = batch_size, verbose = verbose, device = device
  )

  # Denormalize back to the original Y scale.
  mu0_hat <- fit_out$mu0 * y_sd + y_mean
  mu1_hat <- fit_out$mu1 * y_sd + y_mean
  ite <- mu1_hat - mu0_hat
  ate <- mean(ite)

  res_final <- list(
    tau_hat = ite,
    ate = ate,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    ate_ci  = NULL,
    tau_ci  = NULL,
    # X z-score is central (meta$x_zscore from prep); ganite only standardises
    # Y here and inverse-maps it. cf_predict re-applies the central X z-score to
    # newdata.
    x_zscore = NULL,
    y_zscore = list(Y_mean = y_mean, Y_sd = y_sd),
    device_effective = fit_out$device_effective,
    metadata = list(
      ci_method = "unavailable",
      hyperparameters = list(
        h_dim = as.integer(h_dim),
        alpha = as.numeric(alpha),
        beta = as.numeric(beta),
        iterations = as.integer(iterations),
        batch_size = as.integer(batch_size),
        random_state = as.integer(random_state)
      )
    )
  )

  structure(
    list(
      model = list(model_path = model_path, backend = "ganite.Ganite"),
      res   = res_final
    ),
    class = "cf_model_ganite"
  )
}

#' Predict from GANITE model
#'
#' Returns in-sample estimates when \code{newdata} is \code{NULL}. Otherwise
#' reloads the serialized Inference network from \code{object$model$model$model_path}
#' and re-scores \code{newdata} using the same z-score normalization captured at
#' fit time. The model is persisted to the on-disk cache, so prediction still
#' works after a session restart as long as the cache file exists.
#'
#' @param object cf_model object with \code{method = "ganite"}.
#' @param newdata Optional matrix or data.frame of new covariates with the
#'   same column order/dimensionality as the training \code{X}. If
#'   \code{NULL}, returns in-sample estimates.
#' @param type One of \code{"ate"}, \code{"tau"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param ... Unused.
#' @return The requested prediction; see \code{type}.
#' @keywords internal
predict_cf_ganite <- function(object,
                              newdata = NULL,
                              type = c("all", "ate", "ate_ci",
                                       "tau", "tau_ci",
                                       "mu0", "mu1",
                                       "summary"),
                              ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    return(.cf_insample_predict(object$res, type))
  }

  if (type %in% c("ate_ci", "tau_ci")) {
    .abort_ite_ci_unsupported("ganite")
  }

  ystats <- object$meta$y_zscore
  model_path <- object$model$model$model_path
  if (is.null(ystats) || is.null(model_path)) {
    rlang::abort(paste0("GANITE newdata prediction requires a fresh fit (meta Y z-score stats ",
         "or serialised model missing). Refit with cf_fit(method = 'ganite')."), class = "cfomics_ganite_newdata_unsupported")
  }

  # X z-score is central. cf_predict.R `.apply_newdata_preprocessing` already
  # re-applies meta$x_zscore to user newdata (and OOF held-out rows are already
  # standardised), so predict_cf_ganite receives scaled X and only
  # inverse-transforms Y.
  X_new_norm <- as.matrix(newdata)
  storage.mode(X_new_norm) <- "double"
  expected_p <- length(object$meta$covariate_names)
  if (expected_p > 0L && ncol(X_new_norm) != expected_p) {
    rlang::abort(sprintf("newdata has %d columns but fit was trained with %d.",
                 ncol(X_new_norm), expected_p), class = "cfomics_newdata_dim_mismatch")
  }

  pred <- .cf_py_run(
    function(pypath, model_path, X) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      mod <- reticulate::import_from_path("ganite_wrapper",
                                          path = pypath, convert = FALSE)
      m <- mod$load_ganite(model_path)
      p <- reticulate::py_to_r(mod$predict_ganite_pypi(m, X))
      list(mu0 = as.numeric(p$mu0_hat), mu1 = as.numeric(p$mu1_hat))
    },
    pypath = system.file("python", package = "cfomics"),
    model_path = model_path, X = X_new_norm
  )
  y0_new <- pred$mu0 * ystats$Y_sd + ystats$Y_mean
  y1_new <- pred$mu1 * ystats$Y_sd + ystats$Y_mean
  ite_new <- y1_new - y0_new
  ate_new <- mean(ite_new)

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    ate_ci  = NA_real_,
    tau     = ite_new,
    tau_ci  = NA_real_,
    mu0_hat = y0_new,
    mu1_hat = y1_new
  )
}
