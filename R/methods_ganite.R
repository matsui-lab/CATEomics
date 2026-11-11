#' @title GANITE Method Implementation
#' @description Internal functions for the GANITE causal inference method.
#'
#' GANITE (Generative Adversarial Nets for Inference of Individualized
#' Treatment Effects) estimates individual treatment effects with two
#' adversarial blocks, implemented for CATEomics in
#' \code{inst/python/ganite_impl.py} on PyTorch:
#' \itemize{
#'   \item Counterfactual block: a generator \eqn{G} completes the unobserved
#'     outcome and a discriminator \eqn{D_G} tells which of the pair was
#'     observed. \eqn{G} minimises \eqn{V_{CF} + \alpha \cdot L_S^G}, the
#'     supervised term being the squared error on the observed outcome, or a
#'     cross-entropy when the outcome is binary.
#'   \item ITE block (k = 2, binary T): an inference network \eqn{I} is fitted
#'     to the completed data and a discriminator \eqn{D_I} tells \eqn{I}'s pair
#'     from \eqn{G}'s. \eqn{I} minimises \eqn{V_{ITE} + \beta \cdot L_S^I}
#'     with \eqn{L_S^I = ((\bar y_1 - \bar y_0) - (\hat y_1 - \hat y_0))^2}.
#' }
#' Each block runs \code{iterations} steps, and in each step the discriminator
#' is updated \code{d_steps} times before its counterpart once. The losses of
#' both blocks are recorded in \code{meta$estimator$training_trace}, one row
#' per 100 steps and one for the step each block ended on; the columns are
#' \code{block}, \code{iteration}, and \code{d_g} / \code{g} / \code{g_factual} for the
#' counterfactual block and \code{d_i} / \code{i} / \code{l_s_i} for the ITE
#' block, with the other block's columns \code{NA}. \code{d_i} is the sum of
#' two cross-entropies and \code{d_g} a single one, so the two are not on one
#' scale. With a binary outcome the heads are read through a logistic function,
#' so \code{mu0_hat} / \code{mu1_hat} are probabilities in \code{[0, 1]} and
#' the effect is their difference.
#'
#' Reference: Jinsung Yoon, James Jordon, Mihaela van der Schaar,
#' "GANITE: Estimation of Individualized Treatment Effects using Generative
#' Adversarial Nets", ICLR 2018. \url{https://openreview.net/forum?id=ByKWUeWA-}

#' Get GANITE Python module
#'
#' \code{ganite_wrapper.py} is the R-facing entry over \code{ganite_impl.py}.
#'
#' Seeding is local: the model seeds torch and its own mini-batch generator
#' from \code{random_state}, and \code{.cf_py_run()} snapshots and restores the caller's
#' global numpy / random / torch state. The device is passed on every call.
#'
#' @return The Python module object for the ganite backend.
#' @keywords internal
.get_ganite_module <- function() {
  reticulate::import_from_path(
    module = "ganite_wrapper",
    path = system.file("python", package = "CATEomics")
  )
}

#' Fit GANITE model
#'
#' Internal fit used by \code{cf_fit(method = "ganite")}. The model is an
#' independent implementation of the published method, written for this
#' package in \code{inst/python/ganite_impl.py}; no code of the authors' own
#' release is used. \code{Y} is z-score
#' normalized on the R side (\code{X} is standardised centrally in
#' preprocessing) and the data is handed to that module. The
#' trained Inference network is serialized to disk (\code{torch.save}) under the
#' package model cache; the returned object holds only the model path, and
#' \code{predict(newdata = ...)} reloads it to score fresh samples.
#'
#' @param X Numeric matrix of covariates.
#' @param T Binary treatment vector (0/1).
#' @param Y Numeric outcome vector.
#' @param random_state Integer random seed for reproducibility. It seeds the
#'   torch generator the networks are initialised from and the generator the
#'   mini-batches are drawn with; \code{.cf_py_run()} snapshots and restores
#'   the caller's global numpy / \code{random} / torch state around the call,
#'   so the seeding does not leak into the user's session.
#' @param h_dim Integer hidden-layer width (default 100).
#' @param num_layers Integer number of hidden layers (default \code{2}),
#'   in each of the two trunks and in each discriminator. Every outcome head
#'   adds one more hidden layer of the same width.
#' @param alpha Numeric weight of the supervised factual loss
#'   \eqn{L_S^G} in the generator objective (default \code{2.0}).
#' @param beta Numeric weight of the supervised PEHE-type loss
#'   \eqn{L_S^I} in the inference objective, which balances GAN
#'   (\eqn{V_{ITE}}) and supervised terms (default \code{5.0}).
#' @param iterations Integer number of steps per block (default \code{1000}).
#'   The counterfactual block and the ITE block each run this many.
#' @param batch_size Integer mini-batch size (default 256), reduced to the
#'   number of rows when it exceeds them. \code{h_dim}, \code{num_layers},
#'   \code{iterations}, \code{batch_size} and \code{d_steps} must be at
#'   least 1.
#' @param learning_rate Numeric Adam learning rate for all four networks
#'   (default \code{1e-3}).
#' @param weight_decay Numeric Adam weight decay for all four networks
#'   (default \code{1e-4}).
#' @param d_steps Integer discriminator steps per counterpart step
#'   (default \code{2}), applied in both blocks.
#' @param params Named list of additional hyperparameters for this method's
#'   backend, resolved by [cf_fit()] from its \code{params} argument.
#' @param device Character, \code{"auto"} / \code{"cpu"} / \code{"gpu"}.
#' @param binary_outcome Logical; when \code{TRUE} the outcome is treated as
#'   a class label and the effect is a risk difference.
#' @return A \code{cf_model_ganite} object.
#' @references
#' Yoon, Jordon & van der Schaar (2018) GANITE: Estimation of Individualized
#' Treatment Effects using Generative Adversarial Nets, \emph{ICLR}.
#'
#' Paszke et al. (2019) PyTorch: An Imperative Style, High-Performance Deep
#' Learning Library, \emph{NeurIPS} 32:8024-8035 (implementation).
#'
#' @keywords internal
cf_fit_ganite <- function(X, T, Y,
                          random_state = 0L,
                          h_dim = 100L,
                          num_layers = 2L,
                          alpha = 2.0,
                          beta = 5.0,
                          iterations = 1000L,
                          batch_size = 256L,
                          learning_rate = 1e-3,
                          weight_decay = 1e-4,
                          d_steps = 2L,
                          device = "cpu",
                          binary_outcome = FALSE,
                          params = list()) {
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "double"
  Y_vec <- as.numeric(Y)
  T_vec <- as.numeric(T)

  # A zero would not fail, it would fit nothing: the backend floors the layer
  # count, the step count and the discriminator ratio to 1, and a zero
  # mini-batch leaves every network at its initialisation, so the caller sees
  # estimates with no variation and no message. Reject the value instead.
  for (nm in c("h_dim", "num_layers", "iterations", "batch_size", "d_steps")) {
    v <- get(nm)
    if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 1) {
      rlang::abort(
        c(sprintf("`params$ganite$%s` must be a single integer of at least 1.", nm),
          x = sprintf("Received: %s.", paste(format(v), collapse = ", "))),
        class = "cateomics_invalid_ganite_param")
    }
  }

  # X is standardised centrally in `.cf_preprocess_input` (every column,
  # including binary indicators), so feed it through unchanged. Only Y is
  # standardised here (the GANITE backend expects a normalised outcome) and
  # inverse-transformed after fit.
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
  # The outcome heads emit a probability when the outcome is binary, so the
  # supervised factual term is a cross-entropy and the effect is a difference
  # of probabilities.
  binary_y <- isTRUE(binary_outcome)

  model_path <- .cf_new_model_path("ganite", "pt")
  pypath     <- system.file("python", package = "CATEomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X_norm, T, Y_norm, random_state, binary_y,
             h_dim, num_layers, alpha, beta, iterations, batch_size,
             learning_rate, weight_decay, d_steps, device, extra) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      torch <- reticulate::import("torch", convert = FALSE)
      mod   <- reticulate::import_from_path("ganite_wrapper",
                                            path = pypath, convert = FALSE)
      cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
      if (identical(device, "auto")) device <- if (cuda_ok) "cuda" else "cpu"
      else if (identical(device, "gpu")) device <- if (cuda_ok) "cuda" else "cpu"

      ret <- do.call(mod$run_ganite, c(list(
        X = X_norm, T = as.numeric(T), Y = as.numeric(Y_norm),
        random_state = as.integer(random_state), binary_y = isTRUE(binary_y),
        h_dim = as.integer(h_dim), num_layers = as.integer(num_layers),
        alpha = as.numeric(alpha), beta = as.numeric(beta),
        iterations = as.integer(iterations),
        batch_size = as.integer(batch_size),
        learning_rate = as.numeric(learning_rate),
        weight_decay = as.numeric(weight_decay),
        d_steps = as.integer(d_steps), device = device), extra))
      mod$save_ganite(ret$model, model_path)
      list(
        mu0 = as.numeric(reticulate::py_to_r(ret$mu0_hat)),
        mu1 = as.numeric(reticulate::py_to_r(ret$mu1_hat)),
        device_effective = tryCatch(
          as.character(reticulate::py_to_r(ret$device_effective)),
          error = function(e) NA_character_),
        training_trace = tryCatch(
          reticulate::py_to_r(ret$training_trace),
          error = function(e) NULL))
    },
    pypath = pypath, model_path = model_path,
    X_norm = X_norm, T = T_vec, Y_norm = Y_norm,
    random_state = random_state, binary_y = isTRUE(binary_y),
    h_dim = h_dim, num_layers = num_layers, alpha = alpha, beta = beta,
    iterations = iterations, batch_size = batch_size,
    learning_rate = learning_rate, weight_decay = weight_decay,
    d_steps = d_steps, device = device,
    extra = params
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
        num_layers = as.integer(num_layers),
        alpha = as.numeric(alpha),
        beta = as.numeric(beta),
        iterations = as.integer(iterations),
        # The backend caps the mini-batch at the number of rows, so record
        # what it used rather than what was asked for.
        batch_size = min(as.integer(batch_size), nrow(X_mat)),
        learning_rate = as.numeric(learning_rate),
        weight_decay = as.numeric(weight_decay),
        d_steps = as.integer(d_steps),
        random_state = as.integer(random_state)
      )
    ),
    # Moved to meta$estimator$training_trace in cf_fit Step 7.
    training_trace = .cf_trace_to_frame(fit_out$training_trace)
  )

  structure(
    list(
      model = list(model_path = model_path, backend = "CATEomics.ganite_impl.GaniteModel"),
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
#' @param type One of \code{"all"}, \code{"ate"}, \code{"ate_ci"},
#'   \code{"tau"}, \code{"tau_ci"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"ps"}, \code{"summary"}. GANITE estimates no propensity score, so
#'   \code{"ps"} is \code{NA}; the interval types are unavailable with
#'   \code{newdata}.
#' @param ... Unused.
#' @return The requested prediction; see \code{type}.
#' @keywords internal
predict_cf_ganite <- function(object,
                              newdata = NULL,
                              type = c("all", "ate", "ate_ci",
                                       "tau", "tau_ci",
                                       "mu0", "mu1", "ps",
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
         "or serialised model missing). Refit with cf_fit(method = 'ganite')."), class = "cateomics_ganite_newdata_unsupported")
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
                 ncol(X_new_norm), expected_p), class = "cateomics_newdata_dim_mismatch")
  }

  pred <- .cf_py_run(
    function(pypath, model_path, X, device) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      mod <- reticulate::import_from_path("ganite_wrapper",
                                          path = pypath, convert = FALSE)
      m <- mod$load_ganite(model_path, device = device)
      p <- reticulate::py_to_r(mod$predict_ganite(m, X))
      list(mu0 = as.numeric(p$mu0_hat), mu1 = as.numeric(p$mu1_hat))
    },
    pypath = system.file("python", package = "CATEomics"),
    model_path = model_path, X = X_new_norm,
    # Predict where the fit ran. The reload otherwise takes whatever GPU is
    # visible, so a CPU fit answers from the GPU on a GPU host.
    device = .cf_predict_device(object)
  )
  y0_new <- pred$mu0 * ystats$Y_sd + ystats$Y_mean
  y1_new <- pred$mu1 * ystats$Y_sd + ystats$Y_mean
  ite_new <- y1_new - y0_new
  ate_new <- mean(ite_new)

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    # NULL, not NA: an interval that was not computed is absent, and the
    # stored path and every other method's newdata path say it the same way.
    ate_ci  = NULL,
    tau     = ite_new,
    tau_ci  = NULL,
    mu0_hat = y0_new,
    mu1_hat = y1_new
  )
}
