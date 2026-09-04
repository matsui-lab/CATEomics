#' @title CEVAE Method Implementation
#' @description Internal functions for the CEVAE (Causal Effect Variational
#'   Autoencoder) causal inference method. Delegates the generative / inference
#'   model to Pyro's \code{pyro.contrib.cevae} and uses the in-package
#'   \code{CEVAEExtended} subclass (\code{inst/python/cevae_extended.py}) to
#'   obtain per-sample \eqn{\hat y_0}, \eqn{\hat y_1} (which the upstream
#'   \code{ite()} method hides behind a single difference).
#'
#' Reference: Christos Louizos, Uri Shalit, Joris Mooij, David Sontag, Richard
#'   Zemel, Max Welling, "Causal Effect Inference with Deep Latent-Variable
#'   Models", NeurIPS 2017. \url{https://arxiv.org/abs/1705.08821}
#'
#' Normalization: X is z-scored centrally by cfomics preprocessing (binary /
#'   one-hot columns skipped); the Pyro CEVAE PreWhitener is disabled to avoid
#'   double standardization. Y is z-scored internally and predictions are
#'   returned on the original Y scale.
#'
#' Covariate distribution: the generative model treats \strong{every}
#'   covariate as continuous Gaussian, i.e. \eqn{p(x_j \mid z) =
#'   \mathcal{N}(\mu_j(z), \sigma_j(z))} for all \eqn{j} (Pyro's
#'   \code{Model.x_nn = DiagNormalNet}). Binary / one-hot covariates are
#'   therefore modelled as Gaussian rather than Bernoulli, a mild and common
#'   VAE misspecification. If exact binary-covariate likelihoods matter for
#'   your analysis, prefer a method without this assumption (e.g. \code{grf},
#'   the meta-learners, or \code{tarnet}).


#' Fit CEVAE (Causal Effect Variational Autoencoder) model via Pyro
#'
#' @param X Numeric matrix of covariates.
#' @param T Binary treatment vector (0/1).
#' @param Y Numeric outcome vector.
#' @param covariate_names Character vector of covariate names (for provenance).
#' @param random_state Integer random seed. When set, yields bit-identical
#'   results across runs with the same inputs (constructor-time weight
#'   initialization is seeded too, so re-fits match). Seeding runs inside the
#'   ephemeral basilisk worker process, so it does not leak into the user's
#'   session.
#' @param outcome_dist Outcome distribution: \code{"normal"} (default, for
#'   continuous Y) or \code{"bernoulli"}. The default is \code{"normal"} because
#'   the primary target is continuous omics outcomes. For a binary outcome use
#'   \code{cf_fit(outcome_type = "binary")} (or let it auto-detect a 0/1 Y),
#'   which drives the full binary path (Bernoulli likelihood on the raw 0/1
#'   scale with the Y z-score skipped); setting \code{outcome_dist} alone does
#'   not switch off the Y z-score.
#' @param latent_dim Integer latent dimension (default 20).
#' @param hidden_dim Integer hidden-layer dimension (default 200).
#' @param num_layers Integer number of hidden layers (default 3).
#' @param num_samples Integer Monte-Carlo samples of the latent variable used to
#'   estimate y0/y1 (default 100, the Pyro CEVAE default).
#' @param num_epochs Integer number of training epochs (default 100).
#' @param batch_size Integer mini-batch size (default 100).
#' @param learning_rate Numeric learning rate (default 1e-3).
#' @param ... Additional arguments forwarded to Pyro (currently ignored).
#' @return A \code{cf_model_cevae} list (internal use) with components
#'   \code{model} (the Python \code{CEVAEExtended} instance) and \code{res}
#'   (the result payload conforming to the \code{cfomics_result} contract).
#' @keywords internal
cf_fit_cevae <- function(X, T, Y,
                         covariate_names = NULL,
                         random_state = 0L,
                         outcome_dist = "normal",
                         binary_outcome = FALSE,

                         latent_dim = 20L,
                         hidden_dim = 200L,
                         num_layers = 3L,
                         num_samples = 100L,
                         num_epochs = 100L,
                         batch_size = 100L,
                         learning_rate = 1e-3,
                         learning_rate_decay = 0.1,
                         weight_decay = 1e-4,
                         device = "cpu",
                         ...) {
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "double"
  feature_dim <- ncol(X_mat)
  fit_seed <- if (is.null(random_state)) 0L else as.integer(random_state)

  # X is standardised centrally in `.cf_preprocess_input` (binary/one-hot
  # skipped); the pyro PreWhitener is disabled (see `_NoWhiten` in
  # cevae_extended.py), so the model receives the already-scaled X. cfomics
  # applies a Y z-score R-side; the scaled Y is passed to the model and the
  # returned mu/ite/CIs are inverse-transformed below.
  y_train <- as.numeric(Y)
  y_mean  <- mean(y_train)
  y_sd    <- stats::sd(y_train)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1  # constant-Y guard
  y_scaled <- (y_train - y_mean) / y_sd
  # binary outcome -> Bernoulli likelihood on the raw 0/1 scale.
  if (isTRUE(binary_outcome)) {
    outcome_dist <- "bernoulli"
    y_mean <- 0; y_sd <- 1
    y_scaled <- y_train
  }

  model_path <- .cf_new_model_path("cevae", "pt")
  pypath     <- system.file("python", package = "cfomics")

  fit_out <- .cf_py_run(
    function(pypath, model_path, X, T, y_scaled, outcome_dist, feature_dim,
             latent_dim, hidden_dim, num_layers, num_samples, num_epochs,
             batch_size, learning_rate, learning_rate_decay, weight_decay,
             device, fit_seed) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      torch     <- reticulate::import("torch", convert = FALSE)
      np        <- reticulate::import("numpy", convert = FALSE)
      py_random <- reticulate::import("random", convert = FALSE)
      ce_mod    <- reticulate::import_from_path("cevae_extended",
                                                path = pypath, convert = FALSE)
      # Ephemeral worker: seed for reproducibility (no restore needed).
      torch$manual_seed(as.integer(fit_seed))
      np$random$seed(as.integer(fit_seed))
      py_random$seed(as.integer(fit_seed))

      cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
      if (identical(device, "auto")) device <- if (cuda_ok) "gpu" else "cpu"
      use_cuda <- identical(device, "gpu") && cuda_ok
      torch_device <- if (use_cuda) torch$device("cuda") else torch$device("cpu")

      x_tensor <- torch$tensor(X, dtype = torch$float32)$to(torch_device)
      t_tensor <- torch$tensor(as.numeric(T), dtype = torch$float32)$to(torch_device)
      y_tensor <- torch$tensor(as.numeric(y_scaled), dtype = torch$float32)$to(torch_device)

      cevae <- ce_mod$CEVAEExtended(
        feature_dim = as.integer(feature_dim), outcome_dist = outcome_dist,
        latent_dim = as.integer(latent_dim), hidden_dim = as.integer(hidden_dim),
        num_layers = as.integer(num_layers), num_samples = as.integer(num_samples))
      if (isTRUE(use_cuda)) cevae$to(torch_device)

      cevae$fit(x_tensor, t_tensor, y_tensor,
                num_epochs = as.integer(num_epochs),
                batch_size = as.integer(batch_size),
                learning_rate = learning_rate,
                learning_rate_decay = learning_rate_decay,
                weight_decay = weight_decay)

      preds <- reticulate::py_to_r(cevae$predict_y_with_ci(x_tensor, alpha = 0.05))
      ce_mod$save_cevae(cevae, model_path)
      list(
        y0 = as.numeric(preds$y0_hat), y1 = as.numeric(preds$y1_hat),
        ite = as.numeric(preds$ite),
        ate_ci_lower = as.numeric(preds$ate_ci_lower),
        ate_ci_upper = as.numeric(preds$ate_ci_upper),
        ite_ci_lower = as.numeric(preds$ite_ci_lower),
        ite_ci_upper = as.numeric(preds$ite_ci_upper),
        device_effective = tryCatch(
          as.character(reticulate::py_to_r(cevae$get_effective_device())),
          error = function(e) NA_character_))
    },
    pypath = pypath, model_path = model_path,
    X = X_mat, T = T, y_scaled = y_scaled, outcome_dist = outcome_dist,
    feature_dim = feature_dim, latent_dim = latent_dim, hidden_dim = hidden_dim,
    num_layers = num_layers, num_samples = num_samples, num_epochs = num_epochs,
    batch_size = batch_size, learning_rate = learning_rate,
    learning_rate_decay = learning_rate_decay, weight_decay = weight_decay,
    device = device, fit_seed = fit_seed
  )

  # Inverse z-score: mu_t = mu_t_scaled * sd + mean; tau = tau_scaled * sd
  # (the mean cancels in the difference; CI bounds are about tau, so sd-only).
  mu0_hat <- fit_out$y0 * y_sd + y_mean
  mu1_hat <- fit_out$y1 * y_sd + y_mean
  ite     <- fit_out$ite * y_sd
  ate     <- mean(ite)
  ate_ci  <- c(fit_out$ate_ci_lower, fit_out$ate_ci_upper) * y_sd
  ite_ci  <- data.frame(lower = fit_out$ite_ci_lower * y_sd,
                        upper = fit_out$ite_ci_upper * y_sd)

  structure(
    list(
      model = list(model_path = model_path, backend = "pyro.contrib.cevae.CEVAE"),
      res = list(
        tau_hat = ite,
        ate = ate,
        mu0_hat = mu0_hat,
        mu1_hat = mu1_hat,
        ate_ci = ate_ci,
        tau_ci = ite_ci,
        # Y z-score stats travel in meta; X is standardised centrally in
        # preprocessing, so this method contributes no X stats.
        y_zscore = list(Y_mean = y_mean, Y_sd = y_sd),
        x_zscore = NULL,
        device_effective = fit_out$device_effective,
        metadata = list(
          ci_method = "Bayesian posterior",
          hyperparameters = list(
            outcome_dist = outcome_dist,
            latent_dim = as.integer(latent_dim),
            hidden_dim = as.integer(hidden_dim),
            num_layers = as.integer(num_layers),
            num_samples = as.integer(num_samples),
            num_epochs = as.integer(num_epochs),
            batch_size = as.integer(batch_size),
            learning_rate = as.numeric(learning_rate),
            random_state = as.integer(fit_seed)
          )
        )
      )
    ),
    class = "cf_model_cevae"
  )
}

#' Predict from CEVAE model
#'
#' Returns in-sample estimates when \code{newdata} is \code{NULL}. Otherwise
#' re-scores the trained \code{CEVAEExtended} model on \code{newdata} using
#' \code{predict_y()} so that \code{y0}, \code{y1}, \code{ite}, \code{ate},
#' \code{summary} and \code{samples} are all available for new samples.
#'
#' @param object cf_model object with \code{method = "cevae"}.
#' @param newdata Optional matrix or data.frame of new covariates with the
#'   same column order/dimensionality as the training \code{X}. If
#'   \code{NULL}, returns in-sample estimates.
#' @param type One of \code{"ate"}, \code{"tau"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param ... Unused.
#' @return The requested prediction; see \code{type}.
#' @keywords internal
predict_cf_cevae <- function(object,
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

  X_new <- as.matrix(newdata)
  storage.mode(X_new) <- "double"
  feature_dim_train <- length(object$meta$covariate_names)
  if (!is.null(feature_dim_train) && feature_dim_train > 0 &&
      ncol(X_new) != feature_dim_train) {
    rlang::abort(sprintf("newdata has %d columns but fit was trained with %d.",
                 ncol(X_new), feature_dim_train), class = "cfomics_newdata_dim_mismatch")
  }
  hp_seed <- object$res$metadata$hyperparameters$random_state
  if (is.null(hp_seed)) hp_seed <- 0L

  preds <- .cf_py_run(
    function(pypath, model_path, X, seed) {
      sys <- reticulate::import("sys", convert = FALSE)
      sys$path$insert(0L, pypath)
      torch  <- reticulate::import("torch", convert = FALSE)
      ce_mod <- reticulate::import_from_path("cevae_extended",
                                             path = pypath, convert = FALSE)
      cevae <- ce_mod$load_cevae(model_path)
      # Match the input tensor to the model's parameter device (load_cevae
      # keeps the z-prior / cached device in sync with the parameters).
      model_device <- tryCatch(reticulate::iter_next(cevae$parameters())$device,
                               error = function(e) torch$device("cpu"))
      x_tensor <- torch$tensor(X, dtype = torch$float32)$to(model_device)
      p <- reticulate::py_to_r(
        cevae$predict_y_with_ci_seeded(x_tensor, seed = as.integer(seed), alpha = 0.05))
      list(y0 = as.numeric(p$y0_hat), y1 = as.numeric(p$y1_hat),
           ite = as.numeric(p$ite),
           ate_ci_lower = as.numeric(p$ate_ci_lower),
           ate_ci_upper = as.numeric(p$ate_ci_upper))
    },
    pypath = system.file("python", package = "cfomics"),
    model_path = object$model$model$model_path, X = X_new, seed = as.integer(hp_seed)
  )
  y0_new <- preds$y0
  y1_new <- preds$y1
  ite_new <- preds$ite
  # Compute ate on the R side so it stays float64-consistent with mean(ite_new).
  ate_new <- mean(ite_new)
  ate_ci_new <- c(preds$ate_ci_lower, preds$ate_ci_upper)

  # Re-apply the train-fold Y statistics to the newdata predictions (an
  # extension of the data-leakage defence). meta$y_zscore is populated by
  # cf_fit during preprocessing.
  y_zscore <- object$meta$y_zscore
  if (!is.null(y_zscore)) {
    y_mean <- y_zscore$Y_mean
    y_sd   <- y_zscore$Y_sd
    y0_new <- y0_new * y_sd + y_mean
    y1_new <- y1_new * y_sd + y_mean
    ite_new <- ite_new * y_sd
    ate_new <- mean(ite_new)
    ate_ci_new <- ate_ci_new * y_sd
  }

  # CEVAE provides newdata ATE CI but not per-row ITE CI (Pyro Predictive
  # per-sample posterior is computed inside predict_y_with_ci_seeded for the
  # aggregate mean only; ITE CI would require per-sample posterior intervals).
  ite_ci_new <- NA_real_

  .assemble_predict_output(
    type    = type,
    ate     = ate_new,
    ate_ci  = ate_ci_new,
    tau     = ite_new,
    tau_ci  = ite_ci_new,
    mu0_hat = y0_new,
    mu1_hat = y1_new
  )
}
