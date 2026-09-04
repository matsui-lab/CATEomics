#' Fit Bayesian Causal Forests model
#'
#' Heterogeneous treatment-effect estimation via Bayesian Additive Regression
#' Trees (Hahn, Murray & Carvalho 2020). Delegates to \code{bcf::bcf()} with
#' cfomics-side propensity estimation when \code{pihat} is not supplied.
#'
#' Returns the full \code{cfomics_result} contract (\code{tau_hat},
#' \code{mu0_hat}, \code{mu1_hat}, \code{ate} and summary) derived from the
#' MCMC posterior: \code{tau_hat} is the posterior mean of the CATE draws;
#' counterfactual means use the y-hat posterior at \eqn{T=0} / \eqn{T=1}.
#'
#' Binary (0/1) outcomes are modelled on the native continuous scale, so the
#' estimand is the risk difference; the counterfactual means may fall outside
#' \eqn{[0, 1]} and the credible intervals are approximate (bcf has no binary
#' likelihood). For a proper binary likelihood use a meta/DR learner
#' (classifier nuisance) or a neural method.
#'
#' The BART outcome forest can fail to converge when covariates far outnumber
#' samples (n << p); the fit then aborts with a \code{cfomics_bcf_degenerate}
#' error. This is a numerical property of BART on wide data, not all methods
#' behave this way, so \code{\link{cf_compare}} records the failure in
#' \code{metadata$failed_methods} and continues with the remaining methods.
#'
#' @details
#' When \code{pihat = NULL}, propensity scores are estimated internally with
#' a logistic cfomics base learner before the BCF call; the seed passed via
#' \code{random_state} is used for both that submodel and \code{bcf::bcf()}.
#' MCMC tree samples are persisted under \code{save_tree_directory} so that
#' \code{predict(newdata = ...)} can draw posterior samples via
#' \code{bcf:::predict.bcf} without retraining.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param n_burn Integer, MCMC burn-in iterations. Default \code{1000}.
#' @param n_iter Integer, MCMC iterations kept after burn-in.
#'   Default \code{2000}.
#' @param n_chains Integer, number of MCMC chains, passed to
#'   \code{bcf::bcf()} (which runs the chains internally, pools the posterior
#'   draws across them, and persists every chain's tree samples). When
#'   \code{n_chains > 1L} cfomics computes the Gelman-Rubin \eqn{\hat R} on the
#'   per-chain ATE-per-iteration draws; \eqn{\hat R > 1.1} triggers a
#'   \code{cfomics_bcf_mcmc_not_converged} warning. Default \code{4L} (the bcf
#'   default; needed for a between-chain \eqn{\hat R} convergence check).
#' @param pihat Optional numeric vector of pre-computed propensity scores.
#'   If \code{NULL}, estimated internally via a logistic cfomics base
#'   learner.
#' @param save_tree_directory Directory to persist MCMC tree samples.
#'   Required when \code{predict(newdata = ...)} is used. Defaults to a
#'   per-call \code{tempfile()} under \code{tempdir()} to avoid polluting
#'   the working directory. Set to a stable path if the tree samples must
#'   survive beyond the current R session.
#' @param random_state Integer random seed, forwarded to \code{bcf::bcf()}
#'   via its \code{random_seed} argument and used to seed the internal
#'   propensity submodel. When set, yields bit-identical results across
#'   runs with the same inputs. The caller's global \code{.Random.seed}
#'   is saved on entry and restored on exit so seeding does not leak into
#'   the user's session.
#' @param ... Additional arguments forwarded to \code{bcf::bcf()}.
#' @return A list with \code{model} (the \code{bcf::bcf()} fit object,
#'   retaining \code{save_tree_directory} for prediction on \code{newdata})
#'   and \code{res} conforming to the \code{cfomics_result} contract.
#' @keywords internal
cf_fit_bcf <- function(X, T, Y,
                       covariate_names = NULL,
                       n_burn = 1000L,
                       n_iter = 2000L,
                       n_chains = 4L,
                       pihat = NULL,
                       propensity_score = NULL,  # known PS via kfold/dispatch path
                       trim = c(0.025, 0.975),   # unified propensity clip bounds
                       propensity_X = NULL,
                       save_tree_directory = NULL,
                       random_state = NULL,
                       ...) {
  if (!requireNamespace("bcf", quietly = TRUE)) {
    rlang::abort("Package 'bcf' is required for method='bcf'. Please install it.", class = "cfomics_missing_dependency")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  if (is.null(save_tree_directory)) {
    save_tree_directory <- tempfile("cfomics_bcf_trees_")
  }
  dir.create(save_tree_directory, showWarnings = FALSE, recursive = TRUE)

  # Save/restore the caller's global RNG so `random_state` does not leak
  # into `.Random.seed` through the propensity-score fit (cv.glmnet) or
  # bcf::bcf internals.
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
    # Pin the RNG algorithm so bcf's MCMC / PS submodel are reproducible across
    # serial (cf_fit, Mersenne-Twister) and BiocParallel::bplapply (cf_compare,
    # L'Ecuyer-CMRG) contexts. The saved `.Random.seed` restores the kind.
    set.seed(as.integer(random_state), kind = "Mersenne-Twister")
  }

  # Honour a user-supplied known propensity score arriving via the k-fold /
  # cf_compare dispatch path (where it is passed as `propensity_score` rather
  # than the native `pihat`). The top-level cf_fit switch already translates to
  # `pihat`; this covers the OOF / compare path.
  if (is.null(pihat) && !is.null(propensity_score)) {
    pihat <- if (length(propensity_score) == 1L) {
      rep(as.numeric(propensity_score), n)
    } else {
      as.numeric(propensity_score)
    }
  }

  # Estimate propensity score if not provided; persist the PS model so that
  # newdata prediction can auto-compute pihat_new.
  ps_model <- NULL
  if (is.null(pihat)) {
    # propensity_X can subset the X used for PS estimation; the outcome model
    # (bcf BART) still uses the full X_mat, only the PS step is subset.
    # ps_model$colnames stores the subset names so newdata predict can re-apply
    # the same subset.
    .ps_res <- .resolve_propensity_X(X_mat, propensity_X, covariate_names)
    X_ps_mat <- .ps_res$X_ps
    if (requireNamespace("glmnet", quietly = TRUE)) {
      ps_cv <- glmnet::cv.glmnet(X_ps_mat, W, family = "binomial", alpha = 1)
      ps_model <- list(kind = "glmnet", model = ps_cv,
                       colnames = colnames(X_ps_mat), clip = trim)
      pihat <- as.numeric(predict(ps_cv, newx = X_ps_mat, s = "lambda.min",
                                   type = "response"))
    } else {
      df_ps <- data.frame(W = W, X_ps_mat)
      glm_fit <- stats::glm(W ~ ., data = df_ps, family = stats::binomial)
      ps_model <- list(kind = "glm", model = glm_fit,
                       colnames = colnames(X_ps_mat), clip = trim)
      pihat <- stats::predict(glm_fit, type = "response")
    }
    # Unified `trim` API for propensity clipping.
    pihat <- pmax(pmin(pihat, trim[2]), trim[1])
  } else {
    # pihat supplied directly or translated from a known `propensity_score`.
    # For a constant known PS (RCT) persist the value so newdata prediction
    # (incl. the internal OOF test-fold predict) can reuse it; a per-sample PS
    # cannot be mapped to new samples, so it stays "user_supplied".
    const_ps <- if (!is.null(propensity_score) && length(propensity_score) == 1L) {
      as.numeric(propensity_score)[1L]
    } else if (length(unique(round(pihat, 10L))) == 1L) {
      pihat[1L]
    } else NA_real_
    ps_model <- list(kind = if (is.finite(const_ps)) "constant" else "user_supplied",
                     model = NULL, value = const_ps,
                     colnames = colnames(X_mat), clip = trim)
  }

  n_chains <- as.integer(n_chains)
  if (n_chains < 1L) {
    rlang::abort("`n_chains` must be a positive integer.",
                 class = "cfomics_bcf_error")
  }
  extra_bcf <- list(...)
  # Drop meta-args that bcf::bcf() does not accept. cf_compare / cf_fit may
  # forward base-learner unification args (e.g., model_outcome /
  # model_propensity) to every method via `...`; bcf is BART-native and has
  # no base learner, so filter `...` to bcf::bcf() formals (mirrors the
  # cf_fit_grf formals filter).
  .bcf_formals <- names(formals(bcf::bcf))
  if (length(extra_bcf)) {
    extra_bcf <- extra_bcf[names(extra_bcf) %in% .bcf_formals]
  }

  # Helper: run one chain with the given seed.
  run_bcf_chain <- function(seed_c) {
    chain_args <- list(
      y = Y_vec,
      z = W,
      x_control = X_mat,
      x_moderate = X_mat,
      pihat = pihat,
      nburn = as.integer(n_burn),
      nsim = as.integer(n_iter),
      n_chains = n_chains,
      save_tree_directory = save_tree_directory,
      # Give bcf an explicit per-chain log path inside the (temporary) tree
      # directory so it no longer writes bcf_log_*.txt into the cwd. An
      # explicit path is honoured even by parallel chain workers (unlike setwd).
      log_file = tempfile("bcf_log_", tmpdir = save_tree_directory,
                          fileext = ".txt")
    )
    if (!is.null(seed_c) && !"random_seed" %in% names(extra_bcf)) {
      chain_args$random_seed <- as.integer(seed_c)
    }
    chain_args[names(extra_bcf)] <- extra_bcf
    # Capture bcf's MCMC progress printed to stdout so the wrapper stays quiet.
    # The `log_file` set above keeps bcf's log out of the cwd (explicit path
    # honoured by parallel chain workers).
    .chain_fit <- NULL
    .bcf_err <- NULL
    utils::capture.output(
      .chain_fit <- tryCatch(do.call(bcf::bcf, chain_args),
                             error = function(e) { .bcf_err <<- e; NULL })
    )
    if (!is.null(.bcf_err)) {
      .msg <- conditionMessage(.bcf_err)
      # bcf's BART outcome forest degenerates numerically when covariates far
      # outnumber samples; surface it as an actionable error pointing at the
      # comparison workflow rather than bcf's raw "drmu failed".
      if (grepl("drmu", .msg, ignore.case = TRUE)) {
        rlang::abort(
          c(sprintf("bcf's BART outcome model failed to converge (bcf: %s).", .msg),
            i = "This commonly occurs with n << p (far more covariates than samples), where the BART outcome forest degenerates.",
            i = "grf and the meta-learners typically fit the same data; cf_compare() records per-method failures in metadata$failed_methods and continues with the others."),
          class = "cfomics_bcf_degenerate"
        )
      }
      rlang::abort(sprintf("bcf::bcf() failed: %s", .msg),
                   class = "cfomics_bcf_fit_failed")
    }
    .chain_fit
  }

  # A single bcf::bcf() call runs `n_chains` internal MCMC chains: it pools
  # tau / yhat across chains, writes every chain's tree samples into
  # save_tree_directory, and exposes per-chain coda. In-sample summaries and
  # newdata prediction therefore both use all chains and stay consistent.
  # (An earlier design wrapped an R loop that refit bcf n_chains times, each
  # itself multi-chain, into the shared tree directory -- inflating the chain
  # count and leaving only the last refit's trees on disk for newdata.)
  bcf_fit <- run_bcf_chain(random_state)
  tau_samples <- bcf_fit$tau
  yhat_samples <- bcf_fit$yhat
  rhat_val <- NA_real_
  if (n_chains > 1L) {
    # bcf stacks each chain's `n_iter` post-burn draws row-wise; split back into
    # per-chain ATE-per-iteration to compute the Gelman-Rubin Rhat.
    nsim <- nrow(tau_samples) %/% n_chains
    if (nsim >= 2L) {
      ate_by_iter <- vapply(seq_len(n_chains), function(ci) {
        rows <- ((ci - 1L) * nsim + 1L):(ci * nsim)
        rowMeans(tau_samples[rows, , drop = FALSE])
      }, numeric(nsim))
      rhat_val <- .compute_gelman_rubin_rhat(ate_by_iter)
      if (!is.na(rhat_val) && rhat_val > 1.1) {
        rlang::warn(
          sprintf(
            "BCF MCMC convergence not reached (Rhat = %.3f > 1.1). Increase n_iter or n_burn.",
            rhat_val
          ),
          class = "cfomics_bcf_mcmc_not_converged"
        )
      }
    }
  }

  ite <- colMeans(tau_samples)
  ate <- mean(ite)

  ate_samples <- rowMeans(tau_samples)
  ate_ci <- unname(stats::quantile(ate_samples, c(0.025, 0.975),
                                    names = FALSE))
  ate_ci_lower <- ate_ci[1]
  ate_ci_upper <- ate_ci[2]

  mu_samples <- yhat_samples - outer(rep(1, nrow(tau_samples)), W) * tau_samples
  mu0_hat <- colMeans(mu_samples)
  mu1_hat <- mu0_hat + ite

  # Per-sample tau CI (data.frame lower/upper) from the posterior draws.
  # (The canonical res$summary is built centrally by cf_fit / cf_compare via
  # `.cf_build_summary`, so bcf does not assemble its own.)
  ite_ci_lower <- apply(tau_samples, 2, stats::quantile, 0.025,
                        names = FALSE)
  ite_ci_upper <- apply(tau_samples, 2, stats::quantile, 0.975,
                        names = FALSE)
  ite_ci_df <- data.frame(
    lower = as.numeric(ite_ci_lower),
    upper = as.numeric(ite_ci_upper)
  )

  list(
    model    = bcf_fit,
    ps_model = ps_model,   # persist PS model for newdata auto-pihat
    res = list(
      tau_hat = ite,
      ate = ate,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat,
      ate_ci  = c(ate_ci_lower, ate_ci_upper),
      tau_ci  = ite_ci_df,
      propensity = pihat,
      metadata = list(
        ci_method = "Bayesian posterior",
        n_burn = n_burn,
        n_iter = n_iter,
        n_chains = n_chains,
        rhat = rhat_val,
        save_tree_directory = save_tree_directory
      )
    )
  )
}

#' Univariate Gelman-Rubin Rhat statistic
#'
#' Computes \eqn{\hat R = \sqrt{((n-1)/n) W + (1/n) B \over W}} for a
#' samples_matrix with \code{n_iter} rows and \code{n_chains} columns.
#' Returns \code{NA_real_} when \code{n_chains <= 1L}.
#'
#' @keywords internal
#' @noRd
.compute_gelman_rubin_rhat <- function(samples_matrix) {
  if (!is.matrix(samples_matrix) || ncol(samples_matrix) < 2L) {
    return(NA_real_)
  }
  n_iter <- nrow(samples_matrix)
  if (n_iter < 2L) return(NA_real_)
  chain_means <- colMeans(samples_matrix)
  chain_vars <- apply(samples_matrix, 2, stats::var)
  B <- n_iter * stats::var(chain_means)
  W <- mean(chain_vars)
  if (W <= 0) return(NA_real_)
  var_hat <- ((n_iter - 1) / n_iter) * W + (1 / n_iter) * B
  sqrt(var_hat / W)
}

#' Predict from BCF model
#'
#' Supports both in-sample prediction (\code{newdata = NULL}) and
#' out-of-sample prediction via \code{bcf::predict()} on the saved tree
#' samples. Out-of-sample prediction requires that the fit was performed
#' with a persistent \code{save_tree_directory} and that the user supplies
#' \code{pihat_new} for the new sample.
#'
#' @param object cf_model object with method="bcf"
#' @param newdata Optional matrix/data.frame of new covariates. If NULL,
#'   returns in-sample estimates.
#' @param pihat_new Numeric vector of propensity scores for \code{newdata}.
#'   Required when \code{newdata} is given.
#' @param type Type of prediction: "ate", "tau", "mu0", "mu1", "summary"
#' @param ... Additional arguments forwarded to \code{bcf::predict()}
#' @return Predicted values
#' @keywords internal
predict_cf_bcf <- function(object, newdata = NULL, pihat_new = NULL,
                           type = "all", ...) {
  res <- object$res
  if (is.null(newdata)) {
    return(.cf_insample_predict(res, type))
  }

  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.", class = "cfomics_unsupported_predict_type")
  }

  if (!requireNamespace("bcf", quietly = TRUE)) {
    rlang::abort("Package 'bcf' is required for method='bcf'. Please install it.", class = "cfomics_missing_dependency")
  }

  tree_dir <- res$metadata$save_tree_directory
  if (is.null(tree_dir) || !dir.exists(tree_dir)) {
    rlang::abort(
      paste0(
        "BCF newdata prediction requires a persistent save_tree_directory. ",
        "Refit with cf_fit(..., save_tree_directory = '/stable/path') to ",
        "enable out-of-sample prediction."
      ),
      class = "cfomics_bcf_newdata_no_trees"
    )
  }

  X_new <- as.matrix(newdata)
  n_new <- nrow(X_new)

  # Auto-compute pihat_new from the persisted PS model.
  if (is.null(pihat_new)) {
    ps_model <- object$model$ps_model
    if (identical(ps_model$kind, "constant")) {
      # Known constant PS (RCT) reused for newdata / OOF folds.
      cl <- ps_model$clip %||% c(0.025, 0.975)
      pihat_new <- rep(pmax(pmin(as.numeric(ps_model$value), cl[2L]), cl[1L]), n_new)
    } else if (is.null(ps_model) || identical(ps_model$kind, "user_supplied")) {
      rlang::abort(
        paste0(
          "BCF newdata prediction requires the persisted PS model ",
          "or an explicit 'pihat_new'. Refit with cf_fit(method = 'bcf'), ",
          "or pre-compute pihat_new from your preferred PS model."
        ),
        class = "cfomics_bcf_newdata_no_ps"
      )
    } else {
      # When the model was fit with `propensity_X`, the persisted PS model was
      # trained on the column subset (ps_model$colnames carries those subset
      # names). Subset newdata to match BEFORE the glmnet / glm predict call —
      # otherwise the predict path complains about column-count mismatch.
      X_new_for_ps <- X_new
      px <- object$meta$propensity_X
      if (!is.null(px)) {
        px_idx <- match(px, colnames(X_new))
        if (anyNA(px_idx)) {
          rlang::abort(sprintf(
            "newdata is missing propensity_X column(s): %s",
            paste(px[is.na(px_idx)], collapse = ", ")
          ), class = "cfomics_newdata_missing_column")
        }
        X_new_for_ps <- X_new[, px_idx, drop = FALSE]
      }
      if (identical(ps_model$kind, "glmnet")) {
        pihat_new <- as.numeric(predict(
          ps_model$model, newx = X_new_for_ps,
          s = "lambda.min", type = "response"
        ))
      } else if (identical(ps_model$kind, "glm")) {
        df_new <- as.data.frame(X_new_for_ps)
        if (!is.null(ps_model$colnames)) {
          colnames(df_new) <- ps_model$colnames
        }
        pihat_new <- as.numeric(stats::predict(
          ps_model$model, newdata = df_new, type = "response"
        ))
      } else {
        rlang::abort(sprintf("Unsupported PS model kind '%s'.", ps_model$kind), class = "cfomics_unsupported_ps_model")
      }
      # Clip newdata PS with the fit-time bounds (persisted trim); fall back to
      # the default trim, not 0.01/0.99.
      ps_clip <- ps_model$clip %||% c(0.025, 0.975)
      pihat_new <- pmax(pmin(pihat_new, ps_clip[2]), ps_clip[1])
    }
  }

  if (length(pihat_new) != n_new) {
    rlang::abort(sprintf(
      "length(pihat_new) = %d but nrow(newdata) = %d",
      length(pihat_new), n_new
    ), class = "cfomics_newdata_length_mismatch")
  }

  if (n_new < 2L) {
    rlang::abort(paste0("BCF newdata prediction requires at least 2 rows (bcf::predict ",
         "itself requires z_pred with both 0's and 1's)."), class = "cfomics_bcf_newdata_min_rows")
  }

  z_pred <- rep(c(0, 1), length.out = n_new)

  bcf_pred <- stats::predict(
    object = object$model$model,
    x_predict_control  = X_new,
    x_predict_moderate = X_new,
    pi_pred = pihat_new,
    z_pred  = z_pred,
    save_tree_directory = tree_dir,
    # predict.bcf() spawns a PSOCK cluster whose default log_file is
    # `./bcf_log_<timestamp>.txt` (written to the cwd). Redirect it into
    # tree_dir so newdata / OOF-fold prediction does not leak a log into the
    # cwd (and into tests/testthat during the suite).
    log_file = tempfile("bcf_log_", tmpdir = tree_dir, fileext = ".txt"),
    verbose = FALSE,
    ...
  )

  tau_new_samples  <- bcf_pred$tau
  yhat_new_samples <- bcf_pred$yhat
  n_draws <- nrow(tau_new_samples)
  mu_new_samples <- yhat_new_samples -
    outer(rep(1, n_draws), z_pred) * tau_new_samples

  ite_new <- colMeans(tau_new_samples)
  y0_new  <- colMeans(mu_new_samples)
  y1_new  <- y0_new + ite_new
  ate_new <- mean(ite_new)

  # CI on newdata via posterior quantiles (Bayesian credible intervals).
  ite_ci_new <- data.frame(
    lower = apply(tau_new_samples, 2, stats::quantile, 0.025, names = FALSE),
    upper = apply(tau_new_samples, 2, stats::quantile, 0.975, names = FALSE)
  )
  ate_samples_new <- rowMeans(tau_new_samples)
  ate_ci_new <- unname(stats::quantile(ate_samples_new, c(0.025, 0.975),
                                        names = FALSE))

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
