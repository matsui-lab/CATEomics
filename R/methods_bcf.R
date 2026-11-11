#' Fit Bayesian Causal Forests model
#'
#' Heterogeneous treatment-effect estimation via Bayesian Additive Regression
#' Trees (Hahn, Murray & Carvalho 2020). Delegates to \code{bcf::bcf()} with
#' CATEomics-side propensity estimation when \code{pihat} is not supplied.
#'
#' Returns \code{tau_hat}, \code{mu0_hat}, \code{mu1_hat} and \code{ate}
#' derived from the MCMC posterior; \code{res$summary} is assembled centrally
#' by the entry that called this. \code{tau_hat} is the posterior mean of the CATE draws;
#' counterfactual means use the y-hat posterior at \eqn{T=0} / \eqn{T=1}.
#'
#' Binary (0/1) outcomes are modelled on the native continuous scale, so the
#' estimand is the risk difference; the counterfactual means may fall outside
#' \eqn{[0, 1]} and the credible intervals are approximate (bcf has no binary
#' likelihood).
#'
#' The BART outcome forest can fail to converge when the rows reaching it are
#' not comfortably more numerous than the covariates; the fit then aborts with
#' a \code{cateomics_bcf_degenerate} error, and \code{\link{cf_compare}} records
#' it in \code{metadata$failed_methods} and continues with the remaining
#' methods. \code{cf_compare()} fits each fold on
#' \code{n - ceiling(n / n_folds)} rows, so a fit that succeeds on the whole
#' sample can still fail in a fold.
#'
#' @details
#' When \code{pihat = NULL}, propensity scores are estimated internally with
#' the \code{model_propensity} base learner before the BCF call; the seed
#' passed via \code{random_state} is used for both that submodel and
#' \code{bcf::bcf()}. The fitted submodel is persisted so that
#' \code{predict(newdata = ...)} re-applies the same learner.
#' MCMC tree samples are persisted under \code{save_tree_directory} so that
#' \code{predict(newdata = ...)} can draw posterior samples via
#' \code{bcf:::predict.bcf} without retraining.
#'
#' @param X Covariate matrix.
#' @param T Binary treatment vector (0/1).
#' @param Y Outcome vector.
#' @param covariate_names Character vector of covariate names.
#' @param nburn Integer, MCMC burn-in iterations. Default \code{1000}.
#' @param nsim Integer, MCMC iterations kept after burn-in.
#'   Default \code{2000}.
#' @param n_chains Integer, number of MCMC chains, passed to
#'   \code{bcf::bcf()} (which runs the chains internally, pools the posterior
#'   draws across them, and persists every chain's tree samples). When
#'   \code{n_chains > 1L} CATEomics computes the Gelman-Rubin \eqn{\hat R} on the
#'   per-chain ATE-per-iteration draws; \eqn{\hat R > 1.1} triggers a
#'   \code{cateomics_bcf_mcmc_not_converged} warning. Default \code{4L} (the bcf
#'   default).
#' @param pihat Optional numeric vector of pre-computed propensity scores.
#'   If \code{NULL}, estimated internally with the learner named by
#'   \code{model_propensity}.
#' @param model_propensity Base-learner token for the internal propensity
#'   model, used when \code{pihat} and \code{propensity_score} are both
#'   \code{NULL}. Default \code{"logistic"}, the same default as the other
#'   methods that estimate a propensity. \code{bcf::bcf()} requires
#'   \code{pihat} and names no estimator of its own, so this is CATEomics's
#'   choice. The outcome model is the BART forest and is not pluggable.
#' @param propensity_params Named list of hyperparameters for that propensity
#'   model, validated against the named learner.
#' @param save_tree_directory Directory to persist MCMC tree samples.
#'   Required when \code{predict(newdata = ...)} is used. Defaults to a
#'   per-call \code{tempfile()} under \code{tempdir()}; set a stable path to
#'   keep the tree samples beyond the current R session.
#' @param random_state Integer random seed, forwarded to \code{bcf::bcf()}
#'   via its \code{random_seed} argument and used to seed the internal
#'   propensity submodel. It does not by itself make a fit bit-identical
#'   across runs: \code{bcf::bcf()} defaults its \code{n_threads} to two
#'   fewer than the machine reports and combines each iteration's partial
#'   results in whatever order the threads finish, which the seed does not
#'   fix. Measured on 200 rows with a short chain, two fits with the same seed
#'   differ by about 5e-08 in sample, and by as much as 8e-03 on some training
#'   subsets. Pass \code{params = list(n_threads = 1L)} for a bit-identical
#'   fit. The caller's global \code{.Random.seed} is saved on entry and
#'   restored on exit so seeding does not leak into the user's session.
#' @param params Named list of additional arguments forwarded to
#'   \code{bcf::bcf()} (flat; there is no nested block). Names are checked on
#'   the R side against the arguments of this fitter and of
#'   \code{bcf::bcf()}, so an unknown name aborts before sampling starts.
#' @param propensity_score Known propensity score: a scalar in (0, 1) or a
#'   length-n vector. \code{NULL} (default) estimates it internally.
#' @param propensity_X Optional character vector of covariate names to
#'   restrict the propensity model to. Reached through \code{\link{cf_fit}} or
#'   \code{\link{cf_compare}}, leaving it unset fits the propensity on the
#'   clinical covariates alone; called directly, \code{NULL} means every
#'   column of the design matrix.
#' @param clip Length-2 numeric, the bounds the propensity is held within
#'   (default \code{c(0.025, 0.975)}). Values are replaced by the bound; no
#'   sample is dropped.
#' @return A list with \code{model} (the \code{bcf::bcf()} fit object,
#'   retaining \code{save_tree_directory} for prediction on \code{newdata})
#'   and \code{res} conforming to the \code{cateomics_result} contract.
#' @references
#' Hahn, Murray & Carvalho (2020) Bayesian Regression Tree Models for Causal
#' Inference: Regularization, Confounding, and Heterogeneous Effects,
#' \emph{Bayesian Anal.} 15(3):965-1056.
#'
#' @keywords internal
cf_fit_bcf <- function(X, T, Y,
                       covariate_names = NULL,
                       nburn = 1000L,
                       nsim = 2000L,
                       n_chains = 4L,
                       pihat = NULL,
                       propensity_score = NULL,  # known PS via kfold/dispatch path
                       clip = .CF_PS_CLIP,   # unified propensity clip bounds
                       propensity_X = NULL,
                       model_propensity = "logistic",
                       propensity_params = list(),
                       save_tree_directory = NULL,
                       random_state = NULL,
                       params = list()) {
  if (!requireNamespace("bcf", quietly = TRUE)) {
    rlang::abort(
      c("Package 'bcf' is required for method = 'bcf'.",
        i = "Install with install.packages('bcf')."),
      class = "cateomics_missing_dependency")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  if (is.null(save_tree_directory)) {
    save_tree_directory <- tempfile("cateomics_bcf_trees_")
  }
  dir.create(save_tree_directory, showWarnings = FALSE, recursive = TRUE)

  # Save/restore the caller's global RNG so `random_state` does not leak
  # into `.Random.seed` through the propensity submodel or bcf::bcf internals.
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
    # The shared base-learner path, so the propensity is specified the same way
    # here as for every other method that estimates one, and so the checks that
    # path carries -- the fold count a penalised learner selects over, and the
    # rows it needs -- apply without being restated.
    ps_fit <- .fit_base_learner(X_ps_mat, W, type = "classification",
                                learner = model_propensity,
                                random_state = random_state,
                                .params = propensity_params)
    ps_model <- list(kind = "base_learner", learner = model_propensity,
                     model = ps_fit, colnames = colnames(X_ps_mat), clip = clip)
    pihat <- .predict_base_learner(ps_fit, X_ps_mat)
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
                     colnames = colnames(X_mat), clip = clip)
  }
  # One clip for both branches: an estimated and a supplied propensity are
  # bounded the same way.
  pihat <- pmax(pmin(pihat, clip[2]), clip[1])

  n_chains <- as.integer(n_chains)
  if (n_chains < 1L) {
    rlang::abort("`n_chains` must be a positive integer.",
                 class = "cateomics_bcf_error")
  }
  # Remaining hyperparameters go straight to bcf::bcf(); they were validated
  # against its formals upstream.
  extra_bcf <- params

  # Helper: run the fit with the given seed. bcf itself runs `n_chains`
  # chains, so this is called once.
  run_bcf_chain <- function(seed_c) {
    chain_args <- list(
      y = Y_vec,
      z = W,
      x_control = X_mat,
      x_moderate = X_mat,
      pihat = pihat,
      nburn = as.integer(nburn),
      nsim = as.integer(nsim),
      n_chains = n_chains,
      save_tree_directory = save_tree_directory,
      # Give bcf an explicit per-chain log path inside the (temporary) tree
      # directory so it no longer writes bcf_log_*.txt into the cwd. An
      # explicit path is honoured even by parallel chain workers (unlike setwd).
      log_file = tempfile("bcf_log_", tmpdir = save_tree_directory,
                          fileext = ".txt")
    )
    # `random_seed` is refused in `params`, so nothing can reach it from
    # outside; the per-chain seed is always CATEomics's own.
    if (!is.null(seed_c)) {
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
            i = "The draw of the terminal-node means fails when the rows reaching the forest are not comfortably more numerous than the covariates.",
            i = "`cf_compare()` fits each method on `n - ceiling(n / n_folds)` rows per fold, so a fit that succeeds on the whole sample can still fail in a fold; raising `n_folds` or lowering `n_features` leaves more rows per covariate.",
            i = "grf and the meta-learners typically fit the same data; cf_compare() records per-method failures in metadata$failed_methods and continues with the others."),
          class = "cateomics_bcf_degenerate"
        )
      }
      rlang::abort(sprintf("bcf::bcf() failed: %s", .msg),
                   class = "cateomics_bcf_fit_failed")
    }
    .chain_fit
  }

  # A single bcf::bcf() call runs `n_chains` internal MCMC chains: it pools
  # tau / yhat across chains, writes every chain's tree samples into
  # save_tree_directory, and exposes per-chain coda. In-sample summaries and
  # newdata prediction therefore both use all chains and stay consistent.
  bcf_fit <- run_bcf_chain(random_state)
  tau_samples <- bcf_fit$tau
  yhat_samples <- bcf_fit$yhat
  rhat_val <- NA_real_
  if (n_chains > 1L) {
    # bcf stacks each chain's `nsim` post-burn draws row-wise; split back into
    # per-chain ATE-per-iteration to compute the Gelman-Rubin Rhat.
    draws_per_chain <- nrow(tau_samples) %/% n_chains
    if (draws_per_chain >= 2L) {
      ate_by_iter <- vapply(seq_len(n_chains), function(ci) {
        rows <- ((ci - 1L) * draws_per_chain + 1L):(ci * draws_per_chain)
        rowMeans(tau_samples[rows, , drop = FALSE])
      }, numeric(draws_per_chain))
      rhat_val <- .compute_gelman_rubin_rhat(ate_by_iter)
      if (!is.na(rhat_val) && rhat_val > 1.1) {
        rlang::warn(
          sprintf(
            "BCF MCMC convergence not reached (Rhat = %.3f > 1.1). Increase nsim or nburn.",
            rhat_val
          ),
          class = "cateomics_bcf_mcmc_not_converged"
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
      ps      = pihat,
      # bcf fits its PS submodel on the full sample, so an estimated pihat is
      # an in-sample prediction.
      ps_type = if (ps_model$kind %in% c("constant", "user_supplied"))
                  "known" else "in-sample",
      metadata = list(
        ci_method = "Bayesian posterior",
        nburn = nburn,
        nsim = nsim,
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
#' samples_matrix with one row per kept draw and \code{n_chains} columns.
#' Returns \code{NA_real_} when \code{n_chains <= 1L}.
#'
#' @keywords internal
#' @noRd
.compute_gelman_rubin_rhat <- function(samples_matrix) {
  if (!is.matrix(samples_matrix) || ncol(samples_matrix) < 2L) {
    return(NA_real_)
  }
  n_draws <- nrow(samples_matrix)
  if (n_draws < 2L) return(NA_real_)
  chain_means <- colMeans(samples_matrix)
  chain_vars <- apply(samples_matrix, 2, stats::var)
  B <- n_draws * stats::var(chain_means)
  W <- mean(chain_vars)
  if (W <= 0) return(NA_real_)
  var_hat <- ((n_draws - 1) / n_draws) * W + (1 / n_draws) * B
  sqrt(var_hat / W)
}

#' Predict from BCF model
#'
#' Supports both in-sample prediction (\code{newdata = NULL}) and
#' out-of-sample prediction through \code{bcf}'s predict method
#' (\code{bcf:::predict.bcf}, reached by S3 dispatch) on the saved tree
#' samples. Out-of-sample prediction requires that the fit was performed with
#' a persistent \code{save_tree_directory}.
#'
#' @param object cf_model object with method="bcf"
#' @param newdata Optional matrix/data.frame of new covariates. If NULL,
#'   returns in-sample estimates.
#' @param pihat_new Numeric vector of propensity scores for \code{newdata}.
#'   Optional: when omitted it is recomputed from the propensity model the fit
#'   persisted, or reused when the fit was given a constant known score. It is
#'   required only when the fit was given a per-sample known score
#'   (\code{cateomics_bcf_newdata_no_ps}).
#' @param type One of \code{"all"} (the default), \code{"ate"},
#'   \code{"ate_ci"}, \code{"tau"}, \code{"tau_ci"}, \code{"mu0"},
#'   \code{"mu1"}, \code{"ps"}, \code{"summary"}.
#' @param ... Additional arguments forwarded to \code{bcf}'s predict
#'   method
#' @return Predicted values
#' @keywords internal
predict_cf_bcf <- function(object, newdata = NULL, pihat_new = NULL,
                           type = "all", ...) {
  res <- object$res
  if (is.null(newdata)) {
    return(.cf_insample_predict(res, type))
  }

  if (type == "summary") {
    rlang::abort("type 'summary' is not supported with newdata.", class = "cateomics_unsupported_predict_type")
  }

  if (!requireNamespace("bcf", quietly = TRUE)) {
    rlang::abort(
      c("Package 'bcf' is required for method = 'bcf'.",
        i = "Install with install.packages('bcf')."),
      class = "cateomics_missing_dependency")
  }

  tree_dir <- res$metadata$save_tree_directory
  if (is.null(tree_dir) || !dir.exists(tree_dir)) {
    rlang::abort(
      paste0(
        "BCF newdata prediction requires a persistent save_tree_directory. ",
        "Refit with cf_fit(..., save_tree_directory = '/stable/path') to ",
        "enable out-of-sample prediction."
      ),
      class = "cateomics_bcf_newdata_no_trees"
    )
  }

  X_new <- as.matrix(newdata)
  n_new <- nrow(X_new)

  # Auto-compute pihat_new from the persisted PS model. A caller-supplied
  # pihat_new is a known propensity for the new rows.
  ps_type_new <- if (is.null(pihat_new)) NA_character_ else "known"
  if (is.null(pihat_new)) {
    ps_model <- object$model$ps_model
    if (identical(ps_model$kind, "constant")) {
      ps_type_new <- "known"
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
        class = "cateomics_bcf_newdata_no_ps"
      )
    } else {
      # When the model was fit with `propensity_X`, the persisted PS model was
      # trained on the column subset (ps_model$colnames carries those subset
      # names). Subset newdata to match BEFORE the glmnet / glm predict call —
      # otherwise the predict path complains about column-count mismatch.
      ps_type_new  <- "out-of-sample"
      X_new_for_ps <- X_new
      px <- object$meta$propensity_X
      if (!is.null(px)) {
        px_idx <- match(px, colnames(X_new))
        if (anyNA(px_idx)) {
          rlang::abort(sprintf(
            "newdata is missing propensity_X column(s): %s",
            paste(px[is.na(px_idx)], collapse = ", ")
          ), class = "cateomics_newdata_missing_column")
        }
        X_new_for_ps <- X_new[, px_idx, drop = FALSE]
      }
      if (identical(ps_model$kind, "base_learner")) {
        X_pred <- X_new_for_ps
        if (!is.null(ps_model$colnames) &&
            ncol(X_pred) == length(ps_model$colnames)) {
          colnames(X_pred) <- ps_model$colnames
        }
        pihat_new <- as.numeric(.predict_base_learner(ps_model$model, X_pred))
      } else {
        rlang::abort(sprintf("Unsupported PS model kind '%s'.", ps_model$kind),
                     class = "cateomics_unsupported_ps_model")
      }
      # Clip newdata PS with the fit-time bounds (persisted clip), falling
      # back to the default clip.
      ps_clip <- ps_model$clip %||% c(0.025, 0.975)
      pihat_new <- pmax(pmin(pihat_new, ps_clip[2]), ps_clip[1])
    }
  }

  if (length(pihat_new) != n_new) {
    rlang::abort(sprintf(
      "length(pihat_new) = %d but nrow(newdata) = %d",
      length(pihat_new), n_new
    ), class = "cateomics_newdata_length_mismatch")
  }

  # A caller-supplied pihat_new is bounded like every other propensity: the
  # training folds are clipped, so the held-out rows must be too.
  .bcf_clip <- object$model$ps_model$clip %||% c(0.025, 0.975)
  pihat_new <- pmax(pmin(pihat_new, .bcf_clip[2L]), .bcf_clip[1L])

  if (n_new < 2L) {
    rlang::abort(paste0("BCF newdata prediction requires at least 2 rows (bcf::predict ",
         "itself requires z_pred with both 0's and 1's)."), class = "cateomics_bcf_newdata_min_rows")
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
    mu1_hat = y1_new,
    ps      = as.numeric(pihat_new),
    ps_type = ps_type_new
  )
}
