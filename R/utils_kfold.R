#' Stratified K-fold split (base R, no external deps)
#'
#' Produces a list of K integer vectors (test-fold indices) stratified
#' on a binary treatment vector \code{T}. Each treatment group is split
#' into K (approximately) equal-sized folds and combined.
#'
#' Generates the K-fold split once (a shared split with identical fold
#' structure across all methods and the common nuisance).
#'
#' @param T Binary treatment vector (0/1).
#' @param n_folds Integer K (default 5).
#' @param random_state Integer seed or NULL. When provided, the caller's
#'   global RNG state is saved on entry and restored on exit so seeding
#'   does not leak.
#' @return List of length \code{n_folds}; element k is the integer
#'   indices (1-based) of observations in fold k.
#' @keywords internal
#' @noRd
.create_stratified_folds <- function(T, n_folds = 5L, random_state = NULL) {
  n_folds <- as.integer(n_folds)
  stopifnot(n_folds >= 2L, length(T) >= n_folds)

  if (!is.null(random_state)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    # Pin the RNG algorithm so fold assignment is identical whether generated
    # standalone (cf_fit, Mersenne-Twister) or inside BiocParallel::bplapply
    # (cf_compare), which switches the stream to L'Ecuyer-CMRG. The saved
    # `.Random.seed` restores the caller's kind on exit.
    set.seed(as.integer(random_state), kind = "Mersenne-Twister")
  }

  idx_0 <- which(T == 0)
  idx_1 <- which(T == 1)

  shuffle_split <- function(idx) {
    if (length(idx) < n_folds) {
      # When a stratum has fewer observations than folds, fall back to
      # round-robin allocation so every fold receives at least 0 indices.
      sh <- sample(idx)
      out <- vector("list", n_folds)
      for (k in seq_along(sh)) {
        ki <- ((k - 1L) %% n_folds) + 1L
        out[[ki]] <- c(out[[ki]], sh[k])
      }
      lapply(out, function(x) if (is.null(x)) integer(0) else x)
    } else {
      sh <- sample(idx)
      split(sh, rep(seq_len(n_folds), length.out = length(sh)))
    }
  }

  folds_0 <- shuffle_split(idx_0)
  folds_1 <- shuffle_split(idx_1)
  lapply(seq_len(n_folds), function(k) {
    sort(c(folds_0[[k]], folds_1[[k]]))
  })
}

#' Map fold indices to a per-observation fold ID vector (length n).
#' @keywords internal
#' @noRd
.folds_to_fold_id <- function(folds, n) {
  out <- integer(n)
  for (k in seq_along(folds)) {
    out[folds[[k]]] <- as.integer(k)
  }
  out
}


#' Null-coalescing operator (local convenience for the helper)
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# matrix-signature K-fold OOF helper
# ---------------------------------------------------------------------------

#' K-fold OOF predict on preprocessed matrix X
#'
#' Body of the partition='oof' implementation, called after cf_fit Step 1-5
#' (i.e. on the already-preprocessed X). Runs the K-fold loop on a preprocessed
#' numeric matrix X, so no formula/data input is needed.
#'
#' Design:
#' - Calls each method's low-level `cf_fit_<method>(X, T, Y, ...)` directly
#'   (no cf_fit recursion, hence no double preprocessing; operates on the
#'   already-preprocessed X subset)
#' - Wraps the returned inner fit as `"cf_model"` and obtains OOF tau via
#'   `predict()` (X[test_idx, ] matrix newdata)
#' - The K-fold OOF nuisance is computed once via `.compute_common_nuisance_oof()`
#' - The full-data fit is also called directly; the outer cf_fit Step 6-7 wraps
#'   it with meta + class
#'
#' @param X Preprocessed numeric matrix (n x p)
#' @param T_vec Integer vector (n, binary 0/1)
#' @param Y_vec Numeric vector (n)
#' @param method Character method name (CATE-capable only)
#' @param cov_names Character vector of covariate names (incl. post-expansion
#'   dummy names)
#' @param n_folds K (default 5L)
#' @param random_state Integer seed
#' @param keep_kfold_models Logical
#' @param eval_nuisance Base-learner name (rf/lasso/...)
#' @param precomputed_nuisance Optional list to bypass internal nuisance
#'   computation (for cf_compare's shared nuisance path)
#' @param ... Forward to `cf_fit_<method>` (method-specific args)
#'
#' @return list:
#'   * `fit` = method-specific full-data fit object (unwrapped; class attached
#'     by the caller)
#'   * `oof_ite` = numeric(n) K-fold OOF tau
#'   * `nuisance` = list(mu0_oof, mu1_oof, pi_oof) K-fold OOF nuisance
#'   * `fold_cache` = list of K fit objects (when keep_kfold_models = TRUE)
#'
#' @keywords internal
#' @noRd

# Warn on unknown arguments forwarded via `...` to
# `.kfold_oof_predict_matrix` / per-method `cf_fit_<method>()` from cf_compare.
# Otherwise typos such as `base_learner_outcome` (correct name: `model_outcome`)
# are silently swallowed, running the fit with the default learner instead of
# the user-requested one.
#
# Strategy: an arg is treated as a typo if it is not in
#   (1) the dispatched method's own formals (or any formal of another
#       cfomics method, if the dispatched method has `...`)
#   (2) cf_compare's own formals (top-level wrapper passes its named args
#       via `...` for some methods)
#   (3) the `.kfold_oof_predict_matrix` formals (k-fold helper params)
#   (4) the preprocessing meta-args strip list (feature_method, batch, ...)
# Emits an `rlang::warn()` with the closest fuzzy match as a hint.
.cfomics_allowed_meta_args <- c(
  # preprocessing meta-args (consumed by as_cf_data / cf_fit Step 1-3)
  "partition", "n_folds", "feature_method", "n_features",
  "preprocess", "batch", "pathway_method", "gene_sets",
  "assay_name", "experiments", "assay_cols", "assay_pattern",

  "._skip_nuisance",
  "trim",
  # Facade named-args: cf_fit / cf_compare entry args
  "Y", "T", "X_baseline", "X_rnaseq", "X_mirna", "X_methyl",
  "X_protein", "X_mutation", "X_cnv", "X_assays", "propensity_X",
  # cf_compare / k-fold helper params likely to be re-forwarded via ...
  "metric", "eval_nuisance", "eval_propensity_score", "random_state",
  "keep_kfold_models", "verbose", "device",
  "n_bootstrap", "bootstrap_ci",
  "covariate_names", "data"
)

.cfomics_validate_dot_args <- function(dot_args, method) {
  if (length(dot_args) == 0L) return(invisible(NULL))
  nm_dots <- names(dot_args)
  if (is.null(nm_dots)) return(invisible(NULL))
  nm_dots <- nm_dots[nzchar(nm_dots)]
  if (length(nm_dots) == 0L) return(invisible(NULL))

  ns_pkg <- tryCatch(asNamespace("cfomics"), error = function(e) NULL)
  if (is.null(ns_pkg)) return(invisible(NULL))

  meth_fn <- tryCatch(get(paste0("cf_fit_", method), envir = ns_pkg,
                          inherits = FALSE), error = function(e) NULL)
  if (is.null(meth_fn)) return(invisible(NULL))
  meth_formals_all <- names(formals(meth_fn))
  has_dots <- "..." %in% meth_formals_all
  meth_formals <- setdiff(meth_formals_all, "...")

  all_method_names <- c("grf", "bcf",
                         "slearner", "tlearner", "xlearner", "drlearner",
                         "cevae", "ganite", "tarnet", "dragonnet")
  all_method_formals <- unique(unlist(lapply(all_method_names, function(m) {
    fn <- tryCatch(get(paste0("cf_fit_", m), envir = ns_pkg,
                       inherits = FALSE), error = function(e) NULL)
    if (is.null(fn)) character() else setdiff(names(formals(fn)), "...")
  })))

  for (nm in nm_dots) {
    if (nm %in% meth_formals) next
    if (nm %in% .cfomics_allowed_meta_args) next
    if (has_dots && nm %in% all_method_formals) next
    # Candidate typo. Restrict suggestion candidates to "user-meaningful"
    # method args (excluding implementation-only positional args like X/T/Y)
    # so the hint is more useful than e.g. "fluctuation".
    user_meaningful <- setdiff(all_method_formals,
                                c("X", "T", "Y", "covariate_names",
                                  "random_state"))
    nearest <- character(0)
    # Prefer fuzzy match via agrep (Levenshtein-bounded); fall back to
    # adist-based nearest only if no fuzzy hit.
    hit_idx <- tryCatch(
      agrep(nm, user_meaningful, ignore.case = TRUE,
                   max.distance = list(all = 0.5), value = FALSE),
      error = function(e) integer(0)
    )
    if (length(hit_idx) > 0L) {
      nearest <- user_meaningful[hit_idx[1L]]
    } else if (length(user_meaningful) > 0L) {
      d <- utils::adist(nm, user_meaningful, ignore.case = TRUE)
      # Only suggest if the closest is within (length(nm) / 2) edits.
      if (min(d) <= max(3L, floor(nchar(nm) / 2))) {
        nearest <- user_meaningful[which.min(d)]
      }
    }
    msg <- sprintf(
      paste0("argument '%s' is not recognised by ",
              "cf_fit_%s() and is being silently ignored."),
      nm, method
    )
    if (length(nearest) == 1L) {
      msg <- c(msg, i = sprintf("Did you mean '%s'?", nearest))
    }
    rlang::warn(msg, class = "cfomics_unknown_arg")
  }
  invisible(NULL)
}

.kfold_oof_predict_matrix <- function(X, T_vec, Y_vec, method,
                                       cov_names = NULL,
                                       n_folds = 5L,
                                       random_state = NULL,
                                       keep_kfold_models = FALSE,
                                       eval_nuisance = "rf",
                                       eval_propensity_score = NULL,
                                       precomputed_nuisance = NULL,
                                       n_bootstrap = 500L,
                                       bootstrap_ci = FALSE,
                                       device = "cpu",
                                       propensity_X = NULL,
                                       propensity_score = NULL,
                                       trim = c(0.025, 0.975),
                                       binary_outcome = FALSE,
                                       ...) {
  n <- length(Y_vec)
  if (is.null(cov_names) || length(cov_names) == 0L) {
    cov_names <- paste0("X", seq_len(ncol(X)))
  }

  # meta-args flowing in from cf_compare / cf_fit would fire an
  # `unused argument` error in the per-method `cf_fit_<method>()` or its
  # downstream (e.g., grf::causal_forest), so filter them before dispatching via
  # `.dispatch_method_lowlevel`. `device` is handled as an explicit param on a
  # separate path (see `.dispatch_method_lowlevel`); here we strip partition /
  # feature_method / n_features / preprocess / batch / pathway_method /
  # gene_sets / assay_* / ._skip_nuisance. `propensity_score` and `trim` are
  # explicit formals (not stripped) so they reach the PS-using backends.
  .cfomics_meta_args_strip <- c(
    "partition", "n_folds", "feature_method", "n_features",
    "preprocess", "batch", "pathway_method", "gene_sets",
    "assay_name", "experiments", "assay_cols", "assay_pattern",

  "._skip_nuisance",
    # Facade named-args + propensity_X are accepted at the top-level wrapper
    # (cf_fit / cf_compare) but must not leak into per-method backends (e.g.,
    # grf::causal_forest) via `...`. propensity_X is wired into per-method PS
    # fits explicitly.
    "Y", "T", "X_baseline", "X_rnaseq", "X_mirna", "X_methyl",
    "X_protein", "X_mutation", "X_cnv", "X_assays", "propensity_X"
  )
  .dot_args <- list(...)
  .dot_args <- .dot_args[!names(.dot_args) %in% .cfomics_meta_args_strip]

  # Unknown-arg typo validation. Detect args not consumed by any cf_compare
  # formal, the strip list above, or the formals of the per-method
  # `cf_fit_<method>()`, so typos such as `base_learner_outcome` (correct name:
  # `model_outcome`) are surfaced rather than silently swallowed.
  .cfomics_validate_dot_args(.dot_args, method)

  # 1. Stratified K-fold split (on T)
  folds <- .create_stratified_folds(T_vec, n_folds = n_folds,
                                     random_state = random_state)

  # 2. K-fold OOF nuisance — reuse cf_compare's shared nuisance via
  #    `precomputed_nuisance` when supplied.
  nuisance_oof <- if (is.null(precomputed_nuisance)) {
    .compute_common_nuisance_oof(
      X, T_vec, Y_vec, folds,
      eval_nuisance = eval_nuisance,
      eval_propensity_score = eval_propensity_score,
      random_state = random_state,
      binary_outcome = binary_outcome
    )
  } else {
    precomputed_nuisance
  }

  # 3. Per-fold method dispatch + predict (X matrix subset)
  oof_ite <- rep(NA_real_, n)
  # Per-method OOF potential outcomes, used by the factual loss so it reflects
  # this method's own outcome predictions (not the shared evaluation nuisance).
  oof_mu0 <- rep(NA_real_, n)
  oof_mu1 <- rep(NA_real_, n)
  fold_cache <- if (isTRUE(keep_kfold_models)) {
    vector("list", length(folds))
  } else NULL

  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(n), test_idx)
    X_train <- X[train_idx, , drop = FALSE]
    T_train <- T_vec[train_idx]
    Y_train <- Y_vec[train_idx]

    # If `propensity_score` is a length-n vector it must be subset to the
    # current training fold; scalars are broadcast by the per-method backends.
    ps_for_fold <- if (is.null(propensity_score)) NULL else
                   if (length(propensity_score) == 1L) propensity_score else
                   propensity_score[train_idx]
    fit_k_inner <- do.call(.dispatch_method_lowlevel, c(list(
      method = method, X = X_train, T = T_train, Y = Y_train,
      cov_names = cov_names,
      n_bootstrap = n_bootstrap,
      bootstrap_ci = bootstrap_ci,
      random_state = random_state,
      device = device,
      propensity_X = propensity_X,
      propensity_score = ps_for_fold,
      trim = trim,
      binary_outcome = binary_outcome
    ), .dot_args))
    # Wrap with class for predict.cf_model dispatch (res = inner res,
    # model = the backend container).
    fit_k <- structure(
      list(method = method,
            res   = fit_k_inner$res,
            model = fit_k_inner,
            meta = list(
              covariate_names = cov_names,
              treatment_name = "T", outcome_name = "Y",
              # Populate each method's meta z-score / factor_levels on the
              # per-fold predict path (same logic as cf_fit Step 7): when the
              # per-fold inner fit returns `res$x_zscore` / `res$y_zscore`, copy
              # them here. Without these stats, predict_cf_<method> for methods
              # like ganite aborts with `newdata prediction requires a fresh fit`.
              x_zscore = if (is.list(fit_k_inner) && !is.null(fit_k_inner$res$x_zscore))
                          fit_k_inner$res$x_zscore else NULL,
              y_zscore = if (is.list(fit_k_inner) && !is.null(fit_k_inner$res$y_zscore))
                          fit_k_inner$res$y_zscore else NULL,
              # propensity_X must reach the per-method `predict_cf_<method>()`
              # so it can subset newdata when predicting W.hat / pi_hat using the
              # subset-trained PS forest / classifier.
              propensity_X = propensity_X
            )),
      class = "cf_model"
    )

    # bcf's predict requires a propensity score for the held-out rows. When a
    # per-observation (vector) known PS was supplied, its held-out subset is
    # exactly that value, so thread it in as `pihat_new` (the persisted
    # "user_supplied" PS only covers the training fold). Other methods (and
    # scalar / estimated PS, which predict handles internally) ignore it.
    .pred_extra <- if (identical(method, "bcf") &&
                       !is.null(propensity_score) &&
                       length(propensity_score) > 1L) {
      list(pihat_new = as.numeric(propensity_score)[test_idx])
    } else {
      list()
    }
    pred_out <- do.call(predict, c(
      list(fit_k, newdata = X[test_idx, , drop = FALSE], type = "all"),
      .pred_extra))

    if (length(pred_out$tau_hat) == length(test_idx) && !anyNA(pred_out$tau_hat)) {
      oof_ite[test_idx] <- pred_out$tau_hat
    } else {
      rlang::abort(
        sprintf(
          paste0("partition='oof' failed: method '%s' returned NA / ",
                 "wrong-length tau at fold %d (got length %d, expected %d)."),
          method, k, length(pred_out$tau_hat), length(test_idx)
        ),
        class = "cfomics_partition_oof_failed"
      )
    }
    # Held-out potential outcomes for this method (used by the factual loss).
    # Left NA when a method does not surface finite mu on newdata.
    if (!is.null(pred_out$mu0_hat) &&
        length(pred_out$mu0_hat) == length(test_idx)) {
      oof_mu0[test_idx] <- pred_out$mu0_hat
    }
    if (!is.null(pred_out$mu1_hat) &&
        length(pred_out$mu1_hat) == length(test_idx)) {
      oof_mu1[test_idx] <- pred_out$mu1_hat
    }
    if (isTRUE(keep_kfold_models)) fold_cache[[k]] <- fit_k
  }

  # 4. Full-data fit (method-specific, calls cf_fit_<method> directly)
  full_fit_inner <- do.call(.dispatch_method_lowlevel, c(list(
    method = method, X = X, T = T_vec, Y = Y_vec,
    cov_names = cov_names,
    n_bootstrap = n_bootstrap,
    bootstrap_ci = bootstrap_ci,
    random_state = random_state,
    device = device,
    propensity_X = propensity_X,
    propensity_score = propensity_score,
    trim = trim,
    binary_outcome = binary_outcome
  ), .dot_args))

  list(
    fit        = full_fit_inner,
    oof_ite    = oof_ite,
    oof_mu0    = oof_mu0,
    oof_mu1    = oof_mu1,
    nuisance   = nuisance_oof,
    fold_cache = fold_cache
  )
}


#' Direct dispatch to per-method low-level `cf_fit_<method>` (no cf_fit recursion)
#'
#' @keywords internal
#' @noRd
.dispatch_method_lowlevel <- function(method, X, T, Y,
                                       cov_names,
                                       n_bootstrap = 500L,
                                       bootstrap_ci = FALSE,
                                       random_state = NULL,
                                       device = "cpu",
                                       propensity_X = NULL,
                                       binary_outcome = FALSE,
                                       ...) {
  # `device` is only accepted by the 4 Python NN methods
  # (cevae/ganite/tarnet/dragonnet). The other methods (grf/bcf/meta-learners)
  # have no `device` argument, so it must not be forwarded via `...` (which
  # would fire an `unused argument` error).
  rs_int <- as.integer(random_state %||% 0L)
  switch(method,
    "grf" = cf_fit_grf(X = X, T = T, Y = Y, covariate_names = cov_names,
                       random_state = rs_int, propensity_X = propensity_X, ...),
    "bcf" = cf_fit_bcf(X = X, T = T, Y = Y, covariate_names = cov_names,
                       random_state = rs_int, propensity_X = propensity_X, ...),
    "slearner" = cf_fit_slearner(X = X, T = T, Y = Y, covariate_names = cov_names,
                                  bootstrap_ci = isTRUE(bootstrap_ci),
                                  n_bootstrap  = as.integer(n_bootstrap),
                                  random_state = rs_int,
                                  binary_outcome = binary_outcome, ...),
    "tlearner" = cf_fit_tlearner(X = X, T = T, Y = Y, covariate_names = cov_names,
                                  bootstrap_ci = isTRUE(bootstrap_ci),
                                  n_bootstrap  = as.integer(n_bootstrap),
                                  random_state = rs_int,
                                  binary_outcome = binary_outcome, ...),
    "xlearner" = cf_fit_xlearner(X = X, T = T, Y = Y, covariate_names = cov_names,
                                  bootstrap_ci = isTRUE(bootstrap_ci),
                                  n_bootstrap  = as.integer(n_bootstrap),
                                  random_state = rs_int, propensity_X = propensity_X,
                                  binary_outcome = binary_outcome, ...),
    "drlearner" = cf_fit_drlearner(X = X, T = T, Y = Y, covariate_names = cov_names,
                                    bootstrap_ci = isTRUE(bootstrap_ci),
                                    n_bootstrap  = as.integer(n_bootstrap),
                                    random_state = rs_int, propensity_X = propensity_X,
                                    binary_outcome = binary_outcome, ...),
    "cevae" = cf_fit_cevae(X = X, T = T, Y = Y, covariate_names = cov_names,
                           random_state = rs_int, device = device,
                           binary_outcome = binary_outcome, ...),
    "ganite" = cf_fit_ganite(X = X, T = T, Y = Y,
                              covariate_names = cov_names, random_state = rs_int,
                              device = device, binary_outcome = binary_outcome, ...),
    "tarnet" = cf_fit_tarnet(X = X, T = T, Y = Y, covariate_names = cov_names,
                              random_state = rs_int, device = device,
                              binary_outcome = binary_outcome, ...),
    "dragonnet" = cf_fit_dragonnet(X = X, T = T, Y = Y, covariate_names = cov_names,
                                    random_state = rs_int, device = device,
                                    binary_outcome = binary_outcome, ...),
    stop(sprintf("Unknown method: %s", method))
  )
}
