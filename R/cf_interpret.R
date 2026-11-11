#' Interpret a fitted CATE model via surrogate distillation SHAP.
#'
#' \code{cf_interpret()} distils the model's out-of-fold tau hat into two
#' surrogate models and reports per-feature signed SHAP contributions plus
#' a held-out fidelity score. SHAP comes from the surrogate fitted on all
#' rows; fidelity is measured on rows each surrogate fit did not see:
#'
#' \describe{
#'   \item{\code{"ridge"}}{\code{glmnet} ridge regression of \eqn{\hat\tau}
#'     on z-scored \eqn{X}. For a linear model the per-observation SHAP is
#'     \eqn{\phi_j(x_i) = \beta_j (x_{ij} - \bar x_j)}; we report the mean
#'     absolute contribution as the importance and \code{sign(beta_j)} as
#'     the sign.}
#'   \item{\code{"xgb"}}{\code{xgboost} gradient boosting of \eqn{\hat\tau}
#'     on \eqn{X} with TreeSHAP via
#'     \code{predict(..., predcontrib = TRUE)}.}
#' }
#'
#' Fidelity is \eqn{R^2 = 1 - \mathrm{RSS} / \mathrm{TSS}} of the surrogate
#' on the same \eqn{\hat\tau} target: how much of the estimated effect the
#' surrogate reproduces.
#'
#' A positive SHAP contribution is one that raises \eqn{\hat\tau}, which on
#' CATEomics' convention that a larger \code{Y} is the better outcome (see
#' \code{Y} in \code{\link{cf_fit}}) means it raises the benefit of
#' treatment. For an outcome whose smaller values are better, fit the negated
#' outcome, or every sign here reads the wrong way round.
#'
#' @section Surrogate distillation:
#' The interpretation uses a surrogate distillation path. The contract is uniform across all
#' CATE-capable methods (including \code{cf_ensemble} output) and
#' avoids reaching into method-specific internals. The entire
#' interpretation surface is R-native (\code{glmnet} + \code{xgboost}).
#'
#' @param model A \code{cf_model} or \code{cateomics_ensemble} object.
#' @param surrogate Character vector; subset of \code{c("ridge", "xgb")}.
#'   Defaults to both. Order is preserved in the output.
#' @param tau Optional numeric vector of length n to use as the target
#'   instead of \code{model$res$tau_hat_oof}. \code{NULL} (default) reads the
#'   model's OOF tau-hat, falling back to the in-sample tau-hat with a
#'   warning.
#' @param top_n Integer (default \code{50}); how many features to keep in
#'   the per-surrogate
#'   \code{features} table (sorted by \code{importance}, the mean absolute
#'   SHAP value).
#' @param params Hyperparameters for the surrogate fits, as a list named by
#'   surrogate: \code{list(ridge = list(nfolds = 10), xgb = list(nrounds = 200,
#'   max_depth = 6))}. The names are the ones the R package that runs the
#'   surrogate takes -- \code{glmnet::cv.glmnet} (and what it forwards to
#'   \code{glmnet::glmnet}) for \code{"ridge"}, \code{xgboost::xgb.train}
#'   (and the booster hyperparameters of \code{xgboost::xgb.params}) for
#'   \code{"xgb"}. An unknown name is an error. \code{alpha} is fixed at 0,
#'   so the linear surrogate is always a ridge; the XGB objective is fixed to
#'   squared error because the surrogate regresses on tau-hat. Defaults:
#'   \code{nfolds = 5} for the ridge, and \code{nrounds = 100},
#'   \code{max_depth = 4}, \code{eta = 0.1},
#'   \code{min_child_weight = 10} and \code{colsample_bytree = 0.5} for the
#'   XGB. The resolved values are returned in
#'   \code{vi$<surrogate>$params}.
#' @param random_state Seed (default \code{0}) for the surrogate fits and the
#'   fold split used by the held-out fidelity.
#' @param ... Not used; naming an unknown argument is an error.
#'
#' @return An object of class \code{c("cateomics_interpret", "list")} with
#'   elements:
#'   \describe{
#'     \item{\code{vi}}{Named list of one entry per requested surrogate, each
#'       with \code{features} (\code{data.frame} of \code{feature_name},
#'       \code{importance} (mean|SHAP|), \code{direction}
#'       (\code{sign(cor(x_j, phi_j))}), \code{rank}), the full
#'       per-observation \code{shap} matrix (n x p), \code{fidelity}
#'       (K-fold held-out R^2 of the surrogate on \eqn{\hat\tau};
#'       \code{NA} when there are too few rows to hold any out),
#'       \code{params} (the resolved surrogate settings) and
#'       \code{surrogate} (label).}
#'     \item{\code{X}}{The design matrix the surrogates were fitted on, held
#'       once for the object rather than once per surrogate. It is the
#'       preprocessed design, and it is what the \pkg{shapviz} coercion
#'       pairs with each \code{shap} matrix.}
#'     \item{\code{metadata}}{Method, partition, surrogate names,
#'       \code{tau_source}, \code{n}, \code{p}.}
#'   }
#'
#' @seealso \code{\link{cf_plot_interpret}} for the delegated \pkg{shapviz}
#'   plots (importance, beeswarm, dependence).
#'
#' @references
#' Lundberg & Lee (2017) A Unified Approach to Interpreting Model Predictions,
#' \emph{NeurIPS} 30 (arXiv:1705.07874).
#'
#' Lundberg et al. (2020) From local explanations to global understanding with
#' explainable AI for trees (TreeSHAP), \emph{Nature Machine Intelligence}
#' 2:56-67.
#'
#' @examples
#' # Surrogate SHAP importances from a fast R-native causal forest.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' itp <- cf_interpret(fit, surrogate = "xgb")
#' head(itp$vi$xgb$features)
#'
#' @export
cf_interpret <- function(model,
                         surrogate     = c("ridge", "xgb"),
                         tau           = NULL,
                         top_n         = 50L,
                         params        = NULL,
                         random_state  = 0L,
                         ...) {
  # Every penalised fit this call makes reports its thin split through one
  # line, named for the learners it saw.
  .cf_cv_thin_open()
  on.exit(.cf_cv_thin_close(), add = TRUE)
  .cf_reject_dots(..., .fn = "cf_interpret()")
  if (!inherits(model, "cf_model")) {
    rlang::abort(
      "cf_interpret requires a cf_model (cf_fit / cf_compare / cf_ensemble).",
      class = "cateomics_invalid_input"
    )
  }
  surrogate <- .cf_match_several(surrogate, c("ridge", "xgb"),
                                 "surrogate", "cf_interpret()")
  top_n <- as.integer(top_n)[1L]
  sur_params <- .cf_interpret_params(params, surrogate)

  # ----- 1. Resolve tau target (OOF first, full fallback with warning).
  tau_source <- "oof"
  if (is.null(tau)) {
    tau_vec <- model$res$tau_hat_oof
    if (is.null(tau_vec) || all(is.na(tau_vec))) {
      tau_vec <- model$res$tau_hat
      tau_source <- "in_sample"
      if (!is.null(tau_vec) && any(is.finite(tau_vec))) {
        rlang::warn(
          c("cf_interpret: OOF tau hat unavailable; falling back to in-sample tau.",
            i = "Re-fit with partition = 'oof' for honest interpretation."),
          class = "cateomics_interpret_in_sample_warning"
        )
      }
    }
  } else {
    tau_vec    <- as.numeric(tau)
    tau_source <- "user_supplied"
  }
  if (is.null(tau_vec) || !any(is.finite(tau_vec))) {
    rlang::abort(
      "cf_interpret: no usable tau hat (OOF and in-sample both NA).",
      class = "cateomics_interpret_no_tau"
    )
  }

  # ----- 2. Resolve X (the preprocessed design the model was fitted on).
  X_raw <- model$data$X_preprocessed %||% model$data$X
  if (is.null(X_raw)) {
    rlang::abort(
      c("cf_interpret: model$data$X_preprocessed is missing.",
        i = "cf_compare caches X on every model; ensure the model came from cf_compare or cf_fit (not a stale RDS without X)."),
      class = "cateomics_interpret_no_X"
    )
  }
  X_mat <- if (is.matrix(X_raw)) X_raw else as.matrix(X_raw)
  if (!is.numeric(X_mat)) storage.mode(X_mat) <- "double"
  if (length(tau_vec) != nrow(X_mat)) {
    rlang::abort(
      sprintf(paste0("cf_interpret: tau length (%d) does not match nrow(X) (%d).",
                     " Did you pass a different model than the one that built tau?"),
              length(tau_vec), nrow(X_mat)),
      class = "cateomics_interpret_size_mismatch"
    )
  }
  feature_names <- colnames(X_mat)
  if (is.null(feature_names) || any(feature_names == "")) {
    feature_names <- paste0("X", seq_len(ncol(X_mat)))
    colnames(X_mat) <- feature_names
  }

  # Restrict to rows with finite tau and complete X.
  keep_rows <- is.finite(tau_vec) & stats::complete.cases(X_mat)
  tau_use   <- tau_vec[keep_rows]
  X_use     <- X_mat[keep_rows, , drop = FALSE]
  n_use     <- nrow(X_use)
  if (n_use < 10L) {
    rlang::abort(
      sprintf("cf_interpret: only %d usable rows after NA filter; need >= 10.",
              n_use),
      class = "cateomics_interpret_too_few_rows"
    )
  }
  p <- ncol(X_use)

  # ----- 3. Per-surrogate fit + SHAP + fidelity.
  # Built in the order the caller asked for: `cf_plot_interpret()` and
  # `as_shapviz()` take the first surrogate present when none is named, so the
  # order decides which panel a plot shows by default.
  vi <- list()
  for (sg in surrogate) {
    if (identical(sg, "ridge")) {
      vi$ridge <- .cateomics_interpret_ridge(X_use, tau_use, feature_names,
                                             top_n, sur_params$ridge,
                                             random_state)
    } else if (identical(sg, "xgb")) {
      # xgboost is in Imports, so a working install always has it; the guard is
      # kept for the same reason the one in `.fit_xgboost()` is, and cannot be
      # reached from a loaded package.
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        rlang::warn(
          c("xgboost is not installed; cf_interpret(surrogate = 'xgb') is skipped.",
            i = "Install with install.packages('xgboost')."),
          class = "cateomics_interpret_xgb_unavailable"
        )
      } else {
        vi$xgb <- .cateomics_interpret_xgb(X_use, tau_use, feature_names,
                                           top_n, sur_params$xgb, random_state)
      }
    }
  }

  structure(
    list(
      vi       = vi,
      # The design matrix the surrogates were fitted on, once: it is the same
      # matrix for every surrogate (and the same one the fit carries).
      X        = X_use,
      metadata = list(
        method        = if (!is.null(model$method)) model$method
                        else NA_character_,
        partition     = if (!is.null(model$meta$estimator$partition))
                          model$meta$estimator$partition else NA_character_,
        surrogates    = names(vi),
        tau_source    = tau_source,
        n             = n_use,
        p             = p
      )
    ),
    class = c("cateomics_interpret", "list")
  )
}


# -------------------------------------------------------------------------
# Ridge surrogate (linear SHAP = beta_j * (x_j - mean(x_j)))
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Held-out fidelity
# -------------------------------------------------------------------------

#' K-fold held-out R2 of a surrogate against tau hat
#'
#' `fit_fn(X_train, tau_train)` must run the whole fitting procedure,
#' including any tuning, and `predict_fn(fit, X_test)` must score rows the
#' fold never saw.
#'
#' @param X,tau Surrogate design and target (complete rows only).
#' @param fit_fn,predict_fn Fold fitter and scorer.
#' @param random_state Integer seed for the fold split and the fits.
#' @return Scalar R2, or `NA_real_` when a fold cannot be fitted.
#' @keywords internal
#' @noRd
.cateomics_surrogate_fidelity <- function(X, tau, fit_fn, predict_fn,
                                        random_state = NULL) {
  n <- nrow(X)
  k <- max(2L, min(5L, n %/% 5L))
  if (n < 10L || k < 2L) return(NA_real_)

  old <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  # Pin the algorithm as `.create_stratified_folds` does, so the fold split is
  # the same inside a BiocParallel worker.
  if (!is.null(random_state))
    set.seed(as.integer(random_state), kind = "Mersenne-Twister")
  fold_id <- sample(rep_len(seq_len(k), n))

  # One seeded stream for the whole loop, drawn from in turn.
  pred <- rep(NA_real_, n)
  for (f in seq_len(k)) {
    te <- which(fold_id == f); tr <- setdiff(seq_len(n), te)
    if (length(tr) < 5L) return(NA_real_)
    fit <- tryCatch(fit_fn(X[tr, , drop = FALSE], tau[tr]),
                    error = function(e) NULL)
    if (is.null(fit)) {
      rlang::warn(
        "No held-out fidelity: a surrogate fold could not be fitted.",
        class = "cateomics_fidelity_fold_failed")
      return(NA_real_)
    }
    pred[te] <- tryCatch(predict_fn(fit, X[te, , drop = FALSE]),
                         error = function(e) NA_real_)
  }
  if (anyNA(pred)) {
    rlang::warn(
      sprintf("Held-out fidelity computed on %d of %d rows; the rest could not be scored.",
              sum(!is.na(pred)), length(pred)),
      class = "cateomics_fidelity_partial")
  }
  .cateomics_r_squared(tau, pred)
}


#' Resolve and check the surrogate hyperparameters
#'
#' `params` is named by surrogate.
#' The accepted names are derived from the R function that runs the surrogate,
#' so they are whatever that version of glmnet / xgboost accepts; the shared
#' `.cf_learner_params_check()` already excludes the arguments CATEomics fixes
#' itself (`alpha` for the ridge, `objective` for the XGB, and the data
#' arguments).
#'
#' @param params The user's `params`, or NULL.
#' @param surrogate Character vector of requested surrogates.
#' @return A list with a `ridge` and an `xgb` entry, each a named list.
#' @keywords internal
#' @noRd
.cf_interpret_params <- function(params, surrogate) {
  out <- list(ridge = list(), xgb = list())
  if (is.null(params)) return(out)
  if (!is.list(params) || is.null(names(params)) || any(!nzchar(names(params)))) {
    rlang::abort(
      c("`params` must be a named list, named by surrogate.",
        i = 'For example params = list(xgb = list(nrounds = 200)).'),
      class = "cateomics_invalid_params"
    )
  }
  bad <- setdiff(names(params), c("ridge", "xgb"))
  if (length(bad)) {
    rlang::abort(
      c(sprintf("`params` has no surrogate %s.", paste(bad, collapse = " / ")),
        i = "It accepts: ridge, xgb."),
      class = "cateomics_unknown_params_slot"
    )
  }
  unused <- setdiff(names(params), surrogate)
  if (length(unused)) {
    rlang::abort(
      c(sprintf("`params$%s` is set but that surrogate was not requested.",
                paste(unused, collapse = " / ")),
        i = sprintf("`surrogate` is %s.", paste(surrogate, collapse = ", "))),
      class = "cateomics_unused_params"
    )
  }
  for (sg in names(params)) {
    pl <- params[[sg]]
    if (!is.list(pl) || (length(pl) && is.null(names(pl)))) {
      rlang::abort(
        sprintf("`params$%s` must be a named list.", sg),
        class = "cateomics_invalid_params"
      )
    }
    .cf_learner_params_check(if (identical(sg, "xgb")) "xgboost" else sg,
                             pl, sprintf("params$%s", sg))
    out[[sg]] <- pl
  }
  out
}


#' @keywords internal
#' @noRd
.cateomics_interpret_ridge <- function(X, tau, feature_names, top_n,
                                     user_params = list(),
                                     random_state = NULL) {
  # alpha stays 0: the linear surrogate is a ridge by design, and the signed
  # linear SHAP read off it assumes every feature keeps a coefficient.
  gl_args <- utils::modifyList(list(nfolds = .CF_CV_FOLDS), user_params %||% list())
  # Z-score each feature; SHAP becomes beta_j * standardised x_j, the
  # magnitude scale is "1 SD shift in x_j -> tau hat change".
  center <- colMeans(X)
  sd_vec <- apply(X, 2L, stats::sd)
  sd_vec[!is.finite(sd_vec) | sd_vec == 0] <- 1
  Xz <- sweep(sweep(X, 2L, center, "-"), 2L, sd_vec, "/")

  # cv.glmnet draws random fold ids, advancing the caller's
  # global RNG as a side effect. Save/restore .Random.seed (and seed for
  # reproducibility when random_state is supplied).
  old <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(random_state)) set.seed(as.integer(random_state))

  cv_fit <- do.call(glmnet::cv.glmnet,
                    c(list(x = Xz, y = tau, alpha = 0), gl_args))
  cs <- as.numeric(stats::coef(cv_fit, s = "lambda.min"))
  beta <- cs[-1L]
  names(beta) <- feature_names

  # Held-out fidelity: the z-score statistics and the lambda search are redone
  # inside each fold, so nothing about the scored rows enters the fit.
  fidelity <- .cateomics_surrogate_fidelity(
    X, tau,
    fit_fn = function(X_tr, tau_tr) {
      ctr <- colMeans(X_tr)
      sdv <- apply(X_tr, 2L, stats::sd)
      sdv[!is.finite(sdv) | sdv == 0] <- 1
      Xtr_z <- sweep(sweep(X_tr, 2L, ctr, "-"), 2L, sdv, "/")
      fold_args <- gl_args
      fold_args$nfolds <- min(as.integer(gl_args$nfolds),
                              max(3L, nrow(X_tr) %/% 3L))
      list(fit = do.call(glmnet::cv.glmnet,
                         c(list(x = Xtr_z, y = tau_tr, alpha = 0), fold_args)),
           center = ctr, sd = sdv)
    },
    predict_fn = function(f, X_te) {
      X_te_z <- sweep(sweep(X_te, 2L, f$center, "-"), 2L, f$sd, "/")
      as.numeric(stats::predict(f$fit, newx = X_te_z, s = "lambda.min"))
    },
    random_state = random_state
  )

  # retain the full per-observation linear SHAP matrix
  # phi_{ij} = beta_j * Xz_{ij} (Xz has mean 0) so it can be coerced to a
  # shapviz object. importance = mean|phi_j|, direction = sign(cor(X_j, phi_j)).
  shap_mat <- sweep(Xz, 2L, beta, "*")
  colnames(shap_mat) <- feature_names
  features <- .cateomics_vi_summary(feature_names, shap_mat, X, top_n)

  # `coef` is not carried: the SHAP matrix is beta_j (x_ij - mean_j) / sd_j, so
  # the ratio of a SHAP column to its standardised feature is that coefficient
  # for every row, and the sign is already in `features$direction`. The
  # intercept is mean(tau) just as exactly. The design matrix is held once, by
  # the object.
  list(
    features  = features,
    fidelity  = fidelity,
    surrogate = "ridge_linear_shap",
    params    = gl_args,
    shap      = shap_mat
  )
}


# -------------------------------------------------------------------------
# XGB surrogate (TreeSHAP)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cateomics_interpret_xgb <- function(X, tau, feature_names,
                                    top_n, user_params = list(),
                                    random_state = NULL) {
  # objective stays squared error: the surrogate regresses on tau-hat, which
  # is continuous whatever the outcome type was. Arguments xgb.train takes
  # itself go at the top level, the rest are booster hyperparameters.
  # max_depth, eta and min_child_weight follow SuperLearner::SL.xgboost
  # (4 / 0.1 / 10)., which a fixed max_depth does not, and the column
  # subsample bounds how many candidate features each split scans: with p
  # above n, the best training gain among all of them is large even when none
  # carries signal. ranger, the other tree learner here, subsamples columns by
  # default (mtry); xgboost does not.
  up       <- utils::modifyList(list(nrounds          = 100L,
                                     max_depth        = 4L,
                                     eta              = 0.1,
                                     min_child_weight = 10L,
                                     colsample_bytree = 0.5),
                                user_params %||% list())
  top_lvl  <- intersect(names(up), names(formals(xgboost::xgb.train)))
  # The fixed top-level arguments go through modifyList rather than being
  # passed alongside: a caller's `verbose` would otherwise reach do.call twice
  # and R refuses a formal matched by two actual arguments.
  train_as <- utils::modifyList(list(verbose = 0L), up[top_lvl])
  booster  <- utils::modifyList(
    list(objective = "reg:squarederror", verbosity = 0L),
    up[setdiff(names(up), top_lvl)])
  if (!is.null(random_state)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old)) {
        if (exists(".Random.seed", envir = .GlobalEnv))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(random_state))
  }

  dtr <- xgboost::xgb.DMatrix(data = X, label = tau)
  bst <- do.call(xgboost::xgb.train,
                 c(list(params = booster, data = dtr), train_as))

  shap <- stats::predict(bst, newdata = X, predcontrib = TRUE)
  # The last column is the BIAS; drop it for the per-feature SHAP.
  if (ncol(shap) == length(feature_names) + 1L) {
    shap <- shap[, seq_along(feature_names), drop = FALSE]
  }
  colnames(shap) <- feature_names

  fidelity <- .cateomics_surrogate_fidelity(
    X, tau,
    fit_fn = function(X_tr, tau_tr) {
      do.call(xgboost::xgb.train,
              c(list(params = booster,
                     data = xgboost::xgb.DMatrix(data = X_tr, label = tau_tr)),
                train_as))
    },
    predict_fn = function(f, X_te) as.numeric(stats::predict(f, newdata = X_te)),
    random_state = random_state
  )

  # Retain the full per-observation TreeSHAP matrix + feature values.
  features <- .cateomics_vi_summary(feature_names, shap, X, top_n)

  list(
    features  = features,
    fidelity  = fidelity,
    surrogate = "xgb_tree_shap",
    params    = up,
    shap      = shap
  )
}


# -------------------------------------------------------------------------
# Shared helpers
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cateomics_r_squared <- function(y, y_hat) {
  ok <- is.finite(y) & is.finite(y_hat)
  if (sum(ok) < 2L) return(NA_real_)
  rss <- sum((y[ok] - y_hat[ok])^2)
  tss <- sum((y[ok] - mean(y[ok]))^2)
  if (tss <= 0) return(NA_real_)
  1 - rss / tss
}


#' Build the per-surrogate VI summary table from the full SHAP matrix.
#'
#' \code{importance} is the unsigned mean absolute SHAP
#' (\eqn{\mathrm{mean}_i|\phi_{ij}|}); \code{direction} is
#' \eqn{\mathrm{sign}(\mathrm{cor}(x_j, \phi_j))} (0 when undefined). The
#' importance/direction pair carries the per-feature summary, and the full
#' per-observation matrix is retained separately for \code{as_shapviz()}.
#'
#' @keywords internal
#' @noRd
.cateomics_vi_summary <- function(feature_names, shap_matrix, X, top_n) {
  importance <- colMeans(abs(shap_matrix))
  direction <- vapply(seq_along(feature_names), function(j) {
    xj <- X[, j]
    pj <- shap_matrix[, j]
    if (stats::sd(xj) == 0 || stats::sd(pj) == 0) return(0)
    s <- suppressWarnings(stats::cor(xj, pj))
    if (!is.finite(s)) 0 else sign(s)
  }, numeric(1))
  ord <- order(importance, decreasing = TRUE)
  if (length(ord) > top_n) ord <- ord[seq_len(top_n)]
  data.frame(
    feature_name = feature_names[ord],
    importance   = unname(importance[ord]),
    direction    = unname(direction[ord]),
    rank         = seq_along(ord),
    stringsAsFactors = FALSE
  )
}


#' Coerce a \code{cf_interpret} result to a \pkg{shapviz} object.
#'
#' Builds a \code{shapviz} object from the per-observation SHAP matrix
#' and feature values retained by \code{\link{cf_interpret}} for one surrogate,
#' exposing the full \pkg{shapviz} plotting surface (\code{sv_importance},
#' \code{sv_dependence}, beeswarm, ...). Requires the \pkg{shapviz} package.
#'
#' @param object A \code{cateomics_interpret} object from \code{\link{cf_interpret}}.
#' @param surrogate Which surrogate to coerce (\code{"ridge"} / \code{"xgb"});
#'   \code{NULL} (default) picks the first available.
#' @param ... Passed to \code{shapviz::shapviz}.
#' @return A \code{shapviz} object.
#' @keywords internal
#' @noRd
as_shapviz <- function(object, surrogate = NULL, ...) {
  if (!inherits(object, "cateomics_interpret")) {
    rlang::abort(
      "as_shapviz requires a cateomics_interpret object (from cf_interpret).",
      class = "cateomics_invalid_input")
  }
  if (!requireNamespace("shapviz", quietly = TRUE)) {
    rlang::abort(
      c("Package 'shapviz' is required for cf_plot_interpret().",
        i = "Install with install.packages('shapviz')."),
      class = "cateomics_optional_dependency_missing")
  }
  s  <- .cateomics_resolve_surrogate(object, surrogate)
  vi <- object$vi[[s]]
  if (is.null(vi$shap) || is.null(object$X)) {
    rlang::abort(
      "This cf_interpret result has no retained SHAP matrix (stale object?).",
      class = "cateomics_interpret_no_shap")
  }
  shapviz::shapviz(object = vi$shap, X = as.data.frame(object$X), ...)
}


#' @keywords internal
#' @noRd
.cateomics_resolve_surrogate <- function(object, surrogate = NULL) {
  avail <- names(object$vi)
  if (length(avail) == 0L) {
    rlang::abort("cf_interpret result has no surrogates.",
                 class = "cateomics_interpret_empty")
  }
  if (is.null(surrogate)) return(avail[1L])
  surrogate <- match.arg(surrogate, choices = c("ridge", "xgb"))
  if (!(surrogate %in% avail)) {
    rlang::abort(
      sprintf("surrogate '%s' not in this result (available: %s).",
              surrogate, paste(avail, collapse = ", ")),
      class = "cateomics_interpret_bad_surrogate")
  }
  surrogate
}


#' @export
print.cateomics_interpret <- function(x, top_n = 10L, ...) {
  meta <- x$metadata
  method_label <- if (is.null(meta$method) || is.na(meta$method))
    "<unknown>" else meta$method
  cat("<cateomics_interpret>: method = ", method_label,
      ", n = ", meta$n, ", p = ", meta$p,
      ", tau_source = ", meta$tau_source, "\n", sep = "")
  # The surrogates actually in the object, not the index of them: the two are
  # the same thing recorded twice, and reading the index means a print can show
  # less than the object holds.
  for (s in names(x$vi)) {
    cat("  surrogate = ", s,
        " (fidelity R^2 = ", format(x$vi[[s]]$fidelity, digits = 3),
        ", entries = ", nrow(x$vi[[s]]$features), ")\n", sep = "")
    print(utils::head(x$vi[[s]]$features, top_n), row.names = FALSE)
  }
  invisible(x)
}


#' @export
summary.cateomics_interpret <- function(object, ...) {
  res <- list(
    method     = object$metadata$method,
    n          = object$metadata$n,
    p          = object$metadata$p,
    surrogates = object$metadata$surrogates,
    fidelity   = vapply(object$vi, function(s) s$fidelity, numeric(1)),
    top_n      = vapply(object$vi, function(s) nrow(s$features), integer(1))
  )
  class(res) <- "summary.cateomics_interpret"
  res
}


#' @export
print.summary.cateomics_interpret <- function(x, ...) {
  cat("cf_interpret summary\n")
  method_label <- if (is.null(x$method) || is.na(x$method))
    "<unknown>" else x$method
  cat("  method: ", method_label, "\n", sep = "")
  cat("  n: ", x$n, ", p: ", x$p, "\n", sep = "")
  for (s in x$surrogates) {
    cat("  ", s, " fidelity (R^2): ", format(x$fidelity[s], digits = 3),
        ", entries: ", x$top_n[s], "\n", sep = "")
  }
  invisible(x)
}
