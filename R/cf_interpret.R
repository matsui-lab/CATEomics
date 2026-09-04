#' Interpret a fitted CATE model via surrogate distillation SHAP.
#'
#' \code{cf_interpret()} distils the model's out-of-fold tau hat into two
#' surrogate models and reports per-feature signed SHAP contributions plus
#' a fidelity score:
#'
#' \describe{
#'   \item{\code{"ridge"}}{\code{glmnet} ridge regression of \eqn{\hat\tau}
#'     on z-scored \eqn{X}. For a linear model the per-observation SHAP is
#'     \eqn{\phi_j(x_i) = \beta_j (x_{ij} - \bar x_j)}; we report the mean
#'     absolute contribution as the importance and \code{sign(beta_j)} as
#'     the sign. Robust to \eqn{p \gg n}.}
#'   \item{\code{"xgb"}}{\code{xgboost} gradient boosting of \eqn{\hat\tau}
#'     on \eqn{X} with TreeSHAP via
#'     \code{predict(..., predcontrib = TRUE)}. Picks up non-linearities
#'     and interactions.}
#' }
#'
#' Fidelity is \eqn{R^2 = 1 - \mathrm{RSS} / \mathrm{TSS}} of the surrogate
#' on the same \eqn{\hat\tau} target. Low fidelity (< 0.3 by default)
#' triggers an informational message.
#'
#' @section Surrogate distillation:
#' The interpretation uses a surrogate distillation path rather than
#' method-native variable importance. The contract is uniform across all
#' CATE-capable methods (including \code{cf_ensemble} output) and
#' avoids reaching into method-specific internals. The entire
#' interpretation surface is R-native (\code{glmnet} + \code{xgboost}).
#'
#' @param model A \code{cf_model} or \code{cfomics_ensemble} object.
#' @param surrogate Character vector; subset of \code{c("ridge", "xgb")}.
#'   Defaults to both. Order is preserved in the output.
#' @param tau Optional numeric vector of length n to use as the target
#'   instead of \code{model$res$tau_hat_oof}. The intended use is debugging /
#'   reproducibility; production users should leave it \code{NULL}.
#' @param top_n Integer; how many features to keep in the per-surrogate
#'   \code{features} table (sorted by signed importance magnitude).
#' @param min_fidelity Numeric; threshold below which the function emits
#'   a low-fidelity informational message.
#' @param xgb_nrounds Integer; number of boosting rounds for the XGB
#'   surrogate (default \code{100}).
#' @param xgb_max_depth Integer; max tree depth for the XGB surrogate.
#' @param random_state Optional seed for the XGB surrogate.
#' @param ... Currently unused.
#'
#' @return An object of class \code{c("cfomics_interpret", "list")} with
#'   elements:
#'   \describe{
#'     \item{\code{vi}}{Named list of one entry per requested surrogate, each
#'       with \code{features} (\code{data.frame} of \code{feature_name},
#'       \code{importance} (mean|SHAP|), \code{direction}
#'       (\code{sign(cor(x_j, phi_j))}), \code{rank}), the full
#'       per-observation \code{shap} matrix (n x p) and feature-value matrix
#'       \code{X} for the internal shapviz coercion, \code{fidelity}
#'       (R^2 on \eqn{\hat\tau}), and \code{surrogate} (label).}
#'     \item{\code{metadata}}{Method, partition, surrogate names,
#'       \code{tau_source}, \code{n}, \code{p}, fidelity warnings.}
#'   }
#'
#' @seealso \code{\link{cf_plot_interpret}} for the delegated \pkg{shapviz}
#'   plots (importance, beeswarm, dependence).
#'
#' @references
#' Lundberg & Lee (2017) A Unified Approach to Interpreting Model Predictions,
#' \emph{NeurIPS} 30 (arXiv:1705.07874).
#' Lundberg et al. (2020) From local explanations to global understanding with
#' explainable AI for trees (TreeSHAP), \emph{Nature Machine Intelligence}
#' 2:56–67.
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
                         min_fidelity  = 0.3,
                         xgb_nrounds   = 100L,
                         xgb_max_depth = 4L,
                         random_state  = NULL,
                         ...) {
  if (!inherits(model, "cf_model")) {
    rlang::abort(
      "cf_interpret requires a cf_model (cf_fit / cf_compare / cf_ensemble).",
      class = "cfomics_invalid_input"
    )
  }
  surrogate <- match.arg(surrogate,
                         choices = c("ridge", "xgb"),
                         several.ok = TRUE)
  top_n         <- as.integer(top_n)[1L]
  min_fidelity  <- as.numeric(min_fidelity)[1L]
  xgb_nrounds   <- as.integer(xgb_nrounds)[1L]
  xgb_max_depth <- as.integer(xgb_max_depth)[1L]

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
          class = "cfomics_interpret_in_sample_warning"
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
      class = "cfomics_interpret_no_tau"
    )
  }

  # ----- 2. Resolve X (model$data$X is canonical).
  X_raw <- model$data$X
  if (is.null(X_raw)) {
    rlang::abort(
      c("cf_interpret: model$data$X is missing.",
        i = "cf_compare caches X on every model; ensure the model came from cf_compare or cf_fit (not a stale RDS without X)."),
      class = "cfomics_interpret_no_X"
    )
  }
  X_mat <- if (is.matrix(X_raw)) X_raw else as.matrix(X_raw)
  if (!is.numeric(X_mat)) storage.mode(X_mat) <- "double"
  if (length(tau_vec) != nrow(X_mat)) {
    rlang::abort(
      sprintf(paste0("cf_interpret: tau length (%d) does not match nrow(X) (%d).",
                     " Did you pass a different model than the one that built tau?"),
              length(tau_vec), nrow(X_mat)),
      class = "cfomics_interpret_size_mismatch"
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
      class = "cfomics_interpret_too_few_rows"
    )
  }
  p <- ncol(X_use)

  # ----- 3. Per-surrogate fit + SHAP + fidelity.
  vi <- list()
  low_fidelity <- character(0)

  if ("ridge" %in% surrogate) {
    rid <- .cfomics_interpret_ridge(X_use, tau_use, feature_names, top_n,
                                    random_state)
    vi$ridge <- rid
    if (is.finite(rid$fidelity) && rid$fidelity < min_fidelity) {
      low_fidelity <- c(low_fidelity, "ridge")
    }
  }

  if ("xgb" %in% surrogate) {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      rlang::warn(
        c("xgboost is not installed; cf_interpret(surrogate = 'xgb') is skipped.",
          i = "Install with install.packages('xgboost')."),
        class = "cfomics_interpret_xgb_unavailable"
      )
    } else {
      xgb_out <- .cfomics_interpret_xgb(X_use, tau_use, feature_names,
                                        top_n, xgb_nrounds, xgb_max_depth,
                                        random_state)
      vi$xgb <- xgb_out
      if (is.finite(xgb_out$fidelity) && xgb_out$fidelity < min_fidelity) {
        low_fidelity <- c(low_fidelity, "xgb")
      }
    }
  }

  if (length(low_fidelity) > 0L) {
    rlang::inform(
      sprintf(paste0("cf_interpret: low surrogate fidelity (R^2 < %g) for %s; ",
                     "treat the SHAP signs / magnitudes with care."),
              min_fidelity, paste(low_fidelity, collapse = ", ")),
      class = "cfomics_interpret_low_fidelity"
    )
  }

  structure(
    list(
      vi       = vi,
      metadata = list(
        method        = if (!is.null(model$method)) model$method
                        else NA_character_,
        partition     = if (!is.null(model$meta$partition))
                          model$meta$partition else NA_character_,
        surrogates    = names(vi),
        tau_source    = tau_source,
        n             = n_use,
        p             = p,
        low_fidelity  = low_fidelity,
        min_fidelity  = min_fidelity,
        feature_names = feature_names
      )
    ),
    class = c("cfomics_interpret", "list")
  )
}


# -------------------------------------------------------------------------
# Ridge surrogate (linear SHAP = beta_j * (x_j - mean(x_j)))
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_interpret_ridge <- function(X, tau, feature_names, top_n,
                                     random_state = NULL) {
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

  cv_fit <- glmnet::cv.glmnet(Xz, tau, alpha = 0, nfolds = 5L)
  cs <- as.numeric(stats::coef(cv_fit, s = "lambda.min"))
  intercept <- cs[1L]
  beta      <- cs[-1L]
  names(beta) <- feature_names

  pred  <- as.numeric(stats::predict(cv_fit, newx = Xz, s = "lambda.min"))
  fidelity <- .cfomics_r_squared(tau, pred)

  # retain the full per-observation linear SHAP matrix
  # phi_{ij} = beta_j * Xz_{ij} (Xz has mean 0) so it can be coerced to a
  # shapviz object. importance = mean|phi_j|, direction = sign(cor(X_j, phi_j)).
  shap_mat <- sweep(Xz, 2L, beta, "*")
  colnames(shap_mat) <- feature_names
  features <- .cfomics_vi_summary(feature_names, shap_mat, X, top_n)

  list(
    features  = features,
    fidelity  = fidelity,
    surrogate = "ridge_linear_shap",
    shap      = shap_mat,
    X         = X,
    intercept = intercept,
    coef      = beta
  )
}


# -------------------------------------------------------------------------
# XGB surrogate (TreeSHAP)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_interpret_xgb <- function(X, tau, feature_names,
                                    top_n, nrounds, max_depth,
                                    random_state) {
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
  bst <- xgboost::xgb.train(
    params  = list(objective = "reg:squarederror",
                    eta = 0.1, max_depth = max_depth,
                    verbosity = 0L),
    data    = dtr,
    nrounds = nrounds,
    verbose = 0L
  )

  shap <- stats::predict(bst, newdata = X, predcontrib = TRUE)
  # The last column is the BIAS; drop it for the per-feature SHAP.
  if (ncol(shap) == length(feature_names) + 1L) {
    shap <- shap[, seq_along(feature_names), drop = FALSE]
  }
  colnames(shap) <- feature_names

  pred     <- stats::predict(bst, newdata = X)
  fidelity <- .cfomics_r_squared(tau, as.numeric(pred))

  # Retain the full per-observation TreeSHAP matrix + feature values.
  features <- .cfomics_vi_summary(feature_names, shap, X, top_n)

  list(
    features  = features,
    fidelity  = fidelity,
    surrogate = "xgb_tree_shap",
    shap      = shap,
    X         = X
  )
}


# -------------------------------------------------------------------------
# Shared helpers
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_r_squared <- function(y, y_hat) {
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
.cfomics_vi_summary <- function(feature_names, shap_matrix, X, top_n) {
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
#' @param object A \code{cfomics_interpret} object from \code{\link{cf_interpret}}.
#' @param surrogate Which surrogate to coerce (\code{"ridge"} / \code{"xgb"});
#'   \code{NULL} (default) picks the first available.
#' @param ... Passed to \code{shapviz::shapviz}.
#' @return A \code{shapviz} object.
#' @keywords internal
#' @noRd
as_shapviz <- function(object, surrogate = NULL, ...) {
  if (!inherits(object, "cfomics_interpret")) {
    rlang::abort(
      "as_shapviz requires a cfomics_interpret object (from cf_interpret).",
      class = "cfomics_invalid_input")
  }
  if (!requireNamespace("shapviz", quietly = TRUE)) {
    rlang::abort(
      c("Package 'shapviz' is required for cf_plot_interpret().",
        i = "Install with install.packages('shapviz')."),
      class = "cfomics_missing_dependency")
  }
  s  <- .cfomics_resolve_surrogate(object, surrogate)
  vi <- object$vi[[s]]
  if (is.null(vi$shap) || is.null(vi$X)) {
    rlang::abort(
      "This cf_interpret result has no retained SHAP matrix (stale object?).",
      class = "cfomics_interpret_no_shap")
  }
  shapviz::shapviz(object = vi$shap, X = as.data.frame(vi$X), ...)
}


#' @keywords internal
#' @noRd
.cfomics_resolve_surrogate <- function(object, surrogate = NULL) {
  avail <- names(object$vi)
  if (length(avail) == 0L) {
    rlang::abort("cf_interpret result has no surrogates.",
                 class = "cfomics_interpret_empty")
  }
  if (is.null(surrogate)) return(avail[1L])
  surrogate <- match.arg(surrogate, choices = c("ridge", "xgb"))
  if (!(surrogate %in% avail)) {
    rlang::abort(
      sprintf("surrogate '%s' not in this result (available: %s).",
              surrogate, paste(avail, collapse = ", ")),
      class = "cfomics_interpret_bad_surrogate")
  }
  surrogate
}


#' @export
print.cfomics_interpret <- function(x, top_n = 10L, ...) {
  meta <- x$metadata
  method_label <- if (is.null(meta$method) || is.na(meta$method))
    "<unknown>" else meta$method
  cat("<cfomics_interpret>: method = ", method_label,
      ", n = ", meta$n, ", p = ", meta$p,
      ", tau_source = ", meta$tau_source, "\n", sep = "")
  for (s in meta$surrogates) {
    cat("  surrogate = ", s,
        " (fidelity R^2 = ", format(x$vi[[s]]$fidelity, digits = 3),
        ", entries = ", nrow(x$vi[[s]]$features), ")\n", sep = "")
    print(utils::head(x$vi[[s]]$features, top_n), row.names = FALSE)
  }
  if (length(meta$low_fidelity) > 0L) {
    cat("  low fidelity (< ", meta$min_fidelity, "): ",
        paste(meta$low_fidelity, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}


#' @export
summary.cfomics_interpret <- function(object, ...) {
  res <- list(
    method     = object$metadata$method,
    n          = object$metadata$n,
    p          = object$metadata$p,
    surrogates = object$metadata$surrogates,
    fidelity   = vapply(object$vi, function(s) s$fidelity, numeric(1)),
    top_n      = vapply(object$vi, function(s) nrow(s$features), integer(1)),
    low_fidelity = object$metadata$low_fidelity
  )
  class(res) <- "summary.cfomics_interpret"
  res
}


#' @export
print.summary.cfomics_interpret <- function(x, ...) {
  cat("cf_interpret summary\n")
  method_label <- if (is.null(x$method) || is.na(x$method))
    "<unknown>" else x$method
  cat("  method: ", method_label, "\n", sep = "")
  cat("  n: ", x$n, ", p: ", x$p, "\n", sep = "")
  for (s in x$surrogates) {
    cat("  ", s, " fidelity (R^2): ", format(x$fidelity[s], digits = 3),
        ", entries: ", x$top_n[s], "\n", sep = "")
  }
  if (length(x$low_fidelity) > 0L) {
    cat("  low fidelity surrogates: ",
        paste(x$low_fidelity, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
