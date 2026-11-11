# ---------------------------------------------------------------------------
# Super Learner for the evaluation nuisance
#
# A convex combination of base learners whose weights are solved on
# out-of-fold predictions. The members are fitted with the same
# `.fit_base_learner()` the single-learner path uses, so the two differ only in
# that one learner becomes a weighted combination.
# ---------------------------------------------------------------------------

#' The `eval_nuisance` token that selects the Super Learner
#' @keywords internal
#' @noRd
.CF_SL_TOKEN <- "super_learner"

#' Learner names `eval_nuisance` accepts
#'
#' The single learners of `.fit_base_learner()` plus the Super Learner token.
#' `cf_fit()` validates its own argument against this same vector, so the two
#' entry points cannot drift apart.
#'
#' @keywords internal
#' @noRd
.CF_EVAL_LEARNERS <- c(.CF_SL_TOKEN, "linear", "logistic", "lasso",
                       "ridge", "elastic_net", "rf", "xgboost")

#' Library members
#'
#' Each entry maps a member label to the base learner it instantiates, the
#' settings the label itself fixes, and the names of those settings. A label
#' that encodes a setting (`xgboost_d2`) cannot also take it through
#' `eval_params`: the same knob would then have two names.
#'
#' `mean` is handled inside this file rather than added to
#' `.fit_base_learner()`'s pool, because an intercept-only learner is only
#' useful as a library member.
#'
#' @keywords internal
#' @noRd
.CF_SL_MEMBERS <- list(
  mean        = list(learner = "mean",        params = list(),               fixes = character(0)),
  ridge       = list(learner = "ridge",       params = list(),               fixes = character(0)),
  elastic_net = list(learner = "elastic_net", params = list(),               fixes = character(0)),
  lasso       = list(learner = "lasso",       params = list(),               fixes = character(0)),
  rf          = list(learner = "rf",          params = list(),               fixes = character(0)),
  xgboost_d2  = list(learner = "xgboost",     params = list(max_depth = 2L), fixes = "max_depth"),
  xgboost_d4  = list(learner = "xgboost",     params = list(max_depth = 4L), fixes = "max_depth")
)

#' Default library
#'
#' Spans the capacity axis that matters when p exceeds n: an intercept, the
#' three glmnet penalties, a forest, and two boosted fits at depth 2 and 4.
#' The remaining xgboost settings come from `.fit_xgboost()`'s defaults. Three
#' of them follow `SuperLearner::SL.xgboost`; the round count and the column
#' subsample are this package's own.
#'
#' @keywords internal
#' @noRd
.CF_SL_LIBRARY_DEFAULT <- c("mean", "ridge", "elastic_net", "lasso", "rf",
                            "xgboost_d2", "xgboost_d4")

#' Is this `eval_nuisance` value the Super Learner?
#' @keywords internal
#' @noRd
.cf_sl_is_token <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && identical(x, .CF_SL_TOKEN)
}

#' Resolve `eval_library` for one slot
#'
#' Takes the same shapes as `eval_nuisance` and `eval_params`: a character
#' vector applying to both slots, or `list(outcome =, propensity =)` with the
#' missing slot falling back to its partner.
#'
#' @param eval_library NULL, a character vector of member labels, or a named
#'   list with `outcome` / `propensity` entries.
#' @param slot `"outcome"` or `"propensity"`.
#' @return Character vector of member labels.
#' @keywords internal
#' @noRd
.cf_sl_resolve_library <- function(eval_library, slot) {
  if (is.null(eval_library)) return(.CF_SL_LIBRARY_DEFAULT)
  if (is.list(eval_library)) {
    keyed <- intersect(names(eval_library), c("outcome", "propensity"))
    if (!length(keyed)) {
      rlang::abort(
        c("`eval_library` must be a character vector or a list with `outcome` / `propensity`.",
          i = "Received a list with no slot names."),
        class = "cateomics_invalid_eval_library"
      )
    }
    if (length(setdiff(names(eval_library), keyed))) {
      rlang::abort(
        c("`eval_library` mixes slot entries with other names.",
          i = "Use either a character vector or list(outcome = , propensity = )."),
        class = "cateomics_invalid_eval_library"
      )
    }
    lib <- eval_library[[slot]] %||% eval_library[[setdiff(c("outcome", "propensity"), slot)]]
    return(.cf_sl_check_library(lib))
  }
  .cf_sl_check_library(eval_library)
}

#' @keywords internal
#' @noRd
.cf_sl_check_library <- function(lib) {
  if (is.null(lib)) return(.CF_SL_LIBRARY_DEFAULT)
  if (!is.character(lib) || !length(lib) || anyNA(lib)) {
    rlang::abort(
      sprintf("`eval_library` must be a non-empty character vector (got: %s).",
              paste(class(lib), collapse = "/")),
      class = "cateomics_invalid_eval_library"
    )
  }
  if (anyDuplicated(lib)) {
    rlang::abort(
      sprintf("`eval_library` repeats a member: %s.",
              paste(unique(lib[duplicated(lib)]), collapse = ", ")),
      class = "cateomics_invalid_eval_library"
    )
  }
  unknown <- setdiff(lib, names(.CF_SL_MEMBERS))
  if (length(unknown)) {
    rlang::abort(
      c(sprintf("Unknown `eval_library` member(s): %s.", paste(unknown, collapse = ", ")),
        i = sprintf("Available: %s.", paste(names(.CF_SL_MEMBERS), collapse = ", "))),
      class = "cateomics_invalid_eval_library"
    )
  }
  lib
}

# ---------------------------------------------------------------------------
# Member fit / predict
# ---------------------------------------------------------------------------

#' Fit one library member
#'
#' @return A `"cf_base_learner"`, or a `"cf_sl_mean"` object for the
#'   intercept-only member.
#' @keywords internal
#' @noRd
.cf_sl_fit_member <- function(x, y, type, member, random_state = NULL,
                              params = list()) {
  spec <- .CF_SL_MEMBERS[[member]]
  if (identical(spec$learner, "mean")) {
    return(structure(list(value = mean(as.numeric(y))), class = "cf_sl_mean"))
  }
  .fit_base_learner(
    x, y, type = type, learner = spec$learner,
    random_state = random_state,
    .params = utils::modifyList(spec$params, params %||% list())
  )
}

#' @keywords internal
#' @noRd
.cf_sl_predict_member <- function(object, newdata) {
  if (inherits(object, "cf_sl_mean")) {
    return(rep(object$value, NROW(newdata)))
  }
  .predict_base_learner(object, newdata)
}

# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------

#' Solve the library weights
#'
#' Non-negative least squares on the out-of-fold prediction matrix, then
#' divided by the sum. The quadratic program is the one `cf_ensemble()` uses.
#'
#' @param Z Numeric matrix (n x K) of out-of-fold predictions.
#' @param y Numeric target of length n.
#' @return Named numeric of length K summing to 1. Columns that could not be
#'   produced carry weight 0.
#' @keywords internal
#' @noRd
.cf_sl_solve_weights <- function(Z, y, where = "evaluation nuisance") {
  K <- ncol(Z)
  members <- colnames(Z)
  ok <- apply(Z, 2L, function(z) all(is.finite(z)))
  if (!any(ok)) {
    rlang::abort(
      c(sprintf("Every Super Learner member failed to fit (%s).", where),
        i = "Check the outcome / treatment vector and the covariate matrix."),
      class = "cateomics_sl_all_members_failed"
    )
  }
  if (!all(ok)) {
    rlang::warn(
      c(sprintf("Super Learner member(s) dropped (%s): %s.",
                where, paste(members[!ok], collapse = ", ")),
        i = "The weights were solved over the members that did fit."),
      class = "cateomics_sl_member_failed"
    )
  }
  Zk <- Z[, ok, drop = FALSE]
  D  <- crossprod(Zk)
  # Jitter the diagonal so the program stays positive definite when two
  # members predict almost the same thing, as ridge and elastic net can.
  diag(D) <- diag(D) + 1e-8 * max(diag(D), 1)
  d <- as.numeric(crossprod(Zk, y))
  kk <- ncol(Zk)
  qp <- tryCatch(
    quadprog::solve.QP(Dmat = D, dvec = d, Amat = diag(kk),
                       bvec = rep(0, kk), meq = 0L),
    error = function(e) NULL
  )
  wk <- if (is.null(qp)) {
    rlang::warn(
      sprintf("Super Learner weight solver failed (%s); falling back to equal weights.",
              where),
      class = "cateomics_sl_qp_failed"
    )
    rep(1 / kk, kk)
  } else {
    w <- qp$solution
    # Clear the solver's numerical dust before normalising, so a member the
    # program left out does not come back as 1e-17 of the total.
    w[w < 1e-8] <- 0
    if (sum(w) > 0) {
      w / sum(w)
    } else {
      # Nothing predicted the target better than zero: fall back to equal
      # weights over the members that fitted, and say so.
      rlang::warn(
        c(sprintf("Every Super Learner member got zero weight (%s).", where),
          i = "Falling back to equal weights over the members that fitted."),
        class = "cateomics_sl_zero_weights"
      )
      rep(1 / kk, kk)
    }
  }
  out <- stats::setNames(numeric(K), members)
  out[ok] <- wk
  out
}

# ---------------------------------------------------------------------------
# Out-of-fold assembly
# ---------------------------------------------------------------------------

#' Out-of-fold Super Learner predictions for one nuisance quantity
#'
#' Fills an n x K matrix of held-out member predictions over the folds it is
#' given, solves the weights on it, and returns the weighted combination. The
#' fold models are what the weights are read off; the returned prediction for a
#' row comes from the fold that did not train on it, exactly as the
#' single-learner path does.
#'
#' @param x Covariate matrix (n x p).
#' @param y Target: the outcome for mu0 / mu1 / m, the treatment for the
#'   propensity.
#' @param folds List of test-row index vectors.
#' @param members Character vector of member labels.
#' @param type `"regression"` or `"classification"`.
#' @param train_rows Row indices a member may train on (the control arm for
#'   mu0, the treated arm for mu1, every row otherwise).
#' @param weight_rows Row indices the weights are solved on (the arm whose
#'   conditional mean is being estimated).
#' @return List with `pred` (numeric n) and `weights` (named numeric).
#' @keywords internal
#' @noRd
.cf_sl_oof_quantity <- function(x, y, folds, members, type,
                                random_state = NULL, member_params = list(),
                                train_rows = NULL, weight_rows = NULL,
                                where = "evaluation nuisance") {
  n <- length(y)
  K <- length(members)
  train_rows  <- train_rows  %||% seq_len(n)
  weight_rows <- weight_rows %||% seq_len(n)
  Z <- matrix(NA_real_, nrow = n, ncol = K, dimnames = list(NULL, members))

  trained_folds <- 0L
  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- intersect(setdiff(seq_len(n), test_idx), train_rows)
    if (length(train_idx) <= 1L) next
    trained_folds <- trained_folds + 1L
    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_test  <- x[test_idx, , drop = FALSE]
    for (j in seq_len(K)) {
      pred <- tryCatch({
        fit <- .cf_sl_fit_member(x_train, y_train, type = type,
                                 member = members[j],
                                 random_state = random_state,
                                 params = member_params[[members[j]]] %||% list())
        .cf_sl_predict_member(fit, x_test)
      }, error = function(e) NULL)
      if (!is.null(pred) && length(pred) == length(test_idx)) {
        Z[test_idx, j] <- as.numeric(pred)
      }
    }
  }

  # No fold had enough rows of this arm to fit anything. That is not a member
  # failure, and the single-learner path leaves such a quantity at zero rather
  # than stopping, so this one does the same and names the real cause.
  if (trained_folds == 0L) {
    rlang::warn(
      c(sprintf("No fold had enough rows to fit the %s.", where),
        i = sprintf("%d row(s) are available for it across the whole sample.",
                    length(train_rows)),
        i = "That quantity is left at zero, as a single learner would leave it."),
      class = "cateomics_sl_no_training_rows"
    )
    return(list(pred = numeric(n),
                weights = stats::setNames(numeric(K), members)))
  }

  # Non-negative least squares does not determine the weights unless the rows
  # it is solved on outnumber the members. The fold structure stays honest
  # below that -- no row is predicted by a model that saw it -- so nothing
  # else catches it.
  if (length(weight_rows) <= K) {
    rlang::abort(
      c(sprintf("The %s is weighted on %d row(s) over %d library members.",
                where, length(weight_rows), K),
        i = "The weights are not determined unless the rows outnumber the members.",
        i = sprintf("Use at most %d member(s) in `eval_library`, or name a single learner in `eval_nuisance`.",
                    max(1L, length(weight_rows) - 1L))),
      class = "cateomics_sl_library_too_large"
    )
  }

  # A member is usable when every row the weights are read on has a prediction.
  usable <- Z[weight_rows, , drop = FALSE]
  w <- .cf_sl_solve_weights(usable, as.numeric(y[weight_rows]), where = where)
  # Rows no fold could train for have no member prediction. Zeroing them and
  # keeping the full weight would shrink the combination towards zero on those
  # rows without saying so, so the weight is renormalised per row over the
  # members that do predict there. Rows where nothing predicts stay at zero,
  # which is where the single-learner path leaves them through its
  # `numeric(n)` start.
  ok_mat <- is.finite(Z)
  Zf <- Z
  Zf[!ok_mat] <- 0
  wsum <- as.numeric(ok_mat %*% w)
  pred <- as.numeric(Zf %*% w)
  scale_rows <- wsum > 0 & abs(wsum - 1) > 1e-12
  pred[scale_rows] <- pred[scale_rows] / wsum[scale_rows]
  pred[wsum <= 0] <- 0
  list(pred = pred, weights = w)
}

#' Full-data Super Learner for one nuisance quantity
#'
#' Used by the in-sample path (`partition = "none"`), which has no outer folds.
#' The weights come from an internal split, because non-negative least squares
#' on in-sample predictions puts nearly all the weight on whichever member
#' interpolates; the returned prediction is then the weighted combination of
#' members refitted on every row, which is in-sample exactly as that path
#' already is.
#'
#' @inheritParams .cf_sl_oof_quantity
#' @param n_folds Fold count for the internal weight split.
#' @param strata Vector to stratify the internal split on (the treatment).
#' @keywords internal
#' @noRd
.cf_sl_full_quantity <- function(x, y, members, type, strata,
                                 n_folds = 5L, random_state = NULL,
                                 member_params = list(),
                                 train_rows = NULL, weight_rows = NULL,
                                 where = "evaluation nuisance") {
  n <- length(y)
  train_rows  <- train_rows  %||% seq_len(n)
  weight_rows <- weight_rows %||% seq_len(n)
  folds <- .create_stratified_folds(strata, n_folds = n_folds,
                                    random_state = random_state,
                                    context = "sl_weights")
  wz <- .cf_sl_oof_quantity(x, y, folds, members, type,
                            random_state = random_state,
                            member_params = member_params,
                            train_rows = train_rows,
                            weight_rows = weight_rows,
                            where = where)
  w <- wz$weights
  x_train <- x[train_rows, , drop = FALSE]
  y_train <- y[train_rows]
  pred <- numeric(n)
  used  <- 0
  lost  <- character(0)
  for (j in seq_along(members)) {
    if (w[[members[j]]] <= 0) next
    p <- tryCatch({
      fit <- .cf_sl_fit_member(x_train, y_train, type = type,
                               member = members[j],
                               random_state = random_state,
                               params = member_params[[members[j]]] %||% list())
      .cf_sl_predict_member(fit, x)
    }, error = function(e) NULL)
    if (is.null(p) || length(p) != n || !all(is.finite(p))) {
      lost <- c(lost, members[j])
      next
    }
    pred <- pred + w[[members[j]]] * as.numeric(p)
    used <- used + w[[members[j]]]
  }
  # A member can fit on the internal split and then fail on every row. Its
  # share would otherwise stay in the sum as zero and pull the prediction
  # towards it, so the remaining weight is renormalised and the loss is named.
  if (length(lost)) {
    rlang::warn(
      c(sprintf("Super Learner member(s) failed the full-data fit (%s): %s.",
                where, paste(lost, collapse = ", ")),
        i = "Their weight was redistributed over the members that did fit."),
      class = "cateomics_sl_member_failed"
    )
    if (used > 0) pred <- pred / used
  }
  list(pred = pred, weights = w)
}

#' Format library weights for printing
#'
#' Members carrying at least 0.005, largest first.
#' @keywords internal
#' @noRd
.cf_sl_format_weights <- function(w, digits = 2L) {
  if (is.null(w) || !length(w)) return(NA_character_)
  # Below half a point the member is not part of the combination in any
  # readable sense, and listing it crowds out the ones that are.
  w <- w[is.finite(w) & w >= 0.005]
  if (!length(w)) return("none")
  w <- sort(w, decreasing = TRUE)
  # The leading zero stays: dropping it put ".86" and "1.00" in one line, which
  # reads as two different quantities.
  paste(sprintf("%s %s", names(w), formatC(w, format = "f", digits = digits)),
        collapse = "  ")
}

#' Does either slot use the Super Learner?
#' @keywords internal
#' @noRd
.cf_sl_any_slot_is_sl <- function(eval_nuisance) {
  spec <- .cateomics_resolve_eval_nuisance(eval_nuisance)
  .cf_sl_is_token(spec$outcome) || .cf_sl_is_token(spec$propensity)
}

#' The library each slot actually used
#'
#' Recorded in the metadata so a fit says which members its weights refer to.
#' A slot fitted by a single learner has no library.
#'
#' @return List with `outcome` / `propensity`, each a character vector or NULL.
#' @keywords internal
#' @noRd
.cf_sl_libraries_used <- function(eval_nuisance, eval_library) {
  spec <- .cateomics_resolve_eval_nuisance(eval_nuisance)
  list(
    outcome = if (.cf_sl_is_token(spec$outcome)) {
      .cf_sl_resolve_library(eval_library, "outcome")
    } else NULL,
    propensity = if (.cf_sl_is_token(spec$propensity)) {
      .cf_sl_resolve_library(eval_library, "propensity")
    } else NULL
  )
}

#' Label the evaluation nuisance for printing
#'
#' Adds the member count when the Super Learner is in use, so the line says
#' what the weights below it refer to.
#'
#' @keywords internal
#' @noRd
.cf_eval_nuisance_label <- function(metadata) {
  lab <- metadata$eval_nuisance
  if (!is.character(lab) || length(lab) != 1L) {
    lab <- paste(unlist(lab), collapse = " / ")
  }
  lib <- metadata$eval_library
  k <- length(lib$outcome %||% lib$propensity %||% character(0))
  if (k > 0L) sprintf("%s (%d members)", lab, k) else lab
}

#' One line of weights for the first nuisance slot that has any
#'
#' Which slot that is varies with what was fitted, so the name goes in the
#' line: the weights of `m` and of `pi` are different quantities.
#'
#' @param sl_weights The `eval_sl_weights` record.
#' @return A string of the form "m-hat: ridge 0.86  rf 0.14", or NA when no
#'   slot carries weights.
#' @keywords internal
#' @noRd
.cf_sl_weight_line <- function(sl_weights) {
  if (is.null(sl_weights)) return(NA_character_)
  labs <- c(m = "m-hat", mu0 = "mu0-hat", mu1 = "mu1-hat", pi = "pi-hat")
  for (slot in names(labs)) {
    w <- sl_weights[[slot]]
    if (is.null(w)) next
    txt <- .cf_sl_format_weights(w)
    if (is.na(txt)) next
    return(sprintf("%s: %s", labs[[slot]], txt))
  }
  NA_character_
}
