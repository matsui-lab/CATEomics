#' The fold count every penalty selection in this package uses
#'
#' One number for every split that selects a penalty, so a single threshold
#' describes all of them. It is written into `inst/python/learner_factory.py`
#' for `LassoCV` / `ElasticNetCV` and the three `LogisticRegressionCV`
#' branches -- a continuous target's `"ridge"` is a closed-form
#' leave-one-out there and forms no split -- passed to `cv.glmnet` as `nfolds` by
#' `.fit_glmnet`, and used for the ridge surrogate in `cf_interpret`. R cannot
#' read the Python value across the process boundary, so it is restated here
#' and `test-input-validation.R` reads that file to check the two have not
#' drifted apart.
#'
#' @keywords internal
#' @noRd
.CF_CV_FOLDS <- 5L

#' Base-learner tokens whose fit selects a penalty by cross-validation
#'
#' The same three tokens in every role and on both sides of the process
#' boundary: in R they become `cv.glmnet` with `alpha` 1 / 0 / 0.5, in Python
#' `LassoCV` / `RidgeCV` / `ElasticNetCV` for a continuous target and
#' `LogisticRegressionCV` for a class label. `"linear"` / `"logistic"` are
#' plain fits and `"rf"` / `"xgboost"` do not cross-validate, so none of them
#' is affected.
#'
#' @keywords internal
#' @noRd
.CF_CV_TOKENS <- c("lasso", "ridge", "elastic_net")

#' The one token that forms no partition, and where
#'
#' `"ridge"` on a continuous target in Python is `RidgeCV` left at
#' scikit-learn's default, which is the closed-form leave-one-out: ridge has an
#' identity giving the leave-one-out error from the hat matrix, and the other
#' two penalties do not. No number of rows can leave a fold short of it, so the
#' fold rule does not apply there. It applies everywhere else, including:
#'
#' * `"ridge"` on a class label in Python, which is `LogisticRegressionCV` and
#'   does split (stratified).
#' * `"ridge"` in R in either role, which is `cv.glmnet` at `alpha = 0` and
#'   always splits. `cv.glmnet` has no closed-form leave-one-out -- only
#'   `nfolds = n`, by brute force -- so this asymmetry belongs to the two
#'   libraries and is not one this package introduces.
#'
#' @keywords internal
#' @noRd
.CF_CV_NO_SPLIT_PY_ROWS <- "ridge"

#' Rows a learner is left with after the splits that sit above it
#'
#' Folds nest. A per-group model inside an out-of-fold quantity is fitted on
#' the group's rows outside the held-out fold, and when the method also
#' cross-fits internally (only the DR learner does) that part is split again
#' before the learner ever sees it. With a stratified split of `m` rows over
#' `K` folds the held-out share is at most `ceiling(m / K)`, so each nesting
#' level subtracts that.
#'
#' @param m Rows (or rows of one class) entering the outermost split.
#' @param n_folds Outer fold count, or NULL when there is no outer split.
#' @param cv The method's internal cross-fit count, or NULL when it has none.
#' @return The count reaching the learner.
#' @keywords internal
#' @noRd
.cf_cv_rows_left <- function(m, n_folds = NULL, cv = NULL) {
  left <- as.integer(m)
  for (k in list(n_folds, cv)) {
    if (is.null(k)) next
    k <- as.integer(k)
    if (is.na(k) || k < 2L) next
    left <- left - ceiling(left / k)
  }
  as.integer(left)
}

#' Decide on a penalty-selection split from the count reaching it
#'
#' One rule for every penalty selection in the package, stated on the count
#' the learner actually receives:
#'
#' * below `.CF_CV_FOLDS` the split cannot be formed at all -- some fold would
#'   hold no rows -- so the penalty cannot be selected and this stops.
#' * below twice that the split exists but every fold holds fewer than two
#'   rows, so the penalty is selected on almost nothing. That is a fit, so it
#'   is reported and the fit continues.
#'
#' The two implementations behave differently on their own: scikit-learn
#' raises a `ValueError` when the rows cannot fill the folds, while glmnet
#' returns a fit without stopping or warning (measured on a continuous target
#' at 3, 5 and 8 rows). Neither signals the thin-but-formed band. Checking
#' both here with the same number gives one cause one answer.
#'
#' @param count Rows, or rows of the smallest class, reaching the learner.
#' @param learner The resolved base-learner token.
#' @param role Which slot, for the message.
#' @param what `"rows"` or `"class"`, which the message names.
#' @param nested Text naming the splits already subtracted, or NULL.
#' @return Invisibly `TRUE`. Errors or warns otherwise.
#' @keywords internal
#' @noRd
.cf_check_cv_split <- function(count, learner, role = NULL,
                               what = c("rows", "class"), nested = NULL,
                               cv_folds = .CF_CV_FOLDS) {
  what <- match.arg(what)
  if (!is.character(learner) || length(learner) != 1L ||
      !(learner %in% .CF_CV_TOKENS)) {
    return(invisible(TRUE))
  }
  # The count this package writes is a default a caller can lower, so the
  # decision follows the number actually in force rather than the constant:
  # refusing four rows to a three-fold split would refuse a fit glmnet performs.
  cv <- as.integer(cv_folds)
  if (is.na(cv) || cv < 2L) return(invisible(TRUE))
  count <- as.integer(count)
  # A count that could not be formed decides nothing. Every caller derives it
  # from a `sum()` or a `table()`, so this is unreachable from user data; the
  # alternative is `if (NA)` stopping with an unclassed message.
  if (is.na(count)) return(invisible(TRUE))
  # Name the knob the caller set when it is known. Inside `.fit_glmnet` it is
  # not: the token may have arrived from any of the `model_*` slots, from
  # `eval_nuisance`, or from a Super Learner member, so the learner is named
  # instead of inventing an argument the caller never wrote.
  subject <- if (is.null(role)) {
    sprintf("The %s", gsub("_", " ", learner))
  } else {
    sprintf("`%s = \"%s\"`", role, learner)
  }
  held <- if (identical(what, "class")) "the smallest class has" else "it has"
  after <- if (is.null(nested)) "" else sprintf(" after %s", nested)
  alt <- if (identical(what, "class")) "logistic" else "linear"
  remedy <- if (is.null(role)) {
    sprintf("Name a learner that does not cross-validate (`\"%s\"`, `\"rf\"`, `\"xgboost\"`), raise `n_folds`, or collect more rows.",
            alt)
  } else {
    sprintf("Set `%s = \"%s\"`, which does not cross-validate, raise `n_folds`, or collect more rows.",
            role, alt)
  }
  if (count < cv) {
    rlang::abort(
      c(sprintf("%s selects its penalty over %d folds, and %s %d row(s)%s.",
                subject, cv, held, max(count, 0L), after),
        i = "That many folds cannot be formed, so the penalty cannot be selected.",
        i = remedy),
      class = "cateomics_cv_split_too_small"
    )
  }
  if (count < 2L * cv) {
    msg <- sprintf("%s selects its penalty over %d folds, and %s %d row(s)%s.",
                   subject, cv, held, count, after)
    # One call fits a penalised learner many times -- once per fold, once per
    # Super Learner member, once per candidate -- and every one of them is the
    # same fact about the same data. Inside a public entry the facts are
    # collected and reported together; called on its own, this says it at once.
    if (.cf_cv_thin$depth > 0L) {
      .cf_cv_thin$rec[[length(.cf_cv_thin$rec) + 1L]] <-
        list(subject = subject, count = count, cv = cv, msg = msg,
             remedy = remedy)
      return(invisible(TRUE))
    }
    rlang::warn(
      c(msg,
        i = "That is fewer than two rows per fold on average, so the penalty is chosen on almost no information.",
        i = remedy),
      class = "cateomics_cv_split_thin"
    )
  }
  invisible(TRUE)
}

# What the penalised learners of the call now running had to work with. The
# depth counter is what makes a nested entry report once rather than twice.
.cf_cv_thin <- new.env(parent = emptyenv())
.cf_cv_thin$depth <- 0L
.cf_cv_thin$rec   <- list()

#' Start collecting thin-split reports for one public call
#'
#' @keywords internal
#' @noRd
.cf_cv_thin_open <- function() {
  .cf_cv_thin$depth <- .cf_cv_thin$depth + 1L
  if (.cf_cv_thin$depth == 1L) .cf_cv_thin$rec <- list()
  invisible(NULL)
}

#' Report them, once, as the call returns
#'
#' One fit of one learner reads exactly as it did before. Several -- the Super
#' Learner's penalised members, the same learner on each fold, a candidate's
#' own learner -- are one line naming the learners and the smallest count they
#' had, because the remedy is the same for all of them.
#'
#' @keywords internal
#' @noRd
.cf_cv_thin_close <- function() {
  .cf_cv_thin$depth <- max(0L, .cf_cv_thin$depth - 1L)
  if (.cf_cv_thin$depth > 0L) return(invisible(NULL))
  rec <- .cf_cv_thin$rec
  .cf_cv_thin$rec <- list()
  if (!length(rec)) return(invisible(NULL))
  if (length(rec) == 1L) {
    rlang::warn(
      c(rec[[1L]]$msg,
        i = "That is fewer than two rows per fold on average, so the penalty is chosen on almost no information.",
        i = rec[[1L]]$remedy),
      class = "cateomics_cv_split_thin")
    return(invisible(NULL))
  }
  # The subject reads as a sentence opener ("The ridge selects ...") or as the
  # argument the caller wrote; listed, it needs neither the article nor a verb.
  subjects <- unique(sub("^The ", "",
                         vapply(rec, function(z) z$subject, character(1))))
  smallest <- min(vapply(rec, function(z) z$count, numeric(1)))
  cv       <- max(vapply(rec, function(z) z$cv, numeric(1)))
  rlang::warn(
    c(sprintf("%d penalised fit(s) selected a penalty over %d folds with as few as %d row(s).",
              length(rec), as.integer(cv), as.integer(smallest)),
      i = sprintf("The learners: %s.", paste(subjects, collapse = ", ")),
      i = "That is fewer than two rows per fold on average, so the penalty is chosen on almost no information.",
      i = rec[[1L]]$remedy),
    class = "cateomics_cv_split_thin")
  invisible(NULL)
}

#' The learner a slot resolves to, as the method will see it
#'
#' `cf_fit`'s own formals are all NULL: the defaults live on each method's
#' fitter, and the slots a method does not have are absent from it entirely
#' (a propensity learner is taken by bcf, x-learner and dr-learner; a final
#' model by x-learner and dr-learner). A check that
#' reads the public formal therefore sees NULL on the default path and skips
#' silently, which is the path that matters most. The fitter is asked instead.
#'
#' `model_final = NULL` on a method that has the slot means "follow
#' `model_outcome`", so that is what it resolves to.
#'
#' @param fitter The method's fitter, or NULL when it has none.
#' @param slot One of the three slot names.
#' @param user What the caller passed, or NULL.
#' @param outcome The already-resolved outcome learner, for `model_final`.
#' @return A single token, or NULL when the method has no such slot.
#' @keywords internal
#' @noRd
.cf_resolve_learner_slot <- function(fitter, slot, user, outcome = NULL) {
  if (!is.null(user)) return(user)
  if (is.null(fitter) || !(slot %in% names(formals(fitter)))) return(NULL)
  d <- formals(fitter)[[slot]]
  if (is.character(d) && length(d) == 1L) return(d)
  if (identical(slot, "model_final")) return(outcome)
  NULL
}

#' The settings a result records as having run
#'
#' Both entry points record the same four slots, so they resolve them here
#' rather than each writing the rule out. A slot is recorded only when it
#' describes something the fit actually did:
#'
#' * the three learner slots take the fitter's own default when the caller
#'   named none, and are `NA` for a method with no such slot;
#' * the propensity learner is `NA` when a known propensity was supplied,
#'   because then no learner is fitted for it at all;
#' * `clip` is `NA` for the methods that do not estimate a propensity, since
#'   the bounds never reach them.
#'
#' @param method The method name.
#' @param model_outcome,model_final,model_propensity,propensity_score,clip
#'   The arguments as given.
#' @return A list with `model_outcome`, `model_final`, `model_propensity` and
#'   `clip`.
#' @keywords internal
#' @noRd
.cf_spec_slots <- function(method, model_outcome, model_final,
                           model_propensity, propensity_score, clip) {
  fitter <- .cf_method_fitter(method)
  has    <- function(slot) !is.null(fitter) && slot %in% names(formals(fitter))
  mo <- .cf_resolve_learner_slot(fitter, "model_outcome", model_outcome)
  mf <- .cf_resolve_learner_slot(fitter, "model_final", model_final,
                                 outcome = mo)
  mp <- if (!is.null(propensity_score)) NULL else
    .cf_resolve_learner_slot(fitter, "model_propensity", model_propensity)
  list(
    model_outcome    = mo %||% NA_character_,
    model_final      = mf %||% NA_character_,
    model_propensity = mp %||% NA_character_,
    clip             = if (has("clip")) clip %||% .CF_PS_CLIP else NA_real_
  )
}


#' Check a penalty selection whose split is on rows
#'
#' The continuous-target case in Python, where the learner receives a plain row
#' split. `"ridge"` is exempt here and only here: see
#' `.CF_CV_NO_SPLIT_PY_ROWS`.
#'
#' @param rows Rows entering the outermost split.
#' @param n_folds Outer fold count, or NULL when there is no outer split.
#' @param cv The method's internal cross-fit count, or NULL when it has none.
#' @param learner The resolved base-learner token.
#' @param role Which slot, for the message.
#' @keywords internal
#' @noRd
.cf_check_cv_rows <- function(rows, n_folds, cv, learner, role,
                              cv_folds = .CF_CV_FOLDS) {
  if (!is.character(learner) || length(learner) != 1L ||
      !(learner %in% setdiff(.CF_CV_TOKENS, .CF_CV_NO_SPLIT_PY_ROWS))) {
    return(invisible(TRUE))
  }
  .cf_check_cv_split(.cf_cv_rows_left(rows, n_folds, cv), learner, role,
                     what = "rows", nested = .cf_cv_nested_text(n_folds, cv),
                     cv_folds = cv_folds)
}

#' Check a penalty selection whose split is stratified on a class label
#'
#' The propensity always, and a binary outcome within an arm. The smallest
#' class is what has to fill the folds.
#'
#' @param strata The classification target the split is stratified on.
#' @param n_folds Outer fold count, or NULL when there is no outer split.
#' @param cv The method's internal cross-fit count, or NULL when it has none.
#' @param learner The resolved base-learner token.
#' @param role Which slot, for the message.
#' @keywords internal
#' @noRd
.cf_check_cv_classes <- function(strata, n_folds, cv, learner, role,
                                 cv_folds = .CF_CV_FOLDS) {
  if (!is.character(learner) || length(learner) != 1L ||
      !(learner %in% .CF_CV_TOKENS)) {
    return(invisible(TRUE))
  }
  tb <- table(strata[!is.na(strata)])
  if (length(tb) < 2L) return(invisible(TRUE))
  left <- vapply(as.integer(tb), .cf_cv_rows_left, integer(1L),
                 n_folds = n_folds, cv = cv)
  .cf_check_cv_split(min(left), learner, role, what = "class",
                     nested = .cf_cv_nested_text(n_folds, cv),
                     cv_folds = cv_folds)
}

#' Name the splits a count has already been reduced by
#'
#' @keywords internal
#' @noRd
.cf_cv_nested_text <- function(n_folds, cv) {
  parts <- c(if (!is.null(n_folds)) "the outer split",
             if (!is.null(cv)) "the internal cross-fit")
  if (!length(parts)) NULL else paste(parts, collapse = " and ")
}

#' Refuse an internal cross-fit the treatment groups cannot support
#'
#' Distinct from a penalised learner's own cross-validation, which selects a
#' penalty: this is the split that estimates the nuisances out of fold so the
#' final stage is orthogonal. Among the shipped methods only the DR learner
#' has it (`cv`, default 5); the metalearners do not cross-fit, and grf and bcf
#' split on their own terms.
#'
#' econml stratifies that split on the treatment (`check_cv(cv, classifier =)`
#' returns a `StratifiedKFold` for a discrete treatment). scikit-learn does
#' **not** stop when a group is smaller than the split count -- it warns and
#' carries on with folds where a group can be absent from the training part
#' (measured: a minority of 1 leaves one fold with no rows of it at all). That
#' warning is a Python one, so it does not reach an R caller.
#'
#' The split runs on the training part of the outer fold, so the group count
#' reaching it is already reduced by that split.
#'
#' @param T Binary treatment vector.
#' @param n_folds Outer fold count, or NULL when there is no outer split.
#' @param cv The method's internal cross-fit count.
#' @param method Method name, for the message.
#' @return Invisibly `TRUE`. Errors otherwise.
#' @keywords internal
#' @noRd
.cf_check_groups_for_inner_cv <- function(T, n_folds, cv, method) {
  cv <- as.integer(cv)
  if (is.na(cv) || cv < 2L) return(invisible(TRUE))
  tb <- table(as.integer(T))
  if (length(tb) < 2L) return(invisible(TRUE))
  sizes <- stats::setNames(as.integer(tb), names(tb))
  label <- c("0" = "control", "1" = "treated")
  # What the outer split leaves, then what the inner split leaves of that.
  outer_left <- if (is.null(n_folds)) sizes else sizes - ceiling(sizes / as.integer(n_folds))
  inner_left <- outer_left - ceiling(outer_left / cv)
  short <- which(inner_left < 2L)
  if (!length(short)) return(invisible(TRUE))
  nm <- names(sizes)[short]
  txt <- paste(sprintf("%s group: %d row(s)%s, leaving %d to the cross-fit",
                       ifelse(nm %in% names(label), label[nm], nm), sizes[short],
                       if (is.null(n_folds)) "" else sprintf(" (%d after the outer split)", outer_left[short]),
                       pmax(inner_left[short], 0L)),
               collapse = "; ")
  rlang::abort(
    c(sprintf("%s cross-fits internally with cv = %d, and that needs at least 2 rows of each treatment group to train on (%s).",
              method, cv, txt),
      i = "scikit-learn warns rather than stops here, and that warning does not reach R.",
      i = sprintf("Lower it with `params = list(%s = list(cv = 2))`, raise `n_folds`, or collect more of that group.",
                  method)),
    class = "cateomics_group_too_small_for_cv"
  )
}

#' Check every split that sits below one method's fit
#'
#' Called from both public entries so a configuration refused by one is
#' refused by the other. Two different splits live here:
#'
#' * the internal cross-fit, which estimates the nuisances out of fold so the
#'   final stage is orthogonal. Only the DR learner has one. It needs an
#'   arm-specific model to be fittable at all, so it takes the same two-row
#'   floor as the outer fold.
#' * each penalised learner's own cross-validation, which selects the penalty.
#'   A different question with a different threshold, and it nests below both
#'   splits above it.
#'
#' Both are checked here rather than left to upstream because scikit-learn's
#' signals happen inside the Python process and never reach an R caller.
#'
#' @param method Resolved method name.
#' @param T Binary treatment vector (already resolved, not a column name).
#' @param Y Outcome vector.
#' @param binary_outcome Whether the outcome is a class label.
#' @param n_folds Outer fold count, or NULL when there is no outer split.
#' @param method_params This method's resolved `params` block.
#' @param model_outcome,model_propensity,model_final What the caller passed,
#'   or NULL to take the method's own default.
#' @return Invisibly `TRUE`. Errors or warns otherwise.
#' @keywords internal
#' @noRd
.cf_check_method_splits <- function(method, T, Y, binary_outcome, n_folds,
                                    method_params = list(),
                                    model_outcome = NULL,
                                    model_propensity = NULL,
                                    model_final = NULL) {
  fitter <- .cf_method_fitter(method)
  inner_cv <- NULL
  if (!is.null(fitter) && "cv" %in% names(formals(fitter))) {
    inner_cv <- method_params$cv %||% formals(fitter)$cv
    .cf_check_groups_for_inner_cv(T, n_folds = n_folds, cv = inner_cv,
                                  method = method)
  }
  if (!(method %in% c("slearner", "tlearner", "xlearner", "drlearner"))) {
    return(invisible(TRUE))
  }
  # A slot's own fold count, when the caller set one. Distinct from
  # `method_params$cv`, which is the internal cross-fit: `.cf_apply_params()`
  # splits them by the fitter's formals, so `params = list(cv = 3)` is the
  # cross-fit and `params = list(outcome = list(cv = 3))` is the learner's.
  slot_cv <- function(slot) {
    v <- method_params[[slot]]$cv
    if (is.null(v)) .CF_CV_FOLDS else as.integer(v)
  }
  out <- .cf_resolve_learner_slot(fitter, "model_outcome", model_outcome)
  ps  <- .cf_resolve_learner_slot(fitter, "model_propensity", model_propensity)
  fin <- .cf_resolve_learner_slot(fitter, "model_final", model_final,
                                  outcome = out)
  # The propensity's target is always a class label, so its split is
  # stratified on the treatment.
  .cf_check_cv_classes(T, n_folds, inner_cv, ps, "model_propensity",
                       cv_folds = slot_cv("propensity"))
  # The binding count for the outcome slot is the smaller arm. That is what
  # the t- and x-learner actually fit on, arm by arm; the s-learner and the
  # DR-learner fit one interaction-augmented model on every row, so for those
  # two the check is the conservative one. Reported once, for the arm that
  # binds -- checking each arm separately said the same sentence twice on any
  # balanced design, which carries no information the first one did not.
  if (isTRUE(binary_outcome)) {
    worst <- which.min(vapply(c(0L, 1L), function(arm) {
      tb <- table(Y[as.integer(T) == arm])
      if (length(tb) < 2L) NA_integer_ else min(as.integer(tb))
    }, integer(1L)))
    if (length(worst)) {
      .cf_check_cv_classes(Y[as.integer(T) == c(0L, 1L)[worst]], n_folds,
                           inner_cv, out, "model_outcome",
                           cv_folds = slot_cv("outcome"))
    }
  } else {
    .cf_check_cv_rows(min(sum(as.integer(T) == 0L), sum(as.integer(T) == 1L)),
                      n_folds, inner_cv, out, "model_outcome",
                      cv_folds = slot_cv("outcome"))
  }
  # The final stage regresses the pseudo-outcome over every row, after the
  # outer fold but not inside the internal cross-fit.
  .cf_check_cv_rows(length(Y), n_folds, NULL, fin, "model_final",
                    cv_folds = slot_cv("final"))
  invisible(TRUE)
}

#' Refuse a fold structure the treatment groups cannot support
#'
#' Every out-of-fold quantity that is estimated separately per treatment group
#' (the evaluation nuisance's mu0 / mu1, and each per-group model a method
#' fits) is trained on that group's rows **outside** the held-out fold. With a
#' stratified split of `m` rows of a group over `K` folds, the held-out share
#' is at most `ceiling(m / K)`, so the training part has at least
#' `m - ceiling(m / K)` of them. Below two there is nothing to fit.
#'
#' Left unchecked, the evaluation nuisance returns a hard zero for the
#' affected rows, which then enters the doubly-robust pseudo-outcome.
#'
#' Note the direction: *raising* `n_folds` is what relaxes this, because each
#' fold then holds out a thinner slice of a scarce group. A group of two
#' cannot satisfy it at any `n_folds`, since the held-out share is never below
#' one.
#'
#' @param T Binary treatment vector.
#' @param n_folds Number of folds the split will use.
#' @param context Which split is being built; selects the remedy named in the
#'   message.
#' @return Invisibly `TRUE`. Errors otherwise.
#' @keywords internal
#' @noRd
.cf_check_groups_for_folds <- function(T, n_folds,
                                       context = c("folds", "sl_weights")) {
  context <- match.arg(context)
  n_folds <- as.integer(n_folds)
  tb <- table(as.integer(T))
  if (length(tb) < 2L) return(invisible(TRUE))   # single-group data is caught upstream
  sizes <- stats::setNames(as.integer(tb), names(tb))
  label <- c("0" = "control", "1" = "treated")
  train_left <- sizes - ceiling(sizes / n_folds)
  short <- which(train_left < 2L)
  if (!length(short)) return(invisible(TRUE))
  nm  <- names(sizes)[short]
  txt <- paste(sprintf("%s group: %d row(s), leaving %d to train on",
                       ifelse(nm %in% names(label), label[nm], nm),
                       sizes[short], train_left[short]),
               collapse = "; ")
  # A group of two or fewer cannot satisfy this at any number of folds, since
  # the held-out share is never below one. Saying "raise n_folds" there would
  # send the caller after a remedy that does not exist.
  hopeless <- any(sizes[short] <= 2L)
  need_k <- if (hopeless) NA_integer_ else {
    m <- min(sizes[short]); as.integer(ceiling(m / (m - 2L)))
  }
  fold_advice <- if (hopeless) {
    "No number of folds can satisfy this: a fold always holds out at least one row of the group."
  } else {
    sprintf("Raising `n_folds` relaxes this (%d or more would do): each fold then holds out fewer of that group.",
            need_k)
  }
  remedy <- if (identical(context, "sl_weights")) {
    c(i = paste0("The Super Learner solves its weights on held-out predictions, ",
                 "so it needs this even when the nuisance itself is in-sample."),
      i = fold_advice,
      i = "Or name a single learner in `eval_nuisance`, or collect more of that group.")
  } else {
    c(i = fold_advice,
      i = "Or use `partition = \"none\"`, or collect more of that group.")
  }
  rlang::abort(
    c(sprintf("%d-fold cross-fitting needs at least 2 rows of each treatment group in every fold's training part (%s).",
              n_folds, txt),
      remedy),
    class = "cateomics_group_too_small_for_folds"
  )
}

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
#' @param context Which quantity the split is for. Only the wording of the
#'   refusal differs: `"sl_weights"` says the rows are what the Super
#'   Learner weights are solved on.
#' @keywords internal
#' @noRd
.create_stratified_folds <- function(T, n_folds = 5L, random_state = NULL,
                                     context = c("folds", "sl_weights")) {
  n_folds <- as.integer(n_folds)
  stopifnot(n_folds >= 2L, length(T) >= n_folds)
  .cf_check_groups_for_folds(T, n_folds, context = context)

  # A missing seed is read as 0, as the dispatch does: a split drawn from the
  # ambient RNG state would make everything built on these folds depend on
  # what the caller happened to run first.
  random_state <- as.integer(random_state %||% 0L)
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
      # round-robin allocation, so the stratum is spread across the folds
      # instead of filling the first ones.
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


#' Null-coalescing operator, used package-wide
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
#' @param eval_nuisance Learner name, or `"super_learner"` (the default)
#' @param precomputed_nuisance Optional list to bypass internal nuisance
#'   computation (for cf_compare's shared nuisance path)
#' @param params Resolved hyperparameter block for this method
#'
#' @return list:
#'   * `fit` = method-specific full-data fit object (unwrapped; class attached
#'     by the caller)
#'   * `oof_ite` = numeric(n) K-fold OOF tau
#'   * `oof_mu0`, `oof_mu1` = numeric(n) K-fold OOF potential outcomes (NA for
#'     methods that surface none)
#'   * `oof_ps` = numeric(n) K-fold OOF propensity (NA for methods without one)
#'   * `fold_id` = integer(n) the fold that held each row out, so the split
#'     that produced the out-of-fold quantities travels with them
#'   * `nuisance` = list(mu0_oof, mu1_oof, pi_oof, m_oof) K-fold OOF nuisance
#'   * `fold_cache` = list of K fit objects (when keep_kfold_models = TRUE)
#'
#' @keywords internal
#' @noRd


.kfold_oof_predict_matrix <- function(X, T_vec, Y_vec, method,
                                       cov_names = NULL,
                                       n_folds = 5L,
                                       random_state = NULL,
                                       keep_kfold_models = FALSE,
                                       eval_nuisance = "super_learner",
                                       eval_library = NULL,
                                       eval_propensity_score = NULL,
                                       precomputed_nuisance = NULL,
                                       n_bootstrap = 100L,
                                       bootstrap_ci = FALSE,
                                       device = "cpu",
                                       propensity_X = NULL,
                                       propensity_score = NULL,
                                       clip = c(0.025, 0.975),
                                       binary_outcome = FALSE,
                                       model_outcome = NULL,
                                       model_final = NULL,
                                       model_propensity = NULL,
                                       propensity_params = list(),
                                       eval_params = list(),
                                       params = list()) {
  n <- length(Y_vec)
  if (is.null(cov_names) || length(cov_names) == 0L) {
    cov_names <- paste0("X", seq_len(ncol(X)))
  }

  # Arguments CATEomics supplies itself, as opposed to the user hyperparameters
  # in `params`. `.cf_dispatch_fit()` hands each method only the ones its own
  # signature declares.
  .control <- list(
    random_state     = random_state,
    device           = device,
    binary_outcome   = binary_outcome,
    bootstrap_ci     = isTRUE(bootstrap_ci),
    n_bootstrap      = as.integer(n_bootstrap),
    clip             = clip,
    propensity_X     = propensity_X,
    model_outcome    = model_outcome,
    model_final      = model_final,
    model_propensity = model_propensity,
    propensity_params = propensity_params
  )

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
      binary_outcome = binary_outcome,
      eval_params = eval_params,
      eval_library = eval_library,
      propensity_X = propensity_X
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
  # Held-out propensity for the methods that estimate one, assembled the same
  # way as the OOF potential outcomes. Stays NA for methods without a PS.
  oof_ps  <- rep(NA_real_, n)
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
    # A fold fit is read for its held-out predictions only: its interval is
    # never looked at, so resampling it would be `n_bootstrap` refits per fold
    # thrown away. The interval the caller sees comes from the full-data fit
    # below, which keeps the request.
    fit_k_inner <- .cf_dispatch_fit(
      method = method, X = X_train, T = T_train, Y = Y_train,
      cov_names = cov_names,
      control = c(utils::modifyList(.control, list(bootstrap_ci = FALSE)),
                  list(propensity_score = ps_for_fold)),
      params  = params)
    # Wrap with class for predict.cf_model dispatch (res = inner res,
    # model = the backend container).
    fit_k <- structure(
      list(method = method,
            res   = fit_k_inner$res,
            model = fit_k_inner,
            meta = list(
              covariate_names = cov_names,
              # `print` and `summary` on a cf_model report the shape they were
              # fitted on, so the fold objects carry it too: the rows this fold
              # trained on, and the columns of the preprocessed design. Without
              # them the shape line is dropped from `print` (sprintf() returns
              # character(0) for a NULL) and `summary` cannot be built.
              n = nrow(X_train),
              p = ncol(X_train),
              # The fold fit sees the already-preprocessed design, so the only
              # names it needs are the outcome and the treatment, which the
              # SE / MAE newdata branches drop from a container before building
              # X. The full record lives on the user-facing result.
              named_args = list(Y = "Y", T = "T"),
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
        class = "cateomics_partition_oof_failed"
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
    if (!is.null(pred_out$ps) &&
        length(pred_out$ps) == length(test_idx)) {
      oof_ps[test_idx] <- as.numeric(pred_out$ps)
    }
    if (isTRUE(keep_kfold_models)) fold_cache[[k]] <- fit_k
  }

  # 4. Full-data fit (method-specific, calls cf_fit_<method> directly)
  full_fit_inner <- .cf_dispatch_fit(
    method = method, X = X, T = T_vec, Y = Y_vec,
    cov_names = cov_names,
    control = c(.control, list(propensity_score = propensity_score)),
    params  = params)

  list(
    fit        = full_fit_inner,
    oof_ite    = oof_ite,
    oof_mu0    = oof_mu0,
    oof_mu1    = oof_mu1,
    oof_ps     = oof_ps,
    # Which fold held each row out. Every other out-of-fold quantity above is
    # indexed by row, so the split that produced them belongs beside them.
    fold_id    = .folds_to_fold_id(folds, n),
    nuisance   = nuisance_oof,
    fold_cache = fold_cache
  )
}


