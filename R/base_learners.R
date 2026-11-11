#' Fit a base learner
#'
#' Unified interface to seven base learners used across CATEomics methods for
#' propensity score, outcome model, and final stage fits.
#'
#' Supported learners:
#' * `"linear"` / `"logistic"` — aliases. Uses [stats::lm()] for regression
#'   and [stats::glm()] with `family = binomial()` for classification.
#' * `"lasso"` — [glmnet::cv.glmnet()] with `alpha = 1`.
#' * `"ridge"` — [glmnet::cv.glmnet()] with `alpha = 0`.
#' * `"elastic_net"` — [glmnet::cv.glmnet()] with `alpha = 0.5`.
#' * `"rf"` — [ranger::ranger()]; classification mode sets `probability = TRUE`.
#' * `"xgboost"` — [xgboost::xgb.train()].
#'
#' For the `cv.glmnet` variants, the regularisation strength is selected via
#' `s = "lambda.min"` at prediction time.
#'
#' @param x Numeric matrix or data.frame of predictors (`n` rows, `p` columns).
#' @param y Numeric vector of outcomes. For `type = "classification"` must be
#'   a 0/1 vector or a two-level factor.
#' @param type Either `"regression"` (continuous `y`) or `"classification"`
#'   (binary `y`).
#' @param learner One of `"linear"`, `"logistic"`, `"lasso"`, `"ridge"`,
#'   `"elastic_net"`, `"rf"`, `"xgboost"`.
#' @param random_state Optional integer seed. When non-`NULL`, pinned with
#'   [set.seed()] before the fit, which is where the penalty CV fold
#'   assignment, ranger and xgboost all draw from. No value is handed to
#'   ranger, which reads `seed = 0` as "no seed".
#' @param .params Named list of validated user hyperparameters, merged into the
#'   arguments below.
#' @param ... Additional arguments forwarded to the underlying fit function
#'   (`lm`, `glm`, `cv.glmnet`, `ranger` or `xgb.train`).
#'
#' @return An object of class `"cf_base_learner"` — a list with elements:
#'   * `learner`: the learner name.
#'   * `type`: `"regression"` or `"classification"`.
#'   * `model`: the fitted underlying model.
#'   * `feature_names`: column names of `x` (used to align `newdata` at predict
#'     time).
#'
#' @keywords internal
#' @noRd
.fit_base_learner <- function(x, y,
                              type = c("regression", "classification"),
                              learner = "linear",
                              random_state = NULL,
                              .params = list(),
                              ...) {
  type <- match.arg(type)
  .validate_learner_spec(learner)

  x_mat <- .coerce_to_matrix(x)
  y_in <- .coerce_y(y, type)

  # Save/restore the caller's global RNG so `random_state` does not leak
  # into `.Random.seed` via `set.seed()` (used by cv.glmnet fold assignment)
  # or downstream callers.
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
    set.seed(random_state)
  }

  # `.params` carries the user hyperparameters (validated upstream); `...`
  # keeps the internal call sites that pass arguments positionally by name.
  extra <- c(list(...), .params %||% list())
  model <- switch(learner,
    linear       = ,
    logistic     = do.call(.fit_linear, c(list(x_mat, y_in, type), extra)),
    lasso        = do.call(.fit_glmnet, c(list(x_mat, y_in, type, alpha = 1), extra)),
    ridge        = do.call(.fit_glmnet, c(list(x_mat, y_in, type, alpha = 0), extra)),
    elastic_net  = do.call(.fit_glmnet, c(list(x_mat, y_in, type, alpha = 0.5), extra)),
    rf           = do.call(.fit_ranger,
                           c(list(x_mat, y_in, type, random_state = random_state), extra)),
    xgboost      = do.call(.fit_xgboost,
                           c(list(x_mat, y_in, type, random_state = random_state), extra))
  )

  structure(
    list(
      learner       = learner,
      type          = type,
      model         = model,
      feature_names = colnames(x_mat)
    ),
    class = "cf_base_learner"
  )
}

#' Predict from a base learner
#'
#' @param object A `"cf_base_learner"` produced by [.fit_base_learner()].
#' @param newdata Numeric matrix or data.frame with the same columns used at fit
#'   time.
#'
#' @return Numeric vector of length `nrow(newdata)`. For regression, the point
#'   estimate. For classification, the probability of class 1 (`P(y = 1 | x)`).
#'
#' @keywords internal
#' @noRd
.predict_base_learner <- function(object, newdata) {
  if (!inherits(object, "cf_base_learner")) {
    rlang::abort(
      "`object` must be a cf_base_learner produced by .fit_base_learner().",
      class = "cateomics_base_learner_error"
    )
  }
  newdata_mat <- .align_newdata(newdata, object$feature_names)

  switch(object$learner,
    linear       = ,
    logistic     = .predict_linear(object$model, newdata_mat, object$type),
    lasso        = ,
    ridge        = ,
    elastic_net  = .predict_glmnet(object$model, newdata_mat, object$type),
    rf           = .predict_ranger(object$model, newdata_mat, object$type),
    xgboost      = .predict_xgboost(object$model, newdata_mat, object$type)
  )
}

#' The pluggable base learners, named once
#'
#' Read by `.validate_learner_spec()` and by the comparison entry point, so a
#' token is added in one place.
#'
#' @keywords internal
#' @noRd
.CF_BASE_LEARNERS <- c("linear", "logistic", "lasso", "ridge",
                       "elastic_net", "rf", "xgboost")


#' Validate a learner specification
#'
#' @param learner Scalar character, the learner name.
#' @return Invisibly `TRUE`. Errors otherwise.
#' @keywords internal
#' @noRd
.validate_learner_spec <- function(learner) {
  valid <- .CF_BASE_LEARNERS
  if (!is.character(learner) || length(learner) != 1L || !(learner %in% valid)) {
    rlang::abort(
      sprintf(
        "`learner` must be one of: %s (got: %s)",
        paste(shQuote(valid), collapse = ", "),
        paste(deparse(learner), collapse = " ")
      ),
      class = "cateomics_base_learner_error"
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.coerce_to_matrix <- function(x) {
  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      rlang::abort(
        "`x` must be a numeric matrix or a data.frame of numeric columns.",
        class = "cateomics_base_learner_error"
      )
    }
    x_out <- x
  } else if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      rlang::abort(
        "All columns of `x` must be numeric.",
        class = "cateomics_base_learner_error"
      )
    }
    x_out <- as.matrix(x)
  } else {
    rlang::abort(
      "`x` must be a numeric matrix or a data.frame.",
      class = "cateomics_base_learner_error"
    )
  }
  if (is.null(colnames(x_out))) {
    colnames(x_out) <- paste0("V", seq_len(ncol(x_out)))
  }
  x_out
}

.coerce_y <- function(y, type) {
  if (type == "regression") {
    if (!is.numeric(y)) {
      rlang::abort(
        "For regression, `y` must be a numeric vector.",
        class = "cateomics_base_learner_error"
      )
    }
    return(as.numeric(y))
  }
  # classification: accept 0/1 numeric or two-level factor
  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      rlang::abort(
        "For classification, factor `y` must have exactly two levels.",
        class = "cateomics_base_learner_error"
      )
    }
    return(as.numeric(y) - 1L)
  }
  y_num <- suppressWarnings(as.numeric(y))
  if (any(is.na(y_num)) || !all(y_num %in% c(0, 1))) {
    rlang::abort(
      "For classification, `y` must be a 0/1 vector or a two-level factor.",
      class = "cateomics_base_learner_error"
    )
  }
  y_num
}

.align_newdata <- function(newdata, feature_names) {
  nd <- .coerce_to_matrix(newdata)
  if (is.null(feature_names)) {
    return(nd)
  }
  if (!all(feature_names %in% colnames(nd))) {
    missing <- setdiff(feature_names, colnames(nd))
    rlang::abort(
      sprintf("`newdata` is missing columns: %s",
              paste(shQuote(missing), collapse = ", ")),
      class = "cateomics_base_learner_error"
    )
  }
  nd[, feature_names, drop = FALSE]
}

# ---- linear / logistic ----------------------------------------------------

.fit_linear <- function(x, y, type, ...) {
  df <- data.frame(.y = y, as.data.frame(x, check.names = FALSE),
                   check.names = FALSE)
  if (type == "regression") {
    stats::lm(.y ~ ., data = df, ...)
  } else {
    stats::glm(.y ~ ., data = df, family = stats::binomial(), ...)
  }
}

.predict_linear <- function(model, newdata, type) {
  nd <- as.data.frame(newdata, check.names = FALSE)
  if (inherits(model, "glm")) {
    as.numeric(stats::predict(model, newdata = nd, type = "response"))
  } else {
    as.numeric(stats::predict(model, newdata = nd))
  }
}

# ---- cv.glmnet (lasso / ridge / elastic_net) ------------------------------

.fit_glmnet <- function(x, y, type, alpha, ...) {
  # cv.glmnet requires at least 2 columns. Fall back to linear/logistic for
  # single-covariate inputs (typically toy / contract tests), so a penalised
  # token still runs on a degenerate design.
  if (ncol(x) < 2L) {
    return(.fit_linear(x, y, type, ...))
  }
  family <- if (type == "regression") "gaussian" else "binomial"
  # The fold count is this package's, not glmnet's default of 10, so that
  # every penalty selection that forms a split -- here, in the Python factory,
  # and for the ridge surrogate -- is split the same number of ways and one
  # threshold describes all of them. The exception is a continuous target in
  # the Python factory, where `"ridge"` is a closed-form leave-one-out and
  # forms no split. A caller's own `nfolds` still wins.
  dots <- list(...)
  nfolds <- as.integer(dots$nfolds %||% .CF_CV_FOLDS)
  dots$nfolds <- nfolds
  # Checked on the rows actually in hand, so no nesting arithmetic is needed:
  # whatever splits sit above this call have already been applied. glmnet
  # returns a fit on a continuous target without stopping or warning however
  # few rows it has, so the rule the Python side is pre-flighted against is
  # applied here too.
  role <- switch(as.character(alpha), "1" = "lasso", "0" = "ridge", "elastic_net")
  if (identical(family, "binomial")) {
    .cf_check_cv_split(min(as.integer(table(y[!is.na(y)]))), role,
                       what = "class", cv_folds = nfolds)
  } else {
    .cf_check_cv_split(nrow(x), role, what = "rows", cv_folds = nfolds)
  }
  # X is standardised centrally in `.cf_preprocess_input` (every column,
  # binary indicators included). Disable glmnet's internal standardisation so the
  # scale the caller sees is the one the penalty acts on.
  # cv.glmnet can still stop for reasons this package does not model -- a
  # degenerate design, a class that vanishes inside its own split -- with an
  # unclassed message that does not say what the caller has to change, so the
  # counts CATEomics is holding are reported alongside the original text.
  tryCatch(
    do.call(glmnet::cv.glmnet,
            c(list(x = x, y = y, alpha = alpha, family = family,
                   standardize = FALSE), dots)),
    error = function(e) {
      cnt <- if (identical(family, "binomial")) {
        tb <- table(as.integer(y))
        sprintf("%d rows, class counts %s", length(y),
                paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "/"))
      } else {
        sprintf("%d rows", length(y))
      }
      rlang::abort(
        c(sprintf("glmnet could not cross-validate the penalty on this fit (%s).", cnt),
          i = sprintf("It splits internally into %d folds; that many rows of each class are not there.",
                      nfolds),
          i = "Raise the rows available, or name a learner that does not cross-validate.",
          x = sprintf("glmnet: %s", sub("\n.*", "", conditionMessage(e)))),
        class = "cateomics_glmnet_cv_failed"
      )
    }
  )
}

.predict_glmnet <- function(model, newdata, type) {
  # Fallback path (single-covariate): model is an lm/glm from .fit_linear.
  if (inherits(model, c("lm", "glm"))) {
    return(.predict_linear(model, newdata, type))
  }
  ptype <- if (type == "regression") "link" else "response"
  pred <- stats::predict(model, newx = newdata, s = "lambda.min", type = ptype)
  as.numeric(pred)
}

# ---- ranger ---------------------------------------------------------------

.fit_ranger <- function(x, y, type, random_state = NULL, ...) {
  args <- list(x = x, ...)
  if (type == "regression") {
    args$y <- y
  } else {
    args$y <- factor(y, levels = c(0, 1))
    args$probability <- TRUE
  }
  # ranger reads `seed = 0` as "no seed" and draws one from R's RNG, so the
  # package-wide default of 0 left this learner -- and with it the default
  # evaluation nuisance -- following whatever state the caller was in. grf and
  # bcf are seeded the other way round: CATEomics pins R's RNG and lets the
  # backend draw its own seed from it. `.fit_base_learner()` has already
  # pinned it, so nothing is passed down and 0 becomes an ordinary value,
  # while a user's own `seed` still wins.
  do.call(ranger::ranger, args)
}

.predict_ranger <- function(model, newdata, type) {
  pred <- stats::predict(model, data = newdata)$predictions
  if (type == "regression") {
    return(as.numeric(pred))
  }
  # probability forest: matrix with columns "0" and "1"
  if (is.matrix(pred) || is.data.frame(pred)) {
    if ("1" %in% colnames(pred)) {
      return(as.numeric(pred[, "1"]))
    }
    return(as.numeric(pred[, ncol(pred)]))
  }
  as.numeric(pred)
}

# ---- xgboost --------------------------------------------------------------

.fit_xgboost <- function(x, y, type, random_state = NULL, ...) {
  # xgboost is an Import since the default evaluation nuisance library carries
  # two boosted members, so this cannot fire through the package's own entry
  # points. Kept as a guard for a caller that reaches the fitter with the
  # namespace unavailable.
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    rlang::abort(
      c("Package 'xgboost' is required for learner = 'xgboost'.",
        i = "Install with: install.packages('xgboost')"),
      class = "cateomics_xgboost_missing"
    )
  }
  args <- list(...)
  # CATEomics defaults. The sklearn-wrapper spellings (`learning_rate` and so
  # on) are accepted alongside the native ones.
  # nrounds 100 / max_depth 4 / eta 0.1 / min_child_weight 10 /
  # colsample_bytree 0.5. The three middle values are the ones
  # SuperLearner::SL.xgboost uses.
  # `eta` is xgboost's own spelling of `learning_rate`, and the allow-list
  # accepts it because xgb.params() lists it. Resolved here, before the
  # extraction below: left in `args` it would reach the same booster list as
  # the `learning_rate` written there, where `learning_rate` wins whichever
  # comes first, so a validated hyperparameter was silently replaced by the
  # default.
  if ("eta" %in% names(args)) {
    if ("learning_rate" %in% names(args)) {
      rlang::abort(
        c("`eta` and `learning_rate` are the same booster setting.",
          i = "Pass one of them."),
        class = "cateomics_duplicate_param"
      )
    }
    args$learning_rate <- args$eta
    args$eta <- NULL
  }
  nrounds       <- if (is.null(args$nrounds)) 100L else args$nrounds
  max_depth     <- if (is.null(args$max_depth)) 4L else args$max_depth
  learning_rate <- if (is.null(args$learning_rate)) 0.1 else args$learning_rate
  args$nrounds       <- NULL
  args$max_depth     <- NULL
  args$learning_rate <- NULL

  if (!is.null(random_state)) {
    set.seed(as.integer(random_state))
  }
  if (type == "regression") {
    objective <- "reg:squarederror"
  } else {
    objective <- "binary:logistic"
  }
  # Split what the caller passed: xgb.train takes a handful of arguments
  # itself and reads every booster hyperparameter out of `params`. Spreading
  # them all at the top level made xgboost warn ("These should be passed as a
  # list to argument 'params'") and is documented upstream to become an error.
  top_lvl  <- intersect(names(args), names(formals(xgboost::xgb.train)))
  train_as <- args[top_lvl]
  booster  <- utils::modifyList(
    list(objective        = objective,
         max_depth        = as.integer(max_depth),
         learning_rate    = as.numeric(learning_rate),
         min_child_weight = 10L,
         colsample_bytree = 0.5,
         verbosity        = 0L),
    args[setdiff(names(args), top_lvl)])
  do.call(xgboost::xgb.train, c(
    list(
      params  = booster,
      data    = xgboost::xgb.DMatrix(data = as.matrix(x), label = as.numeric(y)),
      nrounds = as.integer(nrounds)
    ),
    train_as
  ))
}

.predict_xgboost <- function(model, newdata, type) {
  pred <- stats::predict(model, newdata = xgboost::xgb.DMatrix(data = as.matrix(newdata)))
  as.numeric(pred)
}

# ---- central covariate z-score helpers ------------------------------------
#
# CATEomics ships omics covariates on heterogeneous scales (gene expression vs.
# clinical variables). `.cf_preprocess_input` standardises the assembled X once
# (before fold splitting / dispatch) via these helpers so every scale-sensitive
# backend (base learners, catenets, ganite, cevae) sees a common representation
# and regularisation is not biased toward large-magnitude features. Every
# column is standardised, binary indicators included; glmnet is called with
# `standardize = FALSE` since standardisation is handled here.


#' Compute central covariate z-score statistics (column mean / sd).
#'
#' Every column is standardised, binary indicators included: a penalised
#' learner shrinks a narrow column harder than a wide one, and 0/1 columns are
#' narrower than z-scored continuous ones. Neither the R learners
#' (glmnet runs with standardize = FALSE) nor the Python ones standardise
#' internally, so this is the only place it happens. Constant columns (sd 0 or
#' non-finite) fall back to sd 1. Uses only X, so it is leakage-free and is
#' computed once on the full training matrix.
#' @keywords internal
#' @noRd
.cf_compute_x_zscore <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X_mean <- colMeans(X, na.rm = TRUE)
  X_sd   <- apply(X, 2L, stats::sd, na.rm = TRUE)
  X_sd[!is.finite(X_sd) | X_sd == 0] <- 1
  list(X_mean = X_mean, X_sd = X_sd)
}


#' Apply stored central z-score statistics to a covariate matrix.
#'
#' Column-wise `(x - mean) / sd`. Statistics come from [[.cf_compute_x_zscore]]
#' on the training data; re-applied to newdata for a leakage-safe transform.
#' @keywords internal
#' @noRd
.cf_apply_x_zscore <- function(X, stats) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X <- sweep(X, 2L, stats$X_mean, "-")
  sweep(X, 2L, stats$X_sd, "/")
}



