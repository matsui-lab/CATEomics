#' Fit a base learner
#'
#' Unified interface to seven base learners used across cfomics methods for
#' propensity score, outcome model, and final stage fits.
#'
#' Supported learners:
#' * `"linear"` / `"logistic"` — aliases. Uses [stats::lm()] for regression
#'   and [stats::glm()] with `family = binomial()` for classification.
#' * `"lasso"` — [glmnet::cv.glmnet()] with `alpha = 1`. Default for
#'   high-dim omics nuisance models (Chernozhukov et al. 2018 DML).
#' * `"ridge"` — [glmnet::cv.glmnet()] with `alpha = 0`.
#' * `"elastic_net"` — [glmnet::cv.glmnet()] with `alpha = 0.5`.
#' * `"rf"` — [ranger::ranger()]; classification mode sets `probability = TRUE`.
#' * `"xgboost"` — [xgboost::xgboost()]; requires
#'   the `xgboost` package as Suggests dependency.
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
#' @param random_state Optional integer seed. When non-`NULL`, applied via
#'   [set.seed()] before `cv.glmnet` fold assignment and passed as `seed =` to
#'   `ranger`.
#' @param ... Additional arguments forwarded to the underlying fit function
#'   (`lm`, `glm`, `cv.glmnet`, or `ranger`).
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

  model <- switch(learner,
    linear       = ,
    logistic     = .fit_linear(x_mat, y_in, type, ...),
    lasso        = .fit_glmnet(x_mat, y_in, type, alpha = 1,   ...),
    ridge        = .fit_glmnet(x_mat, y_in, type, alpha = 0,   ...),
    elastic_net  = .fit_glmnet(x_mat, y_in, type, alpha = 0.5, ...),
    rf           = .fit_ranger(x_mat, y_in, type, random_state = random_state, ...),
    xgboost      = .fit_xgboost(x_mat, y_in, type, random_state = random_state, ...)
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
      class = "cfomics_base_learner_error"
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

#' Validate a learner specification
#'
#' @param learner Scalar character, the learner name.
#' @return Invisibly `TRUE`. Errors otherwise.
#' @keywords internal
#' @noRd
.validate_learner_spec <- function(learner) {
  valid <- c("linear", "logistic", "lasso", "ridge", "elastic_net", "rf", "xgboost")
  if (!is.character(learner) || length(learner) != 1L || !(learner %in% valid)) {
    rlang::abort(
      sprintf(
        "`learner` must be one of: %s (got: %s)",
        paste(shQuote(valid), collapse = ", "),
        paste(deparse(learner), collapse = " ")
      ),
      class = "cfomics_base_learner_error"
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
        class = "cfomics_base_learner_error"
      )
    }
    x_out <- x
  } else if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      rlang::abort(
        "All columns of `x` must be numeric.",
        class = "cfomics_base_learner_error"
      )
    }
    x_out <- as.matrix(x)
  } else {
    rlang::abort(
      "`x` must be a numeric matrix or a data.frame.",
      class = "cfomics_base_learner_error"
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
        class = "cfomics_base_learner_error"
      )
    }
    return(as.numeric(y))
  }
  # classification: accept 0/1 numeric or two-level factor
  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      rlang::abort(
        "For classification, factor `y` must have exactly two levels.",
        class = "cfomics_base_learner_error"
      )
    }
    return(as.numeric(y) - 1L)
  }
  y_num <- suppressWarnings(as.numeric(y))
  if (any(is.na(y_num)) || !all(y_num %in% c(0, 1))) {
    rlang::abort(
      "For classification, `y` must be a 0/1 vector or a two-level factor.",
      class = "cfomics_base_learner_error"
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
      class = "cfomics_base_learner_error"
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
  # single-covariate inputs (typically toy / contract tests). This preserves
  # the user-facing default `"lasso"` for high-dim omics while still working
  # on degenerate-dim test data.
  if (ncol(x) < 2L) {
    return(.fit_linear(x, y, type, ...))
  }
  family <- if (type == "regression") "gaussian" else "binomial"
  # X is standardised centrally in `.cf_preprocess_input` (continuous columns
  # z-scored, binary/one-hot kept native). Disable glmnet's internal
  # standardisation so it does not re-standardise the binary indicators;
  # continuous columns already have unit scale.
  glmnet::cv.glmnet(x = x, y = y, alpha = alpha, family = family,
                    standardize = FALSE, ...)
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
  if (!is.null(random_state) && is.null(args$seed)) {
    args$seed <- as.integer(random_state)
  }
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
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    rlang::abort(
      c("Package 'xgboost' is required for learner = 'xgboost'.",
        i = "Install with: install.packages('xgboost')"),
      class = "cfomics_xgboost_missing"
    )
  }
  args <- list(...)
  # Sensible cfomics defaults (light to avoid GPU/long-running fits in CI).
  # Use xgboost >= 3.x parameter names (`x`, `y`, `learning_rate`).
  nrounds       <- if (is.null(args$nrounds)) 100L else args$nrounds
  max_depth     <- if (is.null(args$max_depth)) 6L else args$max_depth
  learning_rate <- if (is.null(args$learning_rate)) 0.3 else args$learning_rate
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
  do.call(xgboost::xgb.train, c(
    list(
      params = list(objective = objective,
                    max_depth = as.integer(max_depth),
                    learning_rate = as.numeric(learning_rate),
                    verbosity = 0L),
      data = xgboost::xgb.DMatrix(data = as.matrix(x), label = as.numeric(y)),
      nrounds = as.integer(nrounds)
    ),
    args
  ))
}

.predict_xgboost <- function(model, newdata, type) {
  pred <- stats::predict(model, newdata = xgboost::xgb.DMatrix(data = as.matrix(newdata)))
  as.numeric(pred)
}

# ---- central covariate z-score helpers ------------------------------------
#
# cfomics ships omics covariates on heterogeneous scales (gene expression vs.
# clinical variables). `.cf_preprocess_input` standardises the assembled X once
# (before fold splitting / dispatch) via these helpers so every scale-sensitive
# backend (base learners, catenets, ganite, cevae) sees a common representation
# and regularisation is not biased toward large-magnitude features. Binary /
# one-hot columns are kept on their native scale; glmnet is called with
# `standardize = FALSE` since standardisation is handled here.


#' Detect binary / one-hot columns of a numeric matrix (<= 2 distinct finite
#' values). Used by the central X z-score to keep such columns native (mean 0,
#' sd 1 skipped).
#' @keywords internal
#' @noRd
.cfomics_binary_cols <- function(x) {
  apply(x, 2L, function(col) {
    length(unique(col[is.finite(col)])) <= 2L
  })
}


#' Compute central covariate z-score statistics (column mean / sd).
#'
#' Binary / one-hot columns (<= 2 distinct finite values, [[.cfomics_binary_cols]])
#' are left on their native scale (mean 0, sd 1) so 0/1 indicators keep their
#' encoding instead of being standardised. Constant columns (sd 0 / non-finite)
#' also fall back to sd 1. Uses only X (no outcome), so it is a leakage-free
#' covariate transform and can be computed once on the full training matrix.
#' @keywords internal
#' @noRd
.cf_compute_x_zscore <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X_mean <- colMeans(X, na.rm = TRUE)
  X_sd   <- apply(X, 2L, stats::sd, na.rm = TRUE)
  is_bin <- if (ncol(X) > 0L) .cfomics_binary_cols(X) else logical(0)
  X_mean[is_bin] <- 0
  X_sd[is_bin]   <- 1
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



