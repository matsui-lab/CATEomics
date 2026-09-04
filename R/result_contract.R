#' @title Result Contract: cf_model Object Structure
#'
#' @description
#' This file documents the formal structure (contract) for \code{cf_model}
#' objects returned by \code{\link{cf_fit}}. All causal inference methods in
#' cfomics MUST return objects conforming to this structure. (The historical
#' \code{cfomics_result} class was dropped; objects now carry the single class
#' \code{cf_model}.)
#'
#' \strong{Four-pillar layout.} A \code{cf_model} is a flat list with four
#' pillars plus \code{method}:
#'
#' \describe{
#'   \item{method}{Character string identifying the method (e.g. "grf").}
#'   \item{res}{Numeric results (see Res Structure).}
#'   \item{model}{Method-specific fitted backend container (structure varies
#'     by method).}
#'   \item{meta}{Metadata (see Meta Structure).}
#'   \item{data}{List with \code{X_preprocessed}, \code{Y}, \code{T} (the
#'     model-visible training data; used by cf_interpret / cf_compare /
#'     cf_evaluate / cf_rule for nuisance recomputation).}
#'   \item{fold_cache}{Optional: K per-fold fits when
#'     \code{keep_kfold_models = TRUE}.}
#' }
#'
#' @section Res Structure:
#'
#' The \code{res} pillar MUST include:
#'
#' \describe{
#'   \item{tau_hat}{Numeric vector of Conditional Average Treatment Effects
#'     (individual \eqn{\hat\tau_i}), length \code{n}.}
#'   \item{ate}{Numeric scalar, Average Treatment Effect.}
#'   \item{mu0_hat, mu1_hat}{Numeric vectors of \emph{modelled} counterfactual
#'     outcomes under \eqn{do(T=0)} / \eqn{do(T=1)}, length \code{n}.}
#'   \item{summary}{List of summary statistics (see Summary Structure).}
#' }
#'
#' Optional \code{res} fields (populated when available):
#' \code{tau_hat_oof} (K-fold OOF \eqn{\hat\tau}), \code{ate_ci} (length-2
#' numeric), \code{tau_ci} (data.frame \code{lower}/\code{upper}),
#' \code{ci_method} (character label), \code{nuisance} (list
#' \code{mu0_oof}/\code{mu1_oof}/\code{pi_oof}), and method-specific numeric
#' extras (e.g. \code{propensity}, \code{variance}).
#'
#' @section Summary Structure:
#'
#' The \code{res$summary} list MUST contain (IQR form):
#'
#' \describe{
#'   \item{ate}{Numeric, point estimate of ATE.}
#'   \item{ate_ci_lower, ate_ci_upper}{Numeric CI bounds (\code{NA_real_}
#'     when no CI is available).}
#'   \item{tau_mean, tau_sd}{Mean / SD of \eqn{\hat\tau}.}
#'   \item{tau_iqr}{List with \code{q0}, \code{q25}, \code{q50}, \code{q75},
#'     \code{q100} quantiles of \eqn{\hat\tau} (min / Q1 / median / Q3 / max).}
#' }
#'
#' @section Meta Structure:
#'
#' The \code{meta} pillar contains identification fields (\code{formula},
#' \code{n}, \code{p}, \code{outcome_name}, \code{treatment_name},
#' \code{covariate_names}, \code{method}) plus:
#'
#' \describe{
#'   \item{estimand}{List with \code{type} (one of ATE / ATT / ATC / ITE /
#'     CATE) and \code{target_population}.}
#'   \item{provenance}{List with \code{backend}, \code{rng_seed},
#'     \code{software_versions}, \code{device}.}
#'   \item{estimator}{List with \code{bootstrap},
#'     \code{preprocess}, and (when known) \code{ci_method}.}
#' }
#'
#' Preprocessing slots consumed by \code{\link{predict.cf_model}} on newdata
#' (\code{x_zscore}, \code{y_zscore}, \code{factor_levels}, \code{pca_*},
#' \code{input_type}, \code{assay_spec}, \code{covariate_names_clinical_raw},
#' \code{named_args}, \code{propensity_X}, \code{partition}) are retained at
#' \code{meta} top level.
#'
#' @section Estimands:
#'
#' \describe{
#'   \item{ATE}{\eqn{E[Y(1) - Y(0)]}.}
#'   \item{ATT}{\eqn{E[Y(1) - Y(0) | T=1]}.}
#'   \item{ATC}{\eqn{E[Y(1) - Y(0) | T=0]}.}
#'   \item{ITE}{\eqn{Y_i(1) - Y_i(0)}.}
#'   \item{CATE}{\eqn{E[Y(1) - Y(0) | X=x]}.}
#' }
#'
#' @section Causal Assumptions:
#'
#' All methods rely on SUTVA (no interference, single treatment version),
#' Ignorability (\eqn{(Y(0), Y(1)) \perp T | X}), and Positivity
#' (\eqn{0 < P(T=1|X) < 1}). Use \code{\link{cf_diagnostics}} to assess these.
#'
#' @section S3 Methods:
#'
#' \describe{
#'   \item{\code{predict(object, type, ...)}}{Extract predictions (ate, tau,
#'     mu0, mu1).}
#'   \item{\code{summary(object, ...)}}{Display summary with CIs.}
#'   \item{\code{print(object, ...)}}{Print overview of results.}
#'   \item{\code{\link{cf_plot_cate}(object, ...)}}{Visualize the
#'     \eqn{\hat\tau} distribution.}
#' }
#'
#' @section Class Hierarchy:
#'
#' Objects carry class \code{"cf_model"}.
#'
#' @name cfomics_result
#' @aliases cf_model
#' @return This is a documentation-only page describing the structure of the
#'   \code{cfomics_result} (\code{cf_model}) object returned by
#'   \code{\link{cf_fit}}; it does not define a function and has no return value.
#' @seealso \code{\link{cf_fit}}, \code{\link{predict.cf_model}},
#'   \code{\link{summary.cf_model}}
NULL

#' Validate cfomics_result Structure
#'
#' Internal function to validate that a result object conforms to the
#' four-pillar cfomics_result contract.
#'
#' @param x Object to validate
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_cfomics_result <- function(x) {
  if (!inherits(x, "cf_model")) {
    rlang::abort(
      "Object must have class 'cf_model'",
      class = "cfomics_invalid_result"
    )
  }

  # Four-pillar top level. `method` + `res` + `meta` are mandatory; `model` /
  # `data` are validated when present (cf_fit always sets them, but internal
  # K-fold wrappers may build lighter objects).
  required_toplevel <- c("method", "res", "meta")
  missing_toplevel <- setdiff(required_toplevel, names(x))
  if (length(missing_toplevel) > 0) {
    rlang::abort(
      sprintf("Missing required components: %s", paste(missing_toplevel, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }

  if (!is.character(x$method) || length(x$method) != 1) {
    rlang::abort(
      "'method' must be a single character string",
      class = "cfomics_invalid_result"
    )
  }

  # meta is a single flat metadata pillar with `estimand` / `provenance` /
  # `estimator` buckets (diagnostics are on-demand via cf_diagnostics).
  required_meta <- c("formula", "n", "p", "outcome_name", "treatment_name",
                     "covariate_names", "estimand",
                     "provenance", "estimator")
  missing_meta <- setdiff(required_meta, names(x$meta))
  if (length(missing_meta) > 0) {
    rlang::abort(
      sprintf("Missing required meta fields: %s", paste(missing_meta, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }

  if (!is.list(x$meta$estimand) || is.null(x$meta$estimand$type)) {
    rlang::abort(
      "`meta$estimand` must be a list with field `type`.",
      class = "cfomics_invalid_result"
    )
  }
  valid_estimands <- c("ATE", "ATT", "ATC", "ITE", "CATE")
  if (!x$meta$estimand$type %in% valid_estimands) {
    rlang::abort(
      sprintf("`meta$estimand$type` must be one of: %s (got: %s)",
              paste(valid_estimands, collapse = ", "), x$meta$estimand$type),
      class = "cfomics_invalid_result"
    )
  }

  # Results live in the flat `res` pillar.
  res <- x$res
  if (is.null(res) || !is.list(res)) {
    rlang::abort(
      "`res` is NULL or not a list; all cfomics methods must return a res pillar.",
      class = "cfomics_invalid_result"
    )
  }

  required_res <- c("tau_hat", "ate", "mu0_hat", "mu1_hat", "summary")
  missing_res <- setdiff(required_res, names(res))
  if (length(missing_res) > 0) {
    rlang::abort(
      sprintf("Missing required res fields: %s",
              paste(missing_res, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }

  n <- x$meta$n
  .is_num_or_all_na <- function(v, label, expected_len = NULL) {
    if (!(is.numeric(v) || all(is.na(v)))) {
      rlang::abort(
        sprintf("`res$%s` must be numeric or all NA.", label),
        class = "cfomics_invalid_result"
      )
    }
    if (!is.null(expected_len) && length(v) != expected_len) {
      rlang::abort(
        sprintf("`res$%s` must have length n = %d (got %d).",
                label, expected_len, length(v)),
        class = "cfomics_invalid_result"
      )
    }
  }
  .is_num_or_all_na(res$tau_hat, "tau_hat", expected_len = n)
  .is_num_or_all_na(res$mu0_hat, "mu0_hat", expected_len = n)
  .is_num_or_all_na(res$mu1_hat, "mu1_hat", expected_len = n)

  if (!((is.numeric(res$ate) && length(res$ate) == 1L) || is.na(res$ate))) {
    rlang::abort(
      "`res$ate` must be a single numeric value or NA.",
      class = "cfomics_invalid_result"
    )
  }

  # Estimand <-> tau_hat consistency.
  if (x$meta$estimand$type %in% c("ITE", "CATE")) {
    if (all(is.na(res$tau_hat))) {
      rlang::abort(
        sprintf("Method '%s' reports estimand = '%s' but returned all-NA tau_hat.",
                x$method, x$meta$estimand$type),
        class = "cfomics_invalid_result"
      )
    }
  }

  required_summary <- c(
    "ate", "ate_ci_lower", "ate_ci_upper",
    "tau_mean", "tau_sd", "tau_iqr"
  )
  missing_summary <- setdiff(required_summary, names(res$summary))
  if (length(missing_summary) > 0) {
    rlang::abort(
      sprintf("Missing required summary fields: %s",
              paste(missing_summary, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }
  if (!is.list(res$summary$tau_iqr) ||
      !all(c("q0", "q25", "q50", "q75", "q100") %in% names(res$summary$tau_iqr))) {
    rlang::abort(
      "`res$summary$tau_iqr` must be a list with q0/q25/q50/q75/q100.",
      class = "cfomics_invalid_result"
    )
  }

  # ate_ci / tau_ci unified location (both optional).
  if (!is.null(res$ate_ci)) {
    if (!is.numeric(res$ate_ci) || length(res$ate_ci) != 2L ||
        anyNA(res$ate_ci) || res$ate_ci[1] > res$ate_ci[2]) {
      rlang::abort(
        "`res$ate_ci` must be numeric length-2 with ate_ci[1] <= ate_ci[2].",
        class = "cfomics_invalid_result"
      )
    }
  }
  if (!is.null(res$tau_ci)) {
    if (!is.data.frame(res$tau_ci) ||
        !all(c("lower", "upper") %in% colnames(res$tau_ci)) ||
        nrow(res$tau_ci) != length(res$tau_hat)) {
      rlang::abort(
        "`res$tau_ci` must be a data.frame with columns 'lower'/'upper' and nrow == length(tau_hat).",
        class = "cfomics_invalid_result"
      )
    }
  }

  invisible(TRUE)
}

#' Build the standard res$summary list (IQR form)
#'
#' Centralised so every method + cf_fit produce an identical summary shape.
#' Called by cf_fit after lifting the per-method res into the flat pillar.
#'
#' @param ate Numeric scalar ATE.
#' @param ate_ci Length-2 numeric \code{c(lower, upper)} or NULL.
#' @param tau_hat Numeric vector of \eqn{\hat\tau}.
#' @return A list with ate / ate_ci_lower / ate_ci_upper / tau_mean / tau_sd /
#'   tau_iqr(q0,q25,q50,q75,q100).
#' @keywords internal
#' @noRd
.cf_build_summary <- function(ate, ate_ci = NULL, tau_hat = NA_real_) {
  ci_lower <- if (!is.null(ate_ci) && length(ate_ci) == 2L) ate_ci[1] else NA_real_
  ci_upper <- if (!is.null(ate_ci) && length(ate_ci) == 2L) ate_ci[2] else NA_real_

  has_tau <- !all(is.na(tau_hat))
  if (has_tau) {
    tv <- as.numeric(tau_hat)
    qs <- stats::quantile(tv, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
    tau_iqr <- list(q0 = qs[1], q25 = qs[2], q50 = qs[3], q75 = qs[4], q100 = qs[5])
    tau_mean <- mean(tv)
    tau_sd   <- stats::sd(tv)
  } else {
    tau_iqr <- list(q0 = NA_real_, q25 = NA_real_, q50 = NA_real_,
                    q75 = NA_real_, q100 = NA_real_)
    tau_mean <- NA_real_
    tau_sd   <- NA_real_
  }

  # Standard-error of the ATE, recovered from the 95% CI when available
  # (used by print/summary). NA when no CI is available.
  ate_se <- if (is.finite(ci_lower) && is.finite(ci_upper)) {
    (ci_upper - ci_lower) / (2 * 1.959964)
  } else {
    NA_real_
  }

  list(
    ate          = as.numeric(ate),
    ate_ci_lower = as.numeric(ci_lower),
    ate_ci_upper = as.numeric(ci_upper),
    ate_se       = ate_se,
    tau_mean     = tau_mean,
    tau_sd       = tau_sd,
    tau_iqr      = tau_iqr
  )
}

#' Create Standard Meta Object
#'
#' Internal helper to create a standardized \code{meta} pillar for
#' \code{cfomics_result} (a single flat metadata pillar).
#'
#' @param formula The formula used for fitting
#' @param n Number of observations
#' @param p Number of covariates
#' @param outcome_name Name of outcome variable
#' @param treatment_name Name of treatment variable
#' @param covariate_names Character vector of covariate names
#' @param method Character string of method name
#' @param estimand Character, the estimand being estimated (default "ATE")
#' @param target_population Character, target population (default "all")
#' @param bootstrap Logical, whether bootstrap was used
#' @param random_state Random seed used
#' @param preprocess List describing preprocessing steps
#' @param pca_loadings PCA gene x n_PC rotation matrix (feature_method = "pca")
#' @param input_type "data.frame" / "SE" / "MAE"
#' @param assay_spec list of input-path arguments
#' @param factor_levels Named list of factor column train levels (newdata predict)
#' @param x_zscore list(X_mean, X_sd) for X z-score inverse
#' @param y_zscore list(Y_mean, Y_sd) for Y inverse transform (NN methods)
#' @param covariate_names_clinical_raw raw clinical column names (pre one-hot)
#' @param named_args structured named-args bag (NULL for legacy formula path)
#' @param propensity_X PS-specific covariate subset
#' @param pca_center,pca_scale,pca_gene_names PCA newdata projection stats
#' @return A list with the standardized meta structure
#' @keywords internal
create_cf_meta <- function(formula,
                           n,
                           p,
                           outcome_name,
                           treatment_name,
                           covariate_names,
                           method,
                           estimand = "ATE",
                           target_population = "all",
                           bootstrap = FALSE,
                           random_state = NULL,
                           preprocess = NULL,
                           pca_loadings = NULL,
                           input_type = NULL,
                           assay_spec = NULL,
                           factor_levels = NULL,
                           x_zscore = NULL,
                           y_zscore = NULL,
                           covariate_names_clinical_raw = NULL,
                           named_args = NULL,
                           propensity_X = NULL,
                           pca_center = NULL,
                           pca_scale = NULL,
                           pca_gene_names = NULL) {

  # Capture software versions
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  pkg_versions <- list(
    cfomics = tryCatch(as.character(utils::packageVersion("cfomics")),
                       error = function(e) NA_character_),
    reticulate = tryCatch(as.character(utils::packageVersion("reticulate")),
                          error = function(e) NA_character_)
  )

  method_pkg <- switch(method,
    grf       = "grf",
    bcf       = "bcf",
    NULL
  )
  if (!is.null(method_pkg)) {
    pkg_versions[[method_pkg]] <- tryCatch(
      as.character(utils::packageVersion(method_pkg)),
      error = function(e) NA_character_
    )
  }

  list(
    formula = formula,
    n = as.integer(n),
    p = as.integer(p),
    outcome_name = outcome_name,
    treatment_name = treatment_name,
    covariate_names = covariate_names,
    method = method,

    # ----- Bucket 1: estimand (causal target). Assumptions are documented in
    # roxygen rather than stored on every object. --
    estimand = list(
      type              = estimand,
      target_population = target_population
    ),

    # ----- Bucket 2: provenance (reproducibility metadata) -----
    provenance = list(
      backend           = .infer_backend(method),
      rng_seed          = random_state,
      software_versions = list(
        R         = r_version,
        cfomics   = pkg_versions$cfomics,
        timestamp = Sys.time(),
        packages  = pkg_versions
      )
    ),

    # ----- Bucket 3: estimator (method-specific configuration) -----
    estimator = list(
      bootstrap  = bootstrap,
      preprocess = preprocess
    ),

    # ----- named-args facade record (top-level for easy access) -----
    named_args   = named_args,
    propensity_X = propensity_X,

    # ----- preprocessing slots consumed by cf_predict newdata pipeline -----
    # (meta$rowdata is intentionally not stored: the full per-gene rowData
    # table was write-only and inflated every model by ~20k rows.)
    pca_loadings   = pca_loadings,
    pca_center     = pca_center,
    pca_scale      = pca_scale,
    pca_gene_names = pca_gene_names,
    input_type     = input_type,
    assay_spec     = assay_spec,
    factor_levels  = factor_levels,
    x_zscore       = x_zscore,
    y_zscore       = y_zscore,
    covariate_names_clinical_raw = covariate_names_clinical_raw
  )
}

#' Infer the implementation backend used by a method (provenance bucket).
#' @keywords internal
#' @noRd
.infer_backend <- function(method) {
  switch(method,
    grf       = "grf (R)",
    bcf       = "bcf (R)",
    slearner  = ,
    tlearner  = ,
    xlearner  = "econml.metalearners (Python)",
    drlearner = "econml.dr (Python)",
    ganite    = "PyPI ganite 0.1.2 (PyTorch, Yoon 2018 official)",
    cevae     = "pyro-ppl (Python)",
    tarnet    = ,
    dragonnet = "catenets (Python)",
    "unknown"
  )
}
