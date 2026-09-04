#' Predict treatment effects from a cf_model
#'
#' @description
#' Unified \code{predict} API.
#'
#' When \code{newdata} is supplied, the cf_model is dispatched to its
#' method-specific predictor with a uniform set of types. \code{type = "all"}
#' returns a flat named list of \code{ate / ate_ci / tau_hat / tau_ci /
#' mu0_hat / mu1_hat}; individual elements may be \code{NA_real_} when not
#' supported (e.g. CIs on methods without inference). When \code{newdata =
#' NULL}, the stored training-sample fields in \code{object$res} are returned.
#'
#' Ensemble objects from \code{\link{cf_ensemble}} hold no single backend, so
#' \code{newdata} prediction is unsupported; \code{predict()} with no
#' \code{newdata} returns the stored ensemble (OOF) tau.
#'
#' @section Method-class behaviour with newdata:
#' \describe{
#'   \item{Class A (direct potential-outcome output, 7 methods)}{
#'     slearner / tlearner / bcf / cevae / ganite / tarnet / dragonnet —
#'     all of \code{tau / mu0 / mu1 / ate} are produced.
#'   }
#'   \item{Class B' (intermediate nuisance, 2 methods)}{
#'     drlearner / xlearner — produce all four fields but
#'     \code{mu1 - mu0 != tau} by design (Neyman-orthogonal CATE final
#'     stage / propensity-weighted Stage-3 blend).
#'   }
#'   \item{Class C (Robinson residualisation, 1 method)}{
#'     grf reconstructs \code{mu0 = Y.hat - W.hat * tau, mu1 = Y.hat +
#'     (1 - W.hat) * tau}.
#'   }
#' }
#'
#' @param object \code{cf_model} object returned by \code{\link{cf_fit}}.
#' @param newdata Optional new covariate matrix / data.frame. When
#'   \code{NULL} (default), the training-sample predictions stored in
#'   \code{object$res} are returned.
#' @param type One of \code{"all"} (default), \code{"ate"}, \code{"ate_ci"},
#'   \code{"tau"}, \code{"tau_ci"}, \code{"mu0"}, \code{"mu1"},
#'   \code{"summary"}.
#' @param ... Additional arguments forwarded to the per-method predict handler.
#' @return Depending on \code{type}:
#' \describe{
#'   \item{\code{"all"}}{named list \code{$ate / $ate_ci / $tau_hat / $tau_ci /
#'     $mu0_hat / $mu1_hat}, with \code{NA_real_} for unsupported fields.}
#'   \item{\code{"ate"}}{numeric scalar.}
#'   \item{\code{"ate_ci"}}{\code{c(lower, upper)} numeric or
#'     \code{NA_real_}.}
#'   \item{\code{"tau"}}{numeric vector of \eqn{\hat\tau} of length
#'     \code{n_new} or \code{n}.}
#'   \item{\code{"tau_ci"}}{\code{data.frame(lower, upper)} of \code{n_new}
#'     rows or \code{NA_real_}.}
#'   \item{\code{"mu0"}, \code{"mu1"}}{numeric vectors of potential outcomes.}
#'   \item{\code{"summary"}}{stored summary list (newdata not supported).}
#' }
#' @export
#' @examples
#' set.seed(42)
#' n <- 50
#' demo_data <- data.frame(
#'   X = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' demo_data$Y <- 1.5 * demo_data$T + 0.3 * demo_data$X + rnorm(n, sd = 0.5)
#' fit <- cf_fit(data = demo_data, Y = "Y", T = "T", X_baseline = "X",
#'               method = "grf", partition = "none")
#' predict(fit, type = "ate")
predict.cf_model <- function(object,
                             newdata = NULL,
                             type = c("all", "ate", "ate_ci",
                                      "tau", "tau_ci",
                                      "mu0", "mu1", "summary"),
                             ...) {
  if (!inherits(object, "cf_model")) {
    rlang::abort(
      message = "object must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  if (is.null(object$method) || is.null(object$res)) {
    rlang::abort(
      message = "Invalid cf_model structure: missing 'method' or 'res' component",
      class = "cfomics_invalid_result"
    )
  }

  type <- match.arg(type)

  # Ensembles (cf_ensemble) combine per-method OOF tau and hold no single
  # backend, so there is nothing to score newdata against. With newdata = NULL
  # the stored ensemble fields in res are returned via the shared in-sample
  # assembler, exactly as every predict_cf_<method>() does.
  if (identical(object$method, "ensemble")) {
    if (!is.null(newdata)) {
      rlang::abort(
        c("predict(): newdata is not supported for ensemble models.",
          i = "An ensemble has no backend to score new observations.",
          i = "Call predict(object) with no newdata for the stored ensemble tau."),
        class = "cfomics_ensemble_no_newdata"
      )
    }
    return(.cf_insample_predict(object$res, type))
  }

  # Accept SummarizedExperiment / MultiAssayExperiment newdata by extracting the
  # X matrix via `.extract_X_from_bioc_newdata()` before falling into the
  # preprocessing pipeline. data.frame / matrix inputs pass through unchanged.
  if (!is.null(newdata) &&
      inherits(newdata, c("SummarizedExperiment", "MultiAssayExperiment"))) {
    newdata <- .extract_X_from_bioc_newdata(object, newdata)
  }

  # When newdata is supplied, apply the newdata preprocessing pipeline derived
  # from meta. This is a no-op when none of preprocess / factor_levels /
  # x_zscore are populated.
  if (!is.null(newdata)) {
    newdata <- .apply_newdata_preprocessing(object, newdata)
  }

  switch(
    object$method,
    "ganite"    = predict_cf_ganite(object, newdata = newdata, type = type, ...),
    "drlearner" = predict_cf_drlearner(object, newdata = newdata, type = type, ...),
    "cevae"     = predict_cf_cevae(object, newdata = newdata, type = type, ...),
    "grf"       = predict_cf_grf(object, newdata = newdata, type = type, ...),
    "bcf"       = predict_cf_bcf(object, newdata = newdata, type = type, ...),
    "slearner"  = predict_cf_slearner(object, newdata = newdata, type = type, ...),
    "tlearner"  = predict_cf_tlearner(object, newdata = newdata, type = type, ...),
    "xlearner"  = predict_cf_xlearner(object, newdata = newdata, type = type, ...),
    "tarnet"    = predict_cf_tarnet(object, newdata = newdata, type = type, ...),
    "dragonnet" = predict_cf_dragonnet(object, newdata = newdata, type = type, ...),
    stop(sprintf("Unknown method: %s", object$method))
  )
}

# -------------------------------------------------------------------------
# shared helpers (flat assembly + newdata preprocessing)
# -------------------------------------------------------------------------

#' Assemble the canonical predict() return value for a single method
#'
#' \code{type = "all"} returns the flat named list; other types extract
#' the corresponding element. \code{NA_real_} placeholders are used for
#' unsupported fields (e.g. CI when a method does not provide one).
#'
#' @param type Character one of
#'   \code{c("all","ate","ate_ci","tau","tau_ci","mu0","mu1")}.
#' @param ate Scalar.
#' @param ate_ci Length-2 numeric \code{c(lower, upper)} or \code{NA_real_}.
#' @param tau Numeric vector of \eqn{\hat\tau} or \code{NA_real_}.
#' @param tau_ci \code{data.frame(lower, upper)} or \code{NA_real_}.
#' @param mu0_hat,mu1_hat Numeric vectors of counterfactual outcomes or
#'   \code{NA_real_}.
#' @keywords internal
#' @noRd
.assemble_predict_output <- function(type,
                                     ate     = NA_real_,
                                     ate_ci  = NA_real_,
                                     tau     = NA_real_,
                                     tau_ci  = NA_real_,
                                     mu0_hat = NA_real_,
                                     mu1_hat = NA_real_) {
  if (type == "all") {
    return(list(
      ate     = ate,
      ate_ci  = ate_ci,
      tau_hat = tau,
      tau_ci  = tau_ci,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat
    ))
  }
  switch(
    type,
    "ate"    = ate,
    "ate_ci" = ate_ci,
    "tau"    = tau,
    "tau_ci" = tau_ci,
    "mu0"    = mu0_hat,
    "mu1"    = mu1_hat,
    stop(sprintf("Unsupported type '%s' for newdata prediction.", type))
  )
}

#' Assemble the canonical in-sample (newdata = NULL) predict() return value.
#'
#' Single source of truth for the stored-prediction return
#' shape, shared by every \code{predict_cf_<method>()}. Reads the flat
#' \code{res} pillar and maps to the public token set
#' (\code{tau} / \code{tau_ci} / \code{mu0} / \code{mu1}).
#'
#' @param res The \code{res} pillar of a cf_model.
#' @param type One of \code{c("all","ate","ate_ci","tau","tau_ci","mu0","mu1","summary")}.
#' @keywords internal
#' @noRd
.cf_insample_predict <- function(res, type) {
  switch(
    type,
    "all"     = list(
      ate     = res$ate,
      ate_ci  = res$ate_ci %||% NA_real_,
      tau_hat = res$tau_hat,
      tau_ci  = res$tau_ci %||% NA_real_,
      mu0_hat = res$mu0_hat,
      mu1_hat = res$mu1_hat
    ),
    "ate"     = res$ate,
    "ate_ci"  = res$ate_ci %||% NA_real_,
    "tau"     = res$tau_hat,
    "tau_ci"  = res$tau_ci %||% NA_real_,
    "mu0"     = res$mu0_hat,
    "mu1"     = res$mu1_hat,
    "summary" = res$summary,
    stop(sprintf("Unsupported type '%s'.", type))
  )
}

#' Extract X matrix from SummarizedExperiment / MultiAssayExperiment newdata
#'
#' Class-aware extraction of covariates from Bioconductor newdata:
#'
#' * **SE**: pull `colData(newdata)` as a data.frame, materialize the assay
#'   rownames from `assay(newdata, assay_spec$assay_name)` for the requested
#'   layer (resolved via `.identify_omics_columns(object$meta)`), cbind, and
#'   return.
#' * **MAE**: align experiments via `sampleMap()` (same logic as `.extract_mae_assay`),
#'   pull the primary experiment's assay + colData restricted to complete
#'   primary samples, materialize assay rownames, cbind, and return.
#'
#' Layer identification reads `meta$named_args` (named-args path) preferentially
#' and falls back to `meta$assay_spec$assay_cols / assay_pattern` (legacy
#' formula path). When neither is set, only `colData` (no assay) is
#' returned; this matches a fit done with `feature_method = "none"` on
#' clinical-only df input that the user is now trying to apply to an SE.
#'
#' @param object cf_model.
#' @param newdata `SummarizedExperiment` or `MultiAssayExperiment`.
#' @return data.frame ready for `.apply_newdata_preprocessing()`.
#' @keywords internal
#' @noRd
.extract_X_from_bioc_newdata <- function(object, newdata) {
  meta <- object$meta
  layers <- .identify_omics_columns(meta)
  fit_itype <- meta$input_type

  if (inherits(newdata, "SummarizedExperiment") &&
      !inherits(newdata, "MultiAssayExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'SummarizedExperiment' is required for SE newdata.",
          i = "BiocManager::install('SummarizedExperiment')"),
        class = "cfomics_se_package_missing"
      )
    }
    if (!is.null(fit_itype) && !identical(fit_itype, "SE")) {
      rlang::warn(
        c(
          sprintf("Fit input_type was '%s' but newdata is SummarizedExperiment.",
                  fit_itype),
          i = "Continuing with column-name-based extraction; verify the layout matches."
        ),
        class = "cfomics_newdata_input_type_mismatch"
      )
    }
    assay_name <- meta$assay_spec$assay_name %||% 1
    assay_mat  <- SummarizedExperiment::assay(newdata, assay_name)
    col_data   <- as.data.frame(SummarizedExperiment::colData(newdata))
    rownames(col_data) <- colnames(newdata)

    # Materialize assay rows from each layer that has columns (excluding baseline)
    assay_layers <- c("rnaseq", "mirna", "methyl", "protein", "mutation", "cnv")
    for (lyr in assay_layers) {
      cols <- layers[[lyr]]
      if (length(cols) == 0L) next
      missing <- setdiff(cols, rownames(assay_mat))
      if (length(missing) > 0L) {
        rlang::abort(
          c(
            sprintf("newdata SE assay missing %d row(s) used by fit (layer %s):",
                    length(missing), lyr),
            sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))
          ),
          class = "cfomics_newdata_assay_missing_rows"
        )
      }
      for (gene in cols) {
        col_data[[gene]] <- as.numeric(assay_mat[gene, ])
      }
    }
    # Subset / reorder logic depends on whether downstream
    # `.apply_newdata_preprocessing()` will still transform the assay columns.
    # When feature_method = "variance" / "pca" / "pathway", the trained
    # covariate_names are PCs or top-N selections that are NOT yet present
    # in SE; the preprocessing step builds them. So keep the raw layout
    # (baseline + raw assay rownames) and let downstream re-project. For
    # `feature_method = "none"` (or NULL), enforce strict subset against
    # cov_names so the matrix shape immediately matches.
    fm <- (meta$estimator$preprocess %||% meta$preprocess)$feature_method
    needs_refeature <- !is.null(fm) && fm %in% c("variance", "pca", "pathway")
    # Keep batch column visible to .apply_newdata_preprocessing so ComBat
    # newdata refit can read newdata[[batch_column]].
    batch_col <- (meta$estimator$preprocess %||% meta$preprocess)$batch_column
    if (!is.null(batch_col) && batch_col %in% colnames(col_data)) {
      batch_keep <- batch_col
    } else {
      batch_keep <- character(0)
    }
    if (!needs_refeature) {
      cov_names <- meta$covariate_names
      if (!is.null(cov_names)) {
        missing <- setdiff(cov_names, colnames(col_data))
        if (length(missing) > 0L) {
          rlang::abort(
            c(sprintf("newdata SE is missing %d training column(s):",
                      length(missing)),
              sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))),
            class = "cfomics_newdata_missing_columns"
          )
        }
        keep <- unique(c(cov_names, batch_keep))
        col_data <- col_data[, keep, drop = FALSE]
      }
    } else {
      # Strip the outcome / treatment columns so they don't survive as
      # covariates; baseline + raw assay rownames + batch remain.
      drop_cols <- c(meta$outcome_name, meta$treatment_name)
      col_data <- col_data[, setdiff(colnames(col_data), drop_cols), drop = FALSE]
    }
    return(col_data)
  }

  if (inherits(newdata, "MultiAssayExperiment")) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'MultiAssayExperiment' is required for MAE newdata.",
          i = "BiocManager::install('MultiAssayExperiment')"),
        class = "cfomics_mae_package_missing"
      )
    }
    experiments_used <- meta$assay_spec$experiments
    if (is.null(experiments_used)) {
      experiments_used <- names(MultiAssayExperiment::experiments(newdata))
    }
    sample_map <- MultiAssayExperiment::sampleMap(newdata)
    fmap <- sample_map[as.character(sample_map$assay) %in% experiments_used, ,
                       drop = FALSE]
    primary_counts <- table(as.character(fmap$primary))
    complete_primaries <- names(primary_counts)[
      primary_counts == length(experiments_used)
    ]
    if (length(complete_primaries) == 0L) {
      rlang::abort(
        c(
          "0 common primary samples across newdata experiments after sampleMap join.",
          i = "Check MultiAssayExperiment::sampleMap(newdata)."
        ),
        class = "cfomics_newdata_mae_no_overlap"
      )
    }
    col_data <- as.data.frame(
      MultiAssayExperiment::colData(newdata)
    )[complete_primaries, , drop = FALSE]

    # For each layer with assay columns, find the experiment that contains them
    # and materialize. We try in order: rnaseq → primary, others → first matching.
    assay_layers <- c("rnaseq", "mirna", "methyl", "protein", "mutation", "cnv")
    for (lyr in assay_layers) {
      cols <- layers[[lyr]]
      if (length(cols) == 0L) next
      # Try each experiment; first one with all requested rows wins
      matched_exp <- NULL
      for (exp_name in experiments_used) {
        exp_obj <- MultiAssayExperiment::experiments(newdata)[[exp_name]]
        if (inherits(exp_obj, "SummarizedExperiment")) {
          exp_mat <- SummarizedExperiment::assay(exp_obj, 1)
        } else if (is.matrix(exp_obj)) {
          exp_mat <- exp_obj
        } else {
          next
        }
        if (all(cols %in% rownames(exp_mat))) {
          matched_exp <- list(name = exp_name, mat = exp_mat)
          break
        }
      }
      if (is.null(matched_exp)) {
        rlang::abort(
          sprintf("Cannot locate layer '%s' columns in any experiment of newdata MAE.",
                  lyr),
          class = "cfomics_newdata_mae_layer_missing"
        )
      }
      # Align matched_exp$mat columns to complete_primaries via sampleMap
      exp_map <- fmap[as.character(fmap$assay) == matched_exp$name &
                       as.character(fmap$primary) %in% complete_primaries, ,
                     drop = FALSE]
      exp_colnames <- as.character(exp_map$colname)[
        match(complete_primaries, as.character(exp_map$primary))
      ]
      mat_aligned <- matched_exp$mat[cols, exp_colnames, drop = FALSE]
      colnames(mat_aligned) <- complete_primaries
      for (gene in cols) {
        col_data[[gene]] <- as.numeric(mat_aligned[gene, ])
      }
    }

    # X_assays generic list (named list of assay → cols)
    x_assays <- layers$assays
    if (length(x_assays) > 0L) {
      for (a_name in names(x_assays)) {
        cols <- as.character(x_assays[[a_name]])
        if (length(cols) == 0L) next
        if (!(a_name %in% experiments_used)) {
          rlang::abort(
            sprintf("X_assays experiment '%s' not in newdata MAE experiments.",
                    a_name),
            class = "cfomics_newdata_mae_layer_missing"
          )
        }
        exp_obj <- MultiAssayExperiment::experiments(newdata)[[a_name]]
        exp_mat <- if (inherits(exp_obj, "SummarizedExperiment")) {
          SummarizedExperiment::assay(exp_obj, 1)
        } else {
          exp_obj
        }
        exp_map <- fmap[as.character(fmap$assay) == a_name &
                         as.character(fmap$primary) %in% complete_primaries, ,
                       drop = FALSE]
        exp_colnames <- as.character(exp_map$colname)[
          match(complete_primaries, as.character(exp_map$primary))
        ]
        mat_aligned <- exp_mat[cols, exp_colnames, drop = FALSE]
        colnames(mat_aligned) <- complete_primaries
        for (gene in cols) {
          col_data[[gene]] <- as.numeric(mat_aligned[gene, ])
        }
      }
    }
    # Same feature_method-aware subset as SE branch; keep batch column visible
    # for ComBat newdata re-application (parity with SE branch).
    fm_mae <- (meta$estimator$preprocess %||% meta$preprocess)$feature_method
    needs_refeature_mae <- !is.null(fm_mae) && fm_mae %in% c("variance", "pca", "pathway")
    batch_col_mae <- (meta$estimator$preprocess %||% meta$preprocess)$batch_column
    if (!is.null(batch_col_mae) && batch_col_mae %in% colnames(col_data)) {
      batch_keep_mae <- batch_col_mae
    } else {
      batch_keep_mae <- character(0)
    }
    if (!needs_refeature_mae) {
      cov_names <- meta$covariate_names
      if (!is.null(cov_names)) {
        missing <- setdiff(cov_names, colnames(col_data))
        if (length(missing) > 0L) {
          rlang::abort(
            c(sprintf("newdata MAE is missing %d training column(s):",
                      length(missing)),
              sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))),
            class = "cfomics_newdata_missing_columns"
          )
        }
        keep <- unique(c(cov_names, batch_keep_mae))
        col_data <- col_data[, keep, drop = FALSE]
      }
    } else {
      drop_cols <- c(meta$outcome_name, meta$treatment_name)
      col_data <- col_data[, setdiff(colnames(col_data), drop_cols), drop = FALSE]
    }
    return(col_data)
  }

  rlang::abort(
    c(sprintf("Could not extract covariates from Bioconductor newdata of class: %s",
              paste(class(newdata), collapse = "/")),
      i = "newdata must be a SummarizedExperiment or MultiAssayExperiment (data.frame / matrix newdata bypass this path)."),
    class = "cfomics_newdata_unsupported_class"
  )
}


#' Apply fit-time newdata preprocessing (factor levels, z-score, assay
#' extraction, log-CPM, batch, feature_method).
#'
#' When all of \code{meta$preprocess / factor_levels / x_zscore / y_zscore}
#' are \code{NULL} (a clinical-only fit with no preprocessing), this function
#' is an identity transform: the caller's \code{newdata} is returned
#' unchanged. Otherwise it applies the full preprocessing pipeline (assay
#' extraction / log-CPM / ComBat / feature_method / X z-score), driven by the
#' \code{meta} slots saved at fit time.
#'
#' @param object \code{cf_model}.
#' @param newdata Raw user-supplied newdata (matrix / data.frame / SE / MAE).
#' @return Preprocessed newdata, ready for the method-specific predictor.
#' @keywords internal
#' @noRd
.apply_newdata_preprocessing <- function(object, newdata) {
  # preprocess lives under meta$estimator$preprocess; keep a flat-meta fallback.
  prep    <- object$meta$estimator$preprocess %||% object$meta$preprocess
  flevels <- object$meta$factor_levels
  xstats  <- object$meta$x_zscore
  cov_raw <- object$meta$covariate_names_clinical_raw

  # No preprocessing meta -> identity.
  if (is.null(prep) && is.null(flevels) && is.null(xstats) && is.null(cov_raw)) {
    return(newdata)
  }

  X_new <- newdata

  # Step 1.5: factor levels re-apply.
  if (!is.null(flevels) && is.data.frame(X_new)) {
    for (var in names(flevels)) {
      if (var %in% colnames(X_new)) {
        orig_chars <- as.character(X_new[[var]])
        X_new[[var]] <- factor(orig_chars, levels = flevels[[var]])
        if (anyNA(X_new[[var]])) {
          unknown <- setdiff(unique(orig_chars), flevels[[var]])
          rlang::abort(
            sprintf(
              "New factor level in newdata for '%s': %s. Train factor levels: %s.",
              var,
              paste(unknown, collapse = ", "),
              paste(flevels[[var]], collapse = ", ")
            ),
            class = "cfomics_new_factor_level"
          )
        }
      } else {
        rlang::abort(
          sprintf("Missing factor column '%s' in newdata.", var),
          class = "cfomics_missing_factor_column"
        )
      }
    }
  }

  # Factor -> one-hot expansion: expand the train clinical formula RHS
  # (including factor_levels) with model.matrix into the same column structure
  # as training, avoiding `sexF`/`sexM` dummy column-name mismatch errors on the
  # predict path.
  if (!is.null(flevels) && !is.null(cov_raw) && is.data.frame(X_new) &&
      length(cov_raw) > 0L) {
    # Extract only the raw clinical columns present in newdata (factor_levels
    # already applied).
    clinical_cols <- intersect(cov_raw, colnames(X_new))
    if (length(clinical_cols) > 0L) {
      clinical_df <- X_new[, clinical_cols, drop = FALSE]
      # Expand with model.matrix including non-factor columns (numeric passes through).
      X_clinical_dummies <- stats::model.matrix(~ . - 1, data = clinical_df)
      attr(X_clinical_dummies, "assign") <- NULL
      attr(X_clinical_dummies, "contrasts") <- NULL
      # Drop the raw clinical columns from X_new, cbind the dummies.
      other_cols <- setdiff(colnames(X_new), clinical_cols)
      X_new_other <- X_new[, other_cols, drop = FALSE]
      if (ncol(X_new_other) == 0L) {
        X_new <- as.data.frame(X_clinical_dummies, check.names = FALSE)
      } else {
        X_new <- cbind(X_new_other, as.data.frame(X_clinical_dummies, check.names = FALSE))
      }
    }
  }

  # log-CPM + ComBat newdata re-application.
  if (!is.null(prep)) {
    layers <- .identify_omics_columns(object$meta)
    rnaseq_cols_all <- layers$rnaseq
    if (length(rnaseq_cols_all) > 0L && is.data.frame(X_new)) {
      avail <- intersect(rnaseq_cols_all, colnames(X_new))
      if (length(avail) > 0L &&
          (isTRUE(prep$enabled) || isTRUE(prep$batch_corrected))) {
        # gene × sample
        assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))
        rownames(assay_new) <- avail

        # newdata preprocessing mirrors the fit order --
        # log-CPM (Step 2) THEN ComBat on the log values (Step 3).

        # Step 2: newdata log-CPM (per-sample independent; no train stats).
        if (isTRUE(prep$enabled) &&
            identical(prep$preprocess_method, "logcpm")) {
          assay_new <- .cf_newdata_logcpm(assay_new)
        }

        # Step 3: ComBat newdata via a train + new joint refit on the
        # log-CPM values (parity with the fit-time order). ComBat reads
        # batch labels only, not Y, so the joint refit exposes the test
        # batch labels (supplied by the caller) but never the test outcome
        # -- leakage-safe. Applying train-fold gamma/delta one-sidedly would
        # need ComBat internal state that sva does not expose, so the joint
        # refit is the only practical implementation.
        if (isTRUE(prep$batch_corrected) && !is.null(prep$train_assay) &&
            !is.null(prep$train_batch)) {
          if (!requireNamespace("sva", quietly = TRUE)) {
            rlang::abort(
              c("Package 'sva' required for ComBat newdata reapplication.",
                i = "BiocManager::install('sva')"),
              class = "cfomics_sva_missing"
            )
          }
          batch_col_name <- prep$batch_column
          batch_new <- if (!is.null(batch_col_name) &&
                            batch_col_name %in% colnames(X_new)) {
            X_new[[batch_col_name]]
          } else {
            rlang::abort(
              c(sprintf("Batch column '%s' missing in newdata.",
                        batch_col_name %||% "<unknown>"),
                i = "Add the batch column to newdata before predict()."),
              class = "cfomics_newdata_batch_missing"
            )
          }
          # Align train_assay (log-CPM) rows to assay_new rows (intersection).
          common_rows <- intersect(rownames(prep$train_assay),
                                   rownames(assay_new))
          ta <- prep$train_assay[common_rows, , drop = FALSE]
          an <- assay_new[common_rows, , drop = FALSE]
          n_train <- ncol(ta); n_new <- ncol(an)
          batch_joint <- c(as.character(prep$train_batch),
                           as.character(batch_new))
          joint_mat <- cbind(ta, an)
          corrected <- sva::ComBat(dat = as.matrix(joint_mat),
                                   batch = batch_joint, mod = NULL,
                                   par.prior = TRUE, prior.plots = FALSE)
          # Keep only the newdata columns.
          assay_new <- corrected[, (n_train + 1L):(n_train + n_new),
                                 drop = FALSE]
        }

        # Write back into X_new column-by-column
        for (g in rownames(assay_new)) {
          X_new[[g]] <- as.numeric(assay_new[g, ])
        }
      }
    }
  }

  # GSVA pathway newdata re-application. GSVA is per-sample independent
  # (gene-rank based), so we re-run gsva on newdata using the same gene_sets
  # the fit used. The
  # selected pathways (when n_pathways < total_pathways) are taken from
  # `selected_features` so newdata pathway columns match training.
  prep_fm <- prep$feature_method
  if (!is.null(prep_fm) && prep_fm == "pathway") {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      rlang::abort(
        c("Package 'GSVA' required for pathway newdata reapplication.",
          i = "BiocManager::install('GSVA')"),
        class = "cfomics_gsva_missing"
      )
    }
    gene_sets <- prep$gene_sets
    if (is.null(gene_sets)) {
      rlang::abort(
        c("feature_method='pathway' but meta$estimator$preprocess$gene_sets is missing.",
          i = "Refit with cf_fit() to enable pathway newdata."),
        class = "cfomics_corrupt_meta"
      )
    }
    layers <- .identify_omics_columns(object$meta)
    rnaseq_cols <- layers$rnaseq
    avail <- intersect(rnaseq_cols, colnames(X_new))
    if (length(avail) == 0L) {
      rlang::abort(
        c("Pathway newdata: no rnaseq columns found in newdata.",
          i = "Pass newdata with the same assay rownames as the fit."),
        class = "cfomics_newdata_no_rnaseq"
      )
    }
    assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))
    rownames(assay_new) <- avail
    baseline_cols <- if (length(layers$baseline) > 0L) {
      intersect(layers$baseline, colnames(X_new))
    } else {
      setdiff(colnames(X_new), rnaseq_cols)
    }
    X_baseline_new <- X_new[, baseline_cols, drop = FALSE]

    pathway_method <- prep$pathway_method %||% "gsva"
    gsva_out <- if (utils::packageVersion("GSVA") >= "1.50.0") {
      param <- GSVA::gsvaParam(
        exprData = as.matrix(assay_new),
        geneSets = gene_sets,
        kcdf     = "auto"
      )
      GSVA::gsva(param)
    } else {
      GSVA::gsva(as.matrix(assay_new), gene_sets, method = pathway_method)
    }
    sel <- prep$selected_features
    if (!is.null(sel)) {
      have <- intersect(sel, rownames(gsva_out))
      if (length(have) != length(sel)) {
        rlang::abort(
          c(sprintf("Pathway newdata: %d training pathway(s) missing from GSVA output:",
                    length(sel) - length(have)),
            sprintf("  %s", paste(utils::head(setdiff(sel, have), 5L),
                                  collapse = ", "))),
          class = "cfomics_newdata_pathway_missing"
        )
      }
      gsva_out <- gsva_out[sel, , drop = FALSE]
    }
    pathway_mat_new <- as.data.frame(t(gsva_out), check.names = FALSE)
    if (ncol(X_baseline_new) == 0L) {
      X_new <- pathway_mat_new
    } else {
      X_new <- cbind(X_baseline_new, pathway_mat_new)
    }
  }

  # feature_method newdata re-application.
  if (!is.null(prep_fm) && prep_fm %in% c("variance", "pca")) {
    layers <- .identify_omics_columns(object$meta)
    # Only a single rnaseq layer is supported here.
    rnaseq_cols <- layers$rnaseq
    if (length(rnaseq_cols) > 0L && is.data.frame(X_new)) {
      avail <- intersect(rnaseq_cols, colnames(X_new))
      if (length(avail) == 0L) {
        rlang::abort(
          c("feature_method newdata: no rnaseq columns found in newdata.",
            i = "Pass newdata with the same assay rownames as the fit.")
          , class = "cfomics_newdata_no_rnaseq"
        )
      }
      assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))  # gene × sample
      ## Use meta$named_args$X_baseline / meta$assay_spec when available; this
      ## prevents extra colData columns from being treated as baseline (which
      ## would inflate the X_new column count and break grf::predict's strict
      ## column-count check). Fall back to setdiff when no named-args info.
      baseline_cols <- if (length(layers$baseline) > 0L) {
        intersect(layers$baseline, colnames(X_new))
      } else {
        setdiff(colnames(X_new), rnaseq_cols)
      }
      X_baseline_new <- X_new[, baseline_cols, drop = FALSE]

      if (prep_fm == "variance") {
        # Use train-fold selected_names; "variance" means top-N by train
        # variance, so newdata must use the same N rows.
        sel <- prep$selected_features
        if (is.null(sel)) {
          rlang::abort(
            "feature_method='variance' but meta$preprocess$selected_features missing.",
            class = "cfomics_corrupt_meta"
          )
        }
        missing_sel <- setdiff(sel, rownames(assay_new))
        if (length(missing_sel) > 0L) {
          rlang::abort(
            c(sprintf("newdata missing %d selected features:", length(missing_sel)),
              sprintf("  %s", paste(utils::head(missing_sel, 5L), collapse = ", "))),
            class = "cfomics_newdata_selected_features_missing"
          )
        }
        assay_sel <- assay_new[sel, , drop = FALSE]
        X_new_assay <- t(assay_sel)
      } else if (prep_fm == "pca") {
        loadings <- object$meta$pca_loadings
        center <- object$meta$pca_center
        scale_v <- object$meta$pca_scale
        gene_order <- object$meta$pca_gene_names
        if (is.null(loadings) || is.null(gene_order)) {
          rlang::abort(
            "feature_method='pca' but meta$pca_loadings / pca_gene_names missing.",
            class = "cfomics_corrupt_meta"
          )
        }
        missing_g <- setdiff(gene_order, rownames(assay_new))
        if (length(missing_g) > 0L) {
          rlang::abort(
            c(sprintf("newdata PCA: %d training gene(s) missing:", length(missing_g)),
              sprintf("  %s", paste(utils::head(missing_g, 5L), collapse = ", "))),
            class = "cfomics_newdata_pca_genes_missing"
          )
        }
        assay_aligned <- assay_new[gene_order, , drop = FALSE]  # gene × sample
        # Apply train center / scale to per-gene rows (the prcomp(t(X), center=TRUE,
        # scale.=TRUE) standardised samples; here we mirror it column-wise on the
        # transposed sample × gene matrix below).
        X_gene_sample <- t(assay_aligned)  # sample × gene
        if (!is.null(center)) {
          X_gene_sample <- sweep(X_gene_sample, 2, center, "-")
        }
        if (!is.null(scale_v) && !isFALSE(scale_v)) {
          X_gene_sample <- sweep(X_gene_sample, 2, scale_v, "/")
        }
        pc_new <- X_gene_sample %*% loadings  # sample × n_PC
        colnames(pc_new) <- colnames(loadings)
        X_new_assay <- as.data.frame(pc_new, check.names = FALSE)
      }

      # Concatenate baseline + new assay features (matching cf_fit Step 5)
      if (ncol(X_baseline_new) == 0L) {
        X_new <- X_new_assay
      } else {
        X_new <- cbind(X_baseline_new, X_new_assay)
      }
    } else if (prep_fm == "variance" && !is.null(prep$selected_features) &&
               is.matrix(X_new) && !is.null(colnames(X_new))) {
      # df + matrix path (no SE/MAE extraction): simply subset by selected_features
      sel <- prep$selected_features
      have <- intersect(sel, colnames(X_new))
      if (length(have) == length(sel)) {
        baseline_cols <- setdiff(colnames(X_new),
                                 c(sel,
                                   object$meta$assay_spec$assay_cols %||% character(0)))
        X_new <- cbind(X_new[, baseline_cols, drop = FALSE],
                       X_new[, sel,           drop = FALSE])
      }
    }
  }

  # Drop the batch column once ComBat has been re-applied; it
  # was only kept around to feed ComBat and must not survive to the
  # downstream predict_cf_<method> call (otherwise the model sees an extra
  # column that wasn't in training).
  if (!is.null(prep) && !is.null(prep$batch_column) &&
      is.data.frame(X_new) && prep$batch_column %in% colnames(X_new)) {
    X_new[[prep$batch_column]] <- NULL
  }

  # Re-apply the central X z-score to newdata. Align columns to the training
  # covariate names first (drops a treatment column or any extra columns present
  # in newdata and restores training order) so the statistics land on the right
  # columns.
  if (!is.null(xstats)) {
    X_mat <- if (is.data.frame(X_new)) as.matrix(X_new) else X_new
    if (is.matrix(X_mat) && is.numeric(X_mat)) {
      cov_names <- object$meta$covariate_names
      if (!is.null(cov_names) && !is.null(colnames(X_mat)) &&
          all(cov_names %in% colnames(X_mat))) {
        X_mat <- X_mat[, cov_names, drop = FALSE]
      }
      # Standardise only when the column count matches the stored stats length;
      # otherwise leave X_new untouched and let the per-method column check
      # raise the dimension-mismatch error (avoids sweep()'s recycling warning
      # on malformed newdata). No-op on every valid path.
      if (ncol(X_mat) == length(xstats$X_mean)) {
        X_mat <- sweep(X_mat, 2, xstats$X_mean, "-")
        X_mat <- sweep(X_mat, 2, xstats$X_sd,   "/")
        X_new <- X_mat
      }
    }
  }

  X_new
}


#' Standard error for newdata CI requested on a method without newdata CI support
#'
#' Used by s/t/x-learner, cevae, ganite, tarnet, dragonnet, and
#' drlearner with \code{model_final != "linear"}.
#' @keywords internal
#' @noRd
.abort_ite_ci_unsupported <- function(method) {
  rlang::abort(
    c(
      sprintf("CI for newdata is not implemented for method '%s'.", method),
      i = "Supported: grf (IJ variance), bcf (posterior), drlearner (linear final via econml StatsModelsInference)."
    ),
    class = "cfomics_ite_ci_unsupported"
  )
}

# -------------------------------------------------------------------------
# summary / print
# -------------------------------------------------------------------------

#' Get summary statistics from a cf_model
#'
#' @param object cf_model object
#' @param ... Additional arguments (unused)
#' @return A list containing summary statistics including ATE, confidence
#'   intervals (if bootstrap was enabled), and ITE distribution statistics.
#' @export
#' @examples
#' set.seed(123)
#' demo_data <- data.frame(
#'   X = rnorm(50),
#'   T = rbinom(50, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X + rnorm(50)
#' fit <- cf_fit(data = demo_data, Y = "Y", T = "T", X_baseline = "X",
#'               method = "grf", partition = "none")
#' summary(fit)  # ATE + CATE (tau) distribution summary.
summary.cf_model <- function(object, ...) {
  if (!inherits(object, "cf_model")) {
    rlang::abort(
      message = "object must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  if (is.null(object$res) || is.null(object$res$summary)) {
    rlang::abort(
      message = "Invalid cf_model: missing res / summary.",
      class = "cfomics_invalid_result"
    )
  }

  res <- object$res
  s   <- res$summary

  ci_method_label <- res$ci_method %||% object$meta$estimator$ci_method
  if (is.null(ci_method_label)) ci_method_label <- "unavailable"

  out <- list(
    method    = object$method,
    n         = object$meta$n,
    p         = object$meta$p,
    estimand  = object$meta$estimand$type,
    ate       = res$ate,
    ate_ci    = c(s$ate_ci_lower, s$ate_ci_upper),
    ci_method = ci_method_label,
    tau_summary = list(
      mean = s$tau_mean,
      sd   = s$tau_sd,
      iqr  = s$tau_iqr
    )
  )
  class(out) <- "summary.cf_model"
  out
}

#' @export
print.summary.cf_model <- function(x, ...) {
  cat("Causal Inference Model Summary\n")
  cat("==============================\n")
  cat(sprintf("Method:    %s\n", x$method))
  cat(sprintf("Estimand:  %s\n", x$estimand))
  cat(sprintf("N samples: %d\n", x$n))
  cat(sprintf("N covars:  %d\n", x$p))
  cat("\nTreatment Effect Estimates:\n")
  cat(sprintf("  ATE:        %.4f\n", x$ate))
  ci_label <- if (is.null(x$ci_method)) "unavailable" else x$ci_method
  if (anyNA(x$ate_ci)) {
    cat(sprintf("  ATE 95%% CI: NA (%s)\n", ci_label))
  } else {
    cat(sprintf("  ATE 95%% CI: [%.4f, %.4f] (%s)\n",
                x$ate_ci[1], x$ate_ci[2], ci_label))
  }
  cat("\nCATE (tau) Distribution:\n")
  cat(sprintf("  Mean: %.4f\n", x$tau_summary$mean))
  cat(sprintf("  SD:   %.4f\n", x$tau_summary$sd))
  q <- x$tau_summary$iqr
  cat("  Min / Q1 / Median / Q3 / Max:\n")
  cat(sprintf("  %.4f / %.4f / %.4f / %.4f / %.4f\n",
              q$q0, q$q25, q$q50, q$q75, q$q100))
  invisible(x)
}

#' @export
print.cf_model <- function(x, ...) {
  if (!inherits(x, "cf_model")) {
    rlang::abort(
      message = "x must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }
  res <- x$res
  cat(sprintf("<cf_model: %s>\n", x$method))
  cat(sprintf("  N = %d, p = %d\n", x$meta$n, x$meta$p))
  if (!is.null(res$summary)) {
    cat(sprintf("  ATE = %.4f (95%% CI: [%.4f, %.4f])\n",
                res$ate,
                res$summary$ate_ci_lower,
                res$summary$ate_ci_upper))
    cat(sprintf("  CATE (tau): mean = %.4f, sd = %.4f\n",
                res$summary$tau_mean,
                res$summary$tau_sd))
  } else {
    cat(sprintf("  ATE = %.4f\n", res$ate))
  }
  invisible(x)
}

# The package-level `%||%` is defined once in utils_kfold.R.
