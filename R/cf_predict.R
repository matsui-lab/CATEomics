#' Predict treatment effects from a cf_model
#'
#' @description
#' Unified \code{predict} API.
#'
#' When \code{newdata} is supplied, the cf_model is dispatched to its
#' method-specific predictor with a uniform set of types. \code{type = "all"}
#' returns a flat named list of \code{ate / ate_ci / tau_hat / tau_ci /
#' mu0_hat / mu1_hat / ps / ps_type}. A per-row quantity a method does not
#' produce is \code{NA_real_}; an interval it does not compute is
#' \code{NULL}. When \code{newdata = NULL}, the stored training-sample fields
#' in \code{object$res} are returned as they are, so an absent interval is
#' \code{NULL} there too.
#'
#' Ensemble objects from \code{\link{cf_ensemble}} hold no backend of their
#' own. With \code{newdata} they score each weighted member and combine the
#' member effects with the stored weights; a member that cannot score
#' \code{newdata} is dropped with a warning and its weight is spread over the
#' rest. \code{mu0_hat} / \code{mu1_hat} / \code{ps} are \code{NA} on
#' \code{newdata}: in-sample the counterfactual means come from the shared
#' evaluation nuisance, which is stored as out-of-fold prediction vectors and
#' so has no value for rows it was not fitted on. The propensity is \code{NA}
#' either way, because an ensemble's \code{res} carries none; the shared
#' nuisance's propensity is in \code{res$nuisance$pi_oof}, which is what
#' \code{\link{cf_diagnostics}} reads. \code{predict()} with no
#' \code{newdata} returns the stored ensemble (OOF) tau.
#'
#' @section Method-class behaviour with newdata:
#' \describe{
#'   \item{Class A (direct potential-outcome output, 7 methods)}{
#'     slearner / tlearner / bcf / cevae / ganite / tarnet / dragonnet —
#'     all of \code{tau / mu0 / mu1 / ate} are produced.
#'   }
#'   \item{Class B' (intermediate nuisance, 2 methods)}{
#'     drlearner / xlearner — produce all four fields, and
#'     \code{mu1 - mu0} does not equal \code{tau}: tau comes from the
#'     orthogonal final stage (drlearner) or the propensity-weighted
#'     Stage-3 blend (xlearner), not from the difference of the two means.
#'   }
#'   \item{Class C (Robinson residualisation, 1 method)}{
#'     grf reconstructs \code{mu0 = Y.hat - W.hat * tau, mu1 = Y.hat +
#'     (1 - W.hat) * tau}.
#'   }
#' }
#'
#' @param object \code{cf_model} object returned by \code{\link{cf_fit}}.
#' @section The estimate a fitted model carries:
#' \code{res$tau_hat} is what the fit predicts for the rows it was given, and
#' \code{res$tau_hat_oof} is what models fitted without each row predict for
#' it; the comparison, the selection and every downstream function read the
#' second. For \code{grf}, \code{res$tau_hat} holds out-of-bag predictions,
#' made with the trees that did not use that row, so handing the same rows
#' back as \code{newdata} returns the whole forest's prediction instead and
#' the two differ. For every other method they are the same numbers.
#'
#' @param newdata Optional new covariate matrix / data.frame. When
#'   \code{NULL} (default), the training-sample predictions stored in
#'   \code{object$res} are returned. It must carry every omics feature the
#'   fit used, not only the selected ones: the library size, the batch means
#'   and the pathway ranks are taken over the features present, so a missing
#'   one moves the rest. A missing feature is an error.
#' @param type One of \code{"all"} (default), \code{"ate"}, \code{"ate_ci"},
#'   \code{"tau"}, \code{"tau_ci"}, \code{"mu0"}, \code{"mu1"}, \code{"ps"},
#'   \code{"summary"}.
#' @param ... Additional arguments forwarded to the per-method predict handler.
#' @return Depending on \code{type}:
#' \describe{
#'   \item{\code{"all"}}{named list \code{$ate / $ate_ci / $tau_hat / $tau_ci /
#'     $mu0_hat / $mu1_hat / $ps / $ps_type}. A per-row quantity the fitted
#'     method does not produce is \code{NA_real_}; an interval that was not
#'     computed is \code{NULL}, so the element is absent from the list.}
#'   \item{\code{"ate"}}{numeric scalar.}
#'   \item{\code{"ate_ci"}}{\code{c(lower, upper)} numeric, or \code{NULL}
#'     when no interval was computed.}
#'   \item{\code{"tau"}}{numeric vector of \eqn{\hat\tau} of length
#'     \code{n_new} or \code{n}.}
#'   \item{\code{"tau_ci"}}{\code{data.frame(lower, upper)} of \code{n_new}
#'     rows, or \code{NULL} when no interval was computed.}
#'   \item{\code{"mu0"}, \code{"mu1"}}{numeric vectors of potential outcomes.}
#'   \item{\code{"ps"}}{numeric vector of propensity scores, or
#'     \code{NA_real_} for methods that estimate none (slearner, tlearner,
#'     tarnet, dragonnet, cevae, ganite). \code{$ps_type} in the \code{"all"} list records
#'     where the scores come from: \code{"oof"} (assembled from the held-out
#'     fold predictions), \code{"oob"}, \code{"cross-fit (approx)"}
#'     (averaged over the cross-fit propensity models),
#'     \code{"in-sample"}, \code{"out-of-sample"} or \code{"known"}
#'     (supplied via \code{propensity_score}).}
#'   \item{\code{"summary"}}{stored summary list; naming \code{newdata}
#'     together with this type is an error.}
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
                                      "mu0", "mu1", "ps", "summary"),
                             ...) {
  if (!inherits(object, "cf_model")) {
    rlang::abort(
      message = "object must be a cf_model from cf_fit()",
      class = "cateomics_invalid_input"
    )
  }

  if (is.null(object$method) || is.null(object$res)) {
    rlang::abort(
      message = "Invalid cf_model structure: missing 'method' or 'res' component",
      class = "cateomics_invalid_result"
    )
  }

  type <- match.arg(type)

  # Ensembles (cf_ensemble) combine per-method OOF tau and hold no backend of
  # their own. With newdata they delegate to the members they weight, each of
  # which reconstructs its own design matrix -- so this returns before the
  # container extraction and preprocessing below, and `newdata` reaches the
  # members untouched. With newdata = NULL the stored ensemble fields in res
  # are returned via the shared in-sample assembler, exactly as every
  # predict_cf_<method>() does.
  if (identical(object$method, "ensemble")) {
    if (!is.null(newdata)) {
      return(.cf_ensemble_predict_newdata(object, newdata, type))
    }
    return(.cf_insample_predict(object$res, type))
  }

  # Accept SummarizedExperiment / MultiAssayExperiment newdata by extracting the
  # X matrix via `.extract_X_from_bioc_newdata()` before falling into the
  # preprocessing pipeline. data.frame / matrix inputs pass through unchanged.
  if (!is.null(newdata) &&
      inherits(newdata, c("SummarizedExperiment", "MultiAssayExperiment"))) {
    newdata <- .extract_X_from_bioc_newdata(object, newdata)
  } else if (!is.null(newdata)) {
    newdata <- .cf_newdata_training_columns(object, newdata)
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
#' the corresponding element. \code{NA_real_} placeholders stand for the
#' fields the fitted method does not produce (e.g. a CI for a method that
#' does not return one).
#'
#' @param type Character one of
#'   \code{c("all","ate","ate_ci","tau","tau_ci","mu0","mu1","ps")}.
#' @param ate Scalar.
#' @param ate_ci Length-2 numeric \code{c(lower, upper)}, or \code{NULL}
#'   when no interval was computed.
#' @param tau Numeric vector of \eqn{\hat\tau} or \code{NA_real_}.
#' @param tau_ci \code{data.frame(lower, upper)}, or \code{NULL} when no
#'   interval was computed.
#' @param mu0_hat,mu1_hat Numeric vectors of counterfactual outcomes or
#'   \code{NA_real_}.
#' @param ps Numeric vector of propensity scores or \code{NA_real_} when the
#'   method estimates none. The field is always present so that the return
#'   shape does not depend on the method.
#' @param ps_type Character tag for the provenance of \code{ps}
#'   (\code{"oof"} / \code{"oob"} / \code{"cross-fit (approx)"} /
#'   \code{"in-sample"} / \code{"out-of-sample"} / \code{"known"}), or
#'   \code{NA_character_}.
#' @keywords internal
#' @noRd
.assemble_predict_output <- function(type,
                                     ate     = NA_real_,
                                     ate_ci  = NULL,
                                     tau     = NA_real_,
                                     tau_ci  = NULL,
                                     mu0_hat = NA_real_,
                                     mu1_hat = NA_real_,
                                     ps      = NA_real_,
                                     ps_type = NA_character_) {
  if (type == "all") {
    return(list(
      ate     = ate,
      ate_ci  = ate_ci,
      tau_hat = tau,
      tau_ci  = tau_ci,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat,
      ps      = ps,
      ps_type = ps_type
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
    "ps"     = ps,
    stop(sprintf("Unsupported type '%s' for newdata prediction.", type))
  )
}

#' The stored propensity, padded to the length of tau-hat
#'
#' A method that estimates none must still return one value per training
#' observation, so the stored-prediction path gives `predict(fit, type = "ps")`
#' the same shape whatever the method.
#'
#' @param res The `res` pillar of a cf_model.
#' @return Numeric vector as long as `res$tau_hat_oof` when that is present,
#'   otherwise `res$tau_hat` (at least 1).
#' @keywords internal
#' @noRd
.cf_stored_ps <- function(res) {
  n <- length(res$tau_hat_oof %||% res$tau_hat)
  ps <- res$ps
  if (is.null(ps) || length(ps) != n) return(rep(NA_real_, max(n, 1L)))
  as.numeric(ps)
}

#' Check an extracted propensity vector before it is published
#'
#' Reading a propensity out of a backend means reaching into its fitted
#' attributes, so a wrong length, an NA or a value outside (0, 1) means the
#' attribute was not what the extractor assumed. Return NULL in that case: a
#' missing propensity is recoverable downstream, a wrong one is not.
#'
#' @param ps Candidate numeric vector (or NULL).
#' @param n Expected length.
#' @return The validated numeric vector, or NULL.
#' @keywords internal
#' @noRd
.cf_validate_ps <- function(ps, n) {
  if (is.null(ps)) return(NULL)
  ps <- suppressWarnings(as.numeric(ps))
  if (length(ps) != n || anyNA(ps)) return(NULL)
  if (any(ps <= 0 | ps >= 1)) return(NULL)
  ps
}

#' Assemble the canonical in-sample (newdata = NULL) predict() return value.
#'
#' Single source of truth for the stored-prediction return
#' shape, shared by every \code{predict_cf_<method>()}. Reads the flat
#' \code{res} pillar and maps to the public token set
#' (\code{tau} / \code{tau_ci} / \code{mu0} / \code{mu1}).
#'
#' @param res The \code{res} pillar of a cf_model.
#' @param type One of
#'   \code{c("all","ate","ate_ci","tau","tau_ci","mu0","mu1","ps","summary")}.
#' @keywords internal
#' @noRd
.cf_insample_predict <- function(res, type) {
  switch(
    type,
    "all"     = list(
      ate     = res$ate,
      ate_ci  = res$ate_ci,
      tau_hat = res$tau_hat %||% res$tau_hat_oof,
      tau_ci  = res$tau_ci,
      mu0_hat = res$mu0_hat,
      mu1_hat = res$mu1_hat,
      ps      = .cf_stored_ps(res),
      ps_type = res$ps_type %||% NA_character_
    ),
    "ate"     = res$ate,
    "ate_ci"  = res$ate_ci,
    "tau"     = res$tau_hat %||% res$tau_hat_oof,
    "tau_ci"  = res$tau_ci,
    "mu0"     = res$mu0_hat,
    "mu1"     = res$mu1_hat,
    "ps"      = .cf_stored_ps(res),
    "summary" = res$summary,
    stop(sprintf("Unsupported type '%s'.", type))
  )
}

#' The assay a layer's name refers to
#'
#' A layer taken from a SummarizedExperiment is named after its assay. When the
#' assays have no names the label is a position, and `assay()` has to be given
#' that as a number -- its character method would look for an assay literally
#' called "1".
#'
#' @param lab The layer's name, or NULL when the fit recorded no layer.
#' @param se The container to read the assay from.
#' @return A name or a position, usable as `assay(se, <this>)`.
#' @keywords internal
#' @noRd
.cf_layer_assay_name <- function(lab, se) {
  if (is.null(lab) || !nzchar(lab)) return(1L)
  real <- tryCatch(SummarizedExperiment::assayNames(se), error = function(e) NULL)
  if (!is.null(real) && lab %in% real) return(lab)
  pos <- suppressWarnings(as.integer(lab))
  if (is.na(pos)) lab else pos
}

#' Extract X matrix from SummarizedExperiment / MultiAssayExperiment newdata
#'
#' Class-aware extraction of covariates from Bioconductor newdata:
#'
#' * **SE**: pull `colData(newdata)` as a data.frame, materialize the assay
#'   rownames of the layer the fit used (resolved via
#'   `.identify_omics_columns(object$meta)`, whose name is the assay), cbind,
#'   and return.
#' * **MAE**: align experiments via `sampleMap()` (same logic as
#'   `.extract_mae_assay`), then rebuild every experiment the fit used, layer
#'   by layer -- each with its own normalisation, ComBat and feature selection
#'   or PCA projection -- under the `<experiment>_` prefixed names, add
#'   colData restricted to complete primary samples, cbind, and return.
#'
#' Layer identification reads the per-layer records the fit stored
#' (`preprocess$layers`). With no layer recorded, only `colData` (no assay) is
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
        class = "cateomics_se_package_missing"
      )
    }
    if (!is.null(fit_itype) && !identical(fit_itype, "SE")) {
      # Clinical columns carry over by name, so a fit with no omics layer reads
      # a SummarizedExperiment unchanged. A layer is named after the assay it
      # came from, though, and a data.frame fit names its layer "omics": there
      # is no assay of that name to read, and going on would fail inside
      # `assay()` with a message about a subscript.
      lyr_name  <- names(layers$layers)[1L]
      se_assays <- tryCatch(SummarizedExperiment::assayNames(newdata),
                            error = function(e) character(0))
      if (length(layers$layers) > 0L &&
          !is.null(lyr_name) && nzchar(lyr_name) &&
          !(lyr_name %in% se_assays) && is.na(suppressWarnings(as.integer(lyr_name)))) {
        rlang::abort(
          c(sprintf("The fit read a '%s' and its omics layer is called '%s', which newdata has no assay for.",
                    fit_itype, lyr_name),
            x = sprintf("newdata has assay(s): %s.",
                        if (length(se_assays)) paste(se_assays, collapse = ", ") else "none named"),
            i = "Score the same kind of container the fit was given, or supply newdata as a data.frame with the same column names."),
          class = "cateomics_newdata_layer_not_in_assays"
        )
      }
      rlang::warn(
        c(
          sprintf("Fit input_type was '%s' but newdata is SummarizedExperiment.",
                  fit_itype),
          i = "Continuing with column-name-based extraction; verify the layout matches."
        ),
        class = "cateomics_newdata_input_type_mismatch"
      )
    }
    # The layer's name is the assay it was taken from. An unnamed assay is
    # labelled by its position, and `assay()` needs that as a number: the
    # character method would look for an assay literally called "1".
    assay_name <- .cf_layer_assay_name(names(layers$layers)[1L], newdata)
    assay_mat  <- SummarizedExperiment::assay(newdata, assay_name)
    used_rows  <- intersect(unlist(layers$layers, use.names = FALSE),
                            rownames(assay_mat))
    if (length(used_rows) > 0L) {
      .cf_check_assay_complete(assay_mat[used_rows, , drop = FALSE],
                               what = "The newdata assay features")
    }
    .cf_check_unique_samples(colnames(newdata),
                             "The newdata SummarizedExperiment")
    # `optional = TRUE` so a non-syntactic column keeps the name the fit
    # recorded; the default would rewrite it and lose the covariate.
    col_data   <- as.data.frame(SummarizedExperiment::colData(newdata),
                                optional = TRUE)
    rownames(col_data) <- colnames(newdata)

    # Materialize assay rows for each omics layer that named its features.
    for (lyr in names(layers$layers)) {
      cols <- layers$layers[[lyr]]
      if (length(cols) == 0L) next
      missing <- setdiff(cols, rownames(assay_mat))
      if (length(missing) > 0L) {
        rlang::abort(
          c(
            sprintf("newdata SE assay missing %d row(s) used by fit (layer %s):",
                    length(missing), lyr),
            sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))
          ),
          class = "cateomics_newdata_assay_missing_rows"
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
    # Re-featuring is only something the preprocessing can do when there is a
    # layer to re-feature. A fit with clinical columns alone records a
    # `feature_method` all the same, and reading it on its own would keep every
    # colData column -- including the outcome, the treatment and any factor the
    # fit never saw -- and hand them to the backend.
    needs_refeature <- !is.null(fm) && fm %in% c("variance", "pca", "pathway") &&
      length(layers$layers) > 0L
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
        # Clinical columns are kept under their raw names, because a factor
        # is one column here and becomes its training dummies only in
        # `.apply_newdata_preprocessing()`.
        raw_clinical <- meta$named_args$X_baseline
        want <- if (is.null(raw_clinical)) cov_names else
          unique(c(raw_clinical, intersect(cov_names, colnames(col_data))))
        missing <- setdiff(want, colnames(col_data))
        if (length(missing) > 0L) {
          rlang::abort(
            c(sprintf("newdata SE is missing %d training column(s):",
                      length(missing)),
              sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))),
            class = "cateomics_newdata_missing_columns"
          )
        }
        keep <- unique(c(want, batch_keep))
        col_data <- col_data[, keep, drop = FALSE]
      }
    } else {
      # Strip the outcome / treatment columns so they don't survive as
      # covariates; baseline + raw assay rownames + batch remain.
      drop_cols <- c(meta$named_args$Y, meta$named_args$T)
      col_data <- col_data[, setdiff(colnames(col_data), drop_cols), drop = FALSE]
    }
    return(col_data)
  }

  if (inherits(newdata, "MultiAssayExperiment")) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'MultiAssayExperiment' is required for MAE newdata.",
          i = "BiocManager::install('MultiAssayExperiment')"),
        class = "cateomics_mae_package_missing"
      )
    }
    experiments_used <- names((meta$estimator$preprocess %||% meta$preprocess)$layers)
    if (is.null(experiments_used) || length(experiments_used) == 0L) {
      experiments_used <- names(MultiAssayExperiment::experiments(newdata))
    }
    sample_map <- MultiAssayExperiment::sampleMap(newdata)
    fmap <- sample_map[as.character(sample_map$assay) %in% experiments_used, ,
                       drop = FALSE]
    primary_counts <- table(as.character(fmap$primary))
    complete_primaries <- names(primary_counts)[
      primary_counts == length(experiments_used)
    ]
    .cf_check_unique_samples(
      rownames(MultiAssayExperiment::colData(newdata)),
      "The newdata MultiAssayExperiment colData")
    # Same as the fit path: `table()` sorts, and the prediction has to come
    # back in the order the caller handed the samples in.
    complete_primaries <- .cf_mae_primary_order(newdata, complete_primaries)
    if (length(complete_primaries) == 0L) {
      rlang::abort(
        c(
          "0 common primary samples across newdata experiments after sampleMap join.",
          i = "Check MultiAssayExperiment::sampleMap(newdata)."
        ),
        class = "cateomics_newdata_mae_no_overlap"
      )
    }
    # `optional = TRUE` so a non-syntactic column keeps the name the fit
    # recorded; the default would rewrite it and lose the covariate.
    col_data <- as.data.frame(
      MultiAssayExperiment::colData(newdata), optional = TRUE
    )[complete_primaries, , drop = FALSE]

    # Rebuild each layer the way the fit did: the rows that entered
    # preprocessing, aligned through sampleMap, then that layer's own
    # normalisation, ComBat and feature selection, then the `<experiment>_`
    # prefix the training covariates carry.
    prep_mae <- meta$estimator$preprocess %||% meta$preprocess
    layer_recs <- prep_mae$layers
    layer_rows <- lapply(layer_recs, `[[`, "features")
    if (!is.null(layer_rows) && length(layer_rows) > 0L) {
      for (lyr in names(layer_rows)) {
        cols <- as.character(layer_rows[[lyr]])
        if (length(cols) == 0L) next
        exp_obj <- MultiAssayExperiment::experiments(newdata)[[lyr]]
        exp_mat <- if (inherits(exp_obj, "SummarizedExperiment")) {
          SummarizedExperiment::assay(exp_obj, 1)
        } else if (is.matrix(exp_obj)) {
          exp_obj
        } else {
          NULL
        }
        if (is.null(exp_mat)) {
          rlang::abort(
            sprintf("newdata MAE has no usable experiment '%s'.", lyr),
            class = "cateomics_newdata_mae_layer_missing"
          )
        }
        missing_rows <- setdiff(cols, rownames(exp_mat))
        if (length(missing_rows) > 0L) {
          rlang::abort(
            c(sprintf("newdata experiment '%s' is missing %d training feature(s):",
                      lyr, length(missing_rows)),
              sprintf("  %s", paste(utils::head(missing_rows, 5L), collapse = ", "))),
            class = "cateomics_newdata_mae_layer_missing"
          )
        }
        exp_map <- fmap[as.character(fmap$assay) == lyr &
                          as.character(fmap$primary) %in% complete_primaries, ,
                        drop = FALSE]
        exp_colnames <- as.character(exp_map$colname)[
          match(complete_primaries, as.character(exp_map$primary))
        ]
        mat <- exp_mat[cols, exp_colnames, drop = FALSE]
        colnames(mat) <- complete_primaries
        .cf_check_assay_complete(mat, what = sprintf("Layer '%s' of newdata", lyr))
        block <- .cf_mae_newdata_layer(
          mat, prep = prep_mae, layer = lyr, col_data = col_data
        )
        for (g in rownames(block)) col_data[[g]] <- as.numeric(block[g, ])
      }
    }

    cov_names <- meta$covariate_names
    if (!is.null(cov_names)) {
      # Keep the clinical columns under their raw names: a factor is still one
      # column here and `.apply_newdata_preprocessing()` expands it into the
      # training dummies afterwards.
      layer_features <- unlist(lapply(layer_recs, `[[`, "selected_features"),
                               use.names = FALSE)
      if (is.null(layer_features)) {
        layer_features <- intersect(cov_names, colnames(col_data))
      }
      keep <- c(intersect(meta$named_args$X_baseline, colnames(col_data)),
                layer_features)
      missing <- setdiff(c(meta$named_args$X_baseline, layer_features),
                         colnames(col_data))
      if (length(missing) > 0L) {
        rlang::abort(
          c(sprintf("newdata MAE is missing %d training column(s):",
                    length(missing)),
            sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))),
          class = "cateomics_newdata_missing_columns"
        )
      }
      col_data <- col_data[, unique(keep), drop = FALSE]
    } else {
      drop_cols <- c(meta$named_args$Y, meta$named_args$T)
      col_data <- col_data[, setdiff(colnames(col_data), drop_cols),
                           drop = FALSE]
    }
    return(col_data)
  }

  rlang::abort(
    c(sprintf("Could not extract covariates from Bioconductor newdata of class: %s",
              paste(class(newdata), collapse = "/")),
      i = "newdata must be a SummarizedExperiment or MultiAssayExperiment (data.frame / matrix newdata bypass this path)."),
    class = "cateomics_newdata_unsupported_class"
  )
}

#' Keep the columns the fit used, for a data.frame or matrix newdata
#'
#' A container newdata is assembled column by column from what the fit
#' recorded, so anything else the caller's object carried never travels. A
#' data.frame or a matrix is handed over as it is, so the same selection has to
#' happen here -- otherwise what reaches the estimator depends on which extra
#' columns the caller's table happened to hold. The outcome and the treatment
#' always go; the covariates stay under their raw names, because a factor is
#' one column here and becomes its training dummies later. When features were
#' selected or projected the raw omics block has to stay whole, because the
#' selection or the rotation is re-applied downstream.
#'
#' @param object The cf_model.
#' @param newdata A data.frame or a matrix with column names.
#' @return `newdata` with the fit's columns (plus the batch column, which
#'   `.apply_newdata_preprocessing()` consumes and then drops).
#' @keywords internal
#' @noRd
.cf_newdata_training_columns <- function(object, newdata) {
  has_names <- (is.data.frame(newdata) || is.matrix(newdata)) &&
    !is.null(colnames(newdata))
  if (!has_names) return(newdata)
  meta <- object$meta
  prep <- meta$estimator$preprocess %||% meta$preprocess
  keep_in_order <- function(x, want) x[, colnames(x) %in% want, drop = FALSE]

  out <- keep_in_order(newdata, setdiff(colnames(newdata),
                                        c(meta$named_args$Y, meta$named_args$T)))
  # A selecting or projecting feature method is only safe to skip when there is
  # an omics layer: the final assembly then rebuilds X from `X_baseline` plus
  # that layer, so a stray column in newdata never reaches the backend. With no
  # layer -- covariates only -- nothing rebuilds, and the columns have to be
  # picked by name here exactly as they are for `feature_method = "none"`.
  fm <- prep$feature_method
  has_layer <- length(prep$layers %||% list()) > 0L
  if (has_layer && !is.null(fm) && fm %in% c("variance", "pca", "pathway"))
    return(out)

  cov_names <- meta$covariate_names
  if (is.null(cov_names)) return(out)
  raw_clinical <- meta$named_args$X_baseline
  want <- if (is.null(raw_clinical)) cov_names else
    unique(c(raw_clinical, intersect(cov_names, colnames(out))))
  missing <- setdiff(want, colnames(out))
  if (length(missing) > 0L) {
    rlang::abort(
      c(sprintf("newdata is missing %d training column(s):", length(missing)),
        sprintf("  %s", paste(utils::head(missing, 5L), collapse = ", "))),
      class = "cateomics_newdata_missing_columns"
    )
  }
  batch_keep <- if (!is.null(prep$batch_column) &&
                    prep$batch_column %in% colnames(out))
    prep$batch_column else character(0)
  keep_in_order(out, c(want, batch_keep))
}

#' Replay one MAE layer's preprocessing on newdata
#'
#' Mirrors the fit-time order for a single experiment: log-CPM (RNA-seq
#' layers), ComBat against the stored training block, then the layer's own
#' feature selection or PCA projection. Returns a features x samples matrix
#' whose rownames already carry the `<experiment>_` prefix.
#'
#' @param mat Features x samples matrix for this layer, columns aligned.
#' @param prep The preprocessing record (`meta$...$preprocess`); this layer's
#'   record, the feature method and the batch column are read off it.
#' @param layer Experiment name (the prefix).
#' @param col_data Sample table, for the batch labels.
#' @keywords internal
#' @noRd
.cf_mae_newdata_layer <- function(mat, prep, layer, col_data) {
  # The caller iterates the names of `prep$layers`, so the record is there.
  rec <- prep$layers[[layer]]
  feature_method <- prep$feature_method
  batch_column   <- prep$batch_column
  if (isTRUE(rec$enabled) || isTRUE(rec$batch_corrected)) {
    batch_new <- NULL
    if (isTRUE(rec$batch_corrected)) {
      if (is.null(batch_column) || !batch_column %in% colnames(col_data)) {
        rlang::abort(
          c(sprintf("Batch column '%s' missing in newdata.",
                    batch_column %||% "<unknown>"),
            i = "Add the batch column to newdata before predict()."),
          class = "cateomics_newdata_batch_missing"
        )
      }
      batch_new <- col_data[[batch_column]]
    }
    mat <- .cf_newdata_joint_block(
      mat, .cf_layer_train_omics(prep, layer),
      normalise = rec$enabled, method = rec$preprocess_method,
      batch_corrected = rec$batch_corrected,
      train_batch = rec$train_batch, batch_new = batch_new)$new
  }
  if (identical(feature_method, "variance") && !is.null(rec$selected)) {
    mat <- mat[rec$selected, , drop = FALSE]
  } else if (identical(feature_method, "pca") && !is.null(rec$pca_loadings)) {
    gene_order <- rec$pca_gene_names
    missing_g <- setdiff(gene_order, rownames(mat))
    if (length(missing_g) > 0L) {
      rlang::abort(
        c(sprintf("newdata PCA: %d training feature(s) missing from layer '%s':",
                  length(missing_g), layer),
          sprintf("  %s", paste(utils::head(missing_g, 5L), collapse = ", "))),
        class = "cateomics_newdata_pca_genes_missing"
      )
    }
    xs <- t(mat[gene_order, , drop = FALSE])
    if (!is.null(rec$pca_center)) xs <- sweep(xs, 2, rec$pca_center, "-")
    if (!is.null(rec$pca_scale) && !isFALSE(rec$pca_scale)) {
      xs <- sweep(xs, 2, rec$pca_scale, "/")
    }
    mat <- t(xs %*% rec$pca_loadings)
  }
  # The rows were taken by name (or are the PCs of the stored rotation), so
  # prefixing reproduces the training labels; compare rather than assume.
  rownames(mat) <- paste0(layer, "_", rownames(mat))
  if (!is.null(rec$selected_features) &&
      !identical(rownames(mat), rec$selected_features)) {
    rlang::abort(
      c(sprintf("Rebuilt layer '%s' does not match the training features.", layer),
        i = sprintf("Expected %d, got %d.",
                    length(rec$selected_features), nrow(mat))),
      class = "cateomics_newdata_mae_layer_mismatch"
    )
  }
  mat
}



#' Put newdata through the fit's normalisation chain, jointly with the
#' training block
#'
#' log-CPM, ComBat and GSVA all read the other samples in the matrix they are
#' handed: the TMM reference, the batch means, the gene ranks. Running them on the new samples alone therefore puts the result on
#' a different scale from the training values the model was fitted on. The
#' whole chain runs on the training block and the new samples together, and
#' both halves come back so a later step (pathway scoring) can use the same
#' pair.
#'
#' @param assay_new Features x samples matrix of new samples.
#' @param train_omics The stored training block, as it entered normalisation.
#' @param normalise,method Whether the fit normalised, and how.
#' @param batch_corrected,train_batch,batch_new ComBat inputs, or NULL.
#' @return list(train = , new = ), or train = NULL when nothing was stored.
#' @keywords internal
#' @noRd
.cf_newdata_joint_block <- function(assay_new, train_omics,
                                    normalise = FALSE, method = NULL,
                                    batch_corrected = FALSE,
                                    train_batch = NULL, batch_new = NULL) {
  do_lognorm <- isTRUE(normalise) && identical(method, "logcpm")
  if (is.null(train_omics)) {
    out <- if (do_lognorm) .cf_newdata_logcpm(assay_new) else assay_new
    return(list(train = NULL, new = out))
  }
  rows <- intersect(rownames(train_omics), rownames(assay_new))
  if (length(rows) == 0L) {
    rlang::abort(
      c("newdata shares no omics feature with the training block.",
        i = "Pass newdata with the same feature names as the fit."),
      class = "cateomics_newdata_omics_missing"
    )
  }
  n_train <- ncol(train_omics)
  n_new   <- ncol(assay_new)
  joint <- cbind(train_omics[rows, , drop = FALSE],
                 assay_new[rows, , drop = FALSE])
  # newdata that repeats a training sample (self-prediction) gives the joint
  # matrix a repeated column name, which edgeR reports to the user in its own
  # words. Nothing downstream reads these names -- the split below is
  # positional and the features are the rownames -- so make them unique for
  # the call and put the caller's names back on the way out.
  orig_cn <- colnames(joint)
  if (!is.null(orig_cn) && anyDuplicated(orig_cn)) {
    colnames(joint) <- make.unique(orig_cn)
  }
  if (do_lognorm) {
    joint <- .apply_lognorm(joint, input_type = "SE", preprocess = TRUE)$assay
  }
  if (isTRUE(batch_corrected)) {
    if (!requireNamespace("sva", quietly = TRUE)) {
      rlang::abort(
        c("Package 'sva' required for ComBat newdata reapplication.",
          i = "BiocManager::install('sva')"),
        class = "cateomics_sva_missing"
      )
    }
    joint <- sva::ComBat(
      dat = as.matrix(joint),
      batch = c(as.character(train_batch), as.character(batch_new)),
      mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  }
  if (!is.null(orig_cn)) colnames(joint) <- orig_cn
  list(train = joint[, seq_len(n_train), drop = FALSE],
       new   = joint[, n_train + seq_len(n_new), drop = FALSE])
}

#' Which columns of a preprocessed newdata frame are the clinical block
#'
#' The omics block is rebuilt from scratch on the predict path, so everything
#' else has to be carried over. Matching the raw formula names against
#' `X_new` misses a factor that has already been expanded into its training
#' dummies, which is why the training covariate names decide.
#'
#' @param meta The model's meta.
#' @param prep The preprocess record.
#' @param X_new The newdata frame, after factor expansion.
#' @param omics_cols The raw omics column names.
#' @param layers Output of `.identify_omics_columns()`.
#' @keywords internal
#' @noRd
.cf_newdata_baseline_cols <- function(meta, prep, X_new, omics_cols, layers) {
  train <- meta$covariate_names
  if (!is.null(train)) {
    return(intersect(setdiff(train, prep$selected_features %||% character(0)),
                     colnames(X_new)))
  }
  if (length(layers$baseline) > 0L) {
    return(intersect(layers$baseline, colnames(X_new)))
  }
  setdiff(colnames(X_new), omics_cols)
}

#' Apply fit-time newdata preprocessing (factor levels, z-score, assay
#' extraction, log-CPM, batch, feature_method).
#'
#' When all of \code{meta$estimator$preprocess / factor_levels / x_zscore}
#' and \code{meta$named_args$X_baseline} are \code{NULL} (a clinical-only fit
#' with no preprocessing), this function
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
  cov_raw <- object$meta$named_args$X_baseline

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
            class = "cateomics_new_factor_level"
          )
        }
      } else {
        rlang::abort(
          sprintf("Missing factor column '%s' in newdata.", var),
          class = "cateomics_missing_factor_column"
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

  # log-CPM + ComBat newdata re-application. A MAE was rebuilt layer by layer
  # in `.extract_X_from_bioc_newdata`, so these single-matrix steps would
  # apply a second time.
  mae_input <- identical(object$meta$input_type, "MAE")
  # Every omics feature the fit used has to be here. Normalisation divides by
  # the library size of the features present and GSVA ranks a set against
  # them, so a feature that is merely missing, not even selected, moves the
  # rest. The SummarizedExperiment path checks this at extraction; this is the
  # same check for the data.frame and matrix paths.
  if (!is.null(prep) && !mae_input && is.data.frame(X_new)) {
    train_cols <- unlist(.identify_omics_columns(object$meta)$layers,
                         use.names = FALSE)
    missing_omics <- setdiff(train_cols, colnames(X_new))
    if (length(train_cols) > 0L && length(missing_omics) > 0L) {
      rlang::abort(
        c(sprintf("newdata is missing %d of the %d omics features the fit used:",
                  length(missing_omics), length(train_cols)),
          sprintf("  %s", paste(utils::head(missing_omics, 5L), collapse = ", ")),
          i = "They set the library size and the ranking the scores are built on, so the remaining features would be on a different scale."),
        class = "cateomics_newdata_omics_missing"
      )
    }
  }
  # The omics block of newdata goes through the fit's chain once, jointly
  # with the stored training block; `train_block` is kept for the pathway
  # step below, which needs the same pair.
  train_block <- NULL
  if (!is.null(prep) && !mae_input) {
    layers <- .identify_omics_columns(object$meta)
    # One layer here by construction: the MAE path is handled above, and a
    # SummarizedExperiment or a data.frame carries a single omics block.
    lrec  <- if (length(prep$layers)) prep$layers[[1L]] else NULL
    lname <- if (length(prep$layers)) names(prep$layers)[[1L]] else NULL
    rnaseq_cols_all <- unlist(layers$layers, use.names = FALSE)
    if (length(rnaseq_cols_all) > 0L && is.data.frame(X_new)) {
      avail <- intersect(rnaseq_cols_all, colnames(X_new))
      needs_chain <- isTRUE(lrec$enabled) || isTRUE(lrec$batch_corrected) ||
        identical(prep$feature_method, "pathway")
      if (length(avail) > 0L && needs_chain) {
        assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))
        rownames(assay_new) <- avail

        batch_new <- NULL
        if (isTRUE(lrec$batch_corrected)) {
          bc <- prep$batch_column
          if (is.null(bc) || !bc %in% colnames(X_new)) {
            rlang::abort(
              c(sprintf("Batch column '%s' missing in newdata.",
                        bc %||% "<unknown>"),
                i = "Add the batch column to newdata before predict()."),
              class = "cateomics_newdata_batch_missing"
            )
          }
          batch_new <- X_new[[bc]]
        }

        blocks <- .cf_newdata_joint_block(
          assay_new, .cf_layer_train_omics(prep, lname),
          normalise = lrec$enabled, method = lrec$preprocess_method,
          batch_corrected = lrec$batch_corrected,
          train_batch = lrec$train_batch, batch_new = batch_new)
        train_block <- blocks$train
        assay_new   <- blocks$new

        # Write back into X_new column-by-column
        for (g in rownames(assay_new)) {
          X_new[[g]] <- as.numeric(assay_new[g, ])
        }
      }
    }
  }

  # GSVA pathway newdata re-application: re-run gsva on newdata using the same
  # gene_sets the fit used, jointly with the stored training block (see the
  # comment at the joint scoring below). The
  # selected pathways (when n_pathways < total_pathways) are taken from
  # `selected_features` so newdata pathway columns match training.
  prep_fm <- prep$feature_method
  if (!is.null(prep_fm) && prep_fm == "pathway") {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      rlang::abort(
        c("Package 'GSVA' required for pathway newdata reapplication.",
          i = "BiocManager::install('GSVA')"),
        class = "cateomics_gsva_missing"
      )
    }
    gene_sets <- .cf_record_gene_sets(prep)
    if (is.null(gene_sets)) {
      rlang::abort(
        c("This fit was made with feature_method='pathway' but carries no gene sets.",
          i = "Refit with cf_fit() or cf_compare() to score new samples against the training cohort."),
        class = "cateomics_corrupt_meta"
      )
    }
    layers <- .identify_omics_columns(object$meta)
    rnaseq_cols <- unlist(layers$layers, use.names = FALSE)
    avail <- intersect(rnaseq_cols, colnames(X_new))
    if (length(avail) == 0L) {
      rlang::abort(
        c("Pathway newdata: no rnaseq columns found in newdata.",
          i = "Pass newdata with the same assay rownames as the fit."),
        class = "cateomics_newdata_no_rnaseq"
      )
    }
    assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))
    rownames(assay_new) <- avail
    baseline_cols <- .cf_newdata_baseline_cols(
      object$meta, prep, X_new, rnaseq_cols, layers)
    X_baseline_new <- X_new[, baseline_cols, drop = FALSE]

    pathway_method <- prep$pathway_method %||% "gsva"
    # Every GSVA method reads the other samples: gsva and ssgsea rank genes
    # within the set they are given (and ssgsea divides by the range across
    # it), zscore standardises each gene over them, plage takes their SVD. So
    # scoring newdata alone puts it on a different scale from training. Score
    # the training block and the new samples together and keep the new
    # columns, the same way ComBat is re-applied.
    # `train_block` and `assay_new` came out of one joint pass above, so both
    # halves are already on one scale before they are scored together.
    if (!is.null(train_block)) {
      common_rows <- intersect(rownames(train_block), rownames(assay_new))
      if (length(common_rows) == 0L) {
        rlang::abort(
          c("Pathway newdata: no training gene is present in newdata.",
            i = "Pass newdata with the same assay rownames as the fit."),
          class = "cateomics_newdata_no_rnaseq"
        )
      }
      n_train <- ncol(train_block)
      joint <- cbind(train_block[common_rows, , drop = FALSE],
                     assay_new[common_rows, , drop = FALSE])
      gsva_joint <- .cf_gsva_scores(joint, gene_sets, pathway_method)
      gsva_out <- gsva_joint[, (n_train + 1L):ncol(gsva_joint), drop = FALSE]
    } else {
      gsva_out <- .cf_gsva_scores(assay_new, gene_sets, pathway_method)
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
          class = "cateomics_newdata_pathway_missing"
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

  # feature_method newdata re-application (single-matrix inputs only).
  if (!is.null(prep_fm) && !mae_input && prep_fm %in% c("variance", "pca")) {
    layers <- .identify_omics_columns(object$meta)
    rnaseq_cols <- unlist(layers$layers, use.names = FALSE)
    if (length(rnaseq_cols) > 0L && is.data.frame(X_new)) {
      avail <- intersect(rnaseq_cols, colnames(X_new))
      if (length(avail) == 0L) {
        rlang::abort(
          c("feature_method newdata: no rnaseq columns found in newdata.",
            i = "Pass newdata with the same assay rownames as the fit.")
          , class = "cateomics_newdata_no_rnaseq"
        )
      }
      assay_new <- t(as.matrix(X_new[, avail, drop = FALSE]))  # gene × sample
      ## Use meta$named_args$X_baseline / the layer records when available; this
      ## prevents extra colData columns from being treated as baseline (which
      ## would inflate the X_new column count and break grf::predict's strict
      ## column-count check). Fall back to setdiff when no named-args info.
      baseline_cols <- .cf_newdata_baseline_cols(
        object$meta, prep, X_new, rnaseq_cols, layers)
      X_baseline_new <- X_new[, baseline_cols, drop = FALSE]

      if (prep_fm == "variance") {
        # Use train-fold selected_names; "variance" means top-N by train
        # variance, so newdata must use the same N rows.
        sel <- prep$selected_features
        if (is.null(sel)) {
          rlang::abort(
            "feature_method='variance' but the fit records no selected features (meta$estimator$preprocess$selected_features).",
            class = "cateomics_corrupt_meta"
          )
        }
        missing_sel <- setdiff(sel, rownames(assay_new))
        if (length(missing_sel) > 0L) {
          rlang::abort(
            c(sprintf("newdata missing %d selected features:", length(missing_sel)),
              sprintf("  %s", paste(utils::head(missing_sel, 5L), collapse = ", "))),
            class = "cateomics_newdata_selected_features_missing"
          )
        }
        assay_sel <- assay_new[sel, , drop = FALSE]
        X_new_assay <- t(assay_sel)
      } else if (prep_fm == "pca") {
        # The rotation belongs to the layer it was fitted on, so it is read
        # from that layer's record. One layer here: the MAE path projects each
        # of its own layers further up.
        prec <- if (length(prep$layers)) prep$layers[[1L]] else NULL
        loadings <- prec$pca_loadings
        center <- prec$pca_center
        scale_v <- prec$pca_scale
        gene_order <- prec$pca_gene_names
        if (is.null(loadings) || is.null(gene_order)) {
          rlang::abort(
            "feature_method='pca' but the layer records no rotation.",
            class = "cateomics_corrupt_meta"
          )
        }
        missing_g <- setdiff(gene_order, rownames(assay_new))
        if (length(missing_g) > 0L) {
          rlang::abort(
            c(sprintf("newdata PCA: %d training gene(s) missing:", length(missing_g)),
              sprintf("  %s", paste(utils::head(missing_g, 5L), collapse = ", "))),
            class = "cateomics_newdata_pca_genes_missing"
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
        omics_cols <- unlist(.identify_omics_columns(object$meta)$layers,
                             use.names = FALSE)
        baseline_cols <- setdiff(colnames(X_new),
                                 c(sel, omics_cols %||% character(0)))
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
#' Used by slearner / tlearner / xlearner / ganite / tarnet / dragonnet.
#' cevae returns per-row credible intervals on newdata, and drlearner with a
#' non-linear final model aborts with its own class.
#' @keywords internal
#' @noRd
.abort_ite_ci_unsupported <- function(method) {
  rlang::abort(
    c(
      sprintf("CI for newdata is not implemented for method '%s'.", method),
      i = "Supported: grf (IJ variance), bcf (posterior), cevae (posterior), drlearner (linear final via econml StatsModelsInference)."
    ),
    class = "cateomics_ite_ci_unsupported"
  )
}

# -------------------------------------------------------------------------
# summary / print
# -------------------------------------------------------------------------

#' Get summary statistics from a cf_model
#'
#' @param object cf_model object
#' @param ... Additional arguments (unused)
#' @return A list with the ATE, its confidence interval and the label of the
#'   method that produced it (\code{ci_method}), and the summary of the
#'   \eqn{\hat\tau} distribution.
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
      class = "cateomics_invalid_input"
    )
  }

  if (is.null(object$res) || is.null(object$res$summary)) {
    rlang::abort(
      message = "Invalid cf_model: missing res / summary.",
      class = "cateomics_invalid_result"
    )
  }

  res <- object$res
  s   <- res$summary

  ci_method_label <- res$ci_method %||% object$meta$estimator$ci_method
  if (is.null(ci_method_label)) ci_method_label <- "unavailable"

  out <- list(
    method    = object$method,
    # `[[` not `$`: a field that is not there must come back NULL, and `$`
    # partial-matches a longer sibling (`named_args`) and returns a list.
    n         = object$meta[["n"]],
    p         = object$meta[["p"]],
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
  int <- function(v) if (is.null(v) || !length(v) || !is.finite(v)) NA_integer_ else as.integer(v)
  cat(sprintf("N samples: %d\n", int(x[["n"]])))
  cat(sprintf("N covars:  %d\n", int(x[["p"]])))
  cat("\nTreatment Effect Estimates:\n")
  cat(sprintf("  ATE:        %.4f\n", x$ate))
  ci_label <- if (is.null(x$ci_method)) "unavailable" else x$ci_method
  if (anyNA(x$ate_ci)) {
    cat(sprintf("  ATE 95%% CI: NA (%s)\n", ci_label))
  } else {
    cat(sprintf("  ATE 95%% CI: [%.4f, %.4f] (%s)\n",
                x$ate_ci[1], x$ate_ci[2], ci_label))
  }
  cat("\nCATE (tau) Distribution (in-sample):\n")
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
      class = "cateomics_invalid_input"
    )
  }
  res <- x$res
  cat(sprintf("<cf_model: %s>\n", x$method))
  # `sprintf` fed a NULL drops the line without a word, and `$` would reach a
  # longer sibling instead of saying the field is absent.
  int <- function(v) if (is.null(v) || !length(v) || !is.finite(v)) NA_integer_ else as.integer(v)
  cat(sprintf("  N = %d, p = %d\n", int(x$meta[["n"]]), int(x$meta[["p"]])))
  if (!is.null(res$summary)) {
    cat(sprintf("  ATE = %.4f (95%% CI: [%.4f, %.4f])\n",
                res$ate,
                res$summary$ate_ci_lower,
                res$summary$ate_ci_upper))
    cat(sprintf("  CATE (tau, in-sample): mean = %.4f, sd = %.4f\n",
                res$summary$tau_mean,
                res$summary$tau_sd))
  } else {
    cat(sprintf("  ATE = %.4f\n", res$ate))
  }
  sp <- x$meta$spec
  if (!is.null(sp)) {
    .cf_print_params_note(
      if (length(sp$params)) stats::setNames(list(sp$params), x$method) else NULL,
      sp$propensity_params, sp$eval_params)
  }
  invisible(x)
}

# The package-level `%||%` is defined once in utils_kfold.R.
