# =============================================================================
# bioc_integration.R -- container handling for SummarizedExperiment and
# MultiAssayExperiment inputs.
# =============================================================================
# `cf_fit()` / `cf_compare()` read these containers through the shared
# preprocessing in `preprocessing.R`; the helpers below are what that path
# calls to pull a layer out of a container, to order its samples and to
# refuse an input it cannot align.

#' Put MAE primary sample units back in colData order
#'
#' `table()` returns its names sorted, so the set of complete primaries comes
#' out alphabetically. The fitted object's rows follow this order and results
#' are written back onto colData row by row, so a sorted order would put each
#' sample's estimate on a different sample.
#'
#' @param mae The MultiAssayExperiment.
#' @param ids Character vector of primary sample units.
#' @return `ids`, in the order they appear in colData.
#' @keywords internal
#' @noRd
.cf_mae_primary_order <- function(mae, ids) {
  ref <- rownames(MultiAssayExperiment::colData(mae))
  if (is.null(ref)) return(ids)
  c(ref[ref %in% ids], setdiff(ids, ref))
}

#' Require the container to name each sample once
#'
#' Estimates are matched to samples by name: the design matrix carries the
#' sample identifiers as rownames, results are written back onto the
#' container by name, and `newdata` is aligned to the fit by name. A repeated
#' identifier makes that correspondence undefined -- the two rows can hold
#' different outcomes, and nothing says which one the estimate belongs to.
#'
#' A `data.frame` cannot reach this: R rejects duplicate `row.names`. A
#' `SummarizedExperiment` accepts repeated `colnames` and a
#' `MultiAssayExperiment` accepts a `colData` with repeated rownames, and
#' Left to themselves both fail further downstream -- the SE with a bare R
#' message about row names, the MAE by silently repeating the sample -- so
#' they are refused here, where the message can name the input.
#'
#' @param ids Character vector of sample identifiers.
#' @param what Character, how to name the container in the message.
#' @return `ids`, invisibly.
#' @keywords internal
#' @noRd
.cf_check_unique_samples <- function(ids, what = "the object") {
  if (is.null(ids)) return(invisible(ids))
  dups <- unique(ids[duplicated(ids)])
  if (length(dups) == 0L) return(invisible(ids))
  shown <- paste(utils::head(dups, 5L), collapse = ", ")
  rlang::abort(
    c(sprintf("%s names %d sample(s) more than once.", what, length(dups)),
      i = sprintf("Repeated: %s%s.", shown,
                  if (length(dups) > 5L) ", ..." else ""),
      ">" = "CATEomics matches each estimate to its sample by name, so give every sample its own identifier."),
    class = "cateomics_duplicate_samples"
  )
}

#' Extract one matrix per experiment from a MultiAssayExperiment
#'
#' Helper for cf_fit Step 1 (preprocessing.R
#' `.cf_fit_step1_input_validation`): returns every requested experiment as its
#' own features x samples matrix, aligned to the primary sample units through
#' `MultiAssayExperiment::sampleMap()`. The caller normalises, batch-corrects
#' and feature-selects each layer separately and concatenates the results.
#'
#' @param mae A MultiAssayExperiment.
#' @param formula Y ~ T | X1 + X2 + ... formula.
#' @param experiments Character vector of experiment names to use (default
#'   NULL = all; the first, or the RNA-seq one, is the primary carrier).
#'
#' @return List with components:
#'   * layer_mats: named list of features x samples matrices, one per
#'                 experiment, columns named by primary sample identifier.
#'   * layer_is_rnaseq: named logical, which layers are RNA-seq (log-normalised
#'                 downstream).
#'   * parsed:     output of `parse_cf_formula` against colData(mae) (Y / T /
#'                 X / outcome_name / treatment_name / covariate_names), with
#'                 materialised assay features stripped from X.
#'   * rowdata:    rowData of the primary experiment, or NULL when it is a
#'                 plain matrix.
#'   * complete_primaries: primary sample names retained after sampleMap join.
#'   * primary_exp: name of the primary experiment.
#'
#' @keywords internal
#' @noRd
.extract_mae_assay <- function(mae, formula, experiments = NULL,
                                # Resolved `omics` spec. Layers that name their
                                # features have those rows materialised into
                                # colData so parse_cf_formula can resolve them
                                # at MAE fit time.
                                omics_spec = NULL) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    rlang::abort(
      c("Package 'MultiAssayExperiment' is required for MAE input.",
        i = "BiocManager::install('MultiAssayExperiment')"),
      class = "cateomics_mae_package_missing"
    )
  }
  if (is.null(experiments)) {
    experiments <- names(MultiAssayExperiment::experiments(mae))
  }
  if (length(experiments) == 0L) {
    rlang::abort(
      "MAE contains 0 experiments.",
      class = "cateomics_mae_empty"
    )
  }
  missing_exp <- setdiff(experiments, names(MultiAssayExperiment::experiments(mae)))
  if (length(missing_exp) > 0L) {
    rlang::abort(
      c(
        "Requested experiments are absent from the MAE:",
        sprintf("  %s", paste(missing_exp, collapse = ", "))
      ),
      class = "cateomics_mae_unknown_experiment"
    )
  }
  # Only the RNA-seq assay is internally log-normalised; every other assay
  # passes through as-is (user pre-normalised). Which experiment that is comes
  # from the resolved spec, which has already applied `omics_type` and the
  # layer-name rule; without a spec (an internal direct call) the experiment
  # name itself is read, against `.CF_RNASEQ_LABELS`. When none matches, no
  # assay is normalised (all passthrough) and experiments[1] is used only as
  # the covariate carrier.
  rnaseq_exp <- NULL
  flagged <- if (is.null(omics_spec$types) || length(omics_spec$types) == 0L) {
    experiments[tolower(experiments) %in% .CF_RNASEQ_LABELS]
  } else {
    names(omics_spec$types)[omics_spec$types]
  }
  flagged <- intersect(flagged, experiments)
  if (length(flagged) > 0L) rnaseq_exp <- flagged[1]
  primary_is_rnaseq <- !is.null(rnaseq_exp)
  primary_exp <- if (primary_is_rnaseq) rnaseq_exp else experiments[1]

  # sampleMap join to align primary <-> assay colnames
  sample_map <- MultiAssayExperiment::sampleMap(mae)
  filtered_map <- sample_map[as.character(sample_map$assay) %in% experiments, ,
                             drop = FALSE]
  primary_counts <- table(as.character(filtered_map$primary))
  complete_primaries <- names(primary_counts)[
    primary_counts == length(experiments)
  ]
  # `table()` sorts its names, so take the colData order back: downstream
  # writes results onto the object row by row.
  complete_primaries <- .cf_mae_primary_order(mae, complete_primaries)
  if (length(complete_primaries) == 0L) {
    rlang::abort(
      c(
        sprintf(
          "0 common primary samples across experiments %s after sampleMap join.",
          paste(experiments, collapse = ", ")
        ),
        i = "Check MultiAssayExperiment::sampleMap(mae).",
        ">" = "Reduce `experiments =` to overlapping experiments."
      ),
      class = "cateomics_mae_no_overlap"
    )
  }

  # Every requested experiment contributes a layer. Columns are aligned to the
  # primary sample units through sampleMap; rows are the layer's own features,
  # or the subset `omics` named for that layer.
  align_cols <- function(exp_name) {
    ex_map <- filtered_map[
      as.character(filtered_map$assay) == exp_name &
        as.character(filtered_map$primary) %in% complete_primaries, ,
      drop = FALSE
    ]
    as.character(ex_map$colname)[
      match(complete_primaries, as.character(ex_map$primary))
    ]
  }
  layer_is_rnaseq <- stats::setNames(experiments %in% flagged, experiments)
  feats_spec <- omics_spec$features %||% list()
  unknown_layer <- setdiff(names(feats_spec), experiments)
  if (length(unknown_layer) > 0L) {
    rlang::abort(
      sprintf("`omics` layer '%s' is not an experiment of the MAE.",
              unknown_layer[1]),
      class = "cateomics_mae_layer_resolve_failed"
    )
  }

  layer_mats <- list()
  layer_rowdata <- list()
  for (exp_name in experiments) {
    exp_obj <- MultiAssayExperiment::experiments(mae)[[exp_name]]
    if (inherits(exp_obj, "SummarizedExperiment")) {
      mat <- SummarizedExperiment::assay(exp_obj, 1)
      layer_rowdata[[exp_name]] <- SummarizedExperiment::rowData(exp_obj)
    } else if (is.matrix(exp_obj)) {
      mat <- exp_obj
    } else {
      rlang::abort(
        sprintf("Unsupported experiment class '%s' for cf_fit.",
                class(exp_obj)[1]),
        class = "cateomics_mae_unsupported_primary"
      )
    }
    if (is.null(rownames(mat))) {
      rlang::abort(
        c(sprintf("Experiment '%s' has no feature names (rownames).", exp_name),
          i = "CATEomics labels MAE covariates <experiment>_<feature>, and newdata is matched by those names.",
          ">" = sprintf("Set rownames on experiment '%s'.", exp_name)),
        class = "cateomics_mae_layer_no_rownames"
      )
    }
    want <- as.character(feats_spec[[exp_name]] %||% character(0))
    if (length(want) > 0L) {
      missing_rows <- setdiff(want, rownames(mat))
      if (length(missing_rows) > 0L) {
        rlang::abort(
          c(sprintf("`omics$%s` references %d row(s) absent from that experiment:",
                    exp_name, length(missing_rows)),
            sprintf("  %s", paste(utils::head(missing_rows, 5L), collapse = ", "))),
          class = "cateomics_mae_layer_resolve_failed"
        )
      }
      mat <- mat[want, , drop = FALSE]
    }
    mat <- mat[, align_cols(exp_name), drop = FALSE]
    colnames(mat) <- complete_primaries
    layer_mats[[exp_name]] <- mat
    # Only an RNA-seq layer is log-normalised internally; a counts-like matrix
    # in any other layer would pass through untransformed.
    if (!isTRUE(unname(layer_is_rnaseq[[exp_name]])) && .is_raw_counts_like(mat)) {
      rlang::warn(
        c(
          sprintf("Non-RNA-seq assay '%s' looks like raw counts (integer / non-negative).",
                  exp_name),
          i = "Only an RNA-seq layer is internally log-normalised; every other layer passes through as-is.",
          ">" = sprintf("Normalise '%s' externally before adding it to the MAE.", exp_name)
        ),
        class = "cateomics_mae_raw_counts_detected"
      )
    }
  }

  # Subset colData to primary sample units, then parse the formula.
  # See the note in `.cf_preprocess_input()`: the default makes the names
  # syntactic, which loses the column the caller named.
  col_data <- as.data.frame(MultiAssayExperiment::colData(mae),
                            optional = TRUE)
  col_data <- col_data[complete_primaries, , drop = FALSE]
  if (ncol(col_data) == 0L) {
    rlang::abort(
      "MAE colData has 0 columns; formula cannot be parsed.",
      class = "cateomics_mae_empty_coldata"
    )
  }

  # Layers that named their features have those rows materialised into
  # col_data so parse_cf_formula resolves them, as the SE path does.
  for (lyr in names(feats_spec)) {
    cols <- as.character(feats_spec[[lyr]])
    if (length(cols) == 0L) next
    mat_aligned <- layer_mats[[lyr]]
    for (g in cols) col_data[[g]] <- as.numeric(mat_aligned[g, ])
  }

  # The layers supply the design, so the formula is allowed to name no
  # covariate at all (a colData with only Y and T has none to name).
  parsed <- parse_cf_formula(formula, col_data, allow_no_covariate = TRUE)

  # Same dedup as the SE path, for every layer: rows materialised so the
  # formula could resolve are stripped from parsed$X, because they re-enter
  # through the layer matrices, which is where normalisation and feature
  # selection happen.
  materialised <- unique(unlist(lapply(feats_spec, as.character),
                                use.names = FALSE))
  if (length(materialised) > 0L) {
    keep_cov <- setdiff(parsed$covariate_names, materialised)
    parsed$X <- parsed$X[, keep_cov, drop = FALSE]
    parsed$covariate_names <- keep_cov
  }

  list(
    layer_mats = layer_mats,
    layer_is_rnaseq = layer_is_rnaseq,
    parsed = parsed,
    rowdata = layer_rowdata[[primary_exp]],
    complete_primaries = complete_primaries,
    primary_exp = primary_exp,
    primary_is_rnaseq = primary_is_rnaseq
  )
}


#' Raw-counts heuristic
#'
#' Flags a matrix as raw-counts-like when
#' `is.integer(assay) || (numeric, non-negative, whole-number)`. Proteomics
#' log2 intensities and methylation M-values are excluded.
#'
#' @keywords internal
#' @noRd
.is_raw_counts_like <- function(mat) {
  vals <- as.numeric(mat)
  finite_vals <- vals[is.finite(vals)]
  if (length(finite_vals) == 0L) {
    return(FALSE)
  }
  is.integer(mat) ||
    (all(finite_vals >= 0) && all(finite_vals == floor(finite_vals)))
}


#' Check if object is SummarizedExperiment
#' @noRd
is_se <- function(x) {
  inherits(x, "SummarizedExperiment")
}

#' Check if object is MultiAssayExperiment
#' @noRd
is_mae <- function(x) {
  inherits(x, "MultiAssayExperiment")
}
