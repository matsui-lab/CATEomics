# =============================================================================
# preprocessing.R — cf_fit internal preprocessing-flow helpers
# =============================================================================
#
# Step-to-helper mapping for the cf_fit internal preprocessing flow:
#   Step 1   .cf_fit_step1_input_validation()  (input_type detection + formula
#                                                parsing + Y-leakage guard +
#                                                `omics` spec resolution +
#                                                MAE sampleMap join via
#                                                .extract_mae_assay() in
#                                                bioc_integration.R)
#   Step 2   .apply_lognorm()                  (TMM -> log-CPM)
#   Step 3   .apply_combat()                   (ComBat on log-scale values)
#   Step 3.5 .feature_select_step()            (variance / PCA loadings)
#   Step 4   .apply_pathway_aggregation()      (pathway aggregation)
#   Step 5   .expand_factor_columns()          (factor -> one-hot)
#   Step 6   validate_cateomics_result()
#   Step 7   create_cf_meta() + result$data


# ---------------------------------------------------------------------------
# Step 1: Input validation + input_type detection + assay resolution
# ---------------------------------------------------------------------------

#' Drop omics features that carry missing values, and say which
#'
#' Dropping an all-NA feature here rather than inside normalisation means the same thing happens whether
#' or not a layer is normalised, and the caller is told either way.
#'
#' @param assay_mat Numeric matrix (features x samples), or NULL.
#' @param what Label for the message.
#' @return The matrix without those features.
#' @keywords internal
#' @noRd
.cf_drop_na_features <- function(assay_mat, what = "the omics matrix") {
  if (is.null(assay_mat) || !anyNA(assay_mat)) return(assay_mat)
  bad <- apply(is.na(assay_mat), 1L, any)
  nm  <- rownames(assay_mat)
  show <- if (is.null(nm)) paste0("row ", utils::head(which(bad), 5L)) else
    utils::head(nm[bad], 5L)
  rlang::warn(
    c(sprintf("Dropped %d of the %d input features of %s: they carry missing values.",
              sum(bad), nrow(assay_mat), what),
      i = paste0("First dropped: ", paste(show, collapse = ", "), "."),
      i = sprintf("Feature selection runs on the remaining %d.",
                  nrow(assay_mat) - sum(bad)),
      i = "Impute them beforehand to keep them."),
    class = "cateomics_assay_na_features_dropped"
  )
  out <- assay_mat[!bad, , drop = FALSE]
  if (nrow(out) == 0L) {
    rlang::abort(
      c(sprintf("Every feature of %s carries missing values.", what),
        i = "Impute them, or fit without that layer."),
      class = "cateomics_assay_all_na"
    )
  }
  out
}

#' Reject omics features that carry missing values
#'
#' Runs after the entry-point drop, so what it catches is an NA a transform
#' produced (ComBat, pathway scoring), and PCA drops incomplete rows
#' by design.
#'
#' @param assay_mat Numeric matrix (features x samples), or NULL.
#' @param what What the matrix is, for the message.
#' @return Invisibly TRUE; aborts with `cateomics_assay_has_na` otherwise.
#' @keywords internal
#' @noRd
.cf_check_assay_complete <- function(assay_mat, what = "The omics features") {
  if (is.null(assay_mat) || !anyNA(assay_mat)) return(invisible(TRUE))
  bad  <- which(apply(is.na(assay_mat), 1L, any))
  nm   <- rownames(assay_mat)
  show <- if (is.null(nm)) paste0("row ", utils::head(bad, 5L)) else utils::head(nm[bad], 5L)
  rlang::abort(
    c(sprintf("%s carry missing values in %d of %d features.",
              what, length(bad), nrow(assay_mat)),
      i = paste0("First affected: ", paste(show, collapse = ", "), "."),
      i = "Impute or drop them before fitting."),
    class = "cateomics_assay_has_na"
  )
}

#' Read one layer's stored training block, whichever shape the record uses.
#'
#' A single fit stores each layer's block in that layer's record; a comparison
#' moves them into a shared environment, keyed by layer, so one copy is
#' serialised for all configs rather than one per config.
#'
#' @param prep The preprocessing record (`meta$estimator$preprocess`).
#' @param lyr Name of the layer to read.
#' @keywords internal
#' @noRd
.cf_layer_train_omics <- function(prep, lyr) {
  if (is.null(lyr) || !nzchar(lyr)) return(NULL)
  rec <- prep$layers[[lyr]]
  if (!is.null(rec$train_omics)) return(rec$train_omics)
  if (is.environment(prep$shared)) return(prep$shared$train_omics[[lyr]])
  NULL
}

#' Read the stored gene sets, whichever shape the record uses.
#' @keywords internal
#' @noRd
.cf_record_gene_sets <- function(rec) {
  if (!is.null(rec$gene_sets)) return(rec$gene_sets)
  if (is.environment(rec$shared)) return(rec$shared$gene_sets)
  NULL
}

#' Put the large preprocessing members in one box.
#'
#' Called once per `cf_compare()` so its configs share a single copy.
#'
#' @keywords internal
#' @noRd
.cf_preprocess_share <- function(prep) {
  keep <- isTRUE(prep$step2$enabled) || isTRUE(prep$step3$batch_corrected) ||
          identical(prep$feature_method, "pathway")
  if (!keep) return(NULL)
  box <- new.env(parent = emptyenv())
  # Keyed by layer, so a reader asking about one layer looks in one place
  # whatever the container was. Layers that needed no block hold NULL.
  box$train_omics <- lapply(prep$layer_records %||% list(), `[[`, "train_omics")
  box$gene_sets <- if (identical(prep$feature_method, "pathway"))
    prep$step4$gene_sets else NULL
  box
}

#' The preprocessing record a fit carries for its newdata path.
#'
#' Built from the shared Step 1-5 output by one function for both entries, so
#' `cf_fit` and `cf_compare` cannot describe the same preprocessing
#' differently.
#'
#' @param prep Return value of `.cf_preprocess_input()`.
#' @param share Environment from `.cf_preprocess_share()`, or NULL.
#' @return List describing normalisation, batch correction and feature
#'   selection, plus one record per omics layer for the newdata path to replay.
#' @keywords internal
#' @noRd
.cf_preprocess_record <- function(prep, share = NULL) {
  # The training blocks and `gene_sets` are the only large members. R's
  # serialisation writes a plain matrix once per occurrence even when the
  # occurrences are the same object in memory, so a comparison with K configs
  # would save K copies. An environment is written by reference, so when the
  # caller supplies one the blocks move there and every config holds the same
  # box. Reading goes through `.cf_layer_train_omics()` / `.cf_record_gene_sets()`.
  # Besides `layers` (the per-layer records) and `shared` (the box), two kinds
  # of field live here, and only two:
  #   (1) settings that do not belong to a layer (how features were chosen, the
  #       batch column, the pathway dictionary);
  #   (2) one index across the layers (`selected_features`), because its reader
  #       asks which design columns did *not* come from omics -- a question
  #       about all layers at once.
  # Everything that belongs to a layer is in `layers`, whatever the container,
  # so a reader asking about a layer has one place to look.
  layers <- prep$layer_records %||% list()
  # With a box, each layer's block lives in it under that layer's name; the
  # record then points at the box instead of carrying its own copy.
  if (is.environment(share))
    layers <- lapply(layers, function(l) { l$train_omics <- NULL; l })
  list(
    batch_method       = prep$step3$batch_method,
    batch_column       = prep$step3$batch_column,
    feature_method     = prep$feature_method,
    n_features         = prep$n_features,
    selected_features  = if (prep$feature_method == "pathway") prep$step4$selected_names else prep$step3_5$selected_names,
    pathway_method     = prep$step4$pathway_method,
    gene_sets_source   = prep$step4$gene_sets_source,
    # Persist the gene sets for GSVA newdata re-application: the new samples
    # are scored against the stored training block with this same dictionary.
    # The sets that were actually used: with the default the caller passed
    # NULL and the Hallmark collection was resolved inside step 4.
    gene_sets          = if (is.null(share) && prep$feature_method == "pathway")
                           prep$step4$gene_sets else NULL,
    # The box when the caller supplied one, NULL otherwise.
    shared             = share,

    # One record per omics layer -- always present, with a single entry when
    # the container carries one layer. The newdata path rebuilds each block
    # from these (its features, its normalisation, its ComBat material, its
    # selection or PCA projection, and the block it is re-run against).
    layers             = layers
  )
}

#' Step 1: Input validation + assay resolution
#'
#' Resolves data class (data.frame / SE / MAE), parses the formula, extracts
#' the assay matrix from the resolved `omics` spec, and enforces the Y-leakage
#' hard error on SE inputs.
#'
#' @param formula Y ~ T | X1 + X2 + ... formula.
#' @param data data.frame / SummarizedExperiment / MultiAssayExperiment.
#' @param omics_spec Resolved `omics` spec from `.cf_resolve_omics()`.
#'
#' @return List with components:
#'   * input_type:  "data.frame" / "SE" / "MAE"
#'   * parsed:      output of `parse_cf_formula` (Y / T / X / outcome_name /
#'                  treatment_name / covariate_names)
#'   * assay_mat:   numeric matrix (features × samples) or NULL (no assay)
#'   * rowdata:     rowData(se) DataFrame or NULL
#'
#' @keywords internal
#' @noRd
.cf_fit_step1_input_validation <- function(formula, data,
                                           omics_spec = NULL) {
  assay_name  <- omics_spec$assay_name %||% 1
  experiments <- omics_spec$experiments
  assay_cols  <- omics_spec$assay_cols
  # `omics = "none"`: the covariates are whatever X_baseline names, and no
  # assay matrix is built at all.
  omics_none  <- !is.null(omics_spec) && length(omics_spec$layers) == 0L

  input_type <- if (is_se(data)) {
    "SE"
  } else if (is_mae(data)) {
    "MAE"
  } else if (is.data.frame(data)) {
    "data.frame"
  } else {
    rlang::abort(
      sprintf("Unsupported data class: %s", paste(class(data), collapse = "/")),
      class = "cateomics_unsupported_data_class"
    )
  }
  # Sample identity is what every later join relies on, so reject a repeated
  # one here rather than let it surface as a bare R message (SE) or as a
  # silently repeated sample (MAE).
  if (identical(input_type, "SE")) {
    .cf_check_unique_samples(colnames(data), "The SummarizedExperiment")
  } else if (identical(input_type, "MAE")) {
    .cf_check_unique_samples(
      rownames(MultiAssayExperiment::colData(data)),
      "The MultiAssayExperiment colData")
  }

  rowdata <- NULL

  if (input_type == "SE" && omics_none) {
    # `optional = TRUE` keeps a non-syntactic colData name as it is. The
    # default rewrites it (`HLA-B` -> `HLA.B`), and the formula, which is
    # built from the names the caller gave, then fails to find its column.
    # `check.names` is not an option here: the method ignores `...`.
    col_data <- as.data.frame(SummarizedExperiment::colData(data),
                              optional = TRUE)
    rownames(col_data) <- colnames(data)
    parsed <- parse_cf_formula(formula, col_data)
    assay_mat <- NULL
    rowdata <- SummarizedExperiment::rowData(data)
  } else if (input_type == "MAE" && omics_none) {
    rlang::abort(
      c("`omics = \"none\"` leaves a MultiAssayExperiment with no measurements to use.",
        i = "Name the experiment(s) to use, or pass the colData as a data.frame."),
      class = "cateomics_omics_none_mae"
    )
  } else if (input_type == "SE") {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'SummarizedExperiment' is required for SE input.",
          i = "BiocManager::install('SummarizedExperiment')"),
        class = "cateomics_se_package_missing"
      )
    }
    assay_mat <- SummarizedExperiment::assay(data, assay_name)
    outcome_name <- all.vars(formula)[1]
    # Y in assay rownames is a hard error: Y must live in colData, not the
    # assay.
    if (outcome_name %in% rownames(assay_mat)) {
      rlang::abort(
        c(
          sprintf("Outcome '%s' found in assay rownames.", outcome_name),
          i = "CATEomics requires Y in colData only.",
          ">" = sprintf("Move '%s' to colData(se).", outcome_name)
        ),
        class = "cateomics_y_in_assay"
      )
    }
    # `optional = TRUE` keeps a non-syntactic colData name as it is. The
    # default rewrites it (`HLA-B` -> `HLA.B`), and the formula, which is
    # built from the names the caller gave, then fails to find its column.
    # `check.names` is not an option here: the method ignores `...`.
    col_data <- as.data.frame(SummarizedExperiment::colData(data),
                              optional = TRUE)
    rownames(col_data) <- colnames(data)
    # Features named in `omics` are assay rownames, so materialise them into
    # col_data before parse_cf_formula.
    materialize <- .materialize_se_assay_columns(col_data, assay_mat, omics_spec)
    col_data <- materialize$col_data
    if (length(materialize$assay_cols) > 0L && is.null(assay_cols)) {
      assay_cols <- materialize$assay_cols
    }
    # The assay supplies the design, so the formula is allowed to name no
    # covariate at all (an expression-only colData has none to name).
    parsed <- parse_cf_formula(formula, col_data, allow_no_covariate = TRUE)
    # When named-args materialized assay rows into colData, the formula RHS
    # picked them up as clinical-style covariates. Strip them from
    # parsed$X / parsed$covariate_names — they re-enter the X matrix via the
    # assay_mat path so they get the omics-specific normalisation / batch /
    # feature_select preprocessing and avoid duplicate columns. Also subset
    # assay_mat to the user-requested rows so downstream feature_method =
    # "none" doesn't drag in untouched assay genes.
    if (length(materialize$assay_cols) > 0L) {
      keep_cov <- setdiff(parsed$covariate_names, materialize$assay_cols)
      if (is.data.frame(parsed$X)) {
        parsed$X <- parsed$X[, keep_cov, drop = FALSE]
      } else {
        parsed$X <- parsed$X[, keep_cov, drop = FALSE]
      }
      parsed$covariate_names <- keep_cov
      have_rows <- intersect(materialize$assay_cols, rownames(assay_mat))
      assay_mat <- assay_mat[have_rows, , drop = FALSE]
    }
    rowdata <- SummarizedExperiment::rowData(data)
  } else if (input_type == "MAE") {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'MultiAssayExperiment' is required for MAE input.",
          i = "BiocManager::install('MultiAssayExperiment')"),
        class = "cateomics_mae_package_missing"
      )
    }
    # sampleMap handling lives in bioc_integration.R::.extract_mae_assay
    # (MAE-specific logic is consolidated in bioc_integration.R).
    extracted <- .extract_mae_assay(data, formula, experiments = experiments,
                                    omics_spec = omics_spec)
    # A MAE keeps its layers apart all the way through preprocessing, so there
    # is no single assay matrix here.
    assay_mat <- NULL
    mae_layer_mats      <- extracted$layer_mats
    mae_layer_is_rnaseq <- extracted$layer_is_rnaseq
    parsed    <- extracted$parsed
    rowdata   <- extracted$rowdata
    mae_complete_primaries <- extracted$complete_primaries
    mae_primary_is_rnaseq  <- extracted$primary_is_rnaseq
    mae_primary_exp        <- extracted$primary_exp
  } else {
    # data.frame path. The omics columns of a data.frame are the assay matrix
    # and never enter the formula, so with no clinical covariate the formula
    # names nothing; that is only allowed when those columns do exist.
    parsed <- parse_cf_formula(formula, data,
                               allow_no_covariate = !is.null(assay_cols))
    assay_mat <- .resolve_df_assay(data, assay_cols = assay_cols,
                                   parsed = parsed)
  }

  list(
    input_type = input_type,
    parsed     = parsed,
    assay_mat  = assay_mat,
    # MAE only: one matrix per experiment, and which of them are RNA-seq.
    layer_mats = if (exists("mae_layer_mats", inherits = FALSE))
      mae_layer_mats else NULL,
    layer_is_rnaseq = if (exists("mae_layer_is_rnaseq", inherits = FALSE))
      mae_layer_is_rnaseq else NULL,
    rowdata    = rowdata,
    mae_complete_primaries = if (exists("mae_complete_primaries", inherits = FALSE))
      mae_complete_primaries else NULL,
    # TRUE when the experiment the MAE carries its covariates on is the
    # RNA-seq one. Step 2 does not read it -- normalisation is decided per
    # layer from `layer_is_rnaseq` -- so this records the primary and nothing
    # more.
    mae_primary_is_rnaseq = if (exists("mae_primary_is_rnaseq", inherits = FALSE))
      mae_primary_is_rnaseq else NULL
    # What the omics block resolved to is recorded per layer, in
    # `preprocess$layers`: the layer's name is the assay name on a
    # SummarizedExperiment, the experiment name on a MultiAssayExperiment and a
    # fixed label on a data.frame, and that one name is what the newdata path
    # needs to read the same block back.
  )
}


#' Resolve assay matrix for data.frame inputs
#'
#' For the data.frame path, assay columns are specified via `assay_cols`
#' resolved from `omics`. When NULL the input is treated as clinical-only and
#' `preprocess = TRUE` is a no-op.
#'
#' @keywords internal
#' @noRd
.resolve_df_assay <- function(data, assay_cols = NULL, parsed = NULL) {
  if (!is.null(assay_cols)) {
    if (is.character(assay_cols)) {
      missing_cols <- setdiff(assay_cols, names(data))
      if (length(missing_cols) > 0) {
        rlang::abort(
          c(
            "assay_cols references column(s) absent from data:",
            sprintf("  %s", paste(missing_cols, collapse = ", "))
          ),
          class = "cateomics_assay_cols_missing"
        )
      }
    }
    sub <- data[, assay_cols, drop = FALSE]
  } else {
    # omics = "none": clinical-only path, there is no assay matrix to
    # normalise. The raw-counts assumption warning in Step 2 fires only when
    # the user actually asks for an omics block.
    return(NULL)
  }

  # Conflict check: assay columns must not overlap the formula X covariates.
  if (!is.null(parsed) && length(parsed$covariate_names) > 0) {
    overlap <- intersect(names(sub), parsed$covariate_names)
    if (length(overlap) > 0) {
      rlang::abort(
        c(
          "Assay columns overlap formula covariates:",
          sprintf("  %s", paste(overlap, collapse = ", ")),
          i = "Remove it from one side (X_baseline or omics)."
        ),
        class = "cateomics_column_conflict"
      )
    }
  }

  # Transpose to features x samples matrix form.
  t(as.matrix(sub))
}


# ---------------------------------------------------------------------------
# Step 2: Omics-specific preprocessing (TMM log-CPM)
# ---------------------------------------------------------------------------

#' Step 2: Apply fit-time normalisation.
#'
#' The fit-time normalisation is an explicit
#' \code{zero-variance drop -> TMM (normLibSizes) -> log-CPM} pipeline. TMM =
#' TMM = Robinson & Oshlack (2010), via edgeR (Chen et al. 2025); full
#' citations are in \code{?cf_fit}.
#' \code{filterByExpr} is not applied here.
#'
#' @param assay_mat Numeric matrix (features x samples) or NULL.
#' @param input_type "data.frame" / "SE" / "MAE".
#' @param preprocess Logical, normalisation toggle (default TRUE).
#' @param norm_method TMM-family normalisation passed to
#'   \code{edgeR::normLibSizes} (default "TMM").
#'
#' @return list(assay, enabled, method, weights, filter_kept,
#'   norm_method). When `enabled = FALSE` the input matrix passes
#'   through unchanged.
#'
#' @keywords internal
#' @noRd
.apply_lognorm <- function(assay_mat, input_type = "SE",
                           preprocess = TRUE,
                           norm_method = "TMM") {
  if (!isTRUE(preprocess) || is.null(assay_mat)) {
    return(list(assay = assay_mat, enabled = FALSE,
                method = NULL, weights = NULL,
                filter_kept = NULL, norm_method = NULL))
  }
  if (identical(input_type, "data.frame")) {
    rlang::warn(
      c("`preprocess = TRUE` on data.frame assumes raw counts.",
        i = "Use preprocess = FALSE if the assay is already normalized."),
      class = "cateomics_df_preprocess_warning"
    )
  }
  assay_mat <- as.matrix(assay_mat)

  # Zero-variance features carry no signal and break downstream scaling.
  row_var <- apply(assay_mat, 1L, stats::var)
  keep <- is.finite(row_var) & row_var > 0
  if (any(!keep)) assay_mat <- assay_mat[keep, , drop = FALSE]
  filter_kept <- as.integer(sum(keep))

  # edgeR is in Imports, so a working install always has it and this branch is
  # unreachable from a loaded package. It is kept for the same reason the one in
  # `.fit_xgboost()` is: the alternative is a function that cannot say why it
  # produced a different normalisation.
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    rlang::warn(
      c("`edgeR` not installed; log-CPM computed without TMM normalisation.",
        i = "Install edgeR for TMM library-size normalisation."),
      class = "cateomics_edger_missing"
    )
    lib <- colSums(assay_mat)
    lib[lib == 0] <- 1
    logcpm <- log2(t(t(assay_mat) / lib) * 1e6 + 0.5)
    return(list(assay = logcpm, enabled = TRUE, method = "logcpm",
                weights = NULL, filter_kept = filter_kept,
                norm_method = NULL))
  }

  dge <- edgeR::DGEList(counts = assay_mat)
  # calcNormFactors was renamed to normLibSizes in edgeR 4.x.
  dge <- edgeR::normLibSizes(dge, method = norm_method)
  logcpm <- edgeR::cpm(dge, log = TRUE)
  list(
    assay       = logcpm,
    enabled     = TRUE,
    method      = "logcpm",
    weights     = NULL,
    filter_kept = filter_kept,
    norm_method = norm_method
  )
}


#' Newdata log-CPM.
#'
#' Mirrors \code{.apply_lognorm}'s TMM -> log-CPM on a counts matrix, for
#' newdata at predict time. The library size is the column sum over the
#' features present and the TMM reference comes from the samples handed in,
#' so the caller must pass the same features the fit used.
#'
#' @param assay_mat Numeric counts matrix (features x samples).
#' @return log-CPM matrix (same shape).
#' @keywords internal
#' @noRd
.cf_newdata_logcpm <- function(assay_mat) {
  assay_mat <- as.matrix(assay_mat)
  if (requireNamespace("edgeR", quietly = TRUE)) {
    dge <- edgeR::DGEList(counts = assay_mat)
    dge <- edgeR::normLibSizes(dge, method = "TMM")
    return(edgeR::cpm(dge, log = TRUE))
  }
  lib <- colSums(assay_mat)
  lib[lib == 0] <- 1
  log2(t(t(assay_mat) / lib) * 1e6 + 0.5)
}


# ---------------------------------------------------------------------------
# Step 3: Batch correction (sva::ComBat on log-scale values)
# ---------------------------------------------------------------------------

#' Step 3: Apply sva::ComBat batch correction on log-scale values.
#'
#' Batch correction runs after normalisation, on whatever Step 2 produced
#' (log-CPM for an RNA-seq layer, the input unchanged otherwise), using
#' \code{sva::ComBat} (Johnson 2007).
#'
#' @param assay_mat Numeric matrix (features x samples), log-scale.
#' @param batch_vec Batch label vector (length == ncol(assay_mat)).
#' @param batch_label Character, original column name (recorded in meta).
#'
#' @return List(assay, batch_corrected, batch_method, batch_column).
#'
#' @keywords internal
#' @noRd
.apply_combat <- function(assay_mat, batch_vec, batch_label = NULL) {
  if (is.null(batch_vec) || is.null(assay_mat)) {
    return(list(assay = assay_mat, batch_corrected = FALSE,
                batch_method = NULL, batch_column = batch_label))
  }
  if (!requireNamespace("sva", quietly = TRUE)) {
    rlang::abort(
      c("Package 'sva' is required for batch correction.",
        i = "BiocManager::install('sva')"),
      class = "cateomics_sva_missing"
    )
  }
  if (length(batch_vec) != ncol(assay_mat)) {
    rlang::abort(
      sprintf(
        "Batch vector length (%d) != assay sample count (%d).",
        length(batch_vec), ncol(assay_mat)
      ),
      class = "cateomics_batch_length_mismatch"
    )
  }
  corrected <- sva::ComBat(dat = as.matrix(assay_mat), batch = batch_vec,
                           mod = NULL, par.prior = TRUE,
                           prior.plots = FALSE)
  list(
    assay           = corrected,
    batch_corrected = TRUE,
    batch_method    = "combat",
    batch_column    = batch_label
  )
}


# ---------------------------------------------------------------------------
# Step 3.5: Feature select (variance / pca)
# ---------------------------------------------------------------------------

#' Step 3.5: Feature selection (variance / pca)
#'
#' Saves PCA loadings for downstream gene mapping. `feature_method =
#' "pathway"` skips Step 3.5 and proceeds to Step 4 (pathway aggregation).
#' `"none"` is identity (all features retained).
#'
#' @param assay_mat Numeric matrix (features x samples).
#' @param feature_method One of "none" / "variance" / "pca" / "pathway".
#' @param n_features Integer, top n.
#'
#' @return List(features, pca_loadings, method, n_features, selected_names).
#'   `features` is a matrix (features_selected x samples). For PCA, rows are
#'   PC1, PC2, ..., and `pca_loadings` is the gene x n_PC rotation matrix
#'   (for downstream user-side gene mapping).
#'
#' @keywords internal
#' @noRd
.feature_select_step <- function(assay_mat,
                                 feature_method = c("none", "variance", "pca", "pathway"),
                                 n_features = 100) {
  feature_method <- match.arg(feature_method)
  result <- list(
    features = assay_mat,
    pca_loadings = NULL,
    method = feature_method,
    n_features = n_features,
    selected_names = if (!is.null(assay_mat)) rownames(assay_mat) else NULL
  )
  if (is.null(assay_mat) || feature_method == "none" ||
      feature_method == "pathway") {
    return(result)
  }
  if (feature_method == "variance") {
    if (n_features >= nrow(assay_mat)) {
      return(result)  # all features retained
    }
    vars <- apply(assay_mat, 1, stats::var, na.rm = TRUE)
    top_idx <- order(vars, decreasing = TRUE)[seq_len(n_features)]
    result$features <- assay_mat[top_idx, , drop = FALSE]
    result$selected_names <- rownames(assay_mat)[top_idx]
    return(result)
  }
  # feature_method == "pca"
  mat_clean <- assay_mat[stats::complete.cases(assay_mat), , drop = FALSE]
  # Minimum requirement for prcomp: samples >= 2, features >= 1.
  max_pcs <- min(nrow(mat_clean), ncol(mat_clean) - 1L)
  if (max_pcs < 1L) {
    rlang::abort(
      "Cannot compute PCA: insufficient features x samples.",
      class = "cateomics_pca_insufficient_dim"
    )
  }
  n_pcs <- min(n_features, max_pcs)
  pca <- stats::prcomp(t(mat_clean), center = TRUE, scale. = TRUE)
  pc_features <- t(pca$x[, seq_len(n_pcs), drop = FALSE])
  rownames(pc_features) <- paste0("PC", seq_len(n_pcs))
  result$features <- pc_features
  result$pca_loadings <- pca$rotation[, seq_len(n_pcs), drop = FALSE]
  # Persist train-fold center + scale so newdata PCA re-application is
  # data-leakage-safe (apply train statistics, do not recompute on newdata).
  # prcomp's `$center` / `$scale` are length = n_features (gene-level)
  # vectors; rownames match `rownames(mat_clean)`.
  result$pca_center <- pca$center
  result$pca_scale <- pca$scale
  result$pca_gene_names <- rownames(mat_clean)
  result$selected_names <- rownames(pc_features)
  result
}


# ---------------------------------------------------------------------------
# Step 4: Pathway aggregation (GSVA)
# ---------------------------------------------------------------------------

#' Build the GSVA parameter object for a pathway method
#'
#' GSVA >= 1.50 takes a method-specific parameter object; older versions take
#' `method =`.
#'
#' @param assay_mat Numeric matrix (features x samples).
#' @param gene_sets Named list of gene sets.
#' @param pathway_method One of "gsva", "ssgsea", "plage", "zscore".
#' @return The matrix of pathway scores.
#' @keywords internal
#' @noRd
.cf_gsva_scores <- function(assay_mat, gene_sets, pathway_method) {
  mat <- as.matrix(assay_mat)
  if (utils::packageVersion("GSVA") < "1.50") {
    return(GSVA::gsva(mat, gene_sets, method = pathway_method))
  }
  param <- switch(
    pathway_method,
    "gsva"   = GSVA::gsvaParam(exprData = mat, geneSets = gene_sets,
                               kcdf = "auto"),
    "ssgsea" = GSVA::ssgseaParam(exprData = mat, geneSets = gene_sets),
    "plage"  = GSVA::plageParam(exprData = mat, geneSets = gene_sets),
    "zscore" = GSVA::zscoreParam(exprData = mat, geneSets = gene_sets),
    rlang::abort(sprintf("Unsupported pathway_method '%s'.", pathway_method),
                 class = "cateomics_unsupported_pathway_method")
  )
  GSVA::gsva(param)
}


#' Step 4: Pathway aggregation (GSVA::gsva)
#'
#' `feature_method = "pathway"` maps the all-gene assay to pathway scores
#' (default Hallmark 50 via msigdbr H, Homo sapiens). When the number of
#' pathways exceeds n_features, the top n_features by variance are selected.
#'
#' GSVA 1.50+ uses the `gsvaParam`-based API; 1.49 and earlier use the legacy
#' direct API. `utils::packageVersion` branches to support both.
#'
#' Methods: GSVA = Hänzelmann, Castelo & Guinney (2013), ssGSEA = Barbie et al.
#' (2009), both provided via the GSVA package; full citations are in
#' \code{?cf_fit}.
#'
#' @param assay_mat Numeric matrix (genes x samples) — Step 3 output.
#' @param gene_sets Named list of gene sets, or NULL (default = Hallmark 50).
#' @param pathway_method "gsva" / "ssgsea" / "plage" / "zscore".
#' @param n_features Top-variance pathway count (Inf uses all pathways).
#'
#' @return List(pathway_matrix, gene_sets_source, pathway_method,
#'              n_pathways, selected_names).
#'
#' @keywords internal
#' @noRd
.apply_pathway_aggregation <- function(assay_mat, gene_sets = NULL,
                                       pathway_method = c("gsva", "ssgsea", "plage", "zscore"),
                                       n_features = 100) {
  pathway_method <- match.arg(pathway_method)
  if (is.null(assay_mat)) {
    return(list(
      pathway_matrix = NULL, gene_sets_source = NULL,
      pathway_method = NULL, n_pathways = 0L, selected_names = NULL
    ))
  }
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    rlang::abort(
      c("Package 'GSVA' is required for feature_method = 'pathway'.",
        i = "BiocManager::install('GSVA')"),
      class = "cateomics_gsva_missing"
    )
  }
  gene_sets_source <- "custom"
  if (is.null(gene_sets)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      rlang::abort(
        c("Package 'msigdbr' is required when gene_sets = NULL.",
          i = "install.packages('msigdbr')"),
        class = "cateomics_msigdbr_missing"
      )
    }
    h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    gene_sets <- split(h$gene_symbol, h$gs_name)
    gene_sets_source <- "hallmark_msigdbr"
    rlang::inform(
      "Using default gene_sets = Hallmark 50 (msigdbr H, Homo sapiens).",
      class = "cateomics_default_gene_sets"
    )
  }

  pathway_matrix <- .cf_gsva_scores(assay_mat, gene_sets, pathway_method)

  selected_names <- rownames(pathway_matrix)
  if (is.finite(n_features) && n_features < nrow(pathway_matrix)) {
    vars <- apply(pathway_matrix, 1, stats::var, na.rm = TRUE)
    top_idx <- order(vars, decreasing = TRUE)[seq_len(n_features)]
    pathway_matrix <- pathway_matrix[top_idx, , drop = FALSE]
    selected_names <- rownames(pathway_matrix)
  }

  list(
    pathway_matrix   = pathway_matrix,
    # The resolved sets, not the argument: with the default the caller passed
    # NULL, and newdata has to be scored against the same dictionary.
    gene_sets        = gene_sets,
    gene_sets_source = gene_sets_source,
    pathway_method   = pathway_method,
    n_pathways       = nrow(pathway_matrix),
    selected_names   = selected_names
  )
}


# ---------------------------------------------------------------------------
# Step 5 (pre-dispatch): Factor / character -> one-hot
# ---------------------------------------------------------------------------

#' Step 5 pre-dispatch: Expand factor / character columns to one-hot
#'
#' Uses `model.matrix(~ . - 1, data = X_clinical)` for dummy expansion. Edge
#' cases:
#'   * 1-level factor -> warning + column drop
#'   * > 100 levels   -> warning (high-cardinality risk)
#'   * ordered factor -> as.numeric (order preserved)
#'   * character      -> as.factor -> one-hot
#'
#' When `method = "cevae"` and binary dummies are present, a warning fires
#' because they risk violating the continuous-X distribution assumption
#' (Louizos et al. 2017).
#'
#' @param X_clinical data.frame or matrix of clinical covariates.
#' @param method cf_fit method name (used for CEVAE binary X warning).
#'
#' @return List(X, factor_levels).
#'   * X: numeric matrix (factors expanded, column count = expanded width)
#'   * factor_levels: named list of train-time levels for each factor
#'                    column (re-applied at newdata predict time)
#'
#' @keywords internal
#' @noRd
.expand_factor_columns <- function(X_clinical, method = "grf") {
  if (!is.data.frame(X_clinical)) {
    # Already a matrix means no factors -> return the numeric matrix as-is.
    if (is.matrix(X_clinical) && is.numeric(X_clinical)) {
      return(list(X = X_clinical, factor_levels = NULL))
    }
    X_clinical <- as.data.frame(X_clinical, check.names = FALSE)
  }

  is_factor_like <- vapply(
    X_clinical,
    function(x) is.factor(x) || is.character(x) || is.logical(x),
    logical(1)
  )
  factor_cols <- names(X_clinical)[is_factor_like]

  if (length(factor_cols) == 0) {
    return(list(X = as.matrix(X_clinical), factor_levels = NULL))
  }

  factor_levels <- list()
  cols_to_drop <- character(0)

  for (col in factor_cols) {
    raw <- X_clinical[[col]]
    if (is.logical(raw)) {
      # logical is binary; as.integer gives 0/1 (no factor expansion needed).
      X_clinical[[col]] <- as.integer(raw)
      next
    }
    if (is.ordered(raw)) {
      # ordered factor keeps numeric order
      X_clinical[[col]] <- as.numeric(raw)
      next
    }
    fx <- as.factor(raw)
    n_levels <- nlevels(fx)
    if (n_levels < 2L) {
      rlang::warn(
        sprintf("Factor '%s' has < 2 levels; dropping column.", col),
        class = "cateomics_single_level_factor"
      )
      cols_to_drop <- c(cols_to_drop, col)
      next
    }
    if (n_levels > 100L) {
      rlang::warn(
        sprintf("Factor '%s' has %d levels (> 100, high-dim risk).",
                col, n_levels),
        class = "cateomics_high_cardinality_factor"
      )
    }
    factor_levels[[col]] <- levels(fx)
    X_clinical[[col]] <- fx
  }

  if (length(cols_to_drop) > 0) {
    X_clinical <- X_clinical[, setdiff(names(X_clinical), cols_to_drop),
                             drop = FALSE]
  }

  # one-hot expansion via model.matrix (no intercept). When the remaining
  # columns are all numeric, model.matrix passes them through unchanged and
  # expands only the factors.
  if (ncol(X_clinical) == 0L) {
    X_mat <- matrix(numeric(0), nrow = nrow(X_clinical), ncol = 0L)
  } else {
    X_mat <- stats::model.matrix(~ . - 1, data = X_clinical)
    attr(X_mat, "assign") <- NULL
    attr(X_mat, "contrasts") <- NULL
  }

  # CEVAE binary X warning: count binary dummy columns after expansion.
  if (identical(method, "cevae") && ncol(X_mat) > 0L) {
    is_binary <- apply(X_mat, 2, function(v) {
      uniq <- unique(stats::na.omit(v))
      length(uniq) == 2L && all(uniq %in% c(0, 1))
    })
    if (any(is_binary)) {
      rlang::warn(
        c(
          "This CEVAE implementation models X with a continuous (Normal) likelihood.",
          i = "Binary/dummy features may violate that distribution assumption.",
          i = sprintf("Binary feature(s): %s",
                      paste(colnames(X_mat)[is_binary], collapse = ", ")),
          ">" = "Consider grf / bcf for high-cardinality categorical data."
        ),
        class = "cateomics_cevae_binary_x"
      )
    }
  }

  list(
    X = X_mat,
    factor_levels = if (length(factor_levels) > 0) factor_levels else NULL
  )
}


# ---------------------------------------------------------------------------
# omics layer identification helper
# ---------------------------------------------------------------------------

#' Which columns of a fit came from which omics layer
#'
#' Reads the per-layer records the fit stored (`preprocess$layers`): the name
#' of each is the layer, and `features` is what entered preprocessing from it.
#' `.extract_X_from_bioc_newdata()` uses this to rebuild the same blocks from
#' newdata and apply each layer's own normalisation.
#'
#' @param meta `object$meta` from a fitted cf_model.
#' @return List with `baseline` (the clinical column names the caller gave)
#'   and `layers` (named list of feature-name vectors, one per omics layer).
#' @keywords internal
#' @noRd
.identify_omics_columns <- function(meta) {
  prep <- meta$estimator$preprocess %||% meta$preprocess
  lay  <- prep$layers %||% list()
  list(
    baseline = if (is.null(meta$named_args$X_baseline)) character(0)
               else as.character(meta$named_args$X_baseline),
    layers   = lapply(lay, function(r) as.character(r$features))
  )
}


#' Materialise assay rownames into colData for SE input
#'
#' When the user calls `cf_fit(Y = "Y", T = "T", omics = c("G1", ...),
#' data = se)` the facade builds the formula `Y ~ T | G1 + G2 + ...`,
#' but `parse_cf_formula()` evaluates it against `colData(se)` only — so
#' `G1`/`G2` (assay rownames) are not visible and the parse is refused with
#' `cateomics_unknown_column`. This helper copies the requested assay rows
#' into `col_data` as numeric columns so the formula resolves.
#'
#' A SummarizedExperiment carries one layer at a time, so at most one layer
#' is materialised here. The MAE path is handled separately in
#' `.cf_fit_step1_input_validation()`'s MAE branch via
#' `.extract_mae_assay()`.
#'
#' @param col_data data.frame produced by `as.data.frame(colData(se))`.
#' @param assay_mat features x samples matrix from the selected assay.
#' @param omics_spec resolved `omics` spec from `.cf_resolve_omics()`.
#' @return list(col_data = updated_df, assay_cols = character(0|N) used).
#' @keywords internal
#' @noRd
.materialize_se_assay_columns <- function(col_data, assay_mat, omics_spec) {
  feats <- omics_spec$features %||% list()
  if (is.null(assay_mat) || length(feats) == 0L) {
    return(list(col_data = col_data, assay_cols = character(0)))
  }
  assay_cols_used <- character(0)
  for (lyr in names(feats)) {
    rows <- as.character(feats[[lyr]])
    if (length(rows) == 0L) next
    missing_in_assay <- setdiff(rows, rownames(assay_mat))
    if (length(missing_in_assay) > 0L) {
      rlang::abort(
        c(
          sprintf("`omics` references %d row(s) absent from the assay matrix:",
                  length(missing_in_assay)),
          sprintf("  %s", paste(utils::head(missing_in_assay, 5L), collapse = ", ")),
          i = sprintf("Check rownames(assay(se, \"%s\")).", lyr)
        ),
        class = "cateomics_se_assay_rows_missing"
      )
    }
    for (gene in rows) col_data[[gene]] <- as.numeric(assay_mat[gene, ])
    assay_cols_used <- c(assay_cols_used, rows)
  }
  list(col_data = col_data, assay_cols = assay_cols_used)
}

# ---------------------------------------------------------------------------
# Shared preprocessing pipeline
# ---------------------------------------------------------------------------

#' Resolve the input-type-aware default for `n_features`
#'
#' `NULL` (the cf_fit / cf_compare default) resolves to 50 for a
#' MultiAssayExperiment (per-assay budget, multiple assays) and 100 otherwise
#' (SummarizedExperiment / data.frame). An explicit numeric value passes
#' through untouched.
#'
#' @keywords internal
#' @noRd
.cf_resolve_n_features <- function(n_features, data) {
  if (!is.null(n_features)) return(n_features)
  if (is_mae(data)) 50L else 100L
}

#' Shared Step 1-5 preprocessing pipeline for cf_fit / cf_compare
#'
#' Single source of truth for the SE / MAE / data.frame preprocessing that
#' turns raw input into the covariate matrix `X`, treatment `T` and outcome
#' `Y` that every backend consumes. `cf_fit()` and `cf_compare()` apply the
#' *identical* pipeline to the same input.
#'
#' Pipeline: input validation + assay resolution (Step 1) -> omics
#' log-normalisation (Step 2, TMM -> log-CPM; MAE = RNA-seq assay only)
#' -> batch correction (Step 3, ComBat on log-CPM) -> feature select /
#' aggregation (Step 4, variance / pca / pathway) -> factor expansion + final
#' `X = cbind(X_clinical, t(processed_features))` (Step 5). Binary-outcome
#' resolution and method dispatch stay in the caller.
#'
#' @return A list with `X` / `T` / `Y` / `cov_names` (required by every caller)
#'   plus `parsed` / `step1` / `step2` / `step3` / `step3_5` / `step4` /
#'   `batch_vec` / `input_type` / `rowdata` / `factor_levels` /
#'   `clinical_cov_names` / `feature_method` / `n_features` /
#'   `layer_records` / `x_zscore` / `sample_names_meaningful` (consumed by the
#'   cf_fit meta-populate path).
#'
#' @keywords internal
#' @noRd
.cf_preprocess_input <- function(formula, data,
                                 omics_spec     = NULL,
                                 preprocess     = TRUE,
                                 batch          = NULL,
                                 feature_method = c("variance", "none", "pca", "pathway"),
                                 n_features     = 100,
                                 gene_sets      = NULL,
                                 pathway_method = c("gsva", "ssgsea", "plage", "zscore"),
                                 method         = "grf") {
  feature_method <- match.arg(feature_method)
  pathway_method <- match.arg(pathway_method)

  # Step 1: Input validation + assay resolution (SE/MAE/df each resolve
  # assay_mat separately, parsed$X is clinical-only). Y-leakage hard error,
  # rowData, input_type detection, MAE sampleMap join.
  step1 <- .cf_fit_step1_input_validation(
    formula       = formula,
    data          = data,
    omics_spec    = omics_spec
  )
  parsed     <- step1$parsed
  assay_mat  <- step1$assay_mat
  input_type <- step1$input_type
  rowdata    <- step1$rowdata

  # Missing values are handled once, here, so the outcome does not depend on
  # whether a layer happens to be normalised further down.
  assay_mat <- .cf_drop_na_features(assay_mat)
  if (!is.null(step1$layer_mats)) {
    step1$layer_mats <- stats::setNames(
      lapply(names(step1$layer_mats), function(lyr) {
        .cf_drop_na_features(step1$layer_mats[[lyr]],
                             sprintf("layer '%s'", lyr))
      }),
      names(step1$layer_mats))
  }
  # The per-layer records are assembled below, after this drop, so the
  # features they name are the ones that actually entered preprocessing. A
  # feature dropped here must not survive in them, or prediction would demand
  # a column the fit never used.

  # Pathway aggregation needs a single gene-expression assay; it is not
  # defined for the multi-assay MAE feature matrix.
  if (identical(input_type, "MAE") && identical(feature_method, "pathway")) {
    rlang::abort(
      c("feature_method = 'pathway' is not supported for MultiAssayExperiment input.",
        i = "Pathway aggregation requires a single gene-expression assay; use it with SE or data.frame RNA-seq input."),
      class = "cateomics_mae_pathway_unsupported"
    )
  }

  # Step 2: Omics normalisation (zero-var drop -> TMM -> log-CPM). With
  # assay_mat = NULL (omics = "none") this is a no-op passthrough. Only an
  # RNA-seq layer is normalised; with a single layer that is whatever
  # `preprocess` refers to, with several the layer name (or `omics_type`)
  # decides, so a methylation or mutation matrix is never treated as counts.
  # Single-matrix inputs (SE / data.frame) carry one decision; a MAE decides
  # per layer, in `step1$layer_is_rnaseq`.
  layer_is_rnaseq <- is.null(omics_spec$types) ||
    length(omics_spec$types) == 0L || isTRUE(omics_spec$types[[1L]])
  # A data.frame concatenates its layers into one matrix, so a mixed set cannot
  # be normalised selectively: say so rather than applying one layer's decision
  # to all of them.
  if (isTRUE(preprocess) && !identical(input_type, "MAE") &&
      length(omics_spec$types) > 1L &&
      length(unique(omics_spec$types)) > 1L) {
    rlang::abort(
      c("`preprocess = TRUE` cannot be applied to a mix of RNA-seq and other layers on this input.",
        i = sprintf("RNA-seq: %s / other: %s.",
                    paste(names(omics_spec$types)[omics_spec$types], collapse = ", "),
                    paste(names(omics_spec$types)[!omics_spec$types], collapse = ", ")),
        i = "Normalise outside CATEomics and pass preprocess = FALSE, or fit one layer at a time."),
      class = "cateomics_mixed_layer_preprocess"
    )
  }
  # Step 3's batch vector is resolved once and reused by every layer.
  batch_vec <- NULL
  if (!is.null(batch)) {
    batch_vec <- .resolve_batch_vector(data, batch, input_type, step1)
  }

  # Steps 2-4 for one omics matrix: log-normalise (RNA-seq only) -> ComBat on
  # the log values -> feature select. A MAE runs this once per layer, every
  # other input once in total.
  run_layer <- function(mat, is_rnaseq) {
    s2 <- .apply_lognorm(
      mat, input_type = input_type,
      preprocess = isTRUE(preprocess) && isTRUE(is_rnaseq)
    )
    s3 <- .apply_combat(s2$assay, batch_vec = batch_vec, batch_label = batch)
    s35 <- .feature_select_step(
      s3$assay, feature_method = feature_method, n_features = n_features
    )
    list(step2 = s2, step3 = s3, step3_5 = s35)
  }

  layer_records <- NULL
  if (identical(input_type, "MAE")) {
    # Each experiment is normalised, batch-corrected and feature-selected on
    # its own, then the blocks are stacked. Feature names carry the experiment
    # they came from, so two layers can share a feature name.
    layer_mats <- step1$layer_mats %||% list()
    layer_flags <- step1$layer_is_rnaseq %||%
      stats::setNames(rep(FALSE, length(layer_mats)), names(layer_mats))
    layer_records <- list()
    blocks <- list()
    for (lyr in names(layer_mats)) {
      rec <- run_layer(layer_mats[[lyr]], unname(layer_flags[[lyr]]))
      feats <- rec$step3_5$features
      raw_names <- rec$step3_5$selected_names %||% rownames(feats)
      if (is.null(raw_names)) raw_names <- paste0("feat", seq_len(nrow(feats)))
      prefixed <- paste0(lyr, "_", raw_names)
      rownames(feats) <- prefixed
      blocks[[lyr]] <- feats
      layer_records[[lyr]] <- list(
        # What entered preprocessing for this layer, and what survived
        # selection. The newdata path rebuilds the block from `features` and
        # then re-applies the same selection.
        features          = rownames(layer_mats[[lyr]]),
        selected          = raw_names,
        is_rnaseq         = isTRUE(unname(layer_flags[[lyr]])),
        enabled           = rec$step2$enabled,
        preprocess_method = rec$step2$method,
        batch_corrected   = rec$step3$batch_corrected,
        train_omics       = if (isTRUE(rec$step2$enabled) ||
                                isTRUE(rec$step3$batch_corrected))
          layer_mats[[lyr]] else NULL,
        train_batch       = if (isTRUE(rec$step3$batch_corrected))
          batch_vec else NULL,
        selected_features = prefixed,
        pca_loadings      = rec$step3_5$pca_loadings,
        pca_center        = rec$step3_5$pca_center,
        pca_scale         = rec$step3_5$pca_scale,
        pca_gene_names    = rec$step3_5$pca_gene_names
      )
    }
    final_features <- if (length(blocks) == 0L) NULL else do.call(rbind, blocks)
    feat_names <- unlist(lapply(layer_records, `[[`, "selected_features"),
                         use.names = FALSE)
    # Flat records summarise the layers for the parts of meta that predate the
    # per-layer handling; the per-layer detail lives in `layer_records`.
    step2 <- list(
      assay   = NULL,
      enabled = any(vapply(layer_records, function(r) isTRUE(r$enabled), logical(1))),
      method  = {
        ms <- unique(unlist(lapply(layer_records, `[[`, "preprocess_method")))
        if (length(ms) == 1L) ms else if (length(ms) == 0L) NULL else paste(ms, collapse = "/")
      }
    )
    step3 <- list(
      batch_corrected = any(vapply(layer_records,
                                   function(r) isTRUE(r$batch_corrected),
                                   logical(1))),
      batch_method    = if (is.null(batch_vec)) NULL else "combat",
      batch_column    = batch
    )
    step3_5 <- list(features = final_features, selected_names = feat_names,
                    pca_loadings = NULL, pca_center = NULL, pca_scale = NULL,
                    pca_gene_names = NULL, method = feature_method,
                    n_features = n_features)
    step4 <- list(pathway_matrix = NULL, gene_sets = NULL,
                  gene_sets_source = NULL, pathway_method = NULL,
                  n_pathways = 0L, selected_names = NULL)
  } else {
    one <- run_layer(assay_mat, layer_is_rnaseq)
    step2 <- one$step2
    step3 <- one$step3
    step4_select <- one$step3_5
    step4_pathway <- if (feature_method == "pathway") {
      .apply_pathway_aggregation(
        step3$assay,
        gene_sets      = gene_sets,
        pathway_method = pathway_method,
        n_features     = n_features
      )
    } else {
      list(
        pathway_matrix = NULL, gene_sets = NULL, gene_sets_source = NULL,
        pathway_method = NULL, n_pathways = 0L, selected_names = NULL
      )
    }

    # Step 5 preparation: final features = the selected block (or the pathway
    # scores).
    final_features <- if (feature_method == "pathway") {
      step4_pathway$pathway_matrix
    } else {
      step4_select$features
    }
    feat_names <- if (feature_method == "pathway") {
      step4_pathway$selected_names
    } else {
      step4_select$selected_names
    }
    step3_5 <- step4_select
    step4   <- step4_pathway
    # One layer, recorded the same way as each of an MAE's experiments, so the
    # newdata path reads one shape whatever the container was. The name is the
    # assay name on a SummarizedExperiment and the fixed label on a data.frame.
    layer_records <- if (is.null(assay_mat)) list() else {
      lab <- if (!is.null(omics_spec) && length(omics_spec$layers))
        names(omics_spec$layers)[[1L]] else "omics"
      stats::setNames(list(list(
        features          = rownames(assay_mat),
        selected          = feat_names,
        is_rnaseq         = isTRUE(layer_is_rnaseq),
        enabled           = step2$enabled,
        preprocess_method = step2$method,
        batch_corrected   = step3$batch_corrected,
        # The omics matrix as it entered normalisation. Normalisation, ComBat
        # and pathway scoring all read the other samples (the TMM reference,
        # the batch means, the gene ranks), so newdata is put through the same
        # chain jointly with this block and only its own columns are kept.
        # Stored when any of those steps ran.
        train_omics       = if (isTRUE(step2$enabled) ||
                                isTRUE(step3$batch_corrected) ||
                                identical(feature_method, "pathway"))
          assay_mat else NULL,
        train_batch       = if (isTRUE(step3$batch_corrected)) batch_vec else NULL,
        selected_features = feat_names,
        pca_loadings      = step3_5$pca_loadings,
        pca_center        = step3_5$pca_center,
        pca_scale         = step3_5$pca_scale,
        pca_gene_names    = step3_5$pca_gene_names
      )), lab)
    }
  }

  # Say whether each omics layer was library-size normalised. The decision
  # comes from the layer's name when the caller did not pass `omics_type`
  # (`.CF_RNASEQ_LABELS`), and it changes every number downstream, so report
  # what was done rather than leaving the caller to read it off the record.
  # Only when there is more than one layer: a single layer is whatever
  # `preprocess` refers to. Not when `preprocess = FALSE`, because then nothing
  # is normalised and the caller asked for that.
  if (isTRUE(preprocess) && length(layer_records) > 1L) {
    from_name <- omics_spec$types_from_name %||%
      stats::setNames(rep(TRUE, length(layer_records)), names(layer_records))
    did  <- vapply(layer_records, function(r) isTRUE(r$enabled), logical(1))
    if (any(!did & from_name[names(layer_records)])) {
      rlang::inform(
        c(sprintf("Library-size normalisation (log-CPM): %s.",
                  paste(c(sprintf("applied to %s",
                                  if (any(did)) paste(names(layer_records)[did], collapse = ", ")
                                  else "no layer"),
                          sprintf("not applied to %s",
                                  paste(names(layer_records)[!did], collapse = ", "))),
                        collapse = "; ")),
          i = "The layer name decides this; pass omics_type to change it, e.g. omics_type = c(counts = \"rnaseq\")."),
        class = "cateomics_layer_normalisation"
      )
    }
  }

  # (a) Expand factor / character / logical columns in X_clinical to a one-hot
  # numeric matrix. `.expand_factor_columns` returns a numeric matrix plus
  # `factor_levels` for newdata re-application.
  factor_expand <- .expand_factor_columns(parsed$X, method = method)
  X_clinical_expanded <- factor_expand$X
  factor_levels <- factor_expand$factor_levels
  clinical_cov_names <- if (ncol(X_clinical_expanded) == 0L) {
    character(0)
  } else {
    colnames(X_clinical_expanded)
  }

  # A missing value that reaches the covariate matrix only surfaces much
  # further downstream, where the message no longer points at the input. The
  # check is on the block that actually enters X: a feature that selection
  # left behind, or a row PCA dropped, is not the model's problem.
  .cf_check_assay_complete(final_features)

  # (b) Assemble final X = cbind(X_clinical_expanded, t(processed_features))
  X_clinical <- X_clinical_expanded
  if (is.null(final_features) || (is.matrix(final_features) && nrow(final_features) == 0L)) {
    X <- X_clinical
    cov_names <- clinical_cov_names
  } else {
    if (is.null(feat_names)) feat_names <- rownames(final_features)
    if (is.null(feat_names)) feat_names <- paste0("feat", seq_len(nrow(final_features)))
    X <- cbind(X_clinical, t(final_features))
    cov_names <- c(clinical_cov_names, feat_names)
  }
  dup <- unique(cov_names[duplicated(cov_names)])
  if (length(dup) > 0L) {
    rlang::abort(
      c(sprintf("%d covariate name(s) occur twice in X:", length(dup)),
        sprintf("  %s", paste(utils::head(dup, 5L), collapse = ", ")),
        i = "Rename the clashing feature or column so each covariate is named once."),
      class = "cateomics_duplicate_covariate_names"
    )
  }
  Y <- parsed$Y
  T <- parsed$T

  # (c) Central covariate z-score: standardise the assembled X once, before
  # fold splitting / method dispatch, so every scale-sensitive backend (base
  # learners, catenets, ganite, cevae) sees a common representation, binary
  # indicators included. This is an X-only
  # (outcome-free) transform, so a single global fit is leakage-safe; the train
  # statistics are persisted for newdata via meta$x_zscore. grf / bcf are
  # scale-invariant, so the standardised X leaves their output unchanged.
  x_zscore <- if (ncol(X) > 0L) .cf_compute_x_zscore(X) else NULL
  if (!is.null(x_zscore)) X <- .cf_apply_x_zscore(X, x_zscore)

  list(
    X = X, T = T, Y = Y, cov_names = cov_names,
    parsed = parsed, step1 = step1, step2 = step2, step3 = step3,
    step3_5 = step3_5, step4 = step4, batch_vec = batch_vec,
    # One record per omics layer, whatever the container was: what entered
    # preprocessing, what survived selection, and the block the newdata path
    # replays the chain against.
    layer_records = layer_records,
    input_type = input_type, rowdata = rowdata,
    factor_levels = factor_levels, clinical_cov_names = clinical_cov_names,
    feature_method = feature_method, n_features = n_features,
    x_zscore = x_zscore,
    # Whether the row names of `X` name anything. A container always names its
    # samples; a data.frame may carry R's automatic row numbers, and the two
    # cannot be told apart afterwards because "1" .. "n" is also a legitimate
    # set of integer sample identifiers. R itself knows which it is --
    # `.row_names_info()` is negative for the automatic kind -- so the answer
    # is recorded here instead of being inferred later.
    sample_names_meaningful =
      if (!identical(input_type, "data.frame")) TRUE
      else if (is.data.frame(data)) .row_names_info(data) > 0L
      else FALSE
  )
}
