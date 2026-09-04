# =============================================================================
# preprocessing.R — cf_fit internal preprocessing-flow helpers
# =============================================================================
#
# Step-to-helper mapping for the cf_fit internal preprocessing flow:
#   Step 1   .cf_fit_step1_input_validation()  (input_type detection + formula
#                                                parsing + Y-leakage guard +
#                                                assay_cols / pattern resolution +
#                                                MAE sampleMap join via
#                                                .extract_mae_assay() in
#                                                bioc_integration.R)
#   Step 2   .apply_lognorm()                  (TMM -> log-CPM)
#   Step 3   .apply_combat()                   (ComBat on log-scale values)
#   Step 3.5 .feature_select_step()            (variance / PCA loadings)
#   Step 4   .apply_pathway_aggregation()      (pathway aggregation)
#   Step 5   .expand_factor_columns()          (factor -> one-hot)
#   Step 6   validate_cfomics_result()
#   Step 7   create_cf_meta() + result$data


# ---------------------------------------------------------------------------
# Step 1: Input validation + input_type detection + assay resolution
# ---------------------------------------------------------------------------

#' Step 1: Input validation + assay resolution
#'
#' Resolves data class (data.frame / SE / MAE), parses the formula, extracts
#' the assay matrix (assay_cols / assay_pattern for df, assay_name for SE,
#' experiments for MAE), and enforces the Y-leakage hard error on SE inputs.
#'
#' @param formula Y ~ T | X1 + X2 + ... formula.
#' @param data data.frame / SummarizedExperiment / MultiAssayExperiment.
#' @param assay_name SE: index or name of assay (default 1).
#' @param experiments MAE: experiment names (NULL = all).
#' @param assay_cols df: character / integer / logical vector of column spec.
#' @param assay_pattern df: regex pattern matched against `names(data)`.
#'
#' @return List with components:
#'   * input_type:  "data.frame" / "SE" / "MAE"
#'   * parsed:      output of `parse_cf_formula` (Y / T / X / outcome_name /
#'                  treatment_name / covariate_names)
#'   * assay_mat:   numeric matrix (features × samples) or NULL (no assay)
#'   * rowdata:     rowData(se) DataFrame or NULL
#'   * assay_spec:  list(assay_name, experiments, assay_cols, assay_pattern)
#'
#' @keywords internal
#' @noRd
.cf_fit_step1_input_validation <- function(formula, data,
                                           assay_name = 1,
                                           experiments = NULL,
                                           assay_cols = NULL,
                                           assay_pattern = NULL,
                                           # named-args bag, used to
                                           # materialize SE assay rownames into
                                           # colData so parse_cf_formula can
                                           # resolve assay-gene names at fit time.
                                           named_args = NULL) {
  if (!is.null(assay_cols) && !is.null(assay_pattern)) {
    rlang::abort(
      "Specify either `assay_cols` or `assay_pattern`, not both.",
      class = "cfomics_input_conflict"
    )
  }

  input_type <- if (is_se(data)) {
    "SE"
  } else if (is_mae(data)) {
    "MAE"
  } else if (is.data.frame(data)) {
    "data.frame"
  } else {
    rlang::abort(
      sprintf("Unsupported data class: %s", paste(class(data), collapse = "/")),
      class = "cfomics_unsupported_data_class"
    )
  }

  rowdata <- NULL

  if (input_type == "SE") {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      rlang::abort(
        c("Package 'SummarizedExperiment' is required for SE input.",
          i = "BiocManager::install('SummarizedExperiment')"),
        class = "cfomics_se_package_missing"
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
          i = "cfomics requires Y in colData only.",
          ">" = sprintf("Move '%s' to colData(se).", outcome_name)
        ),
        class = "cfomics_y_in_assay"
      )
    }
    col_data <- as.data.frame(SummarizedExperiment::colData(data))
    rownames(col_data) <- colnames(data)
    # For SE named-args input, X_rnaseq / X_mirna / ... reference assay
    # rownames, so materialize them into col_data before parse_cf_formula.
    # assay_cols is propagated into assay_spec because the newdata path
    # (`.extract_X_from_bioc_newdata`) reads it.
    materialize <- .materialize_se_assay_columns(col_data, assay_mat, named_args)
    col_data <- materialize$col_data
    if (length(materialize$assay_cols) > 0L && is.null(assay_cols)) {
      assay_cols <- materialize$assay_cols
    }
    parsed <- parse_cf_formula(formula, col_data)
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
        class = "cfomics_mae_package_missing"
      )
    }
    # sampleMap handling lives in bioc_integration.R::.extract_mae_assay
    # (MAE-specific logic is consolidated in bioc_integration.R).
    extracted <- .extract_mae_assay(data, formula, experiments = experiments,
                                    named_args = named_args)
    assay_mat <- extracted$assay_mat
    parsed    <- extracted$parsed
    rowdata   <- extracted$rowdata
    mae_complete_primaries <- extracted$complete_primaries
    mae_primary_is_rnaseq  <- extracted$primary_is_rnaseq
  } else {
    # data.frame path
    parsed <- parse_cf_formula(formula, data)
    assay_mat <- .resolve_df_assay(data, assay_cols = assay_cols,
                                   assay_pattern = assay_pattern,
                                   parsed = parsed)
  }

  list(
    input_type = input_type,
    parsed     = parsed,
    assay_mat  = assay_mat,
    rowdata    = rowdata,
    mae_complete_primaries = if (exists("mae_complete_primaries", inherits = FALSE))
      mae_complete_primaries else NULL,
    # TRUE when the MAE's normalised assay is the RNA-seq layer; drives
    # whether cf_fit Step 2 applies log-normalisation for MAE input.
    mae_primary_is_rnaseq = if (exists("mae_primary_is_rnaseq", inherits = FALSE))
      mae_primary_is_rnaseq else NULL,
    # assay_spec carries the **resolved** assay_cols (the materialized
    # rownames for SE named-args, the user-supplied cols/pattern for df). This
    # is what `.extract_X_from_bioc_newdata` reads at predict time to know
    # which columns are assay-derived.
    assay_spec = list(
      assay_name    = assay_name,
      experiments   = experiments,
      assay_cols    = assay_cols,
      assay_pattern = assay_pattern
    )
  )
}


#' Resolve assay matrix for data.frame inputs
#'
#' For the data.frame path, assay columns are specified via `assay_cols`
#' (character / integer / logical) or `assay_pattern` (regex). When both are
#' NULL the input is treated as clinical-only and `preprocess = TRUE` is a
#' no-op (with a warning).
#'
#' @keywords internal
#' @noRd
.resolve_df_assay <- function(data, assay_cols = NULL, assay_pattern = NULL,
                              parsed = NULL) {
  if (!is.null(assay_cols)) {
    if (is.character(assay_cols)) {
      missing_cols <- setdiff(assay_cols, names(data))
      if (length(missing_cols) > 0) {
        rlang::abort(
          c(
            "assay_cols references column(s) absent from data:",
            sprintf("  %s", paste(missing_cols, collapse = ", "))
          ),
          class = "cfomics_assay_cols_missing"
        )
      }
    }
    sub <- data[, assay_cols, drop = FALSE]
  } else if (!is.null(assay_pattern)) {
    cols <- grep(assay_pattern, names(data), value = TRUE)
    if (length(cols) == 0) {
      rlang::abort(
        c(
          sprintf("assay_pattern '%s' matched 0 column(s).", assay_pattern),
          i = "Check the regex against names(data)."
        ),
        class = "cfomics_assay_pattern_no_match"
      )
    }
    sub <- data[, cols, drop = FALSE]
  } else {
    # df with both assay_cols / assay_pattern NULL = clinical-only path; there
    # is no assay matrix to normalise. The raw-counts assumption warning in
    # Step 2 fires only when the user explicitly specifies an assay.
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
          i = "Remove from one side (formula X or assay_cols/pattern)."
        ),
        class = "cfomics_column_conflict"
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
#' Robinson & Oshlack (2010, \emph{Genome Biology} 11:R25) via edgeR (Robinson,
#' McCarthy & Smyth 2010, \emph{Bioinformatics} 26(1):139-140). \code{filterByExpr}
#' is \emph{not} applied here -- low-count filtering is a differential-
#' expression concern, kept out of the causal feature matrix. Microarray
#' inputs are already normalised externally and should be passed with
#' \code{preprocess = FALSE}.
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
      class = "cfomics_df_preprocess_warning"
    )
  }
  assay_mat <- as.matrix(assay_mat)

  # Zero-variance features carry no signal and break downstream scaling.
  row_var <- apply(assay_mat, 1L, stats::var)
  keep <- is.finite(row_var) & row_var > 0
  if (any(!keep)) assay_mat <- assay_mat[keep, , drop = FALSE]
  filter_kept <- as.integer(sum(keep))

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    rlang::warn(
      c("`edgeR` not installed; log-CPM computed without TMM normalisation.",
        i = "Install edgeR for TMM library-size normalisation."),
      class = "cfomics_edger_missing"
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


#' Newdata log-CPM: per-sample independent normalisation.
#'
#' Mirrors \code{.apply_lognorm}'s TMM -> log-CPM on a counts matrix, but
#' for newdata at predict time (no train statistics are reused; log-CPM is
#' a per-sample operation).
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
#' Batch correction runs \emph{after} normalisation on the log-CPM matrix,
#' using \code{sva::ComBat} (Johnson 2007) for both RNA-seq and microarray
#' inputs. Correcting the log values (rather than a negative-binomial batch
#' model on counts followed by log-transformation) keeps normalisation and
#' batch correction in a single, well-defined order.
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
      class = "cfomics_sva_missing"
    )
  }
  if (length(batch_vec) != ncol(assay_mat)) {
    rlang::abort(
      sprintf(
        "Batch vector length (%d) != assay sample count (%d).",
        length(batch_vec), ncol(assay_mat)
      ),
      class = "cfomics_batch_length_mismatch"
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
      class = "cfomics_pca_insufficient_dim"
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

#' Step 4: Pathway aggregation (GSVA::gsva)
#'
#' `feature_method = "pathway"` maps the all-gene assay to pathway scores
#' (default Hallmark 50 via msigdbr H, Homo sapiens). When the number of
#' pathways exceeds n_features, the top n_features by variance are selected.
#'
#' GSVA 1.50+ uses the `gsvaParam`-based API; 1.49 and earlier use the legacy
#' direct API. `utils::packageVersion` branches to support both.
#'
#' Methods: GSVA = Hänzelmann, Castelo & Guinney (2013), \emph{BMC
#' Bioinformatics} 14:7. ssGSEA = Barbie et al. (2009), \emph{Nature}
#' 462(7269):108-112 (both provided via the GSVA package).
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
      class = "cfomics_gsva_missing"
    )
  }
  gene_sets_source <- "custom"
  if (is.null(gene_sets)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      rlang::abort(
        c("Package 'msigdbr' is required when gene_sets = NULL.",
          i = "install.packages('msigdbr')"),
        class = "cfomics_msigdbr_missing"
      )
    }
    h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    gene_sets <- split(h$gene_symbol, h$gs_name)
    gene_sets_source <- "hallmark_msigdbr"
    rlang::inform(
      "Using default gene_sets = Hallmark 50 (msigdbr H, Homo sapiens).",
      class = "cfomics_default_gene_sets"
    )
  }

  pathway_matrix <- if (utils::packageVersion("GSVA") >= "1.50") {
    param <- GSVA::gsvaParam(
      exprData = as.matrix(assay_mat),
      geneSets = gene_sets,
      kcdf = "auto"
    )
    GSVA::gsva(param)
  } else {
    GSVA::gsva(as.matrix(assay_mat), gene_sets, method = pathway_method)
  }

  selected_names <- rownames(pathway_matrix)
  if (is.finite(n_features) && n_features < nrow(pathway_matrix)) {
    vars <- apply(pathway_matrix, 1, stats::var, na.rm = TRUE)
    top_idx <- order(vars, decreasing = TRUE)[seq_len(n_features)]
    pathway_matrix <- pathway_matrix[top_idx, , drop = FALSE]
    selected_names <- rownames(pathway_matrix)
  }

  list(
    pathway_matrix   = pathway_matrix,
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
        class = "cfomics_single_level_factor"
      )
      cols_to_drop <- c(cols_to_drop, col)
      next
    }
    if (n_levels > 100L) {
      rlang::warn(
        sprintf("Factor '%s' has %d levels (> 100, high-dim risk).",
                col, n_levels),
        class = "cfomics_high_cardinality_factor"
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
        class = "cfomics_cevae_binary_x"
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

#' Identify omics layer columns from cf_model meta
#'
#' Returns a named list of column-name vectors per omics layer, sourced from
#' either the named-args path (`meta$named_args$X_<layer>`) or the legacy
#' formula path (`meta$assay_spec$assay_cols / assay_pattern`). Used by
#' `.extract_X_from_bioc_newdata()` to apply layer-specific preprocessing
#' (normalisation only on RNA-seq, etc.) to newdata.
#'
#' Layer keys: baseline / rnaseq / mirna / methyl / protein / mutation / cnv /
#' assays (= X_assays generic list). Empty layers come back as
#' `character(0)`.
#'
#' @param meta `object$meta` from a fitted cf_model.
#' @return Named list with 8 elements (one per layer + assays).
#' @keywords internal
#' @noRd
.identify_omics_columns <- function(meta) {
  empty_result <- list(
    baseline = character(0),
    rnaseq   = character(0), mirna   = character(0),
    methyl   = character(0), protein = character(0),
    mutation = character(0), cnv     = character(0),
    assays   = list()
  )

  # named-args path takes priority
  na <- meta$named_args
  if (!is.null(na)) {
    return(list(
      baseline = if (is.null(na$X_baseline)) character(0) else as.character(na$X_baseline),
      rnaseq   = if (is.null(na$X_rnaseq))   character(0) else as.character(na$X_rnaseq),
      mirna    = if (is.null(na$X_mirna))    character(0) else as.character(na$X_mirna),
      methyl   = if (is.null(na$X_methyl))   character(0) else as.character(na$X_methyl),
      protein  = if (is.null(na$X_protein))  character(0) else as.character(na$X_protein),
      mutation = if (is.null(na$X_mutation)) character(0) else as.character(na$X_mutation),
      cnv      = if (is.null(na$X_cnv))      character(0) else as.character(na$X_cnv),
      assays   = if (is.null(na$X_assays))   list() else na$X_assays
    ))
  }

  # Legacy formula path
  spec <- meta$assay_spec
  cov_names <- meta$covariate_names
  if (!is.null(spec) && !is.null(cov_names)) {
    assay_cols <- spec$assay_cols
    if (is.null(assay_cols) && !is.null(spec$assay_pattern)) {
      assay_cols <- grep(spec$assay_pattern, cov_names, value = TRUE)
    }
    if (!is.null(assay_cols)) {
      assay_cols <- as.character(assay_cols)
      baseline_cols <- setdiff(cov_names, assay_cols)
      out <- empty_result
      out$baseline <- baseline_cols
      out$rnaseq   <- assay_cols  # legacy assumes RNA-seq as primary
      return(out)
    }
  }

  # Fallback: treat everything as baseline
  out <- empty_result
  if (!is.null(cov_names)) out$baseline <- as.character(cov_names)
  out
}


#' Materialize assay rownames into colData for SE named-args input
#'
#' When the user calls `cf_fit(Y = "Y", T = "T", X_rnaseq = c("G1", ...),
#' data = se)` the named-args facade builds the formula `Y ~ T | G1 + G2 + ...`,
#' but `parse_cf_formula()` evaluates it against `colData(se)` only — so
#' `G1`/`G2` (assay rownames) are not visible and the parse fails with
#' "object 'G1' not found". This helper copies the requested assay genes
#' from the assay matrix into `col_data` as numeric columns so the formula
#' resolves; it also fills in `assay_spec$assay_cols` so the newdata path
#' (`.extract_X_from_bioc_newdata`) knows which columns came from the
#' assay.
#'
#' Multi-layer behaviour (X_rnaseq + X_mirna + ... combined): SE holds a
#' single primary assay (at most one X_<omics> dedicated slot is allowed for
#' SE), so only the matching slot is materialized. The MAE path is handled
#' separately in `.cf_fit_step1_input_validation()`'s MAE branch via
#' `.extract_mae_assay()`.
#'
#' @param col_data data.frame produced by `as.data.frame(colData(se))`.
#' @param assay_mat features x samples matrix from `assay(se, assay_name)`.
#' @param named_args structured named-args bag (or NULL = legacy path).
#' @return list(col_data = updated_df, assay_cols = character(0|N) used).
#' @keywords internal
#' @noRd
.materialize_se_assay_columns <- function(col_data, assay_mat, named_args) {
  if (is.null(named_args) || is.null(assay_mat)) {
    return(list(col_data = col_data, assay_cols = character(0)))
  }
  candidate_slots <- c("X_rnaseq", "X_mirna", "X_methyl",
                       "X_protein", "X_mutation", "X_cnv")
  assay_cols_used <- character(0)
  for (slot in candidate_slots) {
    slot_val <- named_args[[slot]]
    if (is.null(slot_val) || length(slot_val) == 0L) next
    slot_val <- as.character(slot_val)
    missing_in_assay <- setdiff(slot_val, rownames(assay_mat))
    if (length(missing_in_assay) > 0L) {
      rlang::abort(
        c(
          sprintf("`%s` references %d row(s) absent from the SE assay matrix:",
                  slot, length(missing_in_assay)),
          sprintf("  %s", paste(utils::head(missing_in_assay, 5L), collapse = ", ")),
          i = "Check rownames(assay(se, assay_name)) against the column names you passed."
        ),
        class = "cfomics_se_assay_rows_missing"
      )
    }
    # Materialize each requested gene as a numeric column in col_data
    for (gene in slot_val) {
      col_data[[gene]] <- as.numeric(assay_mat[gene, ])
    }
    assay_cols_used <- c(assay_cols_used, slot_val)
  }
  list(col_data = col_data, assay_cols = assay_cols_used)
}

# ---------------------------------------------------------------------------
# Shared preprocessing pipeline
# ---------------------------------------------------------------------------

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
#'   `batch_vec` / `input_type` / `assay_spec` / `rowdata` / `factor_levels` /
#'   `clinical_cov_names` / `feature_method` / `n_features` (consumed by the
#'   cf_fit meta-populate path).
#'
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

#' @keywords internal
#' @noRd
.cf_preprocess_input <- function(formula, data,
                                 named_args     = NULL,
                                 assay_name     = 1,
                                 experiments    = NULL,
                                 assay_cols     = NULL,
                                 assay_pattern  = NULL,
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
    assay_name    = assay_name,
    experiments   = experiments,
    assay_cols    = assay_cols,
    assay_pattern = assay_pattern,
    # Forward the named-args bag so SE named-args fit can materialize assay
    # rownames into colData.
    named_args    = named_args
  )
  parsed     <- step1$parsed
  assay_mat  <- step1$assay_mat
  input_type <- step1$input_type
  rowdata    <- step1$rowdata
  assay_spec <- step1$assay_spec

  # Pathway aggregation needs a single gene-expression assay; it is not
  # defined for the multi-assay MAE feature matrix.
  if (identical(input_type, "MAE") && identical(feature_method, "pathway")) {
    rlang::abort(
      c("feature_method = 'pathway' is not supported for MultiAssayExperiment input.",
        i = "Pathway aggregation requires a single gene-expression assay; use it with SE or data.frame RNA-seq input."),
      class = "cfomics_mae_pathway_unsupported"
    )
  }

  # Step 2: Omics normalisation (zero-var drop -> TMM -> log-CPM). When
  # assay_mat = NULL (clinical-only path) this is a no-op passthrough. For MAE
  # only the RNA-seq assay is log-normalised; when the MAE has no RNA-seq
  # assay the covariate-carrier assay passes through untouched.
  assay_preprocess <- preprocess &&
    !(identical(input_type, "MAE") && !isTRUE(step1$mae_primary_is_rnaseq))
  step2 <- .apply_lognorm(
    assay_mat,
    input_type = input_type,
    preprocess = assay_preprocess
  )

  # Step 3: Batch correction on the log-CPM values (sva::ComBat). When
  # `batch != NULL` the batch vector is resolved per input_type.
  batch_vec <- NULL
  if (!is.null(batch)) {
    batch_vec <- .resolve_batch_vector(data, batch, input_type, step1)
  }
  step3 <- .apply_combat(
    step2$assay,
    batch_vec   = batch_vec,
    batch_label = batch
  )

  # Step 4: Feature select / aggregation (feature_method spans 4 values).
  # `.feature_select_step` handles variance/pca; `.apply_pathway_aggregation`
  # handles pathway.
  step4_select <- .feature_select_step(
    step3$assay,
    feature_method = feature_method,
    n_features     = n_features
  )
  step4_pathway <- if (feature_method == "pathway") {
    .apply_pathway_aggregation(
      step3$assay,
      gene_sets      = gene_sets,
      pathway_method = pathway_method,
      n_features     = n_features
    )
  } else {
    list(
      pathway_matrix = NULL, gene_sets_source = NULL,
      pathway_method = NULL, n_pathways = 0L, selected_names = NULL
    )
  }

  # Step 5 preparation: Build final X = cbind(X_clinical, t(processed_features)).
  final_features <- if (feature_method == "pathway") {
    step4_pathway$pathway_matrix
  } else {
    # "none" / "variance" / "pca": already handled by step4_select
    step4_select$features
  }
  feat_names <- if (feature_method == "pathway") {
    step4_pathway$selected_names
  } else {
    step4_select$selected_names
  }
  # Backward-compat aliases for the meta-populate path that references the
  # step3_5 / step4 variable names.
  step3_5 <- step4_select
  step4   <- step4_pathway

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
  Y <- parsed$Y
  T <- parsed$T

  # (c) Central covariate z-score: standardise the assembled X once, before
  # fold splitting / method dispatch, so every scale-sensitive backend (base
  # learners, catenets, ganite, cevae) sees a common representation. Binary /
  # one-hot columns are kept on their native scale. This is an X-only
  # (outcome-free) transform, so a single global fit is leakage-safe; the train
  # statistics are persisted for newdata via meta$x_zscore. grf / bcf are
  # scale-invariant, so the standardised X leaves their output unchanged.
  x_zscore <- if (ncol(X) > 0L) .cf_compute_x_zscore(X) else NULL
  if (!is.null(x_zscore)) X <- .cf_apply_x_zscore(X, x_zscore)

  list(
    X = X, T = T, Y = Y, cov_names = cov_names,
    parsed = parsed, step1 = step1, step2 = step2, step3 = step3,
    step3_5 = step3_5, step4 = step4, batch_vec = batch_vec,
    input_type = input_type, assay_spec = assay_spec, rowdata = rowdata,
    factor_levels = factor_levels, clinical_cov_names = clinical_cov_names,
    feature_method = feature_method, n_features = n_features,
    x_zscore = x_zscore
  )
}
