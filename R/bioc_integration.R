#' Convert Data to cfomics Format
#'
#' Generic function to convert various data structures to a format
#' suitable for cfomics analysis. Supports data.frame, SummarizedExperiment,
#' and MultiAssayExperiment objects.
#'
#' @param data Input data object
#' @param ... Additional arguments passed to methods
#' @return A list with components:
#'   \itemize{
#'     \item Y: numeric vector of outcomes
#'     \item T: integer vector of treatment assignments (0/1)
#'     \item X: numeric matrix of covariates
#'     \item outcome_name: character string
#'     \item treatment_name: character string
#'     \item covariate_names: character vector
#'   }
#' @export
#' @examples
#' # With a data.frame
#' df <- data.frame(Y = rnorm(20), T = rbinom(20, 1, 0.5),
#'                  X1 = rnorm(20), X2 = rnorm(20))
#' cfd <- as_cf_data(df, outcome = "Y", treatment = "T",
#'                   covariates = c("X1", "X2"))
#'
#' \donttest{
#' # With a SummarizedExperiment
#' se <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4)
#' cfd_se <- as_cf_data(se, outcome = "outcome", treatment = "treatment",
#'                      covariates = c("age", "sex"))
#' }
as_cf_data <- function(data, ...) {
  UseMethod("as_cf_data")
}

#' @rdname as_cf_data
#' @param formula Optional formula \code{Y ~ T | X1 + X2 + ...}; an alternative
#'   to the named \code{outcome} / \code{treatment} / \code{covariates}
#'   arguments (which are the recommended API).
#' @export
as_cf_data.data.frame <- function(data, formula = NULL, outcome = NULL,
                                  treatment = NULL, covariates = NULL, ...) {
  if (!is.null(formula)) {
    return(parse_cf_formula(formula, data))
  }
  if (is.null(outcome) || is.null(treatment)) {
    rlang::abort(
      c("Provide `outcome` + `treatment` column names (or a `formula`).",
        i = "e.g. as_cf_data(df, outcome = 'Y', treatment = 'T', covariates = c('X1', 'X2'))"),
      class = "cfomics_invalid_input"
    )
  }
  covs <- if (is.null(covariates)) setdiff(names(data), c(outcome, treatment))
          else covariates
  if (length(covs) == 0L) {
    rlang::abort("At least one covariate is required.",
                 class = "cfomics_no_covariate")
  }
  # Build the canonical two-part formula from the named columns and reuse the
  # parser, so the named-args output is identical to the formula path.
  bt <- function(x) sprintf("`%s`", x)
  f  <- stats::as.formula(sprintf("%s ~ %s | %s", bt(outcome), bt(treatment),
                                  paste(vapply(covs, bt, character(1)),
                                        collapse = " + ")))
  parse_cf_formula(f, data)
}

#' @rdname as_cf_data
#' @export
as_cf_data.default <- function(data, ...) {
  stop("as_cf_data() does not support objects of class: ",
       paste(class(data), collapse = ", "),
       "\nSupported classes: data.frame, SummarizedExperiment, MultiAssayExperiment",
       call. = FALSE)
}

#' Convert SummarizedExperiment to cfomics format
#'
#' @rdname as_cf_data
#' @param outcome Character, name of outcome variable in colData or assay
#' @param treatment Character, name of treatment variable in colData
#' @param covariates Character vector, names of covariate variables in colData
#' @param assay_name Name or index of assay to use (default 1)
#' @details
#' `as_cf_data()` is a \emph{pure converter}: it materialises the colData
#' covariates plus \emph{every} assay feature into a cfomics data list without
#' any log-normalisation or feature selection. Preprocessing (TMM/log-CPM,
#' batch correction, variance/PCA/pathway selection) is applied by
#' [cf_fit()] / [cf_compare()], which run the same preprocessing pipeline
#' internally — do not pre-select here.
#' @export
as_cf_data.SummarizedExperiment <- function(data,
                                            formula = NULL,
                                            outcome = NULL,
                                            treatment = NULL,
                                            covariates = NULL,
                                            assay_name = 1,
                                            ...) {
  # Check for SummarizedExperiment package
 if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required. Install with:\n",
         "  BiocManager::install('SummarizedExperiment')",
         call. = FALSE)
  }

  # Extract colData as data.frame
  col_data <- as.data.frame(SummarizedExperiment::colData(data))
  rownames(col_data) <- colnames(data)

  # If formula is provided, use formula-based extraction
  if (!is.null(formula)) {
    return(extract_from_formula_se(data, formula, col_data, assay_name))
  }

  # Manual specification requires outcome and treatment
  if (is.null(outcome) || is.null(treatment)) {
    stop("Either 'formula' or both 'outcome' and 'treatment' must be specified",
         call. = FALSE)
  }

  extract_manual_se(data, outcome, treatment, covariates, col_data,
                    assay_name)
}

#' Extract data from SummarizedExperiment using formula
#' @noRd
extract_from_formula_se <- function(data, formula, col_data, assay_name) {
  vars <- all.vars(formula)
  outcome_name <- vars[1]

  assay_data <- SummarizedExperiment::assay(data, assay_name)

  # If the outcome is a gene in the assay, materialise it in col_data so the
  # formula parser can find it.
  if (outcome_name %in% rownames(assay_data)) {
    col_data[[outcome_name]] <- as.numeric(assay_data[outcome_name, ])
  } else if (!(outcome_name %in% colnames(col_data))) {
    stop("Outcome '", outcome_name, "' not found in colData or assay rownames",
         call. = FALSE)
  }

  parsed <- parse_cf_formula(formula, col_data)

  # Pure converter: append EVERY assay feature as a covariate (no
  # log-normalisation, no selection — cf_fit / cf_compare do that). Exclude the
  # outcome gene if it lives in the assay so it is not fed back as a covariate.
  assay_rows <- setdiff(rownames(assay_data), outcome_name)
  if (length(assay_rows) > 0) {
    feature_matrix <- assay_data[assay_rows, , drop = FALSE]
    X_assay <- t(feature_matrix)
    parsed$X <- cbind(parsed$X, X_assay)
    parsed$covariate_names <- c(parsed$covariate_names, assay_rows)
  }

  parsed
}

#' Extract data from SummarizedExperiment manually
#' @noRd
extract_manual_se <- function(data, outcome, treatment, covariates, col_data,
                              assay_name) {
  assay_data <- SummarizedExperiment::assay(data, assay_name)

  # Extract outcome (Y)
  if (outcome %in% rownames(assay_data)) {
    Y <- as.numeric(assay_data[outcome, ])
  } else if (outcome %in% colnames(col_data)) {
    Y <- as.numeric(col_data[[outcome]])
  } else {
    stop("Outcome '", outcome, "' not found in assay rownames or colData",
         call. = FALSE)
  }

  # Extract treatment (T)
  if (!(treatment %in% colnames(col_data))) {
    stop("Treatment '", treatment, "' must be in colData", call. = FALSE)
  }
  T_vec <- col_data[[treatment]]

  # Convert to binary
  if (is.logical(T_vec)) T_vec <- as.integer(T_vec)
  if (is.factor(T_vec)) {
    # factor -> 0/1 follows level order (first level = 0 / reference,
    # second = 1 / treated). Surface the coding once so an unexpected level
    # order does not silently flip the treatment sign.
    .lv <- levels(T_vec)
    if (length(.lv) == 2L) {
      rlang::inform(
        sprintf("Treatment factor coded as '%s'=0 (reference), '%s'=1 (treated).",
                .lv[1], .lv[2]),
        class = "cfomics_factor_treatment_coding",
        .frequency = "once",
        .frequency_id = "cfomics_factor_treatment_coding")
    }
    T_vec <- as.integer(T_vec) - 1L
  }
  T_vec <- as.integer(T_vec)

  if (!all(unique(T_vec) %in% c(0L, 1L))) {
    stop("Treatment must be binary (0/1)", call. = FALSE)
  }

  # Extract covariates (X) from colData
  if (is.null(covariates)) {
    # Use all numeric columns except outcome and treatment
    numeric_cols <- vapply(col_data, is.numeric, logical(1))
    covariates <- names(col_data)[numeric_cols]
    covariates <- setdiff(covariates, c(outcome, treatment))
  }

  # Filter to existing columns
  covariates <- intersect(covariates, colnames(col_data))

  # Build covariate matrix
  if (length(covariates) > 0) {
    X_coldata <- as.matrix(col_data[, covariates, drop = FALSE])
  } else {
    X_coldata <- matrix(nrow = ncol(data), ncol = 0)
  }

  # Pure converter: append EVERY assay feature as a covariate (no
  # log-normalisation, no selection). Exclude the outcome gene if it lives in
  # the assay so it is not fed back as a covariate.
  assay_rows <- setdiff(rownames(assay_data), outcome)
  if (length(assay_rows) > 0) {
    feature_matrix <- assay_data[assay_rows, , drop = FALSE]
    X <- cbind(X_coldata, t(feature_matrix))
    all_covariate_names <- c(covariates, assay_rows)
  } else {
    X <- X_coldata
    all_covariate_names <- covariates
  }

  # Ensure X is numeric matrix
  X <- apply(X, 2, as.numeric)
  if (!is.matrix(X)) X <- matrix(X, ncol = 1)

  list(
    Y = Y,
    T = T_vec,
    X = X,
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_names = all_covariate_names,
    source = "SummarizedExperiment",
    n_samples = length(Y),
    n_features_from_assay = length(assay_rows)
  )
}

#' Convert MultiAssayExperiment to cfomics format
#'
#' @rdname as_cf_data
#' @param experiments Character vector of experiment names to include
#' @export
as_cf_data.MultiAssayExperiment <- function(data,
                                            formula = NULL,
                                            outcome = NULL,
                                            treatment = NULL,
                                            covariates = NULL,
                                            experiments = NULL,
                                            ...) {
  # Check for MultiAssayExperiment package
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package required. Install with:\n",
         "  BiocManager::install('MultiAssayExperiment')",
         call. = FALSE)
  }

  # Extract colData
  col_data <- as.data.frame(MultiAssayExperiment::colData(data))

  # If formula provided with colData variables only
  if (!is.null(formula)) {
    vars <- all.vars(formula)
    if (all(vars %in% colnames(col_data))) {
      return(parse_cf_formula(formula, col_data))
    }
  }

  # For multi-omics integration, need outcome and treatment
  if (is.null(outcome) || is.null(treatment)) {
    stop("For MultiAssayExperiment, 'outcome' and 'treatment' must be specified",
         call. = FALSE)
  }

  # Get experiments to use
  if (is.null(experiments)) {
    experiments <- names(MultiAssayExperiment::experiments(data))
  }

  # Extract and integrate features from each experiment
  integrated_features <- list()

  for (exp_name in experiments) {
    exp_data <- MultiAssayExperiment::experiments(data)[[exp_name]]

    # Get assay matrix depending on object type
    if (inherits(exp_data, "SummarizedExperiment")) {
      assay_mat <- SummarizedExperiment::assay(exp_data, 1)
    } else if (inherits(exp_data, "matrix")) {
      assay_mat <- exp_data
    } else {
      warning("Skipping experiment '", exp_name, "': unsupported type ",
              class(exp_data)[1])
      next
    }

    if (nrow(assay_mat) == 0) next

    # Pure converter: keep EVERY feature of each assay (no
    # log-normalisation, no selection — cf_fit / cf_compare handle that).
    selected <- assay_mat

    # Prefix feature names with experiment name
    rownames(selected) <- paste0(exp_name, "_", rownames(selected))
    integrated_features[[exp_name]] <- selected
  }

  if (length(integrated_features) == 0) {
    stop("No features extracted from experiments", call. = FALSE)
  }

  # sampleMap handling. A plain `Reduce(intersect, colnames)` ignores the MAE
  # primary <-> assay-specific colname mapping (sampleMap), so MAEs with
  # differing colnames yield an empty intersection ("No common samples").
  # Use MultiAssayExperiment::sampleMap() to join primary <-> colname and
  # subset to the primary units that are complete across all experiments.
  used_experiments <- names(integrated_features)
  sample_map <- MultiAssayExperiment::sampleMap(data)
  filtered_map <- sample_map[as.character(sample_map$assay) %in% used_experiments, ,
                             drop = FALSE]
  primary_counts <- table(as.character(filtered_map$primary))
  complete_primaries <- names(primary_counts)[
    primary_counts == length(used_experiments)
  ]

  if (length(complete_primaries) == 0) {
    stop(
      sprintf(
        paste0(
          "0 common primary samples across experiments %s after sampleMap join.\n",
          "Hint: check MultiAssayExperiment::sampleMap(mae). Some experiments may ",
          "not cover all primary samples; reduce `experiments =` to overlapping set."
        ),
        paste(used_experiments, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Align each experiment's assay column names to primary names and subset to
  # complete_primaries.
  for (exp_name in used_experiments) {
    map_sub <- filtered_map[
      as.character(filtered_map$assay) == exp_name &
        as.character(filtered_map$primary) %in% complete_primaries,
      , drop = FALSE
    ]
    colname_for_primary <- as.character(map_sub$colname)[
      match(complete_primaries, as.character(map_sub$primary))
    ]
    feat <- integrated_features[[exp_name]]
    feat <- feat[, colname_for_primary, drop = FALSE]
    colnames(feat) <- complete_primaries
    integrated_features[[exp_name]] <- feat
  }

  feature_matrix <- do.call(rbind, integrated_features)

  # Subset colData to complete primary samples (primary identifier = rowname of
  # colData(mae)).
  col_data <- col_data[complete_primaries, , drop = FALSE]

  # Extract Y, T
  if (outcome %in% colnames(col_data)) {
    Y <- as.numeric(col_data[[outcome]])
  } else {
    stop("Outcome '", outcome, "' not found in colData", call. = FALSE)
  }

  T_vec <- col_data[[treatment]]
  if (is.logical(T_vec)) T_vec <- as.integer(T_vec)
  if (is.factor(T_vec)) {
    # factor -> 0/1 follows level order (first level = 0 / reference,
    # second = 1 / treated). Surface the coding once so an unexpected level
    # order does not silently flip the treatment sign.
    .lv <- levels(T_vec)
    if (length(.lv) == 2L) {
      rlang::inform(
        sprintf("Treatment factor coded as '%s'=0 (reference), '%s'=1 (treated).",
                .lv[1], .lv[2]),
        class = "cfomics_factor_treatment_coding",
        .frequency = "once",
        .frequency_id = "cfomics_factor_treatment_coding")
    }
    T_vec <- as.integer(T_vec) - 1L
  }
  T_vec <- as.integer(T_vec)

  if (!all(unique(T_vec) %in% c(0L, 1L))) {
    stop("Treatment must be binary (0/1)", call. = FALSE)
  }

  # Covariates from colData
  if (is.null(covariates)) {
    covariates <- character(0)
  }
  covariates <- intersect(covariates, colnames(col_data))

  if (length(covariates) > 0) {
    covar_data <- as.matrix(col_data[, covariates, drop = FALSE])
    covar_data <- apply(covar_data, 2, as.numeric)
    # For a single-column covariate, apply() returns a vector, which breaks
    # the downstream cbind() shape. Re-wrap as a matrix (same guard as the SE
    # path).
    if (!is.matrix(covar_data)) {
      covar_data <- matrix(covar_data, ncol = 1,
                           dimnames = list(NULL, covariates))
    }
  } else {
    covar_data <- matrix(nrow = nrow(col_data), ncol = 0)
  }

  # Combine covariates and features
  X <- cbind(covar_data, t(feature_matrix))
  all_covariate_names <- c(covariates, rownames(feature_matrix))

  list(
    Y = Y,
    T = T_vec,
    X = X,
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_names = all_covariate_names,
    source = "MultiAssayExperiment",
    n_samples = length(Y),
    experiments_used = names(integrated_features)
  )
}

# Feature selection is not performed in the as_cf_data path. cf_fit /
# cf_compare select features via `.feature_select_step()` (preprocessing.R).

#' Extract primary assay matrix from MultiAssayExperiment for cf_fit 7-step flow
#'
#' Helper for cf_fit Step 1 (preprocessing.R
#' `.cf_fit_step1_input_validation`): extracts only the primary experiment's
#' assay matrix from an MAE and feeds it to the downstream normalisation /
#' batch-correction / feature_method pipeline.
#'
#' Differs from `as_cf_data.MultiAssayExperiment` in scope:
#'   * Returns the primary experiment assay only (experiments[1] or first
#'     experiment when NULL).
#'   * Non-primary experiments are not integrated into X; if any non-primary
#'     experiment matrix looks like raw counts (integer / non-negative numeric
#'     matrix), emit a warning prompting the user to restructure the MAE or
#'     handle externally.
#'   * Uses MultiAssayExperiment::sampleMap() to align primary <-> assay
#'     colnames.
#'
#' @param mae A MultiAssayExperiment.
#' @param formula Y ~ T | X1 + X2 + ... formula.
#' @param experiments Character vector of experiment names to use (default
#'   NULL = all; first becomes primary).
#'
#' @return List with components:
#'   * assay_mat:  features × samples matrix from primary experiment, columns
#'                 aligned to primary sample identifiers.
#'   * parsed:     output of `parse_cf_formula` against colData(mae) (Y / T /
#'                 X / outcome_name / treatment_name / covariate_names).
#'   * rowdata:    rowData of primary SummarizedExperiment-class experiment,
#'                 or NULL when primary is a plain matrix.
#'   * complete_primaries: primary sample names retained after sampleMap join.
#'   * primary_exp: name of the primary experiment used.
#'
#' @keywords internal
#' @noRd
.extract_mae_assay <- function(mae, formula, experiments = NULL,
                                # named-args bag. When supplied, assay
                                # rownames referenced by X_<layer> are
                                # materialized into colData so parse_cf_formula
                                # can resolve them at MAE fit time.
                                named_args = NULL) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    rlang::abort(
      c("Package 'MultiAssayExperiment' is required for MAE input.",
        i = "BiocManager::install('MultiAssayExperiment')"),
      class = "cfomics_mae_package_missing"
    )
  }
  if (is.null(experiments)) {
    experiments <- names(MultiAssayExperiment::experiments(mae))
  }
  if (length(experiments) == 0L) {
    rlang::abort(
      "MAE contains 0 experiments.",
      class = "cfomics_mae_empty"
    )
  }
  missing_exp <- setdiff(experiments, names(MultiAssayExperiment::experiments(mae)))
  if (length(missing_exp) > 0L) {
    rlang::abort(
      c(
        "Requested experiments are absent from the MAE:",
        sprintf("  %s", paste(missing_exp, collapse = ", "))
      ),
      class = "cfomics_mae_unknown_experiment"
    )
  }
  # Only the RNA-seq assay is internally log-normalised; every other assay
  # passes through as-is (user pre-normalised). Identify the RNA-seq
  # experiment by name (rnaseq / rna / rna_seq) or, failing that, by the
  # named-arg X_rnaseq slot whose rows live in a matching experiment. When none
  # matches, no assay is normalised (all passthrough) and experiments[1] is used
  # only as the covariate carrier.
  rnaseq_exp <- NULL
  nm_hit <- experiments[tolower(experiments) %in% c("rnaseq", "rna", "rna_seq")]
  if (length(nm_hit) > 0L) {
    rnaseq_exp <- nm_hit[1]
  } else if (!is.null(named_args) && !is.null(named_args$X_rnaseq) &&
             length(named_args$X_rnaseq) > 0L) {
    sv <- as.character(named_args$X_rnaseq)
    for (cand in experiments) {
      obj <- MultiAssayExperiment::experiments(mae)[[cand]]
      m <- if (inherits(obj, "SummarizedExperiment")) {
        SummarizedExperiment::assay(obj, 1)
      } else if (is.matrix(obj)) obj else next
      if (all(sv %in% rownames(m))) { rnaseq_exp <- cand; break }
    }
  }
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
      class = "cfomics_mae_no_overlap"
    )
  }

  # Primary assay matrix + rowData
  primary_obj <- MultiAssayExperiment::experiments(mae)[[primary_exp]]
  if (inherits(primary_obj, "SummarizedExperiment")) {
    primary_mat <- SummarizedExperiment::assay(primary_obj, 1)
    primary_rowdata <- SummarizedExperiment::rowData(primary_obj)
  } else if (is.matrix(primary_obj)) {
    primary_mat <- primary_obj
    primary_rowdata <- NULL
  } else {
    rlang::abort(
      sprintf("Unsupported primary experiment class '%s' for cf_fit.",
              class(primary_obj)[1]),
      class = "cfomics_mae_unsupported_primary"
    )
  }

  # Align columns to primary names.
  pm_map <- filtered_map[
    as.character(filtered_map$assay) == primary_exp &
      as.character(filtered_map$primary) %in% complete_primaries,
    , drop = FALSE
  ]
  primary_colnames <- as.character(pm_map$colname)[
    match(complete_primaries, as.character(pm_map$primary))
  ]
  primary_mat <- primary_mat[, primary_colnames, drop = FALSE]
  colnames(primary_mat) <- complete_primaries

  # Non-primary experiments: raw-counts auto-detection.
  for (exp_name in setdiff(experiments, primary_exp)) {
    other_obj <- MultiAssayExperiment::experiments(mae)[[exp_name]]
    other_mat <- if (inherits(other_obj, "SummarizedExperiment")) {
      SummarizedExperiment::assay(other_obj, 1)
    } else if (is.matrix(other_obj)) {
      other_obj
    } else {
      next
    }
    if (.is_raw_counts_like(other_mat)) {
      rlang::warn(
        c(
          sprintf(
            "Non-RNA-seq assay '%s' looks like raw counts (integer / non-negative).",
            exp_name
          ),
          i = "Only the RNA-seq assay is internally log-normalised; every other assay passes through as-is.",
          ">" = sprintf("Normalise '%s' externally before adding it to the MAE.", exp_name)
        ),
        class = "cfomics_mae_raw_counts_detected"
      )
    }
  }

  # Subset colData to primary sample units, then parse the formula.
  col_data <- as.data.frame(MultiAssayExperiment::colData(mae))
  col_data <- col_data[complete_primaries, , drop = FALSE]
  if (ncol(col_data) == 0L) {
    rlang::abort(
      "MAE colData has 0 columns; formula cannot be parsed.",
      class = "cfomics_mae_empty_coldata"
    )
  }

  # When called via named-args, X_rnaseq / X_mirna / ... reference assay
  # rownames in the matching experiment. Materialize them into col_data so
  # parse_cf_formula resolves them just like the SE path does in
  # `.materialize_se_assay_columns()`.
  if (!is.null(named_args)) {
    # Map dedicated slots to experiment names by name lookup; fall back to
    # the primary experiment when the layer name (e.g. "rnaseq") matches.
    layer_slot_pairs <- list(
      X_rnaseq   = c("rnaseq",   "rna",      "RNAseq", "rna_seq"),
      X_mirna    = c("mirna",    "miRNA",    "mir"),
      X_methyl   = c("methyl",   "methylation"),
      X_protein  = c("protein",  "proteomics"),
      X_mutation = c("mutation", "mut",      "snv"),
      X_cnv      = c("cnv",      "CNV")
    )
    for (slot_name in names(layer_slot_pairs)) {
      slot_val <- named_args[[slot_name]]
      if (is.null(slot_val) || length(slot_val) == 0L) next
      slot_val <- as.character(slot_val)
      # Find the experiment whose rownames contain all of slot_val
      matched_exp <- NULL
      candidates <- c(layer_slot_pairs[[slot_name]], experiments)
      for (cand in candidates) {
        if (!(cand %in% experiments)) next
        exp_obj <- MultiAssayExperiment::experiments(mae)[[cand]]
        exp_mat <- if (inherits(exp_obj, "SummarizedExperiment")) {
          SummarizedExperiment::assay(exp_obj, 1)
        } else if (is.matrix(exp_obj)) {
          exp_obj
        } else { next }
        if (all(slot_val %in% rownames(exp_mat))) {
          matched_exp <- list(name = cand, mat = exp_mat)
          break
        }
      }
      if (is.null(matched_exp)) {
        rlang::abort(
          c(
            sprintf("Cannot resolve `%s` columns against any MAE experiment.", slot_name),
            sprintf("  Requested rows: %s",
                    paste(utils::head(slot_val, 5L), collapse = ", "))
          ),
          class = "cfomics_mae_layer_resolve_failed"
        )
      }
      # Align matched_exp$mat columns to complete_primaries
      ex_map <- filtered_map[
        as.character(filtered_map$assay) == matched_exp$name &
          as.character(filtered_map$primary) %in% complete_primaries, ,
        drop = FALSE
      ]
      ex_colnames <- as.character(ex_map$colname)[
        match(complete_primaries, as.character(ex_map$primary))
      ]
      mat_aligned <- matched_exp$mat[slot_val, ex_colnames, drop = FALSE]
      colnames(mat_aligned) <- complete_primaries
      for (g in slot_val) {
        col_data[[g]] <- as.numeric(mat_aligned[g, ])
      }
    }
    # X_assays generic list
    if (!is.null(named_args$X_assays) && length(named_args$X_assays) > 0L) {
      for (a_name in names(named_args$X_assays)) {
        cols <- as.character(named_args$X_assays[[a_name]])
        if (length(cols) == 0L) next
        if (!(a_name %in% experiments)) {
          rlang::abort(
            sprintf("X_assays experiment '%s' is not in MAE experiments.", a_name),
            class = "cfomics_mae_layer_resolve_failed"
          )
        }
        exp_obj <- MultiAssayExperiment::experiments(mae)[[a_name]]
        exp_mat <- if (inherits(exp_obj, "SummarizedExperiment")) {
          SummarizedExperiment::assay(exp_obj, 1)
        } else { exp_obj }
        ex_map <- filtered_map[
          as.character(filtered_map$assay) == a_name &
            as.character(filtered_map$primary) %in% complete_primaries, ,
          drop = FALSE
        ]
        ex_colnames <- as.character(ex_map$colname)[
          match(complete_primaries, as.character(ex_map$primary))
        ]
        mat_aligned <- exp_mat[cols, ex_colnames, drop = FALSE]
        colnames(mat_aligned) <- complete_primaries
        for (g in cols) {
          col_data[[g]] <- as.numeric(mat_aligned[g, ])
        }
      }
    }
  }

  parsed <- parse_cf_formula(formula, col_data)

  # Same dedup as the SE path — when named-args materialized assay rows for
  # the primary experiment (only) into col_data, strip those rows from
  # parsed$X and subset primary_mat. Non-primary experiments' columns are
  # materialized in col_data only and are NOT subject to primary_mat
  # preprocessing (they are covariate-stacked as-is).
  if (!is.null(named_args)) {
    primary_layer_cols <- character(0)
    primary_slot_pairs <- list(
      X_rnaseq   = c("rnaseq",   "rna",      "RNAseq", "rna_seq"),
      X_mirna    = c("mirna",    "miRNA",    "mir"),
      X_methyl   = c("methyl",   "methylation"),
      X_protein  = c("protein",  "proteomics"),
      X_mutation = c("mutation", "mut",      "snv"),
      X_cnv      = c("cnv",      "CNV")
    )
    for (slot_name in names(primary_slot_pairs)) {
      slot_val <- named_args[[slot_name]]
      if (is.null(slot_val) || length(slot_val) == 0L) next
      if (primary_exp %in% primary_slot_pairs[[slot_name]] ||
          all(slot_val %in% rownames(primary_mat))) {
        primary_layer_cols <- c(primary_layer_cols, as.character(slot_val))
      }
    }
    if (length(primary_layer_cols) > 0L) {
      keep_cov <- setdiff(parsed$covariate_names, primary_layer_cols)
      if (is.data.frame(parsed$X)) {
        parsed$X <- parsed$X[, keep_cov, drop = FALSE]
      } else {
        parsed$X <- parsed$X[, keep_cov, drop = FALSE]
      }
      parsed$covariate_names <- keep_cov
      have_rows <- intersect(primary_layer_cols, rownames(primary_mat))
      primary_mat <- primary_mat[have_rows, , drop = FALSE]
    }
  }

  list(
    assay_mat = primary_mat,
    parsed = parsed,
    rowdata = primary_rowdata,
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
