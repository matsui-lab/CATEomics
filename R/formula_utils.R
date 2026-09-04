#' Parse cfomics formula `Y ~ T | X1 + X2 + ...`
#'
#' Uses \pkg{Formula} to enforce the two-part RHS syntax explicitly: the
#' treatment goes on the LHS of `|`, the covariates on the RHS. Formulas
#' without `|` (e.g. `Y ~ T + X1`) are rejected at parse time, since
#' `all.vars(formula)[-1]` would silently mis-parse them and produce wrong
#' treatment/covariate assignments.
#'
#' @keywords internal
#' @noRd
parse_cf_formula <- function(formula, data) {
  all_v <- all.vars(formula)
  for (v in all_v) {
    if (v %in% names(data) && any(is.na(data[[v]]))) {
      rlang::abort(sprintf("NA values found in variable '%s'. Please remove or impute missing values.", v), class = "cfomics_na_in_variable")
    }
  }

  f <- Formula::Formula(formula)

  lhs_len <- length(f)[1L]
  rhs_len <- length(f)[2L]
  if (lhs_len != 1L) {
    rlang::abort(
      paste0(
        "cfomics requires a single outcome variable on the LHS.\n",
        "Got ", lhs_len, " LHS parts; expected `Y ~ T | X1 + X2 + ...`."
      ),
      class = "cfomics_formula_lhs"
    )
  }
  if (rhs_len == 1L) {
    rlang::abort(
      paste0(
        "cfomics requires the two-part formula `Y ~ T | X1 + X2 + ...`.\n",
        "Got a single-part RHS; separate the treatment and covariates with `|`.\n",
        "Example: `Y ~ T | X1 + X2` (not `Y ~ T + X1 + X2`)."
      ),
      class = "cfomics_formula_rhs_missing"
    )
  }
  if (rhs_len > 2L) {
    rlang::abort(
      paste0(
        "cfomics requires exactly two RHS parts `T | X`; got ", rhs_len, ".\n",
        "Example: `Y ~ T | X1 + X2`."
      ),
      class = "cfomics_formula_rhs_extra"
    )
  }

  # All formula variables must be columns of `data`; a mistyped outcome /
  # treatment / covariate otherwise surfaces as a bare model.frame error.
  # `.` (dot RHS) is expanded by model.frame, so exclude it from the check.
  missing_cols <- setdiff(all_v, c(names(data), "."))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      sprintf("Column(s) not found in data: %s.",
              paste(missing_cols, collapse = ", ")),
      class = "cfomics_unknown_column"
    )
  }

  mf <- stats::model.frame(f, data = data)

  outcome_name <- all.vars(stats::formula(f, lhs = 1L, rhs = 0L))
  if (length(outcome_name) != 1L) {
    rlang::abort(
      "cfomics requires exactly one outcome variable on the LHS.",
      class = "cfomics_formula_outcome"
    )
  }
  Y <- mf[[outcome_name]]
  if (!is.numeric(Y)) {
    rlang::abort("Outcome Y must be numeric.", class = "cfomics_outcome_not_numeric")
  }

  # Y-type advisory: a 0/1 binary outcome is a supported first-class mode (risk
  # difference) resolved later via `outcome_type`, so it is not warned about
  # here. A count-like outcome has no dedicated mode and is treated as
  # continuous, which is worth flagging.
  Y_ok <- stats::na.omit(Y)
  n_unique_Y <- length(unique(Y_ok))
  binary_Y <- n_unique_Y == 2L && all(Y_ok %in% c(0, 1))
  if (!binary_Y && n_unique_Y > 2L && n_unique_Y <= 20L &&
      all(Y_ok == as.integer(Y_ok)) && all(Y_ok >= 0)) {
    rlang::warn(
      paste0(
        "Y appears count-like (", n_unique_Y, " unique non-negative integers) ",
        "and will be treated as continuous; consider whether that is intended."
      ),
      class = "cfomics_y_type_warning"
    )
  }

  # RHS part 1: treatment. Must be a single variable.
  treatment_vars <- all.vars(stats::formula(f, lhs = 0L, rhs = 1L))
  if (length(treatment_vars) != 1L) {
    rlang::abort(
      paste0(
        "cfomics requires exactly one treatment variable on the LHS of `|`.\n",
        "Got ", length(treatment_vars), ": ",
        paste(treatment_vars, collapse = ", "), ".\n",
        "Example: `Y ~ T | X1 + X2`."
      ),
      class = "cfomics_formula_treatment"
    )
  }
  treatment_name <- treatment_vars
  T <- mf[[treatment_name]]

  if (is.logical(T)) T <- as.integer(T)
  uT <- sort(unique(T))
  if (!all(uT %in% c(0, 1))) {
    rlang::abort("Treatment must be binary (0/1) or logical.", class = "cfomics_treatment_not_binary")
  }

  # RHS part 2: covariates.
  covariate_names <- all.vars(stats::formula(f, lhs = 0L, rhs = 2L))
  if (length(covariate_names) == 0L) {
    rlang::abort(
      "cfomics requires at least one covariate on the RHS of `|`.",
      class = "cfomics_no_covariate"
    )
  }
  if (treatment_name %in% covariate_names) {
    rlang::abort(
      sprintf(
        "Treatment variable '%s' cannot also appear among covariates.",
        treatment_name
      ),
      class = "cfomics_treatment_in_covariates"
    )
  }
  if (outcome_name %in% covariate_names) {
    rlang::abort(
      sprintf(
        "Outcome variable '%s' cannot also appear among covariates.",
        outcome_name
      ),
      class = "cfomics_outcome_in_covariates"
    )
  }
  X <- mf[, covariate_names, drop = FALSE]

  # factor / character / logical covariates are one-hot expanded downstream by
  # `.expand_factor_columns` in cf_fit Step 5, so no hard error here; a
  # non-numeric X passes through as a data.frame and cf_fit expands it. The
  # non-cf_fit paths (`as_cf_data.SE/MAE`) use X directly, so this returns a
  # matrix when all-numeric and a data.frame only when factor/character/logical
  # columns are present (the caller dispatches on type).

  # Allow dots in variable names (R syntactic identifiers). Identifier
  # judgment matches base R: `gene.A`, `gene_A`, `geneA1` -> OK;
  # `gene-A`, `5gene`, `gene[1]`, names with whitespace -> NG.
  # This is only a warning, since the named-args facade backtick-wraps all
  # column names so Formula::Formula() parses fine even with hyphenated
  # identifiers (real gene symbols like `HLA-B`). The warning is still useful
  # to flag legacy-formula users who might be surprised by downstream backtick
  # handling.
  all_names <- c(outcome_name, treatment_name, covariate_names)
  if (!identical(make.names(all_names, unique = FALSE), all_names)) {
    bad <- all_names[make.names(all_names) != all_names]
    rlang::warn(
      c(
        sprintf("Variable name(s) include non-syntactic identifiers: %s",
                paste(utils::head(bad, 5L), collapse = ", ")),
        i = "Formula parsing proceeds via backtick-quoted names (Formula::Formula).",
        i = "If you want pure syntactic names, run: colnames(data) <- make.names(colnames(data), unique = TRUE)."
      ),
      class = "cfomics_nonsyntactic_names",
      .frequency = "regularly",
      .frequency_id = "cfomics_nonsyntactic_names"
    )
  }

  n <- length(Y)
  if (n == 0) {
    rlang::abort("Data has no observations.", class = "cfomics_no_observations")
  }
  if (n < 10) {
    rlang::abort("Insufficient sample size. At least 10 observations required.", class = "cfomics_insufficient_n")
  }

  n_treated <- sum(T == 1)
  n_control <- sum(T == 0)
  if (n_treated == 0 || n_control == 0) {
    rlang::abort("Both treatment and control groups must have observations.", class = "cfomics_group_empty")
  }

  # Return a matrix when X is all-numeric; return a data.frame when it
  # contains factor / character / logical columns (cf_fit Step 5 one-hot
  # expands via `.expand_factor_columns` and records factor_levels in meta).
  is_numeric_X <- vapply(X, is.numeric, logical(1))
  X_out <- if (all(is_numeric_X)) as.matrix(X) else X

  list(
    Y = as.numeric(Y),
    T = as.integer(T),
    X = X_out,
    outcome_name    = outcome_name,
    treatment_name  = treatment_name,
    covariate_names = covariate_names
  )
}


#' Build cfomics formula from named arguments
#'
#' Converts the user-facing named-arguments API
#' (`Y = "outcome", T = "treatment", X_baseline = c(...), X_rnaseq = c(...),
#'   ..., X_assays = list(atac = "ATAC"))`)
#' into the legacy `Y ~ T | X1 + X2 + ...` two-part formula that
#' `parse_cf_formula()` already consumes. This lets the entire downstream
#' pipeline (cf_fit Step 1-7, base learners, dispatch, prediction) stay
#' unchanged: the named-args facade lives at the entry point only.
#'
#' For df input the column names are returned as-is. For SummarizedExperiment /
#' MultiAssayExperiment the column resolution against assays / experiments /
#' colData is the responsibility of `.resolve_X_columns()`; here we just
#' concatenate the supplied character vectors and trust that they refer to
#' resolvable column names downstream.
#'
#' @param Y,T Character scalars naming the outcome / treatment variable.
#' @param X_baseline,X_rnaseq,X_mirna,X_methyl,X_protein,X_mutation,X_cnv
#'   Character vectors of column names per omics dedicated slot. Any may be
#'   `NULL`. At least one X argument (any of the above or `X_assays`) must be
#'   non-empty.
#' @param X_assays Optional named list for generic / future omics
#'   (`list(atac = c("ATAC1", ...))`). The values are flattened into the
#'   formula RHS just like the dedicated slots.
#'
#' @return A list with two elements:
#'   - `formula`: the canonical two-part formula `as.formula("Y ~ T | X1 + ...")`
#'   - `structured`: the structured named-args bag (Y, T, X_*) for downstream
#'     meta-recording; read to drive omics-aware behaviour.
#'
#' @keywords internal
#' @noRd
.build_formula_from_named <- function(Y = NULL, T = NULL,
                                       X_baseline = NULL,
                                       X_rnaseq = NULL, X_mirna = NULL,
                                       X_methyl = NULL, X_protein = NULL,
                                       X_mutation = NULL, X_cnv = NULL,
                                       X_assays = NULL) {
  # ---- Y / T validation ----
  if (is.null(Y) || !is.character(Y) || length(Y) != 1L || is.na(Y) ||
      !nzchar(Y)) {
    rlang::abort(
      "Argument `Y` must be a single non-empty character string naming the outcome variable.",
      class = "cfomics_invalid_Y"
    )
  }
  if (is.null(T) || !is.character(T) || length(T) != 1L || is.na(T) ||
      !nzchar(T)) {
    rlang::abort(
      "Argument `T` must be a single non-empty character string naming the treatment variable.",
      class = "cfomics_invalid_T"
    )
  }

  # ---- X_* validation ----
  X_char_slots <- list(
    X_baseline = X_baseline,
    X_rnaseq   = X_rnaseq,
    X_mirna    = X_mirna,
    X_methyl   = X_methyl,
    X_protein  = X_protein,
    X_mutation = X_mutation,
    X_cnv      = X_cnv
  )
  for (slot_name in names(X_char_slots)) {
    slot_val <- X_char_slots[[slot_name]]
    if (!is.null(slot_val) && !is.character(slot_val)) {
      rlang::abort(
        sprintf("Argument `%s` must be a character vector of column names (got %s).",
                slot_name, class(slot_val)[1L]),
        class = "cfomics_invalid_X_slot"
      )
    }
  }
  if (!is.null(X_assays)) {
    if (!is.list(X_assays) || is.null(names(X_assays)) ||
        any(!nzchar(names(X_assays)))) {
      rlang::abort(
        "Argument `X_assays` must be a named list, e.g., `list(atac = c('ATAC1', 'ATAC2'))`.",
        class = "cfomics_invalid_X_assays"
      )
    }
    for (a_name in names(X_assays)) {
      if (!is.character(X_assays[[a_name]])) {
        rlang::abort(
          sprintf("`X_assays[['%s']]` must be a character vector of column names.",
                  a_name),
          class = "cfomics_invalid_X_assays"
        )
      }
    }
  }

  # ---- Concatenate all X column names ----
  all_X <- c(
    X_baseline, X_rnaseq, X_mirna, X_methyl, X_protein, X_mutation, X_cnv,
    if (!is.null(X_assays)) unlist(X_assays, use.names = FALSE) else NULL
  )
  all_X <- all_X[!is.na(all_X) & nzchar(all_X)]
  if (length(all_X) == 0L) {
    rlang::abort(
      paste0("At least one X column must be supplied via X_baseline / X_rnaseq / ",
             "X_mirna / X_methyl / X_protein / X_mutation / X_cnv / X_assays."),
      class = "cfomics_no_X_supplied"
    )
  }
  if (Y %in% all_X) {
    rlang::abort(
      sprintf("Outcome variable `%s` cannot also appear among X columns.", Y),
      class = "cfomics_Y_in_X"
    )
  }
  if (T %in% all_X) {
    rlang::abort(
      sprintf("Treatment variable `%s` cannot also appear among X columns.", T),
      class = "cfomics_T_in_X"
    )
  }
  dup_X <- duplicated(all_X)
  if (any(dup_X)) {
    rlang::abort(
      sprintf("Duplicate column name(s) in X slots: %s",
              paste(unique(all_X[dup_X]), collapse = ", ")),
      class = "cfomics_duplicate_X"
    )
  }

  # ---- Build formula. Use backticks for syntactic-identifier safety
  # (gene names with dots are already allowed by parse_cf_formula; backticks
  # make non-syntactic names round-trip safely too). ----
  bt <- function(x) sprintf("`%s`", x)
  formula_str <- sprintf(
    "%s ~ %s | %s",
    bt(Y), bt(T),
    paste(vapply(all_X, bt, character(1)), collapse = " + ")
  )
  formula_obj <- stats::as.formula(formula_str, env = parent.frame())

  list(
    formula = formula_obj,
    structured = list(
      Y          = Y,
      T          = T,
      X_baseline = X_baseline,
      X_rnaseq   = X_rnaseq,
      X_mirna    = X_mirna,
      X_methyl   = X_methyl,
      X_protein  = X_protein,
      X_mutation = X_mutation,
      X_cnv      = X_cnv,
      X_assays   = X_assays,
      X_all      = all_X
    )
  )
}


#' Resolve X columns from named arguments based on data class.
#'
#' The facade `.build_formula_from_named()` just concatenates all X_*
#' character vectors into the formula RHS — it does not know whether the
#' caller passed a data.frame / SummarizedExperiment / MultiAssayExperiment.
#' This adds the data-class-aware validation layer so users get a clear error
#' when they combine an X slot with the wrong input type.
#'
#' Rules:
#' - data.frame: any of `X_baseline` / `X_rnaseq` / ... / `X_cnv` may be
#'   used (they are just column-name vectors). `X_assays` is **rejected**
#'   (data.frame has no assay concept).
#' - SummarizedExperiment: at most ONE of the six dedicated `X_<omics>` slots
#'   may be used (SE holds a single primary assay). `X_assays` is rejected.
#' - MultiAssayExperiment: all six dedicated `X_<omics>` slots may be used
#'   in parallel, plus `X_assays` for additional experiments. Most flexible.
#'
#' Column-existence verification (e.g., "`X_baseline = 'foo'` but `foo` is
#' not a column of `data`") is delegated to the downstream
#' `parse_cf_formula()` / `as_cf_data()` path which already runs a full
#' `model.frame` validation; this only enforces the slot-shape rules above.
#'
#' @param data User-supplied data object.
#' @param X_rnaseq,X_mirna,X_methyl,X_protein,X_mutation,X_cnv
#'   Character vectors per omics dedicated slot. Any may be NULL.
#' @param X_assays Optional named list of additional assay slots.
#'
#' @return Invisibly TRUE on success. Errors with a `cfomics_*` class on
#'   validation failure.
#'
#' @keywords internal
#' @noRd
.resolve_X_columns <- function(data,
                                X_rnaseq   = NULL, X_mirna    = NULL,
                                X_methyl   = NULL, X_protein  = NULL,
                                X_mutation = NULL, X_cnv      = NULL,
                                X_assays   = NULL) {
  dedicated <- list(
    X_rnaseq   = X_rnaseq,   X_mirna    = X_mirna,
    X_methyl   = X_methyl,   X_protein  = X_protein,
    X_mutation = X_mutation, X_cnv      = X_cnv
  )
  n_dedicated <- sum(vapply(dedicated, function(x) !is.null(x), logical(1)))
  dedicated_used <- names(dedicated)[vapply(dedicated, function(x) !is.null(x), logical(1))]

  is_df  <- is.data.frame(data)
  is_se  <- !is_df && inherits(data, "SummarizedExperiment")
  is_mae <- !is_df && !is_se && inherits(data, "MultiAssayExperiment")

  if (!is_df && !is_se && !is_mae) {
    rlang::abort(
      sprintf("Unsupported data class: %s. Expected data.frame, SummarizedExperiment, or MultiAssayExperiment.",
              paste(class(data), collapse = ", ")),
      class = "cfomics_invalid_data_type"
    )
  }

  # data.frame: assay concept does not exist
  if (is_df && !is.null(X_assays)) {
    rlang::abort(
      c("Argument `X_assays` is not supported for data.frame input.",
        i = "data.frame has no assay concept; pass column names via X_rnaseq / X_mirna / etc. directly, or switch to SE / MAE input."),
      class = "cfomics_X_assays_unsupported_df"
    )
  }

  # SummarizedExperiment: single primary assay
  if (is_se && n_dedicated > 1L) {
    rlang::abort(
      c(sprintf("SE input supports at most one X_<omics> dedicated argument (got %d: %s).",
                n_dedicated, paste(dedicated_used, collapse = ", ")),
        i = "For multi-assay support, switch to MultiAssayExperiment (MAE) input."),
      class = "cfomics_multiple_omics_unsupported_se"
    )
  }
  if (is_se && !is.null(X_assays)) {
    rlang::abort(
      c("Argument `X_assays` is not supported for SummarizedExperiment input.",
        i = "SE holds a single primary assay; for additional assays switch to MAE input."),
      class = "cfomics_X_assays_unsupported_se"
    )
  }

  invisible(TRUE)
}


#' Resolve a `propensity_X` argument against the covariate matrix.
#'
#' Validates that every name in `propensity_X` exists as a column of `X` (or
#' in `cov_names` when `X` is an unnamed matrix produced by an upstream
#' one-hot / variance-selection step) and returns both an R-style subset
#' matrix (for R-native PS-using methods) and 0-based column indices (for
#' Python ColumnSubsetPropensityModel).
#'
#' `propensity_X = NULL` (default) means "no subset, use full X for PS"; we
#' return the full X and `seq_len(ncol(X)) - 1` as a no-op.
#'
#' @param X Numeric matrix or data.frame.
#' @param propensity_X Character vector of column names. NULL → no-op.
#' @param cov_names Optional column names override (used when `X` is an
#'   unnamed matrix; `colnames(X)` is preferred when available).
#'
#' @return list with:
#'   - `X_ps`        : numeric matrix to use for PS estimation
#'   - `indices_0b`  : 0-based integer column indices (Python convention)
#'   - `subset_names`: character vector of the resolved column names
#'   - `is_subset`   : logical, TRUE only when propensity_X was non-NULL
#'
#' @keywords internal
#' @noRd
.resolve_propensity_X <- function(X, propensity_X = NULL, cov_names = NULL) {
  x_names <- colnames(X)
  if (is.null(x_names) || length(x_names) == 0L) {
    x_names <- cov_names
  }
  if (is.null(propensity_X)) {
    n_col <- if (!is.null(dim(X))) ncol(X) else length(x_names)
    return(list(
      X_ps         = X,
      indices_0b   = seq_len(n_col) - 1L,
      subset_names = x_names %||% paste0("X", seq_len(n_col)),
      is_subset    = FALSE
    ))
  }
  if (!is.character(propensity_X) || length(propensity_X) == 0L) {
    rlang::abort(
      "`propensity_X` must be a non-empty character vector of column names.",
      class = "cfomics_invalid_propensity_X"
    )
  }
  if (is.null(x_names) || any(!nzchar(x_names))) {
    rlang::abort(
      paste0("Cannot resolve `propensity_X`: X has no column names. ",
             "This typically means the upstream preprocessing dropped names; ",
             "consider passing `propensity_X` only with named matrices or ",
             "via the named-args entry to cf_fit."),
      class = "cfomics_propensity_X_no_names"
    )
  }
  missing <- setdiff(propensity_X, x_names)
  if (length(missing) > 0L) {
    rlang::abort(
      sprintf("`propensity_X` columns not present in X: %s",
              paste(missing, collapse = ", ")),
      class = "cfomics_propensity_X_missing"
    )
  }
  idx <- match(propensity_X, x_names)
  X_ps <- if (is.matrix(X)) X[, idx, drop = FALSE] else X[, idx, drop = FALSE]
  list(
    X_ps         = X_ps,
    indices_0b   = as.integer(idx - 1L),
    subset_names = propensity_X,
    is_subset    = TRUE
  )
}


#' Warn when `propensity_X` is supplied to a method that cannot use it.
#'
#' The PS-unused / PS-architecturally-fixed group spans 6 methods:
#' - PS not part of the estimator: `slearner`, `tlearner`
#' - PS architecturally fused into the network: `dragonnet`, `cevae`, `ganite`,
#'   `tarnet`
#'
#' For these methods, propensity_X cannot be honoured; we warn (with the
#' actionable hint that the argument is silently ignored) and return TRUE
#' so callers can branch (skip subsetting or fall through to the full-X path).
#'
#' @param method String method name.
#' @param propensity_X Character vector or NULL.
#' @return Logical: TRUE if a warning fired (i.e., method ignores
#'   propensity_X); FALSE otherwise.
#'
#' @keywords internal
#' @noRd
.warn_propensity_X_ignored <- function(method, propensity_X) {
  if (is.null(propensity_X)) return(FALSE)
  ignored_methods <- c(
    # PS not used by design
    "slearner", "tlearner",
    # PS architecturally fixed (NN with internal PS head sharing all X)
    "dragonnet", "cevae", "ganite", "tarnet"
  )
  if (method %in% ignored_methods) {
    rlang::warn(
      c(sprintf("`propensity_X` ignored for method = '%s'.", method),
        i = "This method does not use a separable propensity score (PS not in the estimator, or NN architecture-fixed).",
        i = "The fit proceeds with the full X covariate set."),
      class = "cfomics_propensity_X_ignored",
      .frequency = "regularly", .frequency_id = paste0("propensity_X_", method)
    )
    return(TRUE)
  }
  FALSE
}


#' Resolve the named-args / formula API dual entry.
#'
#' Single entry point used by `cf_fit` and `cf_compare` to map
#' whichever API the caller supplied (legacy `formula` or named `Y` + `T` +
#' `X_*`) onto the legacy formula that the rest of the pipeline already
#' understands. Returns the canonical formula plus the structured-args bag for
#' downstream meta-recording (the structured bag is `NULL` when the caller went
#' through the legacy formula path).
#'
#' Validation rules:
#' - Both supplied → `cfomics_input_mode_conflict` (cannot guess intent).
#' - Neither supplied → `cfomics_invalid_input` (need at least one).
#' - `formula` supplied → emit a deprecation warning (once per session) and
#'   pass formula through.
#' - Named args supplied → build the formula via `.build_formula_from_named()`.
#'
#' @keywords internal
#' @noRd
.resolve_cf_api <- function(formula = NULL,
                             Y = NULL, T = NULL,
                             X_baseline = NULL,
                             X_rnaseq = NULL, X_mirna = NULL,
                             X_methyl = NULL, X_protein = NULL,
                             X_mutation = NULL, X_cnv = NULL,
                             X_assays = NULL,
                             data = NULL) {
  has_formula <- !is.null(formula)
  has_named   <- !is.null(Y) || !is.null(T) ||
                 !is.null(X_baseline) || !is.null(X_rnaseq) ||
                 !is.null(X_mirna) || !is.null(X_methyl) ||
                 !is.null(X_protein) || !is.null(X_mutation) ||
                 !is.null(X_cnv) || !is.null(X_assays)

  if (has_formula && has_named) {
    rlang::abort(
      c("Specify either `formula` (legacy API) or named arguments (`Y`, `T`, `X_*`), not both.",
        i = "Drop one path. The named-args API is preferred."),
      class = "cfomics_input_mode_conflict"
    )
  }
  if (!has_formula && !has_named) {
    rlang::abort(
      c("Must supply either `formula` (legacy API) or named arguments (`Y`, `T`, `X_*`).",
        i = "Example (named, recommended): cf_fit(Y = 'outcome', T = 'treatment', X_baseline = c('X1','X2'), data = df, method = 'grf')",
        i = "Example (legacy): cf_fit(Y ~ T | X1 + X2, data = df, method = 'grf')"),
      class = "cfomics_invalid_input"
    )
  }

  if (has_formula) {
    rlang::warn(
      c(
        "Direct formula API is deprecated.",
        i = "Use named arguments instead: cf_fit(Y = 'outcome', T = 'treatment', X_baseline = c(...), data = ..., method = ...)",
        "*" = "It will be removed in a future release."
      ),
      class = "cfomics_formula_deprecated",
      .frequency = "once", .frequency_id = "cfomics_formula_deprecation"
    )
    return(list(formula = formula, structured = NULL))
  }

  # Data-class-aware validation of X_* slots (df / SE / MAE rules).
  if (!is.null(data)) {
    .resolve_X_columns(
      data       = data,
      X_rnaseq   = X_rnaseq,   X_mirna    = X_mirna,
      X_methyl   = X_methyl,   X_protein  = X_protein,
      X_mutation = X_mutation, X_cnv      = X_cnv,
      X_assays   = X_assays
    )
  }
  built <- .build_formula_from_named(
    Y = Y, T = T,
    X_baseline = X_baseline,
    X_rnaseq = X_rnaseq, X_mirna = X_mirna, X_methyl = X_methyl,
    X_protein = X_protein, X_mutation = X_mutation, X_cnv = X_cnv,
    X_assays = X_assays
  )
  built
}
