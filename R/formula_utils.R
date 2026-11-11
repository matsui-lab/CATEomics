#' Bring the treatment column to a numeric 0/1 vector.
#'
#' A logical, and a factor or character whose values are exactly "0" and "1",
#' name their own arms, so they are converted. Any other factor or character
#' does not: `as.integer()` on a factor returns the level codes, which are
#' 1-based, so a two-level factor would silently become 1/2 and every arm-wise
#' fit would look for a control group that is not there. Which level is the
#' treated one is the caller's to state, not this function's to guess, so that
#' case is refused with the recoding spelled out.
#'
#' @param T The treatment column as supplied.
#' @param treatment_name Its name, for the message.
#' @return Numeric 0/1 vector.
#' @keywords internal
#' @noRd
.cf_coerce_treatment <- function(T, treatment_name = "T") {
  if (is.logical(T)) return(as.integer(T))
  if (is.factor(T) || is.character(T)) {
    chr <- as.character(T)
    vals <- sort(unique(chr[!is.na(chr)]))
    if (all(vals %in% c("0", "1"))) return(as.integer(chr))
    rlang::abort(
      c(sprintf("`%s` must be a numeric 0/1 vector; %s was supplied.",
                treatment_name, if (is.factor(T)) "a factor" else "a character vector"),
        x = sprintf("Values found: %s.", paste(shQuote(vals), collapse = ", ")),
        i = sprintf("Name the treated arm, e.g. %s = as.integer(%s == %s).",
                    treatment_name, treatment_name,
                    shQuote(if (length(vals)) vals[[length(vals)]] else "treated"))),
      class = "cateomics_treatment_not_binary"
    )
  }
  uT <- sort(unique(T))
  if (!all(uT %in% c(0, 1))) {
    rlang::abort(
      c(sprintf("`%s` must be binary (0/1) or logical.", treatment_name),
        x = sprintf("Values found: %s.",
                    paste(utils::head(format(uT), 5L), collapse = ", "))),
      class = "cateomics_treatment_not_binary"
    )
  }
  as.integer(T)
}

#' Parse CATEomics formula `Y ~ T | X1 + X2 + ...`
#'
#' Uses \pkg{Formula} to enforce the two-part RHS syntax explicitly: the
#' treatment goes on the LHS of `|`, the covariates on the RHS. Formulas
#' without `|` (e.g. `Y ~ T + X1`) are rejected at parse time.
#'
#' @param allow_no_covariate When the covariates arrive as an omics block
#'   rather than as columns of `data`, the formula names none of them and its
#'   covariate part is empty on purpose. The callers that build an assay matrix
#'   pass `TRUE`; everyone else keeps the default, so a formula with nothing
#'   after `|` and no omics block is still refused.
#'
#' @keywords internal
#' @noRd
parse_cf_formula <- function(formula, data, allow_no_covariate = FALSE) {
  all_v <- all.vars(formula)
  for (v in all_v) {
    if (v %in% names(data) && any(is.na(data[[v]]))) {
      rlang::abort(sprintf("NA values found in variable '%s'. Please remove or impute missing values.", v), class = "cateomics_na_in_variable")
    }
  }

  f <- Formula::Formula(formula)

  lhs_len <- length(f)[1L]
  rhs_len <- length(f)[2L]
  if (lhs_len != 1L) {
    rlang::abort(
      paste0(
        "CATEomics requires a single outcome variable on the LHS.\n",
        "Got ", lhs_len, " left-hand-side parts; expected one outcome."
      ),
      class = "cateomics_formula_lhs"
    )
  }
  if (rhs_len == 1L) {
    rlang::abort(
      paste0(
        "Internal: the canonical formula has a single-part right-hand side; ",
        "expected `T | X`."
      ),
      class = "cateomics_formula_rhs_missing"
    )
  }
  if (rhs_len > 2L) {
    rlang::abort(
      paste0(
        "Internal: the canonical formula has ", rhs_len,
        " right-hand-side parts; expected exactly two, `T | X`."
      ),
      class = "cateomics_formula_rhs_extra"
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
      class = "cateomics_unknown_column"
    )
  }

  mf <- stats::model.frame(f, data = data)

  outcome_name <- all.vars(stats::formula(f, lhs = 1L, rhs = 0L))
  if (length(outcome_name) != 1L) {
    rlang::abort(
      "CATEomics requires exactly one outcome variable on the LHS.",
      class = "cateomics_formula_outcome"
    )
  }
  Y <- mf[[outcome_name]]
  if (!is.numeric(Y)) {
    rlang::abort("Outcome Y must be numeric.", class = "cateomics_outcome_not_numeric")
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
      class = "cateomics_y_type_warning"
    )
  }

  # RHS part 1: treatment. Must be a single variable.
  treatment_vars <- all.vars(stats::formula(f, lhs = 0L, rhs = 1L))
  if (length(treatment_vars) != 1L) {
    rlang::abort(
      paste0(
        "CATEomics requires exactly one treatment variable on the LHS of `|`.\n",
        "Got ", length(treatment_vars), ": ",
        paste(treatment_vars, collapse = ", "), ".\n",
        "Example: `Y ~ T | X1 + X2`."
      ),
      class = "cateomics_formula_treatment"
    )
  }
  treatment_name <- treatment_vars
  T <- mf[[treatment_name]]

  T <- .cf_coerce_treatment(T, treatment_name)

  # RHS part 2: covariates.
  covariate_names <- all.vars(stats::formula(f, lhs = 0L, rhs = 2L))
  if (length(covariate_names) == 0L && !allow_no_covariate) {
    rlang::abort(
      "CATEomics requires at least one covariate on the RHS of `|`.",
      class = "cateomics_no_covariate"
    )
  }
  if (treatment_name %in% covariate_names) {
    rlang::abort(
      sprintf(
        "Treatment variable '%s' cannot also appear among covariates.",
        treatment_name
      ),
      class = "cateomics_treatment_in_covariates"
    )
  }
  if (outcome_name %in% covariate_names) {
    rlang::abort(
      sprintf(
        "Outcome variable '%s' cannot also appear among covariates.",
        outcome_name
      ),
      class = "cateomics_outcome_in_covariates"
    )
  }
  X <- mf[, covariate_names, drop = FALSE]

  # factor / character / logical covariates are one-hot expanded downstream by
  # `.expand_factor_columns` in cf_fit Step 5, so no hard error here; a
  # non-numeric X passes through as a data.frame and cf_fit expands it.
  # This returns a matrix when every column is numeric and a data.frame only
  # when factor / character / logical columns are present; the caller
  # dispatches on the type.

  # Allow dots in variable names (R syntactic identifiers). Identifier
  # judgment matches base R: `gene.A`, `gene_A`, `geneA1` -> OK;
  # `gene-A`, `5gene`, `gene[1]`, names with whitespace -> NG.
  # This is only a warning, since the names are backtick-wrapped when the
  # canonical formula is built, so Formula::Formula() parses fine even with
  # hyphenated identifiers (real gene symbols like `HLA-B`). It is still worth
  # saying, because such a name comes back quoted in the design matrix and in
  # every table built from it.
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
      class = "cateomics_nonsyntactic_names"
    )
  }

  n <- length(Y)
  if (n == 0) {
    rlang::abort("Data has no observations.", class = "cateomics_no_observations")
  }
  if (n < 10) {
    rlang::abort("Insufficient sample size. At least 10 observations required.", class = "cateomics_insufficient_n")
  }

  n_treated <- sum(T == 1)
  n_control <- sum(T == 0)
  if (n_treated == 0 || n_control == 0) {
    rlang::abort("Both treatment and control groups must have observations.", class = "cateomics_group_empty")
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


#' Build CATEomics formula from named arguments
#'
#' Converts the user-facing named-arguments API
#' (`Y = "outcome", T = "treatment", X_baseline = c(...), omics = ...`)
#' into the canonical `Y ~ T | X1 + X2 + ...` two-part formula that
#' `parse_cf_formula()` already consumes. This lets the entire downstream
#' pipeline (cf_fit Step 1-7, base learners, dispatch, prediction) stay
#' unchanged: the named-args facade lives at the entry point only.
#'
#' Only explicitly named omics features join the formula RHS: they are
#' materialised into colData upstream so `parse_cf_formula()` can resolve them,
#' and are stripped again once the assay matrix takes over. A whole layer
#' (`"all"`) never enters the formula, and for a data.frame the omics block is
#' the assay matrix rather than formula covariates.
#'
#' @param Y,T Character scalars naming the outcome / treatment variable.
#' @param X_baseline Character vector of non-omics covariate names.
#' @param omics_spec Resolved `omics` spec from `.cf_resolve_omics()`.
#'
#' @return A list with two elements:
#'   - `formula`: the canonical two-part formula `as.formula("Y ~ T | X1 + ...")`
#'   - `structured`: the structured named-args bag (Y, T, X_baseline, omics)
#'     for downstream meta-recording; read to drive omics-aware behaviour.
#'
#' @keywords internal
#' @noRd
.build_formula_from_named <- function(Y = NULL, T = NULL,
                                       X_baseline = NULL,
                                       omics_spec = NULL) {
  # ---- Y / T validation ----
  if (is.null(Y) || !is.character(Y) || length(Y) != 1L || is.na(Y) ||
      !nzchar(Y)) {
    rlang::abort(
      "Argument `Y` must be a single non-empty character string naming the outcome variable.",
      class = "cateomics_invalid_Y"
    )
  }
  if (is.null(T) || !is.character(T) || length(T) != 1L || is.na(T) ||
      !nzchar(T)) {
    rlang::abort(
      "Argument `T` must be a single non-empty character string naming the treatment variable.",
      class = "cateomics_invalid_T"
    )
  }
  if (!is.null(X_baseline) && !is.character(X_baseline)) {
    rlang::abort(
      sprintf("Argument `X_baseline` must be a character vector of column names (got %s).",
              class(X_baseline)[1L]),
      class = "cateomics_invalid_X_slot"
    )
  }

  # Only explicitly listed omics features enter the formula; a whole layer
  # ("all") travels down the assay path instead, where feature selection
  # applies. For a data.frame the omics block never enters the formula: those
  # columns are the assay matrix.
  omics_cols <- if (is.null(omics_spec) ||
                    identical(omics_spec$input_type, "data.frame")) {
    character(0)
  } else {
    unlist(omics_spec$features, use.names = FALSE)
  }
  all_X <- c(X_baseline, omics_cols)
  all_X <- all_X[!is.na(all_X) & nzchar(all_X)]
  # An omics block supplies the design without naming anything in the formula:
  # a whole layer travels down the assay path, and for a data.frame the omics
  # columns are the assay matrix. So the covariate part is empty on purpose
  # whenever a layer is in play, and the formula is written with an empty
  # second part (the covariates join the design matrix later, in Step 5, which
  # already has a branch for no clinical columns at all). Only a call with no
  # omics block left -- `omics = "none"`, or its default on a data.frame -- has
  # genuinely nothing to fit on.
  has_omics_block <- !is.null(omics_spec) && length(omics_spec$layers) > 0L
  if (length(all_X) == 0L && !has_omics_block) {
    rlang::abort(
      c("No covariates: the call leaves the omics block out, so there is nothing to fit on.",
        i = "Name the clinical columns in `X_baseline`, or take an omics layer (`omics = \"all\"`, a layer name, or a feature vector)."),
      class = "cateomics_no_X_supplied"
    )
  }
  if (Y %in% all_X) {
    rlang::abort(
      sprintf("Outcome variable `%s` cannot also appear among X columns.", Y),
      class = "cateomics_Y_in_X"
    )
  }
  if (T %in% all_X) {
    rlang::abort(
      sprintf("Treatment variable `%s` cannot also appear among X columns.", T),
      class = "cateomics_T_in_X"
    )
  }
  dup_X <- duplicated(all_X)
  if (any(dup_X)) {
    rlang::abort(
      sprintf("Duplicate column name(s) across `X_baseline` / `omics`: %s",
              paste(unique(all_X[dup_X]), collapse = ", ")),
      class = "cateomics_duplicate_X"
    )
  }

  # ---- Build formula. Backticks keep non-syntactic feature names safe. ----
  # With nothing to name, the second part is written as `1` so the formula
  # keeps the two-part shape the parser requires while naming no covariate.
  bt <- function(x) sprintf("`%s`", x)
  formula_str <- sprintf(
    "%s ~ %s | %s",
    bt(Y), bt(T),
    if (length(all_X) == 0L) "1"
    else paste(vapply(all_X, bt, character(1)), collapse = " + ")
  )
  # `baseenv()`, not the calling frame: this formula is only ever evaluated
  # against an explicit `data` argument (`parse_cf_formula()` checks every
  # variable against `names(data)` before `model.frame()` sees it), so it needs
  # no enclosing scope. A calling frame would be carried along by every object
  # the user receives, input container and all.
  formula_obj <- stats::as.formula(formula_str, env = baseenv())

  list(
    formula = formula_obj,
    structured = list(
      Y          = Y,
      T          = T,
      X_baseline = X_baseline,
      # Which layers, and the features named in each. Whether a layer is
      # treated as RNA-seq is recorded with that layer, in the preprocessing
      # record, so there is one place to ask.
      omics      = if (is.null(omics_spec)) NULL else omics_spec$layers
    )
  )
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
      class = "cateomics_invalid_propensity_X"
    )
  }
  if (is.null(x_names) || any(!nzchar(x_names))) {
    rlang::abort(
      paste0("Cannot resolve `propensity_X`: X has no column names. ",
             "This typically means the upstream preprocessing dropped names; ",
             "consider passing `propensity_X` only with named matrices or ",
             "via the named-args entry to cf_fit."),
      class = "cateomics_propensity_X_no_names"
    )
  }
  missing <- setdiff(propensity_X, x_names)
  if (length(missing) > 0L) {
    shown <- utils::head(x_names, 10L)
    rlang::abort(
      c(sprintf("`propensity_X` columns not present in X: %s.",
                paste(missing, collapse = ", ")),
        i = sprintf("Available: %s%s.", paste(shown, collapse = ", "),
                    if (length(x_names) > length(shown)) ", ..." else ""),
        i = paste("Preprocessing expands factors, selects features and can",
                  "replace the omics block with components, so the names are",
                  "the ones it produced.")),
      class = "cateomics_propensity_X_missing"
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


#' Refuse a propensity covariate subset that is not in the design matrix.
#'
#' `propensity_X` restricts two propensities -- the estimator's and the
#' evaluation nuisance's -- so a name that preprocessing did not produce (a
#' factor level before expansion, a gene dropped by feature selection, `PC1`
#' when no PCA ran) is refused at the entry, where the caller can see which
#' names exist, rather than deep inside a nuisance fit.
#'
#' @param propensity_X Character vector or NULL.
#' @param x_names Column names of the preprocessed design matrix.
#' @keywords internal
#' @noRd
.cf_check_propensity_X <- function(propensity_X, x_names) {
  if (is.null(propensity_X)) return(invisible(NULL))
  # The resolver already refuses a non-character value and a name the matrix
  # does not have, so it is asked here on an empty matrix carrying those names
  # rather than a second copy of the same two rules.
  invisible(.resolve_propensity_X(
    matrix(numeric(0), nrow = 0L, ncol = length(x_names),
           dimnames = list(NULL, x_names)),
    propensity_X, x_names))
}

#' Warn when `propensity_X` is supplied to a method that cannot use it.
#'
#' The PS-unused / PS-architecturally-fixed group spans 6 methods:
#' - PS not part of the estimator: `slearner`, `tlearner`
#' - PS architecturally fused into the network: `dragonnet`, `cevae`, `ganite`,
#'   `tarnet`
#'
#' For these methods the estimator cannot honour propensity_X, so we warn --
#' naming that it still restricts the evaluation nuisance, which the proxy
#' losses and the evaluation metrics read -- and return TRUE so callers can
#' branch (skip subsetting or fall through to the full-X path).
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
      c(sprintf("`propensity_X` does not reach the estimator for method = '%s'.",
                method),
        i = "This method does not use a separable propensity score (PS not in the estimator, or NN architecture-fixed).",
        i = "This method's fit proceeds with the full X covariate set.",
        i = "It still restricts the evaluation nuisance, which the proxy losses and the evaluation metrics read."),
      class = "cateomics_propensity_X_ignored"
    )
    return(TRUE)
  }
  FALSE
}

#' The class of one condition, for naming it in a report
#'
#' @param cnd A condition.
#' @return The condition's own class, or `""` when it has none.
#' @keywords internal
#' @noRd
.cf_condition_class <- function(cnd) {
  cl <- setdiff(class(cnd),
                c("warning", "message", "condition", "simpleWarning",
                  "simpleMessage", "rlang_warning", "rlang_message",
                  "rlang_condition"))
  if (length(cl)) cl[[1L]] else ""
}

#' A key that collapses a repeat but keeps a different wording
#'
#' A candidate is fitted once per fold, so the same warning arrives several
#' times and only one copy is worth reporting. The class alone is the wrong
#' key: one class covers several messages -- the Super Learner reports a
#' collapsed slot by naming the slot in the text -- and keying on it would
#' keep whichever arrived last. The pair keeps every distinct wording and
#' collapses only exact repeats.
#'
#' @param cnd A condition.
#' @return A single string.
#' @keywords internal
#' @noRd
.cf_condition_key <- function(cnd) {
  paste(.cf_condition_class(cnd), conditionMessage(cnd), sep = "\r")
}

#' Evaluate a call, letting each distinct condition through only once
#'
#' The fold loop raises a fold's warning once per fold. The caller asked for
#' one fit, so they hear about it once; a wording that differs between folds
#' is a different warning and still comes through.
#'
#' @param expr The call to evaluate.
#' @return The value of `expr`.
#' @keywords internal
#' @noRd
.cf_dedup_conditions <- function(expr) {
  seen <- character(0)
  withCallingHandlers(
    expr,
    warning = function(cnd) {
      k <- .cf_condition_key(cnd)
      if (k %in% seen) invokeRestart("muffleWarning")
      seen <<- c(seen, k)
    },
    message = function(cnd) {
      k <- .cf_condition_key(cnd)
      if (k %in% seen) invokeRestart("muffleMessage")
      seen <<- c(seen, k)
    }
  )
}

#' The covariates a propensity is fitted on when none were named
#'
#' A propensity model built on every column would take the omics block with
#' it. The features are measured before treatment, so they are not colliders,
#' but a propensity fitted on thousands of them is unstable where it matters,
#' and treatment is rarely assigned on them, so the adjustment buys little.
#' The default therefore fits the propensity on the clinical covariates and
#' leaves the omics block to the outcome models.
#'
#' With no clinical covariate there is nothing to fall back to but the full
#' set; the caller is told, since
#' the alternative is to refuse a dataset that has no clinical variable at all.
#'
#' @param propensity_X What the caller named, or NULL.
#' @param clinical_cov_names The clinical columns of the design matrix, after
#'   factor expansion.
#' @param propensity_score,eval_propensity_score Known scores, if any.
#' @return A character vector of column names, or NULL for every column.
#' @keywords internal
#' @noRd
.cf_default_propensity_X <- function(propensity_X, clinical_cov_names,
                                     propensity_score = NULL,
                                     eval_propensity_score = NULL) {
  if (!is.null(propensity_X)) return(propensity_X)
  # Both scores known means no propensity is fitted anywhere, so there is
  # nothing to restrict and nothing to report.
  if (!is.null(propensity_score) && !is.null(eval_propensity_score)) return(NULL)
  clin <- as.character(clinical_cov_names %||% character(0))
  if (length(clin)) return(clin)
  rlang::warn(
    c("The propensity is being fitted on every covariate, omics included.",
      i = paste("This fit has no clinical covariate to restrict it to:",
                "every column of the design matrix came from an omics layer."),
      i = "Narrow it with `propensity_X`, or give the known value in `propensity_score` (and `eval_propensity_score`).",
      i = "A randomised trial has a known value: `propensity_score = 0.5`."),
    class = "cateomics_propensity_X_fallback"
  )
  NULL
}


#' Resolve the named arguments into the canonical form.
#'
#' Single entry point used by `cf_fit` and `cf_compare`. It maps the names the
#' caller gave (`Y` + `T` + `X_baseline` + `omics`) onto the two-part formula
#' that the rest of the pipeline parses, and returns that formula together with
#' the structured bag the result `meta` records.
#'
#' The formula is an internal normal form, not an entry point: it is written
#' with backticks and bound to `baseenv()`, and the omics block never appears
#' in it (a whole layer travels down the assay path instead).
#'
#' @keywords internal
#' @noRd
.resolve_cf_api <- function(Y = NULL, T = NULL,
                            X_baseline = NULL,
                            omics = NULL,
                            omics_type = NULL,
                            data = NULL) {
  # `omics` alone cannot stand in for the outcome and the treatment: a call has
  # to name those two whatever else it supplies.
  if (is.null(Y) && is.null(T)) {
    rlang::abort(
      # `omics` is deliberately left out of this list: it accompanies the names
      # and does not replace them, so naming it here would ask again for an
      # argument the caller may already have supplied.
      c("Must supply the names of the outcome and the treatment: `Y` and `T`.",
        i = "Example: cf_fit(data = df, Y = 'outcome', T = 'treatment', X_baseline = c('X1','X2'), method = 'grf')"),
      class = "cateomics_invalid_input"
    )
  }
  spec <- .cf_resolve_omics(omics, omics_type, data = data,
                            Y = Y, T = T, X_baseline = X_baseline)
  built <- .build_formula_from_named(Y = Y, T = T, X_baseline = X_baseline,
                                     omics_spec = spec)
  built$omics_spec <- spec
  built
}
