#' Compare CATE methods on a common K-fold split
#'
#' @description
#' Honest cross-fit comparison of CATE-capable methods on a single
#' stratified K-fold split (Mahajan 2024 / Curth-van der Schaar 2023
#' convention). Computes per-method out-of-fold tau-hat, applies
#' selection metrics (R-loss / DR-loss / factual loss) computed against
#' a shared "evaluation nuisance" (\code{eval_nuisance}, default
#' \code{"rf"}), and returns a wide-format scoreboard plus the full-data
#' fitted models for downstream \code{\link{cf_select}} /
#' \code{\link{cf_evaluate}} / \code{\link{cf_interpret}} /
#' \code{\link{cf_rule}}.
#'
#' @details
#' \strong{Two-level nuisance:} cf_compare's evaluation nuisance
#' (\code{eval_nuisance}, default \code{"rf"}) is fit once with K-fold
#' OOF and is independent of each method's internal nuisance. The
#' evaluation nuisance feeds the selection loss, GATES/RATE/Calibration
#' DR pseudo-outcomes, and BLP-on-X regressions; per-method internal
#' nuisance estimates (e.g., DR-Learner cross-fit) are untouched and
#' remain optimal per-method.
#'
#' \strong{R-squared interpretation caveat:} \code{r_squared_score} is
#' bounded above by \code{Var(tau) / (Var(tau) + irreducible_noise)}
#' and typically sits in 0.1-0.3 for omics data. Use it for relative
#' comparison only; absolute thresholds (e.g., "0.5 = good") are not
#' meaningful.
#'
#' @inheritParams cf_fit
#' @param BPPARAM A \code{BiocParallelParam} instance controlling how the
#'   per-method fits are dispatched (default \code{BiocParallel::SerialParam()},
#'   deterministic and safe alongside the basilisk per-method Python processes;
#'   pass e.g. \code{BiocParallel::MulticoreParam()} to fit methods
#'   concurrently).
#' @param formula Formula \code{Y ~ T | X1 + X2 + ...} (cf_fit syntax).
#' @param data data.frame / SummarizedExperiment / MultiAssayExperiment.
#' @param methods Character vector of method names or \code{"all"}.
#'   \code{"all"} expands to the 10 CATE-capable methods (\code{grf,
#'   bcf, slearner, tlearner, xlearner, drlearner, cevae, ganite,
#'   tarnet, dragonnet}); cf_compare is CATE-only.
#' @param metric Character vector of metric names; one or more of
#'   \code{"r_loss"} (default), \code{"dr_loss"}, \code{"factual_loss"}.
#' @param base_learners Optional character vector of base-learner names
#'   (\code{"linear"}, \code{"logistic"}, \code{"lasso"}, \code{"ridge"},
#'   \code{"elastic_net"}, \code{"rf"}, \code{"xgboost"}; see
#'   \code{.fit_base_learner}). When supplied, the pluggable meta / DR
#'   learners (\code{slearner}, \code{tlearner}, \code{xlearner},
#'   \code{drlearner}) are expanded into one config per base learner, each
#'   setting \code{model_outcome} (and \code{model_final} for xlearner /
#'   drlearner) to that learner; grf / bcf / NN methods stay single (base
#'   ignored). Scoreboard rows are labelled \code{method_base} (e.g.
#'   \code{xlearner_rf}). The propensity model (\code{propensity_score} /
#'   \code{model_propensity}) and the evaluation nuisance
#'   (\code{eval_nuisance}) stay single across all configs. Mutually
#'   exclusive with \code{model_outcome} / \code{model_final}. Default
#'   \code{NULL} reproduces the single-config behaviour.
#' @param n_folds Integer, K for stratified K-fold. Default \code{5L}.
#' @param eval_nuisance Base-learner name for evaluation nuisance
#'   (default \code{"rf"}); see \code{.fit_base_learner}.
#' @param eval_propensity_score Scalar in (0,1) or length-n vector; used
#'   when \code{eval_nuisance$propensity = "known"} (RCT context).
#' @param propensity_score Scalar in (0,1) or length-n vector; forwarded
#'   to each per-method \code{cf_fit_<method>} call as the user-supplied
#'   known propensity score. For grf this
#'   bypasses the internal \code{regression_forest} and pins \code{W.hat};
#'   for xlearner / drlearner it pins the IPW weights.
#' @param trim Length-2 numeric \code{c(lower, upper)} propensity clipping
#'   bounds forwarded to each per-method fit (default \code{c(0.025, 0.975)}).
#' @param random_state Integer seed or NULL.
#' @param keep_kfold_models Logical, retain per-fold partial fit objects
#'   (default \code{FALSE}; memory expensive).
#' @param verbose Logical, emit progress messages.
#' @param device One of \code{"auto"} (default), \code{"cpu"}, \code{"gpu"};
#'   forwarded to each per-method \code{cf_fit} call.
#' @param ... Additional arguments forwarded to per-method fits.
#' @return Object of class \code{c("cfomics_compare", "list")} with
#'   elements \code{scoreboard, models, fold_cache, metadata}.
#' @seealso \code{\link{cf_select}}, \code{\link{cf_evaluate}},
#'   \code{\link{cf_interpret}}, \code{\link{cf_rule}}, \code{\link{cf_fit}}.
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf); known randomised
#' # propensity, and an independent evaluation nuisance so the loss favours
#' # no compared method. Runs without Python.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L,  # 1 chain keeps bcf fast
#'                   propensity_score = 0.5, random_state = 1)
#' cmp
#' @export
cf_compare <- function(formula = NULL, data = NULL, methods,
                       # ---- Facade named-args (mirrors cf_fit) ----
                       Y          = NULL, T          = NULL,
                       X_baseline = NULL, X_rnaseq   = NULL, X_mirna    = NULL,
                       X_methyl   = NULL, X_protein  = NULL, X_mutation = NULL,
                       X_cnv      = NULL, X_assays   = NULL,
                       propensity_X = NULL,
                       # ---- preprocessing controls (mirror cf_fit) ----
                       # cf_compare runs the SAME Step 1-5 preprocessing as
                       # cf_fit via `.cf_preprocess_input()`, so an SE / MAE omics
                       # assay is log-normalised and included as covariates. df
                       # is neutral.
                       assay_name     = 1,
                       experiments    = NULL,
                       assay_cols     = NULL,
                       assay_pattern  = NULL,
                       preprocess     = TRUE,
                       batch          = NULL,
                       feature_method = c("variance", "none", "pca", "pathway"),
                       n_features     = NULL,
                       gene_sets      = NULL,
                       pathway_method = c("gsva", "ssgsea", "plage", "zscore"),
                       metric         = c("r_loss", "dr_loss", "factual_loss"),
                       base_learners  = NULL,  # optional: grid the pluggable
                                               # meta/DR learners over base
                                               # learners (label = method_base)
                       outcome_type   = c("auto", "continuous", "binary"),  # binary
                                                       # Y=RD support, symmetric with cf_fit
                       n_folds        = 5L,
                       eval_nuisance  = "rf",
                       eval_propensity_score = NULL,  # known PS for
                                                       # eval_nuisance$propensity = "known"
                       propensity_score = NULL,       # known PS for the per-method fits
                                                       # (e.g., grf W.hat, xlearner
                                                       # FixedPropensityModel).
                       trim           = c(0.025, 0.975),  # PS clip forwarded to
                                                       # per-method fits
                       random_state   = NULL,
                       keep_kfold_models = FALSE,
                       verbose        = TRUE,
                       # Per-method fits are dispatched through BiocParallel. The
                       # default SerialParam() is deterministic and safe
                       # alongside the basilisk per-method Python processes; pass
                       # e.g. BiocParallel::MulticoreParam(k, RNGseed = seed) to
                       # fit methods concurrently.
                       BPPARAM        = BiocParallel::SerialParam(),
                       device         = c("auto", "cpu", "gpu"),
                       ...) {
  device <- match.arg(device)
  # The default `metric` is the full c("r_loss", "dr_loss", "factual_loss")
  # triplet so cf_select / cf_ensemble can request any of the three losses
  # without an explicit opt-in. rank uses metric[1].
  # Resolve dual entry (formula vs named-args). cf_compare
  # internally dispatches `cf_fit` per method; we resolve to a canonical
  # `formula` here and forward both `formula` + `propensity_X` so each child
  # `cf_fit` does NOT re-trigger the deprecation warning (legacy formula path
  # gets warned once at the cf_compare entry).
  .facade <- .resolve_cf_api(
    formula    = formula,
    Y          = Y, T = T,
    X_baseline = X_baseline,
    X_rnaseq   = X_rnaseq,   X_mirna    = X_mirna,
    X_methyl   = X_methyl,   X_protein  = X_protein,
    X_mutation = X_mutation, X_cnv      = X_cnv,
    X_assays   = X_assays,
    data       = data
  )
  formula <- .facade$formula
  # Named-args bag (NULL under the legacy formula path) drives SE/MAE
  # assay-column materialization inside `.cf_preprocess_input()`.
  .cfomics_named_args <- .facade$structured
  feature_method <- match.arg(feature_method)
  pathway_method <- match.arg(pathway_method)
  if (is.null(data)) {
    rlang::abort("`data` argument is required.", class = "cfomics_invalid_data")
  }
  # Input-type-aware n_features default (MAE = 50, else 100).
  n_features <- .cf_resolve_n_features(n_features, data)
  # Mirror cf_fit's guard: a known propensity (RCT) and an estimated-PS base
  # learner are mutually exclusive. cf_fit aborts on the conflict; cf_compare
  # forwards both via `...`, so enforce the same contract here (otherwise the
  # estimated PS silently wins and the known PS is ignored per method).
  if (!is.null(propensity_score) && "model_propensity" %in% names(list(...))) {
    rlang::abort(
      c("Both `propensity_score` and `model_propensity` are explicitly specified.",
        i = "For RCT (known PS), use `propensity_score = ...` only.",
        i = "For observational study (estimated PS), use `model_propensity = ...` only."),
      class = "cfomics_propensity_args_conflict"
    )
  }
  # `base_learners` (grid over base learners) is mutually exclusive with the
  # single-config `model_outcome` / `model_final` knobs: one selects the
  # comparison space, the other a single fit. Aborting (rather than silently
  # letting one win) mirrors the propensity_score / model_propensity guard.
  if (!is.null(base_learners) &&
      any(c("model_outcome", "model_final") %in% names(list(...)))) {
    rlang::abort(
      c("Both `base_learners` and `model_outcome` / `model_final` are specified.",
        i = "Use `base_learners = c(...)` to grid over base learners, OR",
        i = "`model_outcome` / `model_final` for a single-config comparison."),
      class = "cfomics_base_learners_conflict"
    )
  }
  if (missing(methods) || length(methods) == 0L) {
    rlang::abort(
      c("`methods` is required.",
        i = "Use methods = 'all' to compare all 10 CATE-capable methods, or",
        i = "supply a character vector like c('grf', 'bcf', 'drlearner')."),
      class = "cfomics_methods_required"
    )
  }
  cate_methods_all <- c("grf", "bcf", "slearner", "tlearner", "xlearner",
                        "drlearner", "cevae", "ganite", "tarnet", "dragonnet")
  # Methods whose CATE-driving regression stage(s) take a pluggable base
  # learner (model_outcome / model_final). grf (honest forest), bcf (BART) and
  # the NN methods do not, so `base_learners` leaves them single.
  pluggable_methods <- c("slearner", "tlearner", "xlearner", "drlearner")
  if (length(methods) == 1L && identical(methods, "all")) {
    methods <- cate_methods_all
  }
  bad_unknown <- setdiff(methods, cate_methods_all)
  if (length(bad_unknown) > 0L) {
    rlang::abort(
      sprintf("Unknown CATE methods: %s. Allowed: %s.",
              paste(bad_unknown, collapse = ", "),
              paste(cate_methods_all, collapse = ", ")),
      class = "cfomics_compare_unknown_method"
    )
  }
  # Validate + dedup `base_learners`; inform when no method in the set is
  # pluggable (the argument then has no effect).
  if (!is.null(base_learners)) {
    valid_base <- c("linear", "logistic", "lasso", "ridge",
                    "elastic_net", "rf", "xgboost")
    base_learners <- unique(as.character(base_learners))
    bad_base <- setdiff(base_learners, valid_base)
    if (length(bad_base) > 0L) {
      rlang::abort(
        sprintf("Unknown base_learners: %s. Allowed: %s.",
                paste(bad_base, collapse = ", "),
                paste(valid_base, collapse = ", ")),
        class = "cfomics_base_learners_unknown"
      )
    }
    if (!any(methods %in% pluggable_methods)) {
      rlang::inform(
        c(sprintf("`base_learners` has no effect: none of the methods (%s) use a pluggable base learner.",
                  paste(methods, collapse = ", ")),
          i = sprintf("Pluggable methods: %s.",
                      paste(pluggable_methods, collapse = ", "))),
        class = "cfomics_base_learners_no_pluggable"
      )
    }
  }

  # Symmetry with cf_fit: warn once per PS-free / NN method in the compare set
  # that `propensity_X` is ignored for it (cf_fit emits the same signal per fit,
  # so a mixed compare set gives the user the same notice).
  if (!is.null(propensity_X)) {
    for (m in methods) .warn_propensity_X_ignored(m, propensity_X)
  }

  # Symmetry with cf_fit: validate a known propensity score and reject a
  # per-sample (vector) PS for econml cross-fit methods (drlearner/xlearner),
  # which cannot map row-level PS through internal folds. cf_fit aborts
  # identically (cfomics_propensity_vector_unsupported).
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cfomics_invalid_propensity_score"
      )
    }
    n_obs <- if (inherits(data, c("SummarizedExperiment", "MultiAssayExperiment")))
               ncol(data) else nrow(data)
    if (length(propensity_score) != 1L && length(propensity_score) != n_obs) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match number of samples = %d.",
                length(propensity_score), n_obs),
        class = "cfomics_invalid_propensity_score"
      )
    }
    if (length(propensity_score) > 1L &&
        any(c("drlearner", "xlearner") %in% methods)) {
      rlang::abort(
        c("Vector `propensity_score` is unsupported for drlearner / xlearner.",
          i = "econml cross-fit cannot map per-sample PS through internal folds.",
          i = "Use a scalar (`propensity_score = 0.5`) or drop drlearner / xlearner."),
        class = "cfomics_propensity_vector_unsupported"
      )
    }
  }

  metric <- match.arg(metric, c("r_loss", "dr_loss", "factual_loss"),
                      several.ok = TRUE)

  # Validate `eval_nuisance` structure at entry (parity with cf_fit) instead of
  # letting a malformed value error deep inside the nuisance helpers. Accepts a
  # base-learner name or list(outcome=, propensity=) with the "known" sentinel;
  # `.cfomics_resolve_eval_nuisance` aborts on anything else.
  invisible(.cfomics_resolve_eval_nuisance(eval_nuisance))

  # Typo / unknown-arg validation. Emit a warning here -- BEFORE the per-method
  # suppressWarnings wrap below -- so user-visible typos (e.g.
  # `base_learner_outcome` instead of `model_outcome`) are surfaced rather than
  # silently swallowed.
  for (m in methods) {
    .cfomics_validate_dot_args(list(...), m)
  }

  # Shared Step 1-5 preprocessing — identical to cf_fit via
  # `.cf_preprocess_input()`. SE / MAE omics assays are log-normalised (RNA-seq)
  # and feature-selected (default "variance"), then factor / character columns
  # are one-hot expanded into the final numeric X. A data.frame without an assay
  # is neutral (clinical only). cf_fit(data) and cf_compare(data) build the
  # *same* X from the same input.
  prep <- .cf_preprocess_input(
    formula        = formula,
    data           = data,
    named_args     = .cfomics_named_args,
    assay_name     = assay_name,
    experiments    = experiments,
    assay_cols     = assay_cols,
    assay_pattern  = assay_pattern,
    preprocess     = preprocess,
    batch          = batch,
    feature_method = feature_method,
    n_features     = n_features,
    gene_sets      = gene_sets,
    pathway_method = pathway_method,
    # The preprocessing X is method-independent; `method` only drives the
    # CEVAE binary-X distribution warning, so pass "cevae" when it is in the
    # compare set to keep that warning symmetric with cf_fit.
    method         = if ("cevae" %in% methods) "cevae" else "grf"
  )
  X <- prep$X
  T <- as.integer(prep$T)
  Y <- as.numeric(prep$Y)
  n <- length(Y)

  # Binary Y = risk-difference support, symmetric with cf_fit. Detect (or honour
  # the user override) once, then thread `binary_outcome` to the shared eval
  # nuisance (classifier) and to every per-method fit (RD mode). The loss /
  # DR-metric formulas are outcome-type agnostic.
  outcome_type   <- match.arg(outcome_type)
  binary_outcome <- .resolve_binary_outcome(Y, outcome_type)

  if (!is.null(random_state)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(random_state))
  }

  folds   <- .create_stratified_folds(T, n_folds = n_folds,
                                       random_state = random_state)
  fold_id <- .folds_to_fold_id(folds, n)

  if (isTRUE(verbose)) {
    message("Fitting evaluation nuisance (", eval_nuisance,
            ") OOF on ", n_folds, "-fold split...")
  }
  nuis <- .compute_common_nuisance_oof(X, T, Y, folds,
                                        eval_nuisance = eval_nuisance,
                                        random_state = random_state,
                                        eval_propensity_score = eval_propensity_score,
                                        binary_outcome = binary_outcome)

  oof_tau    <- list()
  oof_mu0    <- list()   # per-method held-out potential outcomes (factual loss)
  oof_mu1    <- list()
  models     <- list()
  fold_cache <- if (isTRUE(keep_kfold_models)) list() else NULL
  failed     <- list()

  ## After factor pre-expansion the matrix may have more columns than
  ## `parsed$covariate_names` (which still reflects the user-facing names).
  ## Use the expanded `colnames(X)` so downstream per-method names reflect the
  ## one-hot dummies (`sexF` / `sexM`).
  cov_names <- colnames(X)
  if (is.null(cov_names) || length(cov_names) == 0L) {
    cov_names <- paste0("X", seq_len(ncol(X)))
  }

  # The per-method K-fold loop uses the matrix-signature helper
  # `.kfold_oof_predict_matrix(X, T, Y, method, ...)`, which operates on the
  # already-preprocessed X (no cf_fit recursion, so no double preprocessing).
  # cf_compare's shared nuisance (`nuis`) is reused via precomputed_nuisance,
  # and folds are generated deterministically from the same random_state.
  # The per-method fit loop runs through BiocParallel::bplapply so parallelism
  # is controlled by the caller's BPPARAM (SerialParam by default). The verbose
  # message is emitted once on the main thread (worker-side messages are lost
  # under parallel backends).
  # Expand the (method x base_learner) grid. When `base_learners` is NULL this
  # is one config per method (label = method, base = NA) = the legacy path.
  # Otherwise the pluggable methods fan out to one config per base learner
  # (label = method_base), while grf / bcf / NN stay single. `label` is the
  # uniform external key used for models / oof_tau / scoreboard rows; each
  # config still fits its BARE method (so predict / dispatch are unchanged).
  configs   <- .cfomics_build_configs(methods, base_learners, pluggable_methods)
  grid_mode <- !is.null(base_learners)

  if (isTRUE(verbose)) {
    message("Fitting ", nrow(configs), " config(s) (BPPARAM = ",
            class(BPPARAM)[1], ")...")
  }

  fit_one_config <- function(i) {
    cfg <- configs[i, , drop = FALSE]
    m   <- cfg$method
    # In grid mode, set this config's CATE-driving regressor(s) to its base
    # learner: model_outcome for all pluggable methods, plus model_final for
    # the ones with a Stage-2 CATE model (xlearner / drlearner). Non-pluggable
    # methods ignore the base. The mutual-exclusion guard guarantees the user's
    # `...` carries no model_outcome / model_final in grid mode.
    base_inject <- list()
    if (grid_mode && m %in% pluggable_methods && !is.na(cfg$base)) {
      base_inject$model_outcome <- cfg$base
      if (m %in% c("xlearner", "drlearner")) {
        base_inject$model_final <- cfg$base
      }
    }
    call_args <- c(
      list(
        X                    = X,
        T_vec                = T,
        Y_vec                = Y,
        method               = m,
        cov_names            = cov_names,
        n_folds              = n_folds,
        random_state         = random_state,
        keep_kfold_models    = isTRUE(keep_kfold_models),
        eval_nuisance        = eval_nuisance,
        precomputed_nuisance = nuis,
        device               = device,
        propensity_score     = propensity_score,
        trim                 = trim,
        propensity_X         = propensity_X,
        eval_propensity_score = eval_propensity_score,
        binary_outcome       = binary_outcome
      ),
      base_inject,
      list(...)
    )
    helper_out <- tryCatch(
      suppressWarnings(suppressMessages(
        do.call(.kfold_oof_predict_matrix, call_args)
      )),
      error = function(e) e
    )
    list(label = cfg$label, method = m, base = cfg$base, out = helper_out)
  }

  per_config <- BiocParallel::bplapply(seq_len(nrow(configs)), fit_one_config,
                                       BPPARAM = BPPARAM)

  for (cret in per_config) {
    lab        <- cret$label
    m          <- cret$method
    helper_out <- cret$out
    if (inherits(helper_out, "error")) {
      failed[[lab]] <- conditionMessage(helper_out)
      oof_tau[[lab]] <- rep(NA_real_, n)
      oof_mu0[[lab]] <- rep(NA_real_, n)
      oof_mu1[[lab]] <- rep(NA_real_, n)
      models[[lab]]  <- NULL
      next
    }
    # Wrap inner full-data fit with cf_model class + attach K-fold artifacts.
    # cf_compare does not go through cf_fit, so class + meta are built here.
    # Four-pillar: res = inner res (+ OOF tau + nuisance), model = the backend
    # container.
    .res <- helper_out$fit$res
    .res$tau_hat_oof <- helper_out$oof_ite
    .res$nuisance    <- helper_out$nuisance
    # The per-method backend returns the raw `ite_*` summary which the main
    # cf_fit path overwrites centrally via `.cf_build_summary`. cf_compare wraps
    # the backend fit directly (not via cf_fit), so rebuild the summary here too
    # — otherwise `summary.cf_model` (which reads tau_mean / tau_iqr) errors on
    # a cf_compare-derived model.
    .res$summary <- .cf_build_summary(
      .res$ate, ate_ci = .res$ate_ci, tau_hat = .res$tau_hat
    )
    full_fit <- structure(
      list(
        method = m,
        res    = .res,
        model  = helper_out$fit,
        meta   = list(
          covariate_names = cov_names,
          treatment_name  = "T",
          outcome_name    = "Y",
          partition       = "oof",
          input_type      = "data.frame",
          n               = n,
          p               = length(cov_names),
          # Carry the central X z-score (from the shared preprocessing) and the
          # NN Y z-score so newdata prediction on a cf_compare model (e.g.
          # cf_ensemble test predict) re-applies the train statistics.
          x_zscore        = prep$x_zscore,
          y_zscore        = .res$y_zscore,
          # cf_compare is CATE-only; carry the estimand so summary.cf_model /
          # print show it (cf_compare wraps the fit directly, not via cf_fit).
          estimand        = list(type = "CATE"),
          # Base-learner config for this row. `spec` lets downstream (cf_select)
          # report / reproduce the winner without parsing the compound label.
          base_learner    = cret$base,
          spec            = list(
            method        = m,
            model_outcome = if (grid_mode && m %in% pluggable_methods &&
                                !is.na(cret$base)) cret$base else NA_character_,
            model_final   = if (grid_mode && m %in% c("xlearner", "drlearner") &&
                                !is.na(cret$base)) cret$base else NA_character_
          )
        )
      ),
      class = "cf_model"
    )
    if (isTRUE(keep_kfold_models) && !is.null(helper_out$fold_cache)) {
      full_fit$fold_cache <- helper_out$fold_cache
      fold_cache[[lab]]   <- helper_out$fold_cache
    }
    oof_tau[[lab]] <- helper_out$oof_ite
    oof_mu0[[lab]] <- helper_out$oof_mu0
    oof_mu1[[lab]] <- helper_out$oof_mu1
    # Cache (Y, T, X) on the model so cf_evaluate / cf_interpret / cf_rule
    # can compute DR pseudo-outcomes without recomputation.
    if (is.null(full_fit$data)) full_fit$data <- list()
    full_fit$data$Y <- Y
    full_fit$data$T <- T
    full_fit$data$X <- X
    models[[lab]] <- full_fit
  }

  # `ids` are the uniform external keys (config labels): equal to `methods` in
  # the legacy path, or `method_base` in grid mode. Everything downstream
  # (scoreboard rows, models names, oof_tau_matrix columns, metadata$methods)
  # is keyed by `ids`; a per-row `method` / `base` decomposition is exposed in
  # the scoreboard (grid mode) and in each model's `meta$spec`.
  ids <- configs$label
  scoreboard <- .build_compare_scoreboard(
    ids        = ids,
    metric     = metric,
    oof_tau    = oof_tau,
    oof_mu0    = oof_mu0,
    oof_mu1    = oof_mu1,
    Y          = Y, T = T, nuis = nuis,
    fold_id    = fold_id,
    failed     = failed,
    config_df  = if (grid_mode) configs else NULL
  )

  oof_tau_matrix <- tryCatch({
    mat <- do.call(cbind, lapply(ids, function(k) {
      v <- oof_tau[[k]]
      if (is.null(v)) rep(NA_real_, n) else as.numeric(v)
    }))
    colnames(mat) <- ids
    mat
  }, error = function(e) NULL)

  metadata <- list(
    partition         = "oof",
    n_folds           = as.integer(n_folds),
    methods           = ids,
    base_learners     = base_learners,
    configs           = configs,
    metric            = metric,
    eval_nuisance     = eval_nuisance,
    timestamp         = Sys.time(),
    random_state      = random_state,
    Y_type            = if (isTRUE(binary_outcome)) "binary" else "continuous",
    n_obs             = n,
    n_covariates      = ncol(X),
    failed_methods    = if (length(failed) > 0L) failed else NULL,
    formula           = formula,
    parsed            = list(covariate_names = cov_names),
    nuisance_oof      = nuis,
    oof_tau_matrix    = oof_tau_matrix
  )

  if (length(failed) > 0L) {
    rlang::warn(
      sprintf(
        "%d/%d config(s) failed: %s. See metadata$failed_methods.",
        length(failed), nrow(configs),
        paste(names(failed), collapse = ", ")
      ),
      class = "cfomics_compare_partial_failure"
    )
  }

  structure(
    list(
      scoreboard = scoreboard,
      models     = models,
      fold_cache = fold_cache,
      metadata   = metadata
    ),
    class = c("cfomics_compare", "list")
  )
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Expand the (method x base_learner) grid into a config table.
#'
#' Returns a data.frame with one row per fit config: `label` (external key),
#' `method` (bare method), `base` (base learner or NA). With
#' `base_learners = NULL` this is one row per method (label == method,
#' base = NA). Otherwise pluggable methods fan out to one row per base learner
#' (label = method_base); non-pluggable methods stay single (base = NA).
#' @keywords internal
#' @noRd
.cfomics_build_configs <- function(methods, base_learners, pluggable_methods) {
  if (is.null(base_learners)) {
    return(data.frame(label = methods, method = methods,
                      base = NA_character_, stringsAsFactors = FALSE))
  }
  rows <- lapply(methods, function(m) {
    if (m %in% pluggable_methods) {
      data.frame(label = paste0(m, "_", base_learners), method = m,
                 base = base_learners, stringsAsFactors = FALSE)
    } else {
      data.frame(label = m, method = m, base = NA_character_,
                 stringsAsFactors = FALSE)
    }
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Build the wide-format scoreboard given per-method OOF tau-hat.
#' @keywords internal
#' @noRd
.build_compare_scoreboard <- function(ids, metric, oof_tau,
                                       Y, T, nuis, fold_id, failed,
                                       oof_mu0 = NULL, oof_mu1 = NULL,
                                       config_df = NULL) {
  base_const <- mean(.compute_dr_pseudo_outcome(Y, T, nuis$mu0_oof,
                                                  nuis$mu1_oof, nuis$pi_oof))
  # R-loss baseline = the R-loss-optimal constant-effect model,
  # theta* = argmin_theta sum((Yres - theta * Tres)^2) = sum(Yres Tres)/sum(Tres^2),
  # so the r_squared_score measures heterogeneity captured beyond the best
  # constant effect (matching econml's RScorer). NA when Tres has no variation.
  Yres_c <- as.numeric(Y) - nuis$m_oof
  Tres_c <- as.numeric(T) - nuis$pi_oof
  sum_TT <- sum(Tres_c^2)
  theta_star <- if (is.finite(sum_TT) && sum_TT > 1e-12) {
    sum(Yres_c * Tres_c) / sum_TT
  } else NA_real_
  # Grid mode exposes `label` (the id) + bare `method` + `base`; the legacy
  # path keeps the single `method` column (= ids) for backward compatibility.
  sb <- if (!is.null(config_df)) {
    data.frame(label = ids, method = config_df$method, base = config_df$base,
               stringsAsFactors = FALSE)
  } else {
    data.frame(method = ids, stringsAsFactors = FALSE)
  }

  for (m_metric in metric) {
    loss_vals     <- numeric(length(ids))
    loss_lo       <- numeric(length(ids))
    loss_hi       <- numeric(length(ids))
    r2_vals       <- numeric(length(ids))

    for (i in seq_along(ids)) {
      lab <- ids[i]                 # config id (bare method, or method_base)
      tau <- oof_tau[[lab]]
      if (is.null(tau) || all(is.na(tau)) || !is.null(failed[[lab]])) {
        loss_vals[i] <- NA_real_; loss_lo[i] <- NA_real_; loss_hi[i] <- NA_real_
        r2_vals[i]   <- NA_real_
        next
      }
      per_obs <- switch(m_metric,
        "r_loss"       = .compute_r_loss_per_obs(Y, T, nuis$m_oof,
                                                  nuis$pi_oof, tau),
        "dr_loss"      = .compute_dr_loss_per_obs(Y, T, nuis$mu0_oof,
                                                   nuis$mu1_oof,
                                                   nuis$pi_oof, tau),
        # Factual loss uses THIS method's held-out potential outcomes
        # (oof_mu0/oof_mu1), so it reflects the method's own outcome fit rather
        # than the shared evaluation nuisance. NA when the method surfaces no mu.
        "factual_loss" = {
          mu0_m <- if (!is.null(oof_mu0)) oof_mu0[[lab]] else NULL
          mu1_m <- if (!is.null(oof_mu1)) oof_mu1[[lab]] else NULL
          if (is.null(mu0_m) || is.null(mu1_m) ||
              all(is.na(mu0_m)) || all(is.na(mu1_m)))
            rep(NA_real_, length(Y))
          else
            .compute_factual_loss_per_obs(Y, T, mu0_m, mu1_m)
        }
      )
      per_fold <- .aggregate_loss_by_fold(per_obs, fold_id)
      ci <- .between_fold_ci(per_fold)
      loss_vals[i] <- ci$value
      loss_lo[i]   <- ci$ci_lower
      loss_hi[i]   <- ci$ci_upper

      base_per_obs <- switch(m_metric,
        "r_loss"  = .compute_r_loss_per_obs(Y, T, nuis$m_oof, nuis$pi_oof,
                                             rep(theta_star, length(Y))),
        "dr_loss" = .compute_dr_loss_per_obs(Y, T, nuis$mu0_oof, nuis$mu1_oof,
                                              nuis$pi_oof,
                                              rep(base_const, length(Y))),
        "factual_loss" = .compute_factual_loss_per_obs(
          Y, T,
          rep(mean(Y[T == 0]), length(Y)),
          rep(mean(Y[T == 1]), length(Y))
        )
      )
      r2_vals[i] <- .compute_r_squared_score(per_obs, base_per_obs)
    }

    sb[[paste0(m_metric)]]              <- loss_vals
    sb[[paste0(m_metric, "_ci_lower")]] <- loss_lo
    sb[[paste0(m_metric, "_ci_upper")]] <- loss_hi
    sb[[sub("_loss$", "_squared_score", m_metric)]] <- r2_vals
  }

  ite_mean <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau))) NA_real_ else mean(tau, na.rm = TRUE)
  }, numeric(1))
  ite_std <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau))) NA_real_ else stats::sd(tau, na.rm = TRUE)
  }, numeric(1))
  ite_q <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau)))
      c(NA_real_, NA_real_, NA_real_)
    else
      unname(stats::quantile(tau, c(0.25, 0.5, 0.75), na.rm = TRUE))
  }, numeric(3))

  sb$ite_mean <- ite_mean
  sb$ite_std  <- ite_std
  sb$ite_q25  <- ite_q[1, ]
  sb$ite_q50  <- ite_q[2, ]
  sb$ite_q75  <- ite_q[3, ]

  primary_metric <- metric[1]
  sb$rank <- rank(sb[[primary_metric]], ties.method = "min", na.last = "keep")

  # Expose a per-metric rank column when multiple metrics are requested. The
  # single `rank` column above ranks by metric[1]; `rank_<metric>` lets
  # cf_select / downstream tools read the rank under any of the requested
  # losses without recomputation.
  if (length(metric) > 1L) {
    for (mm in metric) {
      if (mm %in% names(sb)) {
        sb[[paste0("rank_", mm)]] <- rank(sb[[mm]], ties.method = "min",
                                          na.last = "keep")
      }
    }
  }

  sb
}

#' @export
print.cfomics_compare <- function(x, full = FALSE, ...) {
  cat("<cfomics_compare>: ", nrow(x$scoreboard), " method(s) on ",
      x$metadata$n_folds, "-fold split, n=", x$metadata$n_obs,
      ", p=", x$metadata$n_covariates, "\n", sep = "")
  cat("Metric(s): ", paste(x$metadata$metric, collapse = ", "), "\n", sep = "")
  cat("Eval nuisance: ", x$metadata$eval_nuisance, "\n", sep = "")
  sb <- x$scoreboard
  primary_metric <- x$metadata$metric[1]
  if (isTRUE(full)) {
    print(sb)
  } else {
    id_cols   <- if ("label" %in% colnames(sb)) c("label", "base") else "method"
    keep_cols <- c(id_cols, primary_metric, "ite_mean", "ite_std", "rank")
    print(sb[, intersect(keep_cols, colnames(sb)), drop = FALSE])
  }
  if (!is.null(x$metadata$failed_methods)) {
    cat("Failed methods: ",
        paste(names(x$metadata$failed_methods), collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' Summarise a cf_compare result.
#'
#' Prints the full scoreboard (every metric column) plus the fold/n/p and
#' evaluation-nuisance context, rather than exposing the default list dump.
#'
#' @param object A \code{cfomics_compare} object.
#' @param ... Ignored.
#' @return \code{object}, invisibly.
#' @export
summary.cfomics_compare <- function(object, ...) {
  cat("cf_compare summary\n")
  cat("  ", nrow(object$scoreboard), " method(s), ",
      object$metadata$n_folds, "-fold split, n=", object$metadata$n_obs,
      ", p=", object$metadata$n_covariates, "\n", sep = "")
  cat("  metric(s): ", paste(object$metadata$metric, collapse = ", "),
      " | eval nuisance: ", object$metadata$eval_nuisance, "\n", sep = "")
  cat("  full scoreboard:\n")
  print(object$scoreboard)
  if (!is.null(object$metadata$failed_methods)) {
    cat("  failed methods: ",
        paste(names(object$metadata$failed_methods), collapse = ", "),
        "\n", sep = "")
  }
  invisible(object)
}
