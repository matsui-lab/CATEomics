#' Fit a Causal Inference Model
#'
#' @description
#' Fits a causal inference model using various backends with a unified API.
#' Supports data.frame, SummarizedExperiment, and MultiAssayExperiment objects.
#'
#' @details
#' \strong{Estimand:}
#'
#' All methods estimate the \strong{Average Treatment Effect (ATE)}:
#' \deqn{ATE = E[Y(1) - Y(0)]}
#' where \eqn{Y(1)} and \eqn{Y(0)} are potential outcomes under treatment
#' and control, respectively. A subset of methods also estimates the
#' \strong{Individual Treatment Effect (ITE)} / Conditional Average
#' Treatment Effect:
#' \deqn{ITE_i = E[Y_i(1) - Y_i(0) \mid X_i]}
#'
#' All methods estimate the individual/conditional treatment effect:
#' \code{grf}, \code{bcf}, \code{slearner}, \code{tlearner}, \code{xlearner},
#' \code{drlearner}, \code{ganite}, \code{cevae}, \code{tarnet},
#' \code{dragonnet}. \code{predict()} returns per-sample \eqn{\hat\tau} and
#' potential outcomes for every method.
#'
#' \strong{Counterfactual outcomes \eqn{\hat y_0}, \eqn{\hat y_1}:}
#' \eqn{\hat y_t} is always the \emph{modelled} counterfactual under
#' \eqn{\text{do}(T=t)}, not the observed value for treated / control
#' subsets.
#'
#' \strong{Causal Assumptions:}
#'
#' Valid causal inference requires:
#' \enumerate{
#'   \item \strong{SUTVA}: No interference between units
#'   \item \strong{Ignorability}: \eqn{(Y(0), Y(1)) \perp T | X} (no unmeasured confounding)
#'   \item \strong{Positivity}: \eqn{0 < P(T=1|X) < 1} for all X
#' }
#'
#' Use \code{\link{cf_diagnostics}} on the fitted model to assess these
#' assumptions (overlap, covariate balance, E-value).
#'
#' \strong{Method Selection Guide:}
#' \itemize{
#'   \item \strong{grf}: heterogeneous effects, robust, R-native
#'     (recommended default).
#'   \item \strong{bcf}: Bayesian Additive Regression Trees with
#'     regularized treatment splits.
#'   \item \strong{slearner / tlearner / xlearner}: sklearn-backed
#'     metalearners via econml; see \code{\link{metalearners}}.
#'   \item \strong{drlearner}: Neyman-orthogonal DR-Learner (econml).
#'   \item \strong{ganite}: GAN-based ITE (PyTorch).
#'   \item \strong{cevae}: latent-confounding VAE, requires Pyro.
#'   \item \strong{tarnet / dragonnet}: neural CATE via catenets.
#' }
#'
#' @param formula Formula specifying the model: \code{Y ~ T | X1 + X2 + ...}
#'   \itemize{
#'     \item Y: Outcome variable (numeric)
#'     \item T: Treatment variable (binary: 0/1)
#'     \item X1, X2, ...: Covariates (numeric)
#'   }
#'   For SummarizedExperiment, Y can be a gene name from the assay or a
#'   column from colData.
#' @param data Data source: data.frame, SummarizedExperiment, or
#'   MultiAssayExperiment object.
#' @param Y Outcome column name (character scalar), for the named-args API
#'   (mutually exclusive with \code{formula}).
#' @param T Treatment column name (character scalar, binary 0/1), for the
#'   named-args API (mutually exclusive with \code{formula}).
#' @param X_baseline Character vector of non-omics (clinical / demographic)
#'   covariate column names, for the named-args API.
#' @param X_rnaseq Character vector of feature / column names for the RNA-seq
#'   omics layer (named-args API).
#' @param X_mirna Character vector of feature / column names for the miRNA
#'   omics layer (named-args API).
#' @param X_methyl Character vector of feature / column names for the
#'   methylation omics layer (named-args API).
#' @param X_protein Character vector of feature / column names for the
#'   proteomics omics layer (named-args API).
#' @param X_mutation Character vector of feature / column names for the somatic
#'   mutation omics layer (named-args API).
#' @param X_cnv Character vector of feature / column names for the copy-number
#'   variation omics layer (named-args API).
#' @param X_assays Named list of additional omics assays (each a vector of
#'   feature / column names), for MultiAssayExperiment input.
#' @param method Character, the causal inference method to use:
#'   \describe{
#'     \item{"grf"}{Generalized Random Forest (R-native, recommended default)}
#'     \item{"bcf"}{Bayesian Causal Forests (R-native, requires bcf)}
#'     \item{"slearner" / "tlearner" / "xlearner"}{sklearn-backed metalearners
#'       via econml (Python)}
#'     \item{"drlearner"}{Neyman-orthogonal DR-Learner (econml, Python)}
#'     \item{"cevae"}{Causal Effect VAE (Python/Pyro)}
#'     \item{"ganite"}{GAN-based ITE estimation (Python/PyTorch)}
#'     \item{"tarnet" / "dragonnet"}{representation-learning neural CATE
#'       (catenets, Python)}
#'     \item{"auto"}{Automatically select the best available method}
#'   }
#' @param n_bootstrap Integer, number of bootstrap iterations (default 100).
#' @param random_state Integer random seed for reproducibility (default
#'   \code{0}). When set, yields bit-identical results across runs with
#'   the same inputs. The caller's R \code{.Random.seed} and any
#'   Python-side RNG state touched by the backend (\code{numpy},
#'   \code{torch}, Python \code{random}) are saved on entry and restored
#'   on exit so seeding does not leak into the user's session. See the
#'   per-method documentation for backend-specific scope (e.g.,
#'   \code{ganite} touches TensorFlow's \emph{global} RNG and cannot
#'   save/restore it; see \code{\link{cf_fit_ganite}}).
#' @param assay_name For SummarizedExperiment: name or index of assay to use
#'   (default 1, Bioconductor convention).
#' @param experiments For MultiAssayExperiment: character vector of experiment
#'   names to include (\code{NULL} = all experiments). Only the
#'   RNA-seq assay (identified by experiment name \code{rnaseq} / \code{rna} /
#'   \code{rna_seq} or the \code{X_rnaseq} slot) is internally log-normalised;
#'   every other assay is passed through as supplied (normalise it externally).
#'   Samples are restricted to the sampleMap intersection (complete cases across
#'   the selected experiments).
#' @param assay_cols For data.frame: column specification of assay features
#'   (character vector / integer indices / logical mask). When supplied,
#'   selected columns are treated as omics assay (log-CPM + feature_method
#'   applied). Mutually exclusive with \code{assay_pattern}.
#' @param assay_pattern For data.frame: regular expression matched against
#'   column names to select assay features. Mutually exclusive
#'   with \code{assay_cols}.
#' @param preprocess Logical (default \code{TRUE}). When \code{TRUE} the assay
#'   matrix (SE / MAE primary / df assay_cols/pattern) is preprocessed via
#'   log-CPM normalisation. On data.frame inputs without
#'   \code{assay_cols} / \code{assay_pattern}, \code{preprocess} is silently
#'   disabled (clinical-only path, backward compatible).
#' @param batch Optional character: column name in colData (SE/MAE) or data
#'   (df) containing batch labels for \code{sva::ComBat} correction on the
#'   log-CPM values. \code{NULL} (default) skips batch
#'   correction.
#' @param feature_method One of \code{"variance"} (default), \code{"none"},
#'   \code{"pca"}, \code{"pathway"}.
#'   \code{"variance"}/\code{"pca"} pick top \code{n_features} by row variance
#'   or principal components; \code{"none"} keeps every assay feature.
#'   \code{"pathway"} runs \code{GSVA::gsva} with \code{gene_sets} (default
#'   Hallmark 50 via msigdbr). The default \code{"variance"} means an SE / MAE
#'   omics assay is included as covariates by default (a \code{data.frame}
#'   without an assay is unaffected).
#' @param n_features Integer, number of features / components / pathways to
#'   keep. \code{NULL} (default) resolves to 100 for a SummarizedExperiment /
#'   data.frame and 50 for a MultiAssayExperiment (per-assay budget).
#' @param gene_sets Named list of gene sets, used only when
#'   \code{feature_method = "pathway"}. \code{NULL} (default) loads Hallmark
#'   50 via msigdbr.
#' @param pathway_method One of \code{"gsva"} (default), \code{"ssgsea"},
#'   \code{"plage"}, \code{"zscore"}; passed to \code{GSVA::gsva}.
#' @param partition One of \code{"oof"} (default) or
#'   \code{"none"}. When \code{"oof"}, cf_fit runs K-fold OOF + full-data
#'   fit (\code{n_folds + 1} fits per method) and populates
#'   \code{fit$res$tau_hat_oof} and \code{fit$res$nuisance}. When \code{"none"},
#'   only the single full-data fit is run; \code{fit$res$tau_hat_oof} is
#'   \code{NULL} but \code{fit$res$nuisance} is still populated with full-data
#'   in-sample nuisance.
#' @param n_folds Integer K for the stratified K-fold split when
#'   \code{partition = "oof"}. Default \code{5L}.
#' @param eval_nuisance Character, base-learner name used for the evaluation
#'   nuisance (\code{mu0}, \code{mu1}, \code{pi}) needed by DR-based metrics in
#'   \code{\link{cf_evaluate}} and \code{\link{cf_interpret}}. One of
#'   \code{"rf"} (default), \code{"lasso"}, \code{"ridge"}, \code{"elastic_net"},
#'   \code{"linear"}, \code{"logistic"}, or \code{"xgboost"}.
#' @param eval_propensity_score Optional known propensity score (scalar in
#'   (0, 1) or length-n vector) used for the evaluation nuisance instead of an
#'   estimated one (RCT context). \code{NULL} (default) estimates it. Symmetric
#'   with \code{\link{cf_compare}}.
#' @param keep_kfold_models Logical (default \code{FALSE}). When \code{TRUE}
#'   and \code{partition = "oof"}, the K per-fold fit objects are
#'   retained in \code{fit$fold_cache}.
#' @param bootstrap_ci Logical (default \code{FALSE}); per-method opt-in
#'   bootstrap CI for ATE / ITE. Only honoured
#'   by methods that explicitly accept it.
#' @param device One of \code{"auto"} (default), \code{"cpu"}, \code{"gpu"};
#'   forwarded to the chosen backend. For NN
#'   methods the backend probes \code{torch.cuda.is_available()} and
#'   records the resolved device in \code{meta$provenance$device_effective}.
#' @param trim Numeric length-2 vector \code{c(L, U)} (default
#'   \code{c(0.025, 0.975)}); propensity score clipping bounds applied
#'   uniformly across PS-using methods via the
#'   \code{ClippedPropensityModel}.
#' @param propensity_score Optional scalar or length-\code{n} numeric
#'   vector of known propensity scores (RCT context).
#'   \code{NULL} (default) lets the backend estimate PS internally.
#' @param outcome_type One of \code{"auto"} (default), \code{"continuous"},
#'   or \code{"binary"}. \code{"auto"} detects a 0/1 outcome and
#'   switches every method to a risk-difference / classifier outcome model,
#'   skipping the outcome z-score for NN backends. Force with
#'   \code{"binary"} / \code{"continuous"}.
#' @param propensity_X Optional character vector of column names to use
#'   for PS estimation when \code{propensity_score = NULL}.
#'   PS-using methods (grf/bcf/drlearner/xlearner)
#'   honour the subset; the others emit \code{cfomics_propensity_X_ignored}.
#' @param ._skip_nuisance Internal flag (default \code{FALSE}) used by
#'   \code{cf_compare} / \code{cf_evaluate} to suppress redundant nuisance
#'   refit when shared K-fold nuisance is already available.
#' @param ... Additional method-specific options passed to the backend.
#'
#' @return Object of class \code{"cf_model"} with components:
#'   \describe{
#'     \item{method}{Character string identifying the method used.}
#'     \item{fit}{Method-specific fitted model with \code{res} containing results.}
#'     \item{meta}{Metadata including formula, dimensions, and software versions.}
#'   }
#'   See \code{\link{cfomics_result}} for the full structure specification.
#'
#' @seealso
#'   \code{\link{cfomics_result}} for result object structure,
#'   \code{\link{cf_diagnostics}} for post-fit assumption checks
#' @export
#' @examples
#' # Basic example with a GRF causal forest (R-native).
#' set.seed(123)
#' n <- 100
#' demo_data <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X1 + rnorm(n)
#'
#' # GRF causal forest (R-native, ITE-capable).
#' fit <- cf_fit(data = demo_data, Y = "Y", T = "T",
#'               X_baseline = c("X1", "X2"), method = "grf",
#'               partition = "none")
#' predict(fit, type = "ate")        # point ATE
#' head(predict(fit, type = "tau"))  # per-sample CATE
#' summary(fit)                      # ATE + 95% CI
#'
#' \donttest{
#' # SummarizedExperiment input via the named-args API
#' se <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4)
#' fit_se <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'                  X_baseline = c("age", "sex"), method = "grf",
#'                  n_features = 40)
#' }
cf_fit <- function(formula = NULL, data = NULL,
                   method = c("grf",
                              "bcf",
                              "ganite", "drlearner",
                              "cevae",
                              "slearner", "tlearner", "xlearner",
                              "tarnet", "dragonnet",
                              "auto"),
                   # ---- Facade named arguments ----
                   # Named-args API. Mutually
                   # exclusive with `formula`; one of the two must be supplied.
                   # See `.resolve_cf_api()` for the dispatch logic.
                   Y          = NULL,    # outcome variable name (character scalar)
                   T          = NULL,    # treatment variable name (character scalar)
                   X_baseline = NULL,    # baseline (non-omics): clinical / demographic
                   X_rnaseq   = NULL,    # RNA-seq layer
                   X_mirna    = NULL,    # miRNA layer
                   X_methyl   = NULL,    # methylation layer
                   X_protein  = NULL,    # proteomics layer
                   X_mutation = NULL,    # somatic mutation layer
                   X_cnv      = NULL,    # copy-number variation layer
                   X_assays   = NULL,    # generic named list for additional omics
                   # PS-specific covariate subset. NULL = use full X for PS;
                   # character vector = use that subset (applied per
                   # PS-using method).
                   propensity_X = NULL,
                   n_bootstrap = 100L,
                   bootstrap_ci = FALSE,
                   random_state = 0L,
                   # ---- input path (per data class) ----
                   assay_name    = 1,                    # SE
                   experiments   = NULL,                 # MAE
                   assay_cols    = NULL,                 # df
                   assay_pattern = NULL,                 # df
                   # ---- preprocessing + feature selection ----
                   preprocess     = TRUE,                # TMM -> log-CPM
                   batch          = NULL,                # ComBat batch column (log values)
                   # Default feature selection = "variance" (top-n by
                   # row variance) so SE/MAE omics assays are actually included
                   # as covariates by default. df without an assay is a no-op.
                   feature_method = c("variance", "none", "pca", "pathway"),
                   # NULL resolves to an input-type-aware default —
                   # 100 for SE / data.frame, 50 for MAE (per-assay budget,
                   # multiple assays). An explicit value overrides both.
                   n_features     = NULL,
                   gene_sets      = NULL,                # for pathway; NULL = Hallmark 50
                   pathway_method = c("gsva", "ssgsea", "plage", "zscore"),
                   partition         = c("oof", "none"),
                   n_folds           = 5L,
                   eval_nuisance     = "rf",
                   eval_propensity_score = NULL,  # known eval PS (RCT), symmetric with cf_compare
                   keep_kfold_models = FALSE,
                   # ---- GPU device ----
                   device            = c("auto", "cpu", "gpu"),
                   # ---- PS trim ----
                   # PS clip applied by PS-using methods (grf/bcf/drlearner/xlearner);
                   # PS-free methods ignore it.
                   trim              = c(0.025, 0.975),
                   # ---- known PS specification ----
                   # NULL = estimate (default, observational study) / scalar (RCT 50/50) /
                   # vector length-n (stratified RCT, survey). Mutually exclusive with an
                   # explicit `model_propensity` (cfomics_propensity_args_conflict error).
                   propensity_score  = NULL,
                   # ---- binary outcome support (risk difference / classifier).
                   # "auto" (default) detects a 0/1 outcome and switches every
                   # method to risk-difference / classifier outcome models,
                   # skipping the Y z-score for NN backends. Force with
                   # "binary" / "continuous".
                   outcome_type = c("auto", "continuous", "binary"),
                   # ---- internal (not exposed in roxygen):
                   # skips the inner full-data nuisance computation when a
                   # recursive caller computes K-fold OOF nuisance separately.
                   # Not intended for general users.
                   `._skip_nuisance` = FALSE,
                   ...) {
  # Named-args / formula dual entry. `.resolve_cf_api()` returns the canonical
  # formula + a `structured` bag that we stash on the result `meta` later.
  # When the caller went through the legacy formula API, structured is NULL and
  # we behave as before, with one deprecation warning per session.
  .facade <- .resolve_cf_api(
    formula    = formula,
    Y          = Y,
    T          = T,
    X_baseline = X_baseline,
    X_rnaseq   = X_rnaseq,   X_mirna    = X_mirna,
    X_methyl   = X_methyl,   X_protein  = X_protein,
    X_mutation = X_mutation, X_cnv      = X_cnv,
    X_assays   = X_assays,
    data       = data
  )
  formula <- .facade$formula
  .cfomics_named_args <- .facade$structured  # NULL when legacy formula path
  if (is.null(data)) {
    rlang::abort("`data` argument is required.", class = "cfomics_invalid_data")
  }
  # 8-method propensity_X-ignored warning. We need `method`
  # in canonical form for the lookup but `match.arg(method)` happens further
  # below; use a lightweight pre-resolve here. The warning may fire for the
  # auto / unknown method case too — silently swallow if the lookup fails.
  .method_pre <- tryCatch(match.arg(method), error = function(e) NA_character_)
  if (!is.na(.method_pre)) {
    .warn_propensity_X_ignored(.method_pre, propensity_X)
  }

  # `feature_method` / `pathway_method` are character options; match.arg gives
  # the canonical first-element default when callers omit the argument.
  feature_method <- match.arg(feature_method)
  pathway_method <- match.arg(pathway_method)

  # Resolve the input-type-aware n_features default. NULL (the
  # default) -> 50 for MAE (per-assay budget, multiple assays) / 100 otherwise
  # (SE / data.frame). An explicit numeric value is honoured as-is.
  n_features <- .cf_resolve_n_features(n_features, data)

  partition <- match.arg(partition)
  device <- match.arg(device)  # GPU device validation
  outcome_type <- match.arg(outcome_type)  # binary outcome routing

  # PS trim validation.
  if (!is.numeric(trim) || length(trim) != 2L ||
      trim[1] >= trim[2] || trim[1] <= 0 || trim[2] >= 1) {
    rlang::abort(
      sprintf("`trim` must be length-2 numeric with 0 < trim[1] < trim[2] < 1 (got: %s).",
              paste(deparse(trim), collapse = " ")),
      class = "cfomics_invalid_trim"
    )
  }

  # Known PS validation + mutual-exclusion constraint.
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cfomics_invalid_propensity_score"
      )
    }
    # For SummarizedExperiment / MultiAssayExperiment, samples live in columns
    # (`ncol`) and features in rows; for data.frame the samples are in rows.
    # Use the sample-count interpretation consistently so `propensity_score`
    # vector length matches the number of subjects, not the number of features.
    n_obs <- if (inherits(data, c("SummarizedExperiment", "MultiAssayExperiment"))) {
      ncol(data)
    } else {
      nrow(data)
    }
    if (length(propensity_score) != 1L &&
        length(propensity_score) != n_obs) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match number of samples = %d.",
                length(propensity_score), n_obs),
        class = "cfomics_invalid_propensity_score"
      )
    }
    if ("model_propensity" %in% names(list(...))) {
      rlang::abort(
        c("Both `propensity_score` and `model_propensity` are explicitly specified.",
          i = "For RCT (known PS), use `propensity_score = ...` only.",
          i = "For observational study (estimated PS), use `model_propensity = ...` only."),
        class = "cfomics_propensity_args_conflict"
      )
    }
    # Caveat: aipw / drlearner / xlearner internally K-fold
    # cross-fit via econml; FixedPropensityModel.predict_proba(X)
    # cannot recover the row-level subset for vector PS, so vector PS is
    # rejected with an actionable hint pointing at scalar / bcf.
    if (length(propensity_score) > 1L &&
        isTRUE(.method_pre %in% c("drlearner", "xlearner"))) {
      rlang::abort(
        c(sprintf("Vector `propensity_score` is unsupported for method = '%s'.",
                  .method_pre),
          i = "econml cross-fit pipelines cannot map per-sample PS through internal folds.",
          i = "Use a scalar (`propensity_score = 0.5`) or switch to method = 'bcf' for per-sample PS."),
        class = "cfomics_propensity_vector_unsupported"
      )
    }
  }

  # Device resolution: under "auto", probe GPU availability. GPU-unsupported
  # methods = R-native 2 (grf/bcf) + Python CPU-only 4 (slearner / tlearner /
  # xlearner / drlearner: econml = sklearn CPU); device='gpu' on these warns +
  # continues on CPU. GPU-capable = Python NN 4 (cevae/ganite/tarnet/dragonnet),
  # decided backend-side via `torch.cuda.is_available()`.
  gpu_unsupported_methods <- c("grf", "bcf",
                               "slearner", "tlearner", "xlearner", "drlearner")
  # Use the pre-resolved scalar method (.method_pre); the formal `method` is
  # still the full default vector until match.arg() below, and passing a vector
  # into the `%in%` device checks would error on `if()` (length > 1).
  device_resolved <- .resolve_device(device, .method_pre, gpu_unsupported_methods)
  # Hide CUDA for a CPU fit only for this call's duration, then
  # restore the caller's CUDA_VISIBLE_DEVICES so later GPU fits are not
  # silently disabled.
  if (identical(device_resolved, "cpu")) {
    .old_cuda <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA)
    Sys.setenv(CUDA_VISIBLE_DEVICES = "")
    on.exit({
      if (is.na(.old_cuda)) Sys.unsetenv("CUDA_VISIBLE_DEVICES")
      else Sys.setenv(CUDA_VISIBLE_DEVICES = .old_cuda)
    }, add = TRUE)
  }
  if (identical(device, "gpu") && isTRUE(.method_pre %in% gpu_unsupported_methods)) {
    rlang::warn(
      c(
        sprintf("device='gpu' is unsupported for method '%s' (CPU-only backend).", .method_pre),
        i = "Falling back to CPU. Use Python NN methods (cevae/tarnet/dragonnet/ganite) for GPU acceleration."
      ),
      class = "cfomics_gpu_unsupported_warning"
    )
  }

  valid_eval_nuisance <- c("rf", "lasso", "ridge", "elastic_net",
                           "linear", "logistic", "xgboost")
  if (!is.character(eval_nuisance) || length(eval_nuisance) != 1L ||
      !(eval_nuisance %in% valid_eval_nuisance)) {
    rlang::abort(
      sprintf(
        "`eval_nuisance` must be one of: %s (got: %s).",
        paste(shQuote(valid_eval_nuisance), collapse = ", "),
        paste(deparse(eval_nuisance), collapse = " ")
      ),
      class = "cfomics_invalid_eval_nuisance"
    )
  }

  if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L ||
      n_folds != as.integer(n_folds)) {
    rlang::abort(
      sprintf("`n_folds` must be a single integer >= 2 (got: %s).",
              paste(deparse(n_folds), collapse = " ")),
      class = "cfomics_invalid_n_folds"
    )
  }
  n_folds <- as.integer(n_folds)

  if (!is.logical(keep_kfold_models) || length(keep_kfold_models) != 1L) {
    rlang::abort(
      "`keep_kfold_models` must be a single logical value.",
      class = "cfomics_invalid_keep_kfold_models"
    )
  }

  # Mutual-exclusion check (meaningful on the df path; effectively ignored for
  # SE/MAE, but rejecting the contradictory input early is clearer for the user).
  if (!is.null(assay_cols) && !is.null(assay_pattern)) {
    rlang::abort(
      "Specify either `assay_cols` or `assay_pattern`, not both.",
      class = "cfomics_input_conflict"
    )
  }
  method <- match.arg(method)

  # Handle "auto" method selection
  if (method == "auto") {
    method <- select_best_method()
    message(sprintf("Auto-selected method: %s", method))
  }

  # Validate `...` for unknown / mistyped arguments (symmetric with cf_compare).
  # A typo'd hyper-parameter (e.g. `model_outcom`) would otherwise be silently
  # dropped into a wrapper's ignored `...`. Emits a `cfomics_unknown_arg`
  # warning with the nearest match.
  .cfomics_validate_dot_args(list(...), method)

  # Hybrid CI policy: structurally unsupported `bootstrap_ci`
  # requests error before any fit work is done. `.resolve_ci_method()`
  # raises `cfomics_bootstrap_unsupported_nn` for tarnet/dragonnet/ganite.
  # For grf/bcf/cevae a warning is
  # raised here (native CI is faster); for s/t/x-learner /
  # drlearner (non-linear final) it is silently accepted.
  if (isTRUE(bootstrap_ci)) {
    native_label_map <- list(
      grf       = "asymptotic Wald",
      bcf       = "Bayesian posterior",
      cevae     = "Bayesian posterior",
      drlearner = "asymptotic Wald",  # only when model_final = "linear"
      slearner  = "unavailable",
      tlearner  = "unavailable",
      xlearner  = "unavailable",
      ganite    = "unavailable",
      tarnet    = "unavailable",
      dragonnet = "unavailable"
    )
    native_label <- native_label_map[[method]]
    if (is.null(native_label)) native_label <- "unavailable"
    invisible(.resolve_ci_method(method, bootstrap_ci, native_label))
  }

  # ===================================================================
  # cf_fit 7-step flow
  # ===================================================================

  # partition dispatch sits between Step 5 and Step 6. Step 1-5 (log-CPM + batch +
  # feature_method + factor expansion) always run, and both the K-fold loop and
  # the full-data fit operate on the preprocessed X. This guarantees:
  #   - SE/MAE/df + preprocess=TRUE => Step 2 log-CPM is always applied
  #   - K-fold loop and full-data fit run on the same preprocessed X (BLP-on-X
  #     consistency)
  #   - meta slots (input_type / assay_spec / pca_loadings / preprocess) are
  #     populated in Step 7
  # Step 1-7 run regardless of `partition` (Step 5.5 = the partition='oof'
  # dispatch is the only branch; Step 6-7 are shared by both partitions).
  # `._skip_nuisance` is an explicit signature argument (above) and does not
  # flow into `...`, so method dispatch is unaffected.

  # Steps 1-5: shared preprocessing (input validation + assay resolution ->
  # log-normalisation -> ComBat -> feature select/aggregation -> factor expand
  # -> assemble X = cbind(X_clinical, t(processed_features))). Extracted to
  # `.cf_preprocess_input()` so cf_fit and cf_compare apply the identical
  # Step 1-5 pipeline to the same SE/MAE/df input (single source of truth — no
  # drift). Binary-outcome resolution + method dispatch stay in cf_fit below.
  # `parsed` / `step2` / `step3` / `step3_5`
  # (feature select) / `step4` (pathway) / `batch_vec` / `input_type` /
  # `assay_spec` / `factor_levels` feed the Step 7 meta-populate path.
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
    method         = method
  )
  parsed        <- prep$parsed
  step2         <- prep$step2
  step3         <- prep$step3
  step3_5       <- prep$step3_5
  step4         <- prep$step4
  batch_vec     <- prep$batch_vec
  input_type    <- prep$input_type
  assay_spec    <- prep$assay_spec
  factor_levels <- prep$factor_levels
  X             <- prep$X
  T             <- prep$T
  Y             <- prep$Y
  cov_names     <- prep$cov_names

  # Resolve binary vs continuous outcome ONCE (before dispatch), inform
  # the user on auto-detection, and thread `binary_outcome` to every backend +
  # the evaluation nuisance so binary outcomes use risk-difference / classifier
  # models and NN backends skip the Y z-score.
  binary_outcome <- .resolve_binary_outcome(Y, outcome_type)

  # Wire the known evaluation propensity. A non-NULL `eval_propensity_score`
  # switches the evaluation-nuisance propensity slot to the user-supplied value
  # (RCT context) instead of estimating it; NULL (default) keeps the estimated
  # eval nuisance, so behaviour is unchanged unless a known eval PS is supplied.
  # This mirrors cf_compare, where a numeric `eval_propensity_score` engages the
  # same `propensity = "known"` path. `eval_nuisance` stays the user-facing
  # base-learner name; the list spec is built here only for the nuisance helpers.
  eval_nuisance_spec <- if (!is.null(eval_propensity_score))
    list(outcome = eval_nuisance, propensity = "known")
  else eval_nuisance

  # ===================================================================
  # Step 5.5: partition dispatch
  # ===================================================================
  # partition='oof' runs `.kfold_oof_predict_matrix()` for per-fold method
  # training + a full-data fit on the preprocessed X. The returned
  # full_inner$fit has the same shape as the partition='none' `fit`, so both
  # merge into the shared Step 6-7 path. oof_ite / nuisance / fold_cache are
  # attached below.
  if (identical(partition, "oof")) {
    full_inner <- .kfold_oof_predict_matrix(
      X = X, T_vec = T, Y_vec = Y, method = method,
      cov_names         = cov_names,
      n_folds           = n_folds,
      random_state      = random_state,
      keep_kfold_models = keep_kfold_models,
      eval_nuisance     = eval_nuisance_spec,
      n_bootstrap       = n_bootstrap,
      bootstrap_ci      = bootstrap_ci,
      device            = device_resolved,
      propensity_X      = propensity_X,
      propensity_score  = propensity_score,
      trim              = trim,
      eval_propensity_score = eval_propensity_score,
      binary_outcome    = binary_outcome,
      ...
    )
    fit <- full_inner$fit
    oof_ite_value    <- full_inner$oof_ite
    nuisance_full    <- full_inner$nuisance
    fold_cache_value <- full_inner$fold_cache
  } else {
    # partition = 'none' (Step 5 switch dispatch, in-sample fit)
    oof_ite_value    <- NULL
    fold_cache_value <- NULL
    nuisance_full    <- NULL  # full-data nuisance computed below

  fit <- switch(
    method,
    "ganite" = cf_fit_ganite(X = X, T = T, Y = Y,
                             covariate_names = cov_names,
                             random_state = as.integer(random_state),
                             device = device_resolved,
                             binary_outcome = binary_outcome,
                             ...),
    "drlearner" = cf_fit_drlearner(X = X, T = T, Y = Y,
                                   covariate_names = cov_names,
                                   bootstrap_ci = isTRUE(bootstrap_ci),
                                   n_bootstrap  = as.integer(n_bootstrap),
                                   random_state = as.integer(random_state),
                                   trim = trim,
                                   propensity_score = propensity_score,
                                   propensity_X = propensity_X,
                                   binary_outcome = binary_outcome,
                                   ...),
    "cevae" = cf_fit_cevae(X = X, T = T, Y = Y,
                           covariate_names = cov_names,
                           random_state = as.integer(random_state),
                           device = device_resolved,
                           binary_outcome = binary_outcome,
                           ...),
    "grf" = cf_fit_grf(X = X, T = T, Y = Y,
                       covariate_names = cov_names,
                       random_state = as.integer(random_state),
                       propensity_score = propensity_score,
                       trim = trim,  # symmetric with bcf/drlearner/xlearner
                       propensity_X = propensity_X,
                       ...),
    "bcf" = cf_fit_bcf(X = X, T = T, Y = Y,
                       covariate_names = cov_names,
                       random_state = as.integer(random_state),
                       trim = trim,
                       # known PS: forward to bcf's native `pihat`.
                       pihat = if (is.null(propensity_score)) NULL else
                                  rep(as.numeric(propensity_score), length.out = length(T)),
                       propensity_X = propensity_X,
                       ...),
    "slearner" = cf_fit_slearner(X = X, T = T, Y = Y,
                                 covariate_names = cov_names,
                                 bootstrap_ci = isTRUE(bootstrap_ci),
                                 n_bootstrap  = as.integer(n_bootstrap),
                                 random_state = as.integer(random_state),
                                 binary_outcome = binary_outcome,
                                 ...),
    "tlearner" = cf_fit_tlearner(X = X, T = T, Y = Y,
                                 covariate_names = cov_names,
                                 bootstrap_ci = isTRUE(bootstrap_ci),
                                 n_bootstrap  = as.integer(n_bootstrap),
                                 random_state = as.integer(random_state),
                                 binary_outcome = binary_outcome,
                                 ...),
    "xlearner" = cf_fit_xlearner(X = X, T = T, Y = Y,
                                 covariate_names = cov_names,
                                 bootstrap_ci = isTRUE(bootstrap_ci),
                                 n_bootstrap  = as.integer(n_bootstrap),
                                 random_state = as.integer(random_state),
                                 trim = trim,
                                 propensity_score = propensity_score,
                                 propensity_X = propensity_X,
                                 binary_outcome = binary_outcome,
                                 ...),
    "tarnet"   = cf_fit_tarnet(X = X, T = T, Y = Y,
                               covariate_names = cov_names,
                               random_state = as.integer(random_state),
                               device = device_resolved,
                               binary_outcome = binary_outcome,
                               ...),
    "dragonnet" = cf_fit_dragonnet(X = X, T = T, Y = Y,
                                   covariate_names = cov_names,
                                   random_state = as.integer(random_state),
                                   device = device_resolved,
                                   binary_outcome = binary_outcome,
                                   ...),
    stop(sprintf("Unknown method: %s", method))
  )

    # partition='none' path: fit `eval_nuisance` once on the full data and
    # populate the cf_model top-level `nuisance = list(mu0_oof, mu1_oof,
    # pi_oof)`. `._skip_nuisance = TRUE` skips it.
    if (!isTRUE(`._skip_nuisance`)) {
      nuisance_full <- .compute_full_data_nuisance(
        X, T, Y,
        eval_nuisance         = eval_nuisance_spec,
        random_state          = random_state,
        eval_propensity_score = eval_propensity_score,
        binary_outcome        = binary_outcome
      )
    }
  } # end of partition='none' branch

  # Step 7: Build meta + result. x_zscore / y_zscore are populated per method
  # (cf_fit reads them from fit$res here).
  preprocess_record <- list(
    enabled            = step2$enabled,
    preprocess_method  = step2$method,
    batch_corrected    = step3$batch_corrected,
    batch_method       = step3$batch_method,
    batch_column       = step3$batch_column,
    feature_method     = feature_method,
    n_features         = n_features,
    selected_features  = if (feature_method == "pathway") step4$selected_names else step3_5$selected_names,
    pathway_method     = step4$pathway_method,
    gene_sets_source   = step4$gene_sets_source,
    # Persist the train-fold *log-CPM* assay (pre-ComBat) + batch so newdata
    # ComBat can re-apply via a train+new joint refit. log-CPM is per-sample
    # independent so no train statistics are needed for the normalisation
    # itself; only ComBat needs the joint matrix. Stored only when ComBat was
    # used, keeping the common path compact.
    train_assay        = if (isTRUE(step3$batch_corrected)) step2$assay else NULL,
    train_batch        = if (isTRUE(step3$batch_corrected)) batch_vec else NULL,
    # Persist gene_sets for GSVA newdata re-application. GSVA::gsva is
    # per-sample independent (gene-rank based) so we only need the gene_sets
    # dictionary; the user may have passed a custom list or accepted the
    # default Hallmark 50.
    gene_sets          = if (feature_method == "pathway") gene_sets else NULL
  )

  # Create standardized meta object.
  meta <- create_cf_meta(
    formula = formula,
    n = length(Y),
    p = ncol(X),
    outcome_name = parsed$outcome_name,
    treatment_name = parsed$treatment_name,
    covariate_names = cov_names,
    method = method,
    estimand = "CATE",  # All methods are CATE-capable (produce per-unit tau)
    target_population = "all",
    bootstrap = bootstrap_ci,  # meta records whether a bootstrap CI was requested
    random_state = random_state,
    preprocess     = preprocess_record,
    pca_loadings   = step3_5$pca_loadings,
    # Persist PCA train-fold center / scale / gene order so
    # `.apply_newdata_preprocessing()` can re-apply the same projection to
    # newdata without recomputing on test samples (data-leakage guard).
    pca_center     = step3_5$pca_center,
    pca_scale      = step3_5$pca_scale,
    pca_gene_names = step3_5$pca_gene_names,
    input_type     = input_type,
    assay_spec     = assay_spec,
    factor_levels  = factor_levels,
    # X z-score is applied centrally in `.cf_preprocess_input` (binary/one-hot
    # skipped), so the train stats come from `prep$x_zscore` for all methods and
    # cf_predict.R `.apply_newdata_preprocessing` re-applies them to newdata.
    # `y_zscore` stays NN-only: each NN backend standardises Y on the data it
    # receives and inverse-transforms its predictions, returning `res$y_zscore`.
    x_zscore       = prep$x_zscore,
    y_zscore       = if (is.list(fit) && !is.null(fit$res$y_zscore)) fit$res$y_zscore else NULL,
    # Raw clinical column names from the formula RHS (pre factor expansion),
    # used by `.apply_newdata_preprocessing()` to model.matrix one-hot expand
    # the clinical columns of newdata.
    covariate_names_clinical_raw = parsed$covariate_names,
    # Named-args facade structured info + propensity_X. `.cfomics_named_args`
    # is the `structured` element returned by `.resolve_cf_api()` — NULL when
    # the caller used the legacy formula API.
    named_args   = .cfomics_named_args,
    propensity_X = propensity_X
  )

  # ===================================================================
  # Step 7: assemble the four-pillar cf_model.
  # ===================================================================
  # Flatten the double-fit nesting: each cf_fit_<method>() returns
  # list(model=, [nuisance siblings], res=list(...)); lift `res` to the
  # top-level `res` pillar and keep the remaining backend container (model +
  # siblings such as grf's forest.Y / forest.W) as the `model` pillar.
  res_pillar   <- fit$res
  model_pillar <- fit
  model_pillar$res <- NULL

  # Surface `ci_method` at the res top level (and record it once under
  # meta$estimator) so downstream readers use a single canonical location.
  # The residual res$metadata (small method-specific predict-config such as
  # drlearner model_final / bcf save_tree_directory) is retained for the
  # per-method predictors.
  if (!is.null(res_pillar$metadata)) {
    if (is.null(res_pillar$ci_method)) {
      res_pillar$ci_method <- res_pillar$metadata$ci_method
    }
    meta$estimator$ci_method <- res_pillar$metadata$ci_method
    # ci_method is now canonical at res$ci_method (+ meta$estimator$ci_method);
    # drop the duplicate copy from res$metadata. The remaining res$metadata
    # fields (bcf save_tree_directory / cevae hyperparameters / drlearner
    # model_final) are read by the per-method predictors, so keep them; when
    # nothing else is left, drop the empty container so grf / bcf models carry
    # no res$metadata at all.
    res_pillar$metadata$ci_method <- NULL
    if (length(res_pillar$metadata) == 0L) res_pillar$metadata <- NULL
  } else {
    meta$estimator$ci_method <- res_pillar$ci_method
  }

  # OOF tau + shared nuisance move into the res pillar.
  res_pillar$tau_hat_oof <- oof_ite_value
  if (!is.null(nuisance_full)) res_pillar$nuisance <- nuisance_full

  # Drop the redundant per-sample `samples` block (mu0/mu1/tau were exact
  # duplicates of res$mu0_hat / mu1_hat / tau_hat; posterior / bootstrap raw
  # draws are already summarised into ate_ci / tau_ci). x_zscore / y_zscore
  # are persisted once in meta (below), so drop the res copies too.
  res_pillar$samples  <- NULL
  res_pillar$x_zscore <- NULL
  res_pillar$y_zscore <- NULL

  # Recompute the canonical IQR summary centrally so every method + partition
  # path produces an identical summary shape.
  res_pillar$summary <- .cf_build_summary(
    ate     = res_pillar$ate,
    ate_ci  = res_pillar$ate_ci,
    tau_hat = res_pillar$tau_hat
  )

  meta$partition <- partition
  meta$estimator$outcome_type <- if (isTRUE(binary_outcome)) "binary" else "continuous"

  result <- structure(
    list(
      method = method,
      res    = res_pillar,   # pillar 1: numeric results
      model  = model_pillar, # pillar 2: fitted backend container
      meta   = meta,         # pillar 3: metadata
      # pillar 4: model-visible training data (cf_interpret surrogate SHAP /
      # cf_predict newdata restore / cf_compare-evaluate-rule nuisance
      # recomputation).
      data   = list(
        X_preprocessed = X,
        Y              = Y,
        T              = T
      )
    ),
    class = "cf_model"
  )
  if (!is.null(fold_cache_value)) {
    result$fold_cache <- fold_cache_value
  }

  # Device fidelity record: requested vs effective device (backends that
  # report `res$device_effective` override the resolved default).
  device_eff <- res_pillar$device_effective %||% device_resolved
  result$meta$provenance$device <- device_eff
  result$meta$provenance$device_requested <- device_resolved
  result$meta$provenance$device_effective <- device_eff

  # Memory warning.
  obj_bytes <- as.numeric(utils::object.size(result))
  if (obj_bytes > 100 * 1024^2) {
    rlang::warn(
      c(
        sprintf("cf_model object is large (%.1f MB).", obj_bytes / 1024^2),
        i = "The preprocessed X is retained in result$data for downstream cf_evaluate() / cf_interpret() / cf_rule()."
      ),
      class = "cfomics_large_object_size"
    )
  }

  # Diagnostics are on-demand, not auto-attached. cf_fit does not compute
  # E-value / balance / overlap into the result meta; run cf_diagnostics(model)
  # after fitting for overlap, covariate balance (SMD Love), requested-variable
  # overlap, and the E-value.

  # Validate result structure (contract enforcement). Downstream functions
  # (cf_compare/select/evaluate/interpret) rely on fields like tau_hat_oof /
  # ate_ci / tau_ci, so contract violations abort here rather
  # than surfacing only at NPE time in downstream code.
  validate_cfomics_result(result)

  result
}


#' Resolve batch vector from data based on input_type
#'
#' Internal helper for cf_fit Step 3. SE = colData[[batch]] / MAE =
#' colData[complete_primaries, ][[batch]] / df = data[[batch]].
#'
#' @keywords internal
#' @noRd
.resolve_batch_vector <- function(data, batch, input_type, step1) {
  vec <- switch(
    input_type,
    "SE" = {
      cd <- as.data.frame(SummarizedExperiment::colData(data))
      cd[[batch]]
    },
    "MAE" = {
      cd <- as.data.frame(MultiAssayExperiment::colData(data))
      primaries <- step1$mae_complete_primaries
      if (is.null(primaries)) {
        cd[[batch]]
      } else {
        cd[primaries, , drop = FALSE][[batch]]
      }
    },
    "data.frame" = data[[batch]]
  )
  if (is.null(vec)) {
    rlang::abort(
      sprintf("Batch column '%s' not found in colData / data.", batch),
      class = "cfomics_batch_column_missing"
    )
  }
  vec
}

#' Resolve whether the outcome is binary
#'
#' Detects a 0/1 outcome under \code{outcome_type = "auto"} (exactly two
#' distinct non-NA values, both in \code{{0, 1}}), or honours the explicit
#' \code{"binary"} / \code{"continuous"} override. Emits a one-per-session
#' \code{cfomics_binary_outcome_detected} message on auto-detection.
#'
#' @param Y Numeric outcome vector.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"}, \code{"binary"}.
#' @return Logical scalar: \code{TRUE} when the outcome should be treated as binary.
#' @keywords internal
#' @noRd
.resolve_binary_outcome <- function(Y, outcome_type = "auto") {
  y <- Y[!is.na(Y)]
  detected <- length(unique(y)) == 2L && all(y %in% c(0, 1))
  bin <- switch(outcome_type,
    "auto"       = detected,
    "binary"     = TRUE,
    "continuous" = FALSE,
    detected
  )
  if (identical(outcome_type, "auto") && isTRUE(bin)) {
    rlang::inform(
      c("Binary outcome detected (Y in {0, 1}).",
        i = "Using risk-difference / classifier outcome models; NN backends skip the Y z-score.",
        i = "Pass `outcome_type = 'continuous'` to force regression."),
      class = "cfomics_binary_outcome_detected",
      .frequency = "once", .frequency_id = "cfomics_binary_outcome_detected"
    )
  }
  isTRUE(bin)
}

#' Select Best Available Method
#'
#' Automatically selects the best available causal inference method
#' based on installed packages and Python availability.
#'
#' @param data Data object (used for future data-dependent selection)
#' @return Character string with method name
#' @noRd
select_best_method <- function() {
  # Prefer R-native CATE methods for reliability.
  if (requireNamespace("grf", quietly = TRUE)) {
    return("grf")
  }
  if (requireNamespace("bcf", quietly = TRUE)) {
    return("bcf")
  }
  # Fall back to a Python CATE backend (econml) when available.
  py_available <- tryCatch(
    reticulate::py_available(initialize = FALSE), error = function(e) FALSE)
  if (py_available &&
      tryCatch(reticulate::py_module_available("econml"),
               error = function(e) FALSE)) {
    return("drlearner")
  }
  "grf"
}
