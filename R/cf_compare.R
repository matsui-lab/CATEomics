#' Compare CATE methods on a common K-fold split
#'
#' @description
#' Cross-fit comparison of CATE-capable methods on a single stratified
#' K-fold split. Computes per-method out-of-fold tau-hat, applies
#' selection metrics (R-loss / DR-loss / factual loss) computed against
#' a shared "evaluation nuisance" (\code{eval_nuisance}, default
#' \code{"super_learner"}), and returns a wide-format scoreboard plus the full-data
#' fitted models for downstream \code{\link{cf_select}} /
#' \code{\link{cf_evaluate}} / \code{\link{cf_interpret}} /
#' \code{\link{cf_rule}}.
#'
#' @details
#' \strong{Two-level nuisance:} cf_compare's evaluation nuisance
#' (\code{eval_nuisance}, default \code{"super_learner"}) is fit once with K-fold
#' OOF and is independent of each method's internal nuisance. The
#' evaluation nuisance feeds the selection loss, the DR pseudo-outcome that
#' RATE and GATES read, and the BLP regression of that pseudo-outcome on the
#' centred tau-hat; per-method internal nuisance estimates (e.g. DR-Learner
#' cross-fit) are untouched.
#'
#' @inheritParams cf_fit
#' @param BPPARAM A \code{BiocParallelParam} instance controlling how the
#'   per-method fits are dispatched (default
#'   \code{BiocParallel::SerialParam()}; pass e.g.
#'   \code{BiocParallel::MulticoreParam()} to fit methods concurrently).
#' @param data data.frame / SummarizedExperiment / MultiAssayExperiment.
#' @param methods Character vector of method names or \code{"all"}.
#'   \code{"all"} expands to the 10 CATE-capable methods (\code{grf,
#'   bcf, slearner, tlearner, xlearner, drlearner, cevae, ganite,
#'   tarnet, dragonnet}); cf_compare is CATE-only. Required unless
#'   \code{configs} names the candidate set instead.
#' @param metric Character vector of metric names, any of
#'   \code{"r_loss"}, \code{"dr_loss"}, \code{"factual_loss"}. The default
#'   is all three; the \code{rank} column and \code{cf_select}'s default
#'   follow \code{metric[1]}. Each loss column is the mean of the per-fold
#'   losses (the same per-fold values its interval and \code{cf_select}'s tie
#'   test are built from), while the companion \code{_squared_score} column is
#'   \code{1 - loss / baseline}, both averaged over the pooled observations.
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
#' @param configs Optional named list naming the candidate set outright, for
#'   sets a \code{methods} x \code{base_learners} grid cannot express: a
#'   different base learner per method, or two candidates that differ only in
#'   their hyperparameters. Each element is
#'   \code{list(method = , model_outcome = , model_final = , params = )} with
#'   only \code{method} required, and its name is the config's label -- the key
#'   the scoreboard, \code{\link{cf_select}}, \code{\link{cf_ensemble}} and the
#'   plots use. Hyperparameters live inside the config, so there is no
#'   three-stage resolution: a config's \code{params} is the block for its own
#'   method (\code{list(num.trees = 600)} for grf,
#'   \code{list(outcome = list(max_depth = 2))} for a metalearner), not a
#'   method-keyed list. Mutually exclusive with \code{methods},
#'   \code{base_learners}, \code{model_outcome}, \code{model_final} and
#'   \code{params}; the propensity model and the evaluation nuisance stay
#'   single across configs and remain arguments here. Default \code{NULL} uses
#'   the \code{methods} form.
#' @param n_folds Integer, K for stratified K-fold. Default \code{5L}. It is
#'   also the split the Super Learner weights are solved on, so it reaches the
#'   evaluation nuisance as well as the out-of-fold estimates.
#' @param eval_nuisance Learner for the shared evaluation nuisance.
#'   \code{"super_learner"} (default) combines the members of
#'   \code{eval_library} with non-negative weights solved on the out-of-fold
#'   predictions; a single learner name is also accepted (\code{"rf"},
#'   \code{"lasso"}, \code{"ridge"}, \code{"elastic_net"}, \code{"linear"},
#'   \code{"logistic"}, \code{"xgboost"}). A \code{list(outcome = ,
#'   propensity = )} gives the outcome and propensity slots different
#'   learners, and naming only one slot applies it to both. The weights are
#'   recorded in \code{metadata$eval_sl_weights}.
#' @param eval_library Character vector of Super Learner member names, used
#'   when \code{eval_nuisance = "super_learner"}. \code{NULL} (default) takes
#'   \code{c("mean", "ridge", "elastic_net", "lasso", "rf", "xgboost_d2",
#'   "xgboost_d4")}; \code{"xgboost_d2"} and \code{"xgboost_d4"} are boosted
#'   fits at \code{max_depth} 2 and 4. A \code{list(outcome = ,
#'   propensity = )} gives the two slots different libraries. Naming it
#'   while \code{eval_nuisance} names a single learner has no effect and
#'   warns (\code{cateomics_eval_library_unused}).
#' @param eval_propensity_score Scalar in (0,1) or length-n vector; used
#'   in a randomised-trial context. \code{NULL} (default) estimates the
#'   propensity internally. This value is the only way to say that the
#'   propensity is known, so giving it together with
#'   \code{eval_nuisance$propensity} is an error.
#' @param propensity_score Scalar in (0,1) or length-n vector; forwarded
#'   to each per-method \code{cf_fit_<method>} call as the user-supplied
#'   known propensity score. For grf this
#'   bypasses the internal \code{regression_forest} and pins \code{W.hat};
#'   for xlearner / drlearner it pins the IPW weights.
#' @param clip Length-2 numeric \code{c(lower, upper)} propensity clipping
#'   bounds forwarded to each per-method fit; \code{NULL} (the default)
#'   resolves to \code{c(0.025, 0.975)}. It reaches the four methods that
#'   estimate a propensity; naming it when none of the compared methods does
#'   is reported. Scores outside the bounds are replaced by the bound, which
#'   changes the weights and does not drop samples.
#' @param random_state Integer random seed (default \code{0}) or NULL. It
#'   determines the stratified split, and every loss, interval, squared score
#'   and rank in the scoreboard is computed from out-of-fold estimates, so all
#'   of them depend on which split was drawn. Results are bit-identical
#'   across runs with the same inputs and the same thread count; three
#'   methods also depend on the thread count, which the seed does not reach
#'   (\code{"bcf"}, whose thread count comes from the machine, and the two
#'   catenets methods \code{"tarnet"} and \code{"dragonnet"}). Pass
#'   \code{params = list(bcf = list(n_threads = 1L))} for a bit-identical bcf
#'   fit; the thread count the catenets methods see is set outside R, for
#'   instance through \code{OMP_NUM_THREADS}.
#' @param keep_kfold_models Logical, retain the K per-fold fit objects
#'   (default \code{FALSE}).
#' @param verbose Logical (default \code{TRUE}), emit progress messages.
#' @param device One of \code{"auto"} (default), \code{"cpu"}, \code{"gpu"};
#'   forwarded to each per-method fit. \code{"auto"} and \code{"gpu"} fall
#'   back to the CPU when no GPU is visible.
#' @param model_outcome,model_final,model_propensity Base-learner tokens applied
#'   to every config that takes a pluggable learner. \code{model_outcome} and
#'   \code{model_final} are mutually exclusive with \code{base_learners}, which varies
#'   them; \code{model_propensity} stays single across all configs and can be
#'   combined with it. \code{"lasso"}, \code{"ridge"} and \code{"elastic_net"} select their
#'   penalty by cross-validation over 5 folds of the rows they are given, with
#'   the exception of \code{"ridge"} on a continuous outcome in the Python methods,
#'   which uses a closed-form leave-one-out and forms no folds. Every config is
#'   checked against that before any fitting starts.
#' @param bootstrap_ci,n_bootstrap Forwarded to each per-method fit, for the
#'   methods that accept a bootstrap interval (default 100 replicates):
#'   \code{slearner}, \code{tlearner}, \code{xlearner}, \code{drlearner}.
#'   When only some of the compared methods can resample, a
#'   \code{cateomics_bootstrap_partial} warning names the ones that cannot and
#'   the comparison goes ahead without intervals for them; when none of them
#'   can, the comparison stops rather than running without what was asked
#'   for.
#' @param params Hyperparameters, resolved per config in three stages:
#'   \code{params[[label]]} (e.g. \code{drlearner_rf}), then \code{params[[method]]}, then
#'   the remaining flat entries. Stages are not merged. A method that fits
#'   more than one learner takes a block per slot (\code{outcome}, and
#'   \code{final} for the X- and DR-learner); see \code{\link{cf_fit}}.
#' @param propensity_params Hyperparameters for the propensity model
#'   (bcf / xlearner / drlearner), shared by every config.
#' @param eval_params Hyperparameters for the shared evaluation nuisance.
#'   With the default Super Learner the block is keyed by library member
#'   (\code{list(xgboost_d4 = list(nrounds = 50))}), optionally inside
#'   \code{list(outcome = , propensity = )}; with a single learner the block holds
#'   its hyperparameters directly. Those two slots are the evaluation
#'   nuisance's outcome and propensity models, not a method's own learners.
#' @return Object of class \code{c("cateomics_compare", "list")} with
#'   elements \code{scoreboard, models, fold_cache, metadata}. The scoreboard
#'   has one row per configuration, keyed by \code{label}, with \code{method}
#'   and \code{base} (the base learner that actually ran, \code{NA} for
#'   methods with no pluggable outcome slot), the requested loss columns and
#'   their intervals, the squared-score column of each loss, \code{rank}
#'   (which ranks by \code{metric[1]}), one \code{rank_<metric>} per metric
#'   when more than one was asked for, and the out-of-fold effect distribution as
#'   \code{tau_mean}, \code{tau_sd}, \code{tau_q25}, \code{tau_q50} and
#'   \code{tau_q75}.
#'
#'   Three further columns split that spread by the fold assignment.
#'   \code{between_fold_share} is the share of the variance of the out-of-fold
#'   \eqn{\hat\tau} that the fold assignment explains,
#'   \eqn{\mathrm{var}(\mathrm{ave}(\hat\tau, fold)) / \mathrm{var}(\hat\tau)};
#'   \code{between_fold_sd} and \code{within_fold_sd} are
#'   \code{tau_sd} split by that share, so
#'   \code{within_fold_sd^2 + between_fold_sd^2} equals \code{tau_sd^2}.
#'   Both variances are taken over the rows that have an effect and a fold, the
#'   same rows \code{tau_sd} is taken over, so the share is a proportion.
#'   \strong{All three columns are \code{NA}} when the out-of-fold
#'   \eqn{\hat\tau} is constant or has fewer than two usable rows: there is
#'   no spread to divide, so the share is undefined rather than 0 or 1, and the
#'   two standard deviations follow it. \code{tau_sd} is \code{0} beside them.
#'   \code{print()} shows \code{between_fold_share} next to \code{tau_sd}; all
#'   three are in the returned table and under \code{print(x, full = TRUE)}.
#'
#'   \code{metadata$configs_resolved} is the same
#'   decomposition as a table, and \code{metadata$params} holds the
#'   hyperparameters each label resolved to.
#'
#'   A candidate that failed leaves its error in
#'   \code{metadata$failed_methods}, one that warned leaves the warnings in
#'   \code{metadata$config_warnings}, and one that printed something without
#'   raising a caveat leaves it in \code{metadata$config_messages} (a
#'   backend's own progress output arrives this way). Each is a list with one
#'   named character vector per candidate, the names being the condition
#'   classes and the values the messages; a classless condition has an empty
#'   name. A candidate is fitted once per fold, so the repeats are collapsed
#'   and only a different wording is kept as a separate entry. All three are
#'   \code{NULL} when there is nothing to report.
#' @seealso \code{\link{cf_select}}, \code{\link{cf_evaluate}},
#'   \code{\link{cf_interpret}}, \code{\link{cf_rule}}, \code{\link{cf_fit}}.
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf). The randomised
#' # propensity is given to the estimators and to the evaluation nuisance,
#' # which is fitted independently of the compared methods so the loss favours
#' # none of them. Runs without Python.
#' se  <- cf_example_data(n = 120, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   # 1 chain keeps bcf fast
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)),
#'                   propensity_score = 0.5, eval_propensity_score = 0.5,
#'                   random_state = 1)
#' cmp
#' @references
#' Nie & Wager (2021) Quasi-Oracle Estimation of Heterogeneous Treatment
#' Effects, \emph{Biometrika} 108(2):299-319 (R-loss).
#'
#' Kennedy (2023) Towards optimal doubly robust estimation of heterogeneous
#' causal effects, \emph{Electron. J. Stat.} 17(2):3008-3049 (DR-loss).
#'
#' Mahajan, Mitliagkas, Neal & Syrgkanis (2024) Empirical analysis of model
#' selection for heterogeneous causal effect estimation, \emph{ICLR}
#' (proxy losses for model selection).
#'
#' van der Laan, Polley & Hubbard (2007) Super Learner,
#' \emph{Stat. Appl. Genet. Mol. Biol.} 6(1):Article 25 (evaluation nuisance).
#'
#' Robinson & Oshlack (2010) A scaling normalization method for differential
#' expression analysis of RNA-seq data, \emph{Genome Biol.} 11:R25 (TMM).
#'
#' Chen, Chen, Lun, Baldoni & Smyth (2025) edgeR v4, \emph{Nucleic Acids Res.}
#' 53(2):gkaf018 (implementation).
#'
#' Johnson, Li & Rabinovic (2007) Adjusting batch effects in microarray
#' expression data using empirical Bayes methods, \emph{Biostatistics}
#' 8(1):118-127 (ComBat).
#'
#' Hänzelmann, Castelo & Guinney (2013) GSVA: gene set variation analysis for
#' microarray and RNA-Seq data, \emph{BMC Bioinformatics} 14:7.
#'
#' Barbie et al. (2009) Systematic RNA interference reveals that oncogenic
#' KRAS-driven cancers require TBK1, \emph{Nature} 462(7269):108-112 (ssGSEA).
#'
#' @export
cf_compare <- function(data = NULL, methods = NULL,
                       # ---- Facade named-args (mirrors cf_fit) ----
                       Y          = NULL, T          = NULL,
                       X_baseline = NULL,
                       omics      = NULL, omics_type = NULL,
                       propensity_X = NULL,
                       # ---- preprocessing controls (mirror cf_fit) ----
                       # cf_compare runs the SAME Step 1-5 preprocessing as
                       # cf_fit via `.cf_preprocess_input()`, so an SE / MAE omics
                       # assay is log-normalised and included as covariates. df
                       # is neutral.
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
                       configs        = NULL,  # optional: name the candidate
                                               # set outright, for combinations
                                               # a grid cannot express
                       outcome_type   = c("auto", "continuous", "binary"),  # binary
                                                       # Y=RD support, symmetric with cf_fit
                       n_folds        = 5L,
                       eval_nuisance  = "super_learner",
                       eval_library   = NULL,         # Super Learner members
                       eval_propensity_score = NULL,  # known PS for
                                                       # known evaluation propensity
                       propensity_score = NULL,       # known PS for the per-method fits
                                                       # (e.g., grf W.hat, xlearner
                                                       # FixedPropensityModel).
                       clip           = NULL,  # resolved per fitter; PS clip forwarded to
                                                       # per-method fits
                       random_state   = 0L,
                       bootstrap_ci   = FALSE,
                       n_bootstrap    = 100L,
                       keep_kfold_models = FALSE,
                       verbose        = TRUE,
                       # Per-method fits are dispatched through BiocParallel. The
                       # default SerialParam() is deterministic and safe
                       # alongside the basilisk per-method Python processes; pass
                       # e.g. BiocParallel::MulticoreParam(k, RNGseed = seed) to
                       # fit methods concurrently.
                       BPPARAM        = BiocParallel::SerialParam(),
                       device         = c("auto", "cpu", "gpu"),
                       model_outcome    = NULL,
                       model_final      = NULL,
                       model_propensity = NULL,
                       params            = NULL,
                       propensity_params = NULL,
                       eval_params       = NULL) {
  # Every penalised fit this call makes reports its thin split through one
  # line, named for the learners it saw.
  .cf_cv_thin_open()
  on.exit(.cf_cv_thin_close(), add = TRUE)
  device <- match.arg(device)
  # Resolve "auto" to a concrete backend here, as cf_fit does: the per-method
  # workers compare the string literally, so an unresolved "auto" would keep
  # every Python method on the CPU.
  device_resolved <- .resolve_device(device, method = "python",
                                     r_native_methods = character(0))
  # The default `metric` is the full c("r_loss", "dr_loss", "factual_loss")
  # triplet so cf_select / cf_ensemble can request any of the three losses
  # without an explicit opt-in. rank uses metric[1].
  # Turn the named arguments into the canonical two-part formula the rest of
  # the pipeline consumes. cf_compare does not call `cf_fit`: every config goes
  # through `.kfold_oof_predict_matrix()` -> `.cf_dispatch_fit()`, so this
  # resolution happens once, here.
  .facade <- .resolve_cf_api(
    Y          = Y, T = T,
    X_baseline = X_baseline,
    omics      = omics,
    omics_type = omics_type,
    data       = data
  )
  formula <- .facade$formula
  # The names drive SE / MAE assay-column materialisation inside
  # `.cf_preprocess_input()`.
  .cateomics_named_args <- .facade$structured
  feature_method <- match.arg(feature_method)
  pathway_method <- match.arg(pathway_method)
  if (is.null(data)) {
    rlang::abort("`data` argument is required.", class = "cateomics_invalid_data")
  }
  # Input-type-aware n_features default (MAE = 50, else 100).
  n_features <- .cf_resolve_n_features(n_features, data)
  # Mirror cf_fit's guard: a known propensity (RCT) and an estimated-PS base
  # learner are mutually exclusive (otherwise the estimated PS silently wins
  # and the known PS is ignored per method).
  if (!is.null(propensity_score) && !is.null(model_propensity)) {
    rlang::abort(
      c("Both `propensity_score` and `model_propensity` are explicitly specified.",
        i = "For RCT (known PS), use `propensity_score = ...` only.",
        i = "For observational study (estimated PS), use `model_propensity = ...` only."),
      class = "cateomics_propensity_args_conflict"
    )
  }
  # `base_learners` (grid over base learners) is mutually exclusive with the
  # single-config `model_outcome` / `model_final` knobs: one selects the
  # comparison space, the other a single fit. Aborting (rather than silently
  # letting one win) mirrors the propensity_score / model_propensity guard.
  if (!is.null(base_learners) &&
      (!is.null(model_outcome) || !is.null(model_final))) {
    rlang::abort(
      c("Both `base_learners` and `model_outcome` / `model_final` are specified.",
        i = "Use `base_learners = c(...)` to grid over base learners, OR",
        i = "`model_outcome` / `model_final` for a single-config comparison."),
      class = "cateomics_base_learners_conflict"
    )
  }
  # `configs` names the candidate set outright. It is the representation every
  # path resolves to, so letting it sit beside the arguments that build one
  # would give the same candidate two spellings.
  if (!is.null(configs)) {
    clash <- c(if (length(methods))             "methods",
               if (!is.null(base_learners))     "base_learners",
               if (!is.null(model_outcome))     "model_outcome",
               if (!is.null(model_final))       "model_final",
               if (!is.null(params))            "params")
    if (length(clash)) {
      rlang::abort(
        c(sprintf("`configs` cannot be combined with: %s.",
                  paste(clash, collapse = ", ")),
          i = "`configs` already names each candidate's method, learners and hyperparameters.",
          i = "Drop those arguments, or drop `configs` and use the shortcut form."),
        class = "cateomics_configs_conflict"
      )
    }
    if (!is.list(configs) || !length(configs) ||
        is.null(names(configs)) || anyNA(names(configs)) ||
        any(!nzchar(names(configs))) || anyDuplicated(names(configs))) {
      rlang::abort(
        c("`configs` must be a non-empty list with unique, non-empty names.",
          i = "The names label the configs and are the key the scoreboard, cf_select(), cf_ensemble() and the plots use."),
        class = "cateomics_configs_invalid"
      )
    }
  } else if (length(methods) == 0L) {
    rlang::abort(
      c("`methods` is required.",
        i = "Use methods = 'all' to compare all 10 CATE-capable methods,",
        i = "supply a character vector like c('grf', 'bcf', 'drlearner'), or",
        i = "name the candidate set with `configs`."),
      class = "cateomics_methods_required"
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
      class = "cateomics_compare_unknown_method"
    )
  }
  # Validate + dedup `base_learners`; inform when no method in the set is
  # pluggable (the argument then has no effect).
  if (!is.null(base_learners)) {
    valid_base <- .CF_BASE_LEARNERS
    base_learners <- unique(as.character(base_learners))
    bad_base <- setdiff(base_learners, valid_base)
    if (length(bad_base) > 0L) {
      rlang::abort(
        sprintf("Unknown base_learners: %s. Allowed: %s.",
                paste(bad_base, collapse = ", "),
                paste(valid_base, collapse = ", ")),
        class = "cateomics_base_learners_unknown"
      )
    }
    if (!any(methods %in% pluggable_methods)) {
      rlang::inform(
        c(sprintf("`base_learners` has no effect: none of the methods (%s) use a pluggable base learner.",
                  paste(methods, collapse = ", ")),
          i = sprintf("Pluggable methods: %s.",
                      paste(pluggable_methods, collapse = ", "))),
        class = "cateomics_base_learners_no_pluggable"
      )
    }
  }

  # Resolve the candidate set before the entry checks below, not after: those
  # checks read `methods`, and with `configs` the method set is only known once
  # the list has been resolved. In front of it, a `configs` call would pass
  # an unsupported bootstrap request, a vector propensity and an ignored
  # `propensity_X` without a word.
  if (!is.null(configs)) {
    # `configs` carries its own hyperparameters, so there is no three-stage
    # resolution to validate against labels here.
    resolved  <- .cf_resolve_comparison(NULL, NULL, configs, NULL, NULL, NULL,
                                        pluggable_methods, cate_methods_all)
    cfg_tab   <- resolved$table
    cfg_specs <- resolved$specs
    methods   <- unique(cfg_tab$method)
  } else {
    .cf_params_check_keys(
      params, methods = methods,
      labels = .cateomics_build_configs(methods, base_learners,
                                       pluggable_methods)$label)
    resolved  <- .cf_resolve_comparison(methods, base_learners, NULL, params,
                                        model_outcome, model_final,
                                        pluggable_methods, cate_methods_all)
    cfg_tab   <- resolved$table
    cfg_specs <- resolved$specs
    .cf_params_check_unused(params, methods = cfg_tab$method,
                            labels = cfg_tab$label)
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
  # identically (cateomics_propensity_vector_unsupported).
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cateomics_invalid_propensity_score"
      )
    }
    n_obs <- if (inherits(data, c("SummarizedExperiment", "MultiAssayExperiment")))
               ncol(data) else nrow(data)
    if (length(propensity_score) != 1L && length(propensity_score) != n_obs) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match number of samples = %d.",
                length(propensity_score), n_obs),
        class = "cateomics_invalid_propensity_score"
      )
    }
    if (length(propensity_score) > 1L &&
        any(c("drlearner", "xlearner") %in% methods)) {
      rlang::abort(
        c("Vector `propensity_score` is unsupported for drlearner / xlearner.",
          i = "econml cross-fit cannot map per-sample PS through internal folds.",
          i = "Use a scalar (`propensity_score = 0.5`) or drop drlearner / xlearner."),
        class = "cateomics_propensity_vector_unsupported"
      )
    }
  }

  metric <- .cf_match_several(metric, c("r_loss", "dr_loss", "factual_loss"),
                              "metric", "cf_compare()")

  # Validate `eval_nuisance` at entry with the helper cf_fit uses, so both
  # entries accept the same vocabulary and refuse with the same condition.
  .cf_check_eval_nuisance(eval_nuisance)
  # Resolve the per-slot spec once. A slot served by a supplied constant is
  # `"known"` and fits no learner, so recording the raw argument instead would
  # claim a Super Learner library for it; cf_fit records the resolved spec.
  eval_nuisance_spec <- .cf_eval_nuisance_spec(eval_nuisance,
                                               eval_propensity_score)
  # The estimator-side and evaluation-side known propensities are separate
  # settings, so one of them alone leaves the other estimated.
  .cf_inform_known_ps_asymmetry(propensity_score, eval_propensity_score)
  # `clip` bounds the propensity of every config, so it is validated here with
  # the same helper cf_fit uses rather than forwarded unchecked.
  .cf_check_clip(clip)
  # CI policy for the whole candidate set, before any config is fitted: a grid
  # that falls over halfway through is harder to act on than one refused here.
  # A comparison asks for several methods at once, so a request that only some
  # of them can serve is served where it can be and reported for the rest -- a
  # judgement about the set, which is why the set goes in one call. A single fit
  # has nothing to fall back on, so cf_fit refuses instead. The per-candidate
  # `model_final` still travels with it, because it decides which native
  # interval the bootstrap replaces.
  .cf_check_ci_policy(
    methods      = vapply(cfg_specs, function(.s) .s$method, character(1)),
    bootstrap_ci = bootstrap_ci,
    model_final  = lapply(cfg_specs, function(.s) .s$model_final),
    allow_partial = TRUE,
    labels       = cfg_tab$label
  )
  # `eval_library` only describes the Super Learner. Passing it beside a single
  # named learner is a setting that cannot reach the fit, which is the class of
  # argument this package refuses to leave silent.
  if (!is.null(eval_library) && !.cf_sl_any_slot_is_sl(eval_nuisance)) {
    rlang::warn(
      c("`eval_library` is not used: `eval_nuisance` names a single learner.",
        i = sprintf("Set `eval_nuisance = \"%s\"` for the library to apply.",
                    .CF_SL_TOKEN)),
      class = "cateomics_eval_library_unused"
    )
  }

  propensity_params <- .cf_params_shape(propensity_params, "propensity_params")
  # The propensity block reaches a scikit-learn estimator, whose random_state
  # is in get_params(), so without this the seed is only caught in Python --
  # after the environment boots, and as a ValueError rather than a CATEomics
  # condition.
  .cf_params_reject_rng(propensity_params, "propensity_params")
  eval_params       <- .cf_params_shape(eval_params, "eval_params")
  if (length(propensity_params) &&
      !any(methods %in% c("xlearner", "drlearner", "bcf"))) {
    rlang::abort(
      c("`propensity_params` has no effect: no method in the set has a pluggable propensity model.",
        i = "It applies to xlearner / drlearner / bcf."),
      class = "cateomics_propensity_params_unused"
    )
  }
  # The dispatcher's per-fit report is swallowed here (the fold helper runs
  # with its conditions collected rather than signalled, and a fold failure
  # becoming a `failed_methods` entry),
  # so the set-level report is raised at the entry instead.
  said_at_entry <- .cf_check_unreachable_control(
    methods,
    c(if (!is.null(propensity_score)) "propensity_score",
      if (!is.null(model_outcome))    "model_outcome",
      if (!is.null(model_final))      "model_final",
      if (!is.null(model_propensity)) "model_propensity",
      if (!is.null(clip))             "clip"))
  # Symmetry with cf_fit: a request for the GPU on a backend that has no GPU
  # path is said out loud rather than silently answered with the CPU. Said once
  # for the call, naming every method it applies to.
  if (identical(device, "gpu")) {
    cpu_only <- intersect(methods, .CF_GPU_UNSUPPORTED)
    if (length(cpu_only) > 0L) {
      rlang::warn(
        c(sprintf("device='gpu' is unsupported for %d of the compared method(s): %s (CPU-only backend).",
                  length(cpu_only), paste(cpu_only, collapse = ", ")),
          i = "Those fall back to CPU. Use the Python NN methods (cevae/tarnet/dragonnet/ganite) for GPU acceleration."),
        class = "cateomics_gpu_unsupported_warning"
      )
    }
  }
  if ("bcf" %in% methods && length(propensity_params)) {
    .cf_learner_params_check(model_propensity %||% "logistic",
                             propensity_params, "propensity_params",
                             type = "classification")
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
    omics_spec     = .facade$omics_spec,
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

  # Every configuration's splits are checked before any fitting starts, with
  # the same helper `cf_fit` uses, so a configuration one entry refuses the
  # other refuses too.
  for (.i in seq_along(cfg_specs)) {
    .s <- cfg_specs[[.i]]
    .cf_check_method_splits(
      .s$method, T, Y, binary_outcome, n_folds = n_folds,
      method_params    = .s$params,
      model_outcome    = .s$model_outcome,
      model_propensity = model_propensity,
      model_final      = .s$model_final)
  }

  folds   <- .create_stratified_folds(T, n_folds = n_folds,
                                       random_state = random_state)
  fold_id <- .folds_to_fold_id(folds, n)

  if (isTRUE(verbose)) {
    message("Fitting evaluation nuisance (",
            .cf_eval_nuisance_label(list(eval_nuisance = eval_nuisance)),
            ") OOF on ", n_folds, "-fold split...")
  }
  # A numeric `eval_propensity_score` switches the propensity slot to that value
  # instead of fitting a model, exactly as in cf_fit.
  # Same check as cf_fit, against the design matrix preprocessing produced.
  .cf_check_propensity_X(propensity_X, colnames(X) %||% cov_names)
  # Same resolution as cf_fit, from the same preprocessing record, so the two
  # entries fit the propensity on the same columns.
  propensity_X <- .cf_default_propensity_X(propensity_X,
                                           prep$clinical_cov_names,
                                           propensity_score,
                                           eval_propensity_score)
  # One copy of the large preprocessing members for every configuration: an
  # environment is serialised by reference, a matrix once per occurrence.
  prep_share <- .cf_preprocess_share(prep)
  # The shared nuisance is fitted once per fold and per library member, so a
  # learner that complains about a thin split says the same thing many times
  # for one call. It is fitted outside the per-candidate collector, so the
  # de-duplicator is what keeps it to once -- the conditions still reach the
  # caller, unchanged, just not repeated.
  nuis <- .cf_dedup_conditions(
    .compute_common_nuisance_oof(X, T, Y, folds,
                                 eval_nuisance = .cf_eval_nuisance_spec(
                                   eval_nuisance, eval_propensity_score),
                                 random_state = random_state,
                                 eval_propensity_score = eval_propensity_score,
                                 binary_outcome = binary_outcome,
                                 eval_params = eval_params,
                                 eval_library = eval_library,
                                 propensity_X = propensity_X))
  # The weights describe the nuisance but are not one of its vectors, so they
  # are recorded in the metadata and the list keeps its contract shape.
  eval_sl_weights <- nuis$sl_weights
  nuis$sl_weights <- NULL

  oof_tau    <- list()
  oof_mu0    <- list()   # per-method held-out potential outcomes (factual loss)
  oof_mu1    <- list()
  models     <- list()
  fold_cache <- if (isTRUE(keep_kfold_models)) list() else NULL
  failed     <- list()
  config_warnings <- list()
  config_messages <- list()

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

  if (isTRUE(verbose)) {
    message("Fitting ", nrow(cfg_tab), " config(s) (BPPARAM = ",
            class(BPPARAM)[1], ")...")
  }

  fit_one_config <- function(i) {
    spec <- cfg_specs[[i]]
    m    <- spec$method
    # The learner slots were resolved once, for both entry points, by
    # `.cf_resolve_comparison()`.
    base_inject <- list(model_outcome = spec$model_outcome,
                        model_final = spec$model_final,
                        model_propensity = model_propensity)
    cfg_params <- spec$params
    .cf_params_check(m, cfg_params)
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
        bootstrap_ci         = isTRUE(bootstrap_ci),
        n_bootstrap          = as.integer(n_bootstrap),
        # The spec, as cf_fit passes it. Harmless today because a
        # precomputed nuisance is always supplied, but the helper would
        # otherwise estimate a propensity the caller said was known.
        eval_nuisance        = .cf_eval_nuisance_spec(eval_nuisance,
                                                      eval_propensity_score),
        eval_library         = eval_library,
        precomputed_nuisance = nuis,
        device               = device_resolved,
        propensity_score     = propensity_score,
        clip                 = clip,
        propensity_X         = propensity_X,
        eval_propensity_score = eval_propensity_score,
        binary_outcome       = binary_outcome,
        propensity_params    = propensity_params,
        eval_params          = eval_params,
        params               = cfg_params
      ),
      base_inject
    )
    # A candidate is fitted K + 1 times, so its conditions are collected here
    # rather than signalled: the worker may be a separate process, and one
    # thin fold would otherwise speak once per fold per candidate. They are
    # returned as data and summarised once, after the loop.
    # Warnings and messages are kept apart: a backend's progress line is a
    # message, and reporting it as a warning tells the reader a caveat was
    # raised when none was.
    seen_w <- list(); seen_m <- list()
    helper_out <- withCallingHandlers(
      tryCatch(do.call(.kfold_oof_predict_matrix, call_args),
               error = function(e) e),
      warning = function(cnd) {
        seen_w[[.cf_condition_key(cnd)]] <<- stats::setNames(
          conditionMessage(cnd), .cf_condition_class(cnd))
        invokeRestart("muffleWarning")
      },
      message = function(cnd) {
        seen_m[[.cf_condition_key(cnd)]] <<- stats::setNames(
          conditionMessage(cnd), .cf_condition_class(cnd))
        invokeRestart("muffleMessage")
      }
    )
    list(label = cfg_tab$label[i], method = m, base = cfg_tab$base[i],
         model_outcome = spec$model_outcome, model_final = spec$model_final,
         out = helper_out, params = cfg_params,
         warnings = seen_w, messages = seen_m)
  }

  per_config <- BiocParallel::bplapply(seq_len(nrow(cfg_tab)), fit_one_config,
                                       BPPARAM = BPPARAM)

  for (cret in per_config) {
    lab        <- cret$label
    m          <- cret$method
    helper_out <- cret$out
    # Each entry is one named string; `unlist` would paste the de-duplication
    # key onto the name, so they are concatenated instead.
    if (length(cret$warnings))
      config_warnings[[lab]] <- do.call(c, unname(cret$warnings))
    if (length(cret$messages))
      config_messages[[lab]] <- do.call(c, unname(cret$messages))
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
    # The same lift a single fit does, so a candidate answers res$ci_method.
    .res <- .cf_lift_ci_method(helper_out$fit$res)
    .res$tau_hat_oof <- helper_out$oof_ite
    # Same field as cf_fit stores, so a per-config model and a cf_fit model
    # describe their split the same way. It is the split the helper itself
    # used, not this function's own, so the field is the one that produced
    # `tau_hat_oof` rather than one that merely agrees with it: both come from
    # `.create_stratified_folds()` on the same arguments today, and taking the
    # returned one means a future change to either cannot make the record
    # describe a different split than the estimate. Every candidate shares one
    # split, so this is `metadata$fold_id` repeated on each model, by design.
    stopifnot(identical(helper_out$fold_id, fold_id))
    .res$fold_id     <- helper_out$fold_id
    .res$nuisance    <- helper_out$nuisance
    # Same rule as cf_fit: prefer the held-out propensity when every row has
    # one, otherwise keep the full-model PS the backend already stored. A known
    # propensity keeps its label, since the held-out copy is the same number.
    if (!is.null(helper_out$oof_ps) && !anyNA(helper_out$oof_ps) &&
        !identical(.res$ps_type, "known")) {
      .res$ps      <- helper_out$oof_ps
      .res$ps_type <- "oof"
    }
    # The summary is rebuilt from the out-of-fold effects this row carries,
    # which are not the ones the backend saw. cf_compare wraps the backend fit
    # directly rather than going through cf_fit, so the rebuild has to happen
    # here too — otherwise `summary.cf_model`, which reads tau_mean and
    # tau_iqr, errors on a cf_compare-derived model.
    .res$summary <- .cf_build_summary(
      .res$ate, ate_ci = .res$ate_ci, tau_hat = .res$tau_hat
    )
    full_fit <- structure(
      list(
        method = m,
        res    = .res,
        # The backend container, flattened the same way cf_fit does it: `res`
        # was lifted to the pillar above, so leaving a copy inside would answer
        # the same question from a third field.
        model  = local({ .m <- helper_out$fit; .m$res <- NULL; .m }),
        meta   = .cf_compare_model_meta(
          prep = prep, cov_names = cov_names, n = n,
          method = m, res = .res, random_state = random_state,
          bootstrap_ci = bootstrap_ci, propensity_X = propensity_X,
          named_args = .cateomics_named_args,
          # Per config: the single Python-resolved value is what the backends
          # receive, but recording it for an R-native method would claim a GPU
          # fit that cannot happen.
          device = device,
          device_resolved = .resolve_device(device, m, .CF_GPU_UNSUPPORTED),
          outcome_type = if (isTRUE(binary_outcome)) "binary" else "continuous",
          eval_nuisance_resolved =
            .cateomics_resolve_eval_nuisance(eval_nuisance_spec),
          eval_library = if (.cf_sl_any_slot_is_sl(eval_nuisance_spec))
            .cf_sl_libraries_used(eval_nuisance_spec, eval_library) else NULL,
          eval_sl_weights = eval_sl_weights, share = prep_share,
          # The learners this config actually ran with, not the request: a
          # result has to be reproducible from what is recorded here, and
          # `cf_select()` hands this spec back for exactly that.
          spec = local({
            .s <- .cf_spec_slots(m, cret$model_outcome, cret$model_final,
                                 model_propensity, propensity_score, clip)
            list(
            method        = m,
            model_outcome = .s$model_outcome,
            model_final   = .s$model_final,
            model_propensity  = .s$model_propensity,
            clip              = .s$clip,
            params            = cret$params %||% list(),
            propensity_params = propensity_params,
            eval_params       = eval_params
          )})
        )
      ),
      class = "cf_model"
    )
    # The assembler above copied it into the estimator bucket, where `cf_fit`
    # also keeps it; leaving the result copy in place would make the field
    # appear twice for a comparison row and once for a fit.
    full_fit$res$training_trace <- NULL
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
    full_fit$data$X_preprocessed <- X
    models[[lab]] <- full_fit
    .cf_warn_ci_excludes_estimate(full_fit)
  }

  # `ids` are the uniform external keys (config labels): equal to `methods` in
  # the legacy path, or `method_base` in grid mode. Everything downstream
  # (scoreboard rows, models names, oof_tau_matrix columns, metadata$methods)
  # is keyed by `ids`; a per-row `method` / `base` decomposition is exposed in
  # the scoreboard (grid mode) and in each model's `meta$spec`.
  ids <- cfg_tab$label
  scoreboard <- .build_compare_scoreboard(
    ids        = ids,
    metric     = metric,
    oof_tau    = oof_tau,
    oof_mu0    = oof_mu0,
    oof_mu1    = oof_mu1,
    Y          = Y, T = T, nuis = nuis,
    fold_id    = fold_id,
    failed     = failed,
    # Always decomposed: the scoreboard carries `label` / `method` / `base`
    # for every call, and `base` is the learner that actually ran.
    config_df  = cfg_tab
  )
  # Keep the per-fold losses next to the scoreboard (not on it, so the printed
  # table is unchanged) for cf_select's paired comparison.
  per_fold_loss     <- attr(scoreboard, "per_fold_loss")
  per_fold_baseline <- attr(scoreboard, "per_fold_baseline")
  attr(scoreboard, "per_fold_loss")     <- NULL
  attr(scoreboard, "per_fold_baseline") <- NULL

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
    # The resolved table, named apart from the `configs` request the way
    # `eval_nuisance_resolved` is named apart from `eval_nuisance`.
    configs_resolved  = cfg_tab,
    # Hyperparameters each config resolved to, keyed by label.
    params            = stats::setNames(
      lapply(per_config, function(z) z$params %||% list()),
      vapply(per_config, function(z) z$label, character(1))),
    propensity_params = propensity_params,
    eval_params       = eval_params,
    metric            = metric,
    eval_nuisance     = eval_nuisance,
    eval_library      = if (.cf_sl_any_slot_is_sl(eval_nuisance_spec)) {
      .cf_sl_libraries_used(eval_nuisance_spec, eval_library)
    } else NULL,
    eval_sl_weights   = eval_sl_weights,
    timestamp         = Sys.time(),
    random_state      = random_state,
    n_obs             = n,
    n_covariates      = ncol(X),
    failed_methods    = if (length(failed) > 0L) failed else NULL,
    # One entry per candidate that emitted anything, each a named vector of
    # condition class -> message, de-duplicated across the folds.
    config_warnings   = if (length(config_warnings) > 0L) config_warnings
                        else NULL,
    # The same shape for what a candidate said without raising a caveat: a
    # backend's own progress output lands here, not beside the warnings.
    config_messages   = if (length(config_messages) > 0L) config_messages
                        else NULL,
    covariate_names   = cov_names,
    nuisance_oof      = nuis,
    oof_tau_matrix    = oof_tau_matrix,
    fold_id           = fold_id,
    # methods x folds, per requested metric; cf_select compares two methods on
    # the folds they share.
    per_fold_loss     = per_fold_loss,
    per_fold_baseline = per_fold_baseline
  )

  if (length(failed) > 0L) {
    rlang::warn(
      sprintf(
        "%d/%d config(s) failed: %s. See metadata$failed_methods.",
        length(failed), nrow(cfg_tab),
        paste(names(failed), collapse = ", ")
      ),
      class = "cateomics_compare_partial_failure"
    )
  }
  # The record keeps everything a candidate said, because the candidates can
  # differ. The line does not repeat what the entry already said for the whole
  # set: that is one fact about the call, and it has been said.
  cw_new <- lapply(config_warnings, function(v) v[!(names(v) %in% said_at_entry)])
  cw_new <- cw_new[vapply(cw_new, length, integer(1)) > 0L]
  if (length(cw_new) > 0L) {
    rlang::warn(
      sprintf(
        "%d/%d config(s) raised a warning: %s. See metadata$config_warnings.",
        length(cw_new), nrow(cfg_tab),
        paste(names(cw_new), collapse = ", ")
      ),
      class = "cateomics_compare_config_warnings"
    )
  }
  if (length(config_messages) > 0L) {
    rlang::inform(
      sprintf(
        "%d/%d config(s) printed a message: %s. See metadata$config_messages.",
        length(config_messages), nrow(cfg_tab),
        paste(names(config_messages), collapse = ", ")
      ),
      class = "cateomics_compare_config_messages"
    )
  }

  structure(
    list(
      scoreboard = scoreboard,
      models     = models,
      fold_cache = fold_cache,
      metadata   = metadata
    ),
    class = c("cateomics_compare", "list")
  )
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Resolve the comparison space into one internal representation
#'
#' Both entry points land here: the shortcut arguments (`methods`,
#' `base_learners`, `model_outcome`, `model_final`, `params`) and an explicit
#' `configs` list. Everything downstream reads only the result, so every
#' existing call exercises the same representation an explicit list produces.
#'
#' @param methods,base_learners,params,model_outcome,model_final The shortcut
#'   arguments, or NULL when `configs` is given.
#' @param configs The explicit named list, or NULL.
#' @return A list with `table` (one row per config: `label`, `method`, and
#'   `base` = the base learner that actually runs, `NA` for methods with no
#'   pluggable slot) and `specs` (one `list(method, model_outcome, model_final,
#'   params)` per row, in the same order).
#' @keywords internal
#' @noRd
.cf_resolve_comparison <- function(methods, base_learners, configs, params,
                                   model_outcome, model_final,
                                   pluggable_methods, cate_methods_all) {
  # The base learner each config actually runs with: the explicit value when
  # there is one, otherwise the fitter's own default. Recording the request
  # instead would leave the default path blank.
  effective_base <- function(method, mo) {
    if (!(method %in% pluggable_methods)) return(NA_character_)
    .cf_resolve_learner_slot(.cf_method_fitter(method), "model_outcome",
                             mo) %||% NA_character_
  }
  if (!is.null(configs)) {
    labels <- names(configs)
    allowed <- c("method", "model_outcome", "model_final", "params")
    unreachable <- list()
    specs <- lapply(seq_along(configs), function(i) {
      cfg <- configs[[i]]
      lab <- labels[i]
      if (!is.list(cfg)) {
        rlang::abort(
          sprintf("`configs$%s` must be a list(method = , ...).", lab),
          class = "cateomics_configs_invalid"
        )
      }
      extra <- setdiff(names(cfg), allowed)
      if (length(extra)) {
        rlang::abort(
          c(sprintf("`configs$%s` has unknown entries: %s.", lab,
                    paste(extra, collapse = ", ")),
            i = sprintf("A config takes: %s.", paste(allowed, collapse = ", ")),
            i = "The propensity and evaluation-nuisance settings are single across configs and stay as arguments of cf_compare()."),
          class = "cateomics_configs_invalid"
        )
      }
      m <- cfg$method
      if (!is.character(m) || length(m) != 1L || !(m %in% cate_methods_all)) {
        rlang::abort(
          c(sprintf("`configs$%s$method` must be one of the CATE methods.", lab),
            i = sprintf("Allowed: %s.", paste(cate_methods_all, collapse = ", "))),
          class = "cateomics_configs_invalid"
        )
      }
      for (slot in c("model_outcome", "model_final")) {
        v <- cfg[[slot]]
        if (!is.null(v)) .validate_learner_spec(v)
      }
      # Naming a learner this method has no slot for is collected here and
      # reported once for the call: the same mistake in three configs is one
      # fact about the call, and the line names the configs it applies to.
      lost <- setdiff(
        intersect(c(if (!is.null(cfg$model_outcome)) "model_outcome",
                    if (!is.null(cfg$model_final))   "model_final"),
                  .cf_watched_control_args(m)),
        names(formals(.cf_method_fitter(m))))
      if (length(lost)) {
        unreachable[[lab]] <<- sprintf("%s (%s): %s", lab, m,
                                       paste0("`", lost, "`", collapse = " / "))
      }
      prm <- .cf_params_shape(cfg$params, sprintf("configs$%s$params", lab))
      .cf_params_check(m, prm)
      list(method = m, model_outcome = cfg$model_outcome,
           model_final = cfg$model_final, params = prm)
    })
    if (length(unreachable) > 0L) {
      rlang::warn(
        c(sprintf("%d config(s) name a learner slot their method does not have; the value(s) are not used.",
                  length(unreachable)),
          i = paste(unlist(unreachable, use.names = FALSE), collapse = "; "),
          .cf_control_arg_hints()),
        class = "cateomics_control_arg_unreachable"
      )
    }
    tab <- data.frame(
      label  = labels,
      method = vapply(specs, `[[`, character(1), "method"),
      base   = vapply(specs, function(s) effective_base(s$method,
                                                        s$model_outcome),
                      character(1)),
      stringsAsFactors = FALSE
    )
    return(list(table = tab, specs = specs))
  }
  tab <- .cateomics_build_configs(methods, base_learners, pluggable_methods)
  grid_mode <- !is.null(base_learners)
  specs <- lapply(seq_len(nrow(tab)), function(i) {
    m   <- tab$method[i]
    # In grid mode the config's CATE-driving regressor(s) take its base
    # learner: model_outcome for every pluggable method, plus model_final for
    # the ones with a Stage-2 CATE model. Non-pluggable methods ignore it. The
    # mutual-exclusion guard means model_outcome / model_final are NULL here.
    mo <- model_outcome
    mf <- model_final
    if (grid_mode && m %in% pluggable_methods && !is.na(tab$base[i])) {
      mo <- tab$base[i]
      if (m %in% c("xlearner", "drlearner")) mf <- tab$base[i]
    }
    list(method = m, model_outcome = mo, model_final = mf,
         params = .cf_resolve_params(params, method = m, label = tab$label[i],
                                      methods = methods, labels = tab$label))
  })
  tab$base <- vapply(specs, function(s) effective_base(s$method,
                                                       s$model_outcome),
                     character(1))
  list(table = tab, specs = specs)
}


#' Expand the (method x base_learner) grid into a config table.
#'
#' Returns a data.frame with one row per fit config: \code{label} (external key),
#' \code{method} (bare method), \code{base} (base learner or NA). With
#' \code{base_learners = NULL} this is one row per method (label == method,
#' base = NA). Otherwise pluggable methods fan out to one row per base learner
#' (label = method_base); non-pluggable methods stay single (base = NA).
#' @keywords internal
#' @noRd
.cateomics_build_configs <- function(methods, base_learners, pluggable_methods) {
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
                                       config_df) {
  base_const <- mean(.compute_dr_pseudo_outcome(Y, T, nuis$mu0_oof,
                                                  nuis$mu1_oof, nuis$pi_oof))
  # R-loss baseline = the R-loss-optimal constant-effect model,
  # theta* = argmin_theta sum((Yres - theta * Tres)^2) = sum(Yres Tres)/sum(Tres^2),
  # so the r_squared_score measures heterogeneity captured beyond the best
  # constant effect. NA when Tres has no variation.
  Yres_c <- as.numeric(Y) - nuis$m_oof
  Tres_c <- as.numeric(T) - nuis$pi_oof
  sum_TT <- sum(Tres_c^2)
  theta_star <- if (is.finite(sum_TT) && sum_TT > 1e-12) {
    sum(Yres_c * Tres_c) / sum_TT
  } else NA_real_
  # One shape for every call: `label` is the id, `method` the bare method and
  # `base` the base learner that ran (NA where the method has no such slot).
  sb <- data.frame(label = ids, method = config_df$method,
                   base = config_df$base, stringsAsFactors = FALSE)
  per_fold_loss     <- list()
  per_fold_baseline <- list()

  for (m_metric in metric) {
    loss_vals     <- numeric(length(ids))
    loss_lo       <- numeric(length(ids))
    loss_hi       <- numeric(length(ids))
    r2_vals       <- numeric(length(ids))
    # Per-fold losses, kept so cf_select can compare two methods on the folds
    # they share instead of treating their intervals as independent.
    n_folds_sb    <- length(unique(fold_id))
    per_fold_mat  <- matrix(NA_real_, nrow = length(ids), ncol = n_folds_sb,
                            dimnames = list(ids, NULL))
    per_fold_base <- rep(NA_real_, n_folds_sb)

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
      per_fold_mat[i, ] <- per_fold
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
      if (all(is.na(per_fold_base))) {
        per_fold_base <- .aggregate_loss_by_fold(base_per_obs, fold_id)
      }
    }
    per_fold_loss[[m_metric]] <- per_fold_mat
    per_fold_baseline[[m_metric]] <- per_fold_base

    sb[[paste0(m_metric)]]              <- loss_vals
    sb[[paste0(m_metric, "_ci_lower")]] <- loss_lo
    sb[[paste0(m_metric, "_ci_upper")]] <- loss_hi
    sb[[sub("_loss$", "_squared_score", m_metric)]] <- r2_vals
  }

  tau_mean <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau))) NA_real_ else mean(tau, na.rm = TRUE)
  }, numeric(1))
  tau_sd <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau))) NA_real_ else stats::sd(tau, na.rm = TRUE)
  }, numeric(1))
  tau_q <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau)))
      c(NA_real_, NA_real_, NA_real_)
    else
      unname(stats::quantile(tau, c(0.25, 0.5, 0.75), na.rm = TRUE))
  }, numeric(3))

  sb$tau_mean <- tau_mean
  sb$tau_sd  <- tau_sd
  sb$tau_q25  <- tau_q[1, ]
  sb$tau_q50  <- tau_q[2, ]
  sb$tau_q75  <- tau_q[3, ]

  # How the spread of the out-of-fold tau-hat divides between folds and within
  # them. `between_fold_share` is the share of var(tau_hat_oof) explained by the
  # fold assignment, and the two standard deviations are that split applied to
  # sd(tau_hat_oof), so within^2 + between^2 = tau_sd^2. The share is NA when
  # tau_hat_oof is constant, where var is 0 and the ratio is undefined.
  # Both variances are taken over the same rows, the ones with an effect and a
  # fold. That makes the ratio the analysis-of-variance decomposition
  # SS_between / SS_total exactly, hence in [0, 1], and puts it on the same
  # rows as `tau_sd`. Dropping the NA rows from only one of the two gives a
  # ratio that can leave [0, 1], where the pmax below would quietly report all
  # of the spread as between-fold.
  fold_var <- vapply(ids, function(m) {
    tau <- oof_tau[[m]]
    if (is.null(tau) || all(is.na(tau))) return(NA_real_)
    ok <- !is.na(tau) & !is.na(fold_id)
    t  <- tau[ok]
    if (length(t) < 2L) return(NA_real_)
    v <- stats::var(t)
    if (!is.finite(v) || v <= 0) return(NA_real_)
    stats::var(stats::ave(t, fold_id[ok], FUN = mean)) / v
  }, numeric(1))

  sb$within_fold_sd     <- tau_sd * sqrt(pmax(0, 1 - fold_var))
  sb$between_fold_sd    <- tau_sd * sqrt(pmax(0, fold_var))
  sb$between_fold_share <- fold_var

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

  attr(sb, "per_fold_loss")     <- per_fold_loss
  attr(sb, "per_fold_baseline") <- per_fold_baseline
  sb
}

#' @export
print.cateomics_compare <- function(x, full = FALSE, ...) {
  cat("<cateomics_compare>: ", nrow(x$scoreboard), " method(s) on ",
      x$metadata$n_folds, "-fold split, n=", x$metadata$n_obs,
      ", p=", x$metadata$n_covariates, "\n", sep = "")
  cat("Metric(s): ", paste(x$metadata$metric, collapse = ", "), "\n", sep = "")
  cat("Eval nuisance: ", .cf_eval_nuisance_label(x$metadata), "\n", sep = "")
  # Which members the combination leaned on, directly under the nuisance they
  # describe. The outcome side is shown because every metric reads
  # mu0 / mu1 / m; the propensity weights are in the metadata.
  w_line <- .cf_sl_weight_line(x$metadata$eval_sl_weights)
  if (!is.na(w_line)) cat("  weights for ", w_line, "\n", sep = "")
  # Both prints show statistics of tau-hat, from different estimates: the
  # scoreboard is scored out of fold, while print.cf_model shows the fit's own
  # in-sample values. Saying which is which costs a line and stops the two
  # being read as the same number.
  cat("Losses and tau statistics are out-of-fold.\n")
  sb <- x$scoreboard
  primary_metric <- x$metadata$metric[1]
  if (isTRUE(full)) {
    print(sb)
  } else {
    # `base` says which base learner each row ran with, so it earns its place
    # only when some row has one: a set of methods with no pluggable outcome
    # slot would otherwise print a column of NA.
    id_cols   <- if (any(!is.na(sb$base))) c("label", "base") else "label"
    keep_cols <- c(id_cols, primary_metric, "tau_mean", "tau_sd",
                   "between_fold_share", "rank")
    print(sb[, intersect(keep_cols, colnames(sb)), drop = FALSE])
  }
  if (!is.null(x$metadata$failed_methods)) {
    cat("Failed methods: ",
        paste(names(x$metadata$failed_methods), collapse = ", "),
        "\n", sep = "")
  }
  # Only mention hyperparameters when some config departs from the defaults.
  .cf_print_params_note(x$metadata$params, x$metadata$propensity_params,
                        x$metadata$eval_params, full = isTRUE(full))
  invisible(x)
}


#' One-line note about non-default hyperparameters
#'
#' @keywords internal
#' @noRd
.cf_print_params_note <- function(params, propensity_params, eval_params,
                                  full = FALSE) {
  fmt <- function(x, prefix = "") {
    if (!length(x)) return(character(0))
    unlist(lapply(names(x), function(nm) {
      v <- x[[nm]]
      if (is.list(v)) fmt(v, prefix = paste0(prefix, nm, "$"))
      else sprintf("%s%s=%s", prefix, nm,
                   paste(format(v, trim = TRUE), collapse = ","))
    }), use.names = FALSE)
  }
  parts <- character(0)
  if (length(params)) {
    for (lab in names(params)) {
      f <- fmt(params[[lab]])
      if (length(f)) parts <- c(parts, sprintf("%s: %s", lab,
                                               paste(f, collapse = ", ")))
    }
  }
  pp <- fmt(propensity_params); ep <- fmt(eval_params)
  if (length(pp)) parts <- c(parts, paste("propensity:", paste(pp, collapse = ", ")))
  if (length(ep)) parts <- c(parts, paste("eval:", paste(ep, collapse = ", ")))
  if (!length(parts)) return(invisible(NULL))
  line <- paste(parts, collapse = " | ")
  if (!full && nchar(line) > 100L) line <- paste0(substr(line, 1L, 97L), "...")
  cat("Params: ", line, "\n", sep = "")
  invisible(NULL)
}

#' Summarise a cf_compare result.
#'
#' One screen: how the candidates were scored, which one the selection metric
#' puts first, how far ahead of the next it is, and what its estimates look
#' like. It is deliberately not the scoreboard -- \code{print(x)} shows the
#' columns that fit, and \code{print(x, full = TRUE)} shows every one.
#'
#' @param object A \code{cateomics_compare} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.cateomics_compare} holding the
#'   scoreboard and the comparison's metadata. Printing it gives the one
#'   screen; \code{print(x, full = TRUE)} on the comparison itself gives
#'   every column.
#' @export
summary.cateomics_compare <- function(object, ...) {
  # A summary object of its own class, as `summary` on a fit or an ensemble
  # returns: handing the comparison back meant that printing a stored summary
  # emitted this screen and then dispatched to the comparison's own print,
  # so the compact scoreboard came out a second time.
  structure(list(scoreboard = object$scoreboard,
                 metadata   = object$metadata),
            class = "summary.cateomics_compare")
}


#' @export
print.summary.cateomics_compare <- function(x, ...) {
  sb <- x$scoreboard
  md <- x$metadata
  cat("cf_compare summary\n")
  # A print says what it cannot find rather than stopping part way: an object
  # written by an earlier version, or one a caller has edited, still has to be
  # readable enough to show what is wrong with it.
  if (is.null(sb) || !nrow(sb)) {
    cat("  no scoreboard in this object\n")
    return(invisible(x))
  }
  metric <- if (length(md$metric)) md$metric[[1L]] else NA_character_
  cat("  ", nrow(sb), " candidate(s), ", md$n_folds, "-fold split, n=",
      md$n_obs, ", p=", md$n_covariates, "\n", sep = "")
  dir_txt <- if (metric %in% names(.cf_metric_directions)) {
    if (isTRUE(unname(.cf_metric_directions[[metric]]))) " (higher is better)"
    else " (lower is better)"
  } else ""
  cat("  selection metric: ", metric, dir_txt,
      " | eval nuisance: ", .cf_eval_nuisance_label(md), "\n", sep = "")
  cat("  losses and tau statistics are out-of-fold\n")

  ord <- if ("rank" %in% colnames(sb)) order(sb$rank) else seq_len(nrow(sb))
  sb  <- sb[ord, , drop = FALSE]
  val <- if (metric %in% colnames(sb)) sb[[metric]] else rep(NA_real_, nrow(sb))
  lo  <- sb[[paste0(metric, "_ci_lower")]]
  hi  <- sb[[paste0(metric, "_ci_upper")]]
  # One format for the estimate and for its interval: with different precisions
  # an interval can print as if it excluded the estimate, which is a thing that
  # really happens here and has its own warning.
  num <- function(x) sprintf("%.6g", x)
  ci_txt <- function(i) {
    if (is.null(lo) || is.null(hi) || is.na(lo[i]) || is.na(hi[i])) ""
    else sprintf(" [%s, %s]", num(lo[i]), num(hi[i]))
  }
  # With no score on any candidate there is no best and no ranking: the order
  # would be the order the configs were given in, and `cf_select()` refuses the
  # same x. Say that instead of naming the first row.
  if (all(is.na(val))) {
    cat("\n  no candidate has a ", metric, ": nothing to rank\n", sep = "")
  } else {
  cat("\n  best: ", sb$label[1L], "   ", metric, " = ",
      num(val[1L]), ci_txt(1L), "\n", sep = "")
  if (nrow(sb) >= 2L && !is.na(val[2L])) {
    # The gap to the runner-up is what a selection rests on, so it is given
    # rather than left to be read off two numbers. Whether it separates them is
    # `cf_select()`'s question: it compares the two on the folds they share.
    cat("  next: ", sb$label[2L], "   ", metric, " = ",
        num(val[2L]), ci_txt(2L),
        sprintf("   (gap %+.4g)", val[2L] - val[1L]), "\n", sep = "")
  }
  scored <- sb$label[!is.na(val)]
  if (length(scored) >= 2L) {
    cat("  ranking: ", paste(utils::head(scored, 4L), collapse = " > "),
        if (length(scored) > 4L) sprintf(" > ... (of %d)", length(scored)) else "",
        "\n", sep = "")
  }
  }
  if (!all(is.na(val)) && all(c("tau_mean", "tau_sd") %in% colnames(sb))) {
    bfs <- if ("between_fold_share" %in% colnames(sb)) sb$between_fold_share[1L] else NA_real_
    cat("\n  its estimates: mean ", sprintf("%.4g", sb$tau_mean[1L]),
        ", sd ", sprintf("%.4g", sb$tau_sd[1L]),
        if (is.na(bfs)) "" else sprintf(", between-fold share %.3g", bfs),
        "\n", sep = "")
  }
  if (!is.null(md$failed_methods)) {
    cat("  failed: ", paste(names(md$failed_methods), collapse = ", "),
        "\n", sep = "")
  }
  cat("\n  every column: print(x, full = TRUE)\n")
  invisible(x)
}
