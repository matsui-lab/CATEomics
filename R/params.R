# Hyperparameter plumbing: a single `params` entry point for every method.
#
# User-facing surface:
#   params            method hyperparameters, resolved per fit config
#                     (cf_fit / cf_compare)
#   propensity_params hyperparameters of the method's propensity model
#                     (cf_fit / cf_compare)
#   eval_params       hyperparameters of the evaluation nuisance models
#                     (cf_fit / cf_compare / cf_evaluate)
#
# Resolution is three-stage and does NOT merge across stages: the first stage
# that matches wins, so a user reading `params` can point at exactly one block
# and know what a config received.
#
#   1. params[[label]]   config label in a base-learner grid, e.g. "drlearner_rf"
#   2. params[[method]]  method name, e.g. "drlearner"
#   3. the global block  every entry whose name is neither a label nor a method
#
# Internal fitters keep their own named arguments, which is where the CATEomics
# defaults live; the resolved block is split against those formals and the
# remainder travels in the fitter's `params` argument.

.CF_PARAM_METHODS <- c("grf", "bcf", "slearner", "tlearner", "xlearner",
                       "drlearner", "cevae", "ganite", "tarnet", "dragonnet")

# Nested blocks a method understands inside its own params. Anything else is a
# flat argument of that method's fitter.
.CF_PARAM_SLOTS <- list(
  grf       = "nuisance",
  bcf       = character(0),
  slearner  = "outcome",
  tlearner  = "outcome",
  xlearner  = c("outcome", "final"),
  drlearner = c("outcome", "final"),
  cevae     = character(0),
  ganite    = character(0),
  tarnet    = character(0),
  dragonnet = character(0)
)

# Arguments the fitters receive from CATEomics itself. They are never user
# hyperparameters, so they are excluded when the allow-list is derived from a
# fitter's formals.
.CF_PARAM_CONTROL_ARGS <- c(
  "X", "T", "Y", "covariate_names", "random_state", "device",
  "binary_outcome", "bootstrap_ci", "n_bootstrap", "clip",
  "propensity_score", "propensity_X", "pihat",
  "model_outcome", "model_final", "model_propensity",
  "propensity_params", "binary_y", "params", "..."
)

# The names upstream packages use for their own random seed. `random_state` is
# the single knob that drives every draw CATEomics makes, so a seed set through
# `params` would quietly take over from it for that one fitter while the rest
# of the run still follows `random_state`. grf reads `seed`, bcf reads
# `random_seed`, and the scikit-learn estimators behind the meta-learners read
# `random_state`, so all three spellings are refused.
.CF_PARAM_RNG_ARGS <- c("seed", "random_seed", "random_state")

.CF_PARAM_CONTROL_ARGS <- unique(c(.CF_PARAM_CONTROL_ARGS, .CF_PARAM_RNG_ARGS))

#' Check the shape of a user-supplied params argument
#' @keywords internal
#' @noRd
.cf_params_shape <- function(params, arg = "params") {
  if (is.null(params)) return(list())
  if (!is.list(params) || (length(params) > 0L && is.null(names(params)))) {
    rlang::abort(
      sprintf("`%s` must be a named list (got: %s).", arg, class(params)[1L]),
      class = "cateomics_invalid_params"
    )
  }
  nms <- names(params)
  if (length(params) > 0L && (any(is.na(nms)) || any(!nzchar(nms)))) {
    rlang::abort(sprintf("Every element of `%s` must be named.", arg),
                 class = "cateomics_invalid_params")
  }
  if (anyDuplicated(nms)) {
    rlang::abort(
      sprintf("`%s` has duplicated names: %s.", arg,
              paste(unique(nms[duplicated(nms)]), collapse = ", ")),
      class = "cateomics_invalid_params"
    )
  }
  params
}

#' Reject params keys that are neither a fit config nor a hyperparameter
#'
#' A list-valued entry is read as a per-method / per-config block, so a typo in
#' a method name would otherwise be silently treated as a hyperparameter called
#' "grff" and only surface further down. Catch it at the entry point.
#'
#' @keywords internal
#' @noRd
.cf_params_check_keys <- function(params, methods, labels = NULL,
                                  arg = "params") {
  params <- .cf_params_shape(params, arg)
  if (length(params) == 0L) return(invisible(TRUE))
  slots <- unique(unlist(.CF_PARAM_SLOTS, use.names = FALSE))
  valid <- unique(c(methods, labels, slots))
  for (nm in names(params)) {
    if (nm %in% valid) next
    if (is.list(params[[nm]])) {
      rlang::abort(
        c(sprintf("`%s` has an unknown key: %s.", arg, nm),
          i = sprintf("Expected one of the fitted configs (%s)%s, or a flat hyperparameter.",
                      paste(methods, collapse = ", "),
                      if (length(labels)) paste0(" / labels (",
                                                 paste(labels, collapse = ", "), ")") else "")),
        class = "cateomics_unknown_params_key"
      )
    }
  }
  invisible(TRUE)
}

#' Reject a global params block that no fit config would ever read
#'
#' Stages do not merge, so when every config has its own block the leftover
#' flat entries apply to nothing. Silently dropping them would be the kind of
#' argument that looks set but does not reach the fit.
#'
#' @keywords internal
#' @noRd
.cf_params_check_unused <- function(params, methods, labels = NULL,
                                    arg = "params") {
  params <- .cf_params_shape(params, arg)
  if (length(params) == 0L) return(invisible(TRUE))
  keys <- unique(c(methods, labels))
  global <- params[setdiff(names(params), keys)]
  if (length(global) == 0L) return(invisible(TRUE))
  labels <- labels %||% rep(NA_character_, length(methods))
  reads_global <- vapply(seq_along(methods), function(i) {
    lab <- labels[[i]]
    is.null(params[[methods[[i]]]]) &&
      (is.na(lab) || is.null(params[[lab]]))
  }, logical(1))
  if (any(reads_global)) return(invisible(TRUE))
  rlang::abort(
    c(sprintf("`%s` has entries that no fit would read: %s.",
              arg, paste(names(global), collapse = ", ")),
      i = "Every method already has its own block, and the stages are not merged.",
      i = "Move them into the per-method blocks, or drop the per-method blocks."),
    class = "cateomics_unused_params"
  )
}

#' Resolve the params block that applies to one fit config
#' @keywords internal
#' @noRd
.cf_resolve_params <- function(params, method, label = NULL,
                               methods = NULL, labels = NULL,
                               arg = "params") {
  params <- .cf_params_shape(params, arg)
  if (length(params) == 0L) return(list())
  if (!is.null(label) && !is.null(params[[label]])) {
    return(.cf_params_shape(params[[label]], paste0(arg, "$", label)))
  }
  if (!is.null(params[[method]])) {
    return(.cf_params_shape(params[[method]], paste0(arg, "$", method)))
  }
  keys <- unique(c(methods %||% .CF_PARAM_METHODS, labels))
  params[setdiff(names(params), keys)]
}

#' The internal fitter a method name dispatches to
#'
#' One mapping, so a check that has to look at a method's own formals (the
#' internal cross-fit count, for instance) cannot drift from the allow-list
#' that reads the same formals.
#'
#' @return The fitter function, or NULL for a name with none.
#' @keywords internal
#' @noRd
.cf_method_fitter <- function(method) {
  fitter <- switch(method,
    grf = "cf_fit_grf", bcf = "cf_fit_bcf",
    slearner = "cf_fit_slearner", tlearner = "cf_fit_tlearner",
    xlearner = "cf_fit_xlearner", drlearner = "cf_fit_drlearner",
    cevae = "cf_fit_cevae", ganite = "cf_fit_ganite",
    tarnet = "cf_fit_tarnet", dragonnet = "cf_fit_dragonnet",
    NULL)
  if (is.null(fitter) || !exists(fitter, mode = "function")) return(NULL)
  get(fitter, mode = "function")
}

#' Names a method accepts as flat hyperparameters
#'
#' NULL means "not checked on the R side": the three methods it applies to
#' (`cevae`, `tarnet`, `dragonnet`) forward flat entries to a pyro or catenets
#' constructor that raises on an unknown keyword. The R-side methods need an
#' explicit list because ranger only warns and glmnet ignores silently, and
#' `ganite` needs one because it is implemented in this package.
#'
#' @keywords internal
#' @noRd
.cf_params_allowed <- function(method) {
  fitter <- switch(method,
    grf = "cf_fit_grf", bcf = "cf_fit_bcf",
    slearner = "cf_fit_slearner", tlearner = "cf_fit_tlearner",
    xlearner = "cf_fit_xlearner", drlearner = "cf_fit_drlearner",
    cevae = "cf_fit_cevae", ganite = "cf_fit_ganite",
    tarnet = "cf_fit_tarnet", dragonnet = "cf_fit_dragonnet",
    NULL)
  if (is.null(fitter) || !exists(fitter, mode = "function")) return(NULL)
  own <- setdiff(names(formals(get(fitter, mode = "function"))),
                 .CF_PARAM_CONTROL_ARGS)
  if (method == "grf") {
    if (!requireNamespace("grf", quietly = TRUE)) return(NULL)
    # `Y.hat` / `W.hat` are the orthogonalization predictions CATEomics fits
    # itself and passes by name, so a second copy coming through here
    # would reach the same formal twice. A known propensity has its own
    # argument (`propensity_score`), and how the forests are fitted is open
    # through `params$grf$nuisance`.
    return(setdiff(union(own, names(formals(grf::causal_forest))),
                   c("X", "Y", "W", "Y.hat", "W.hat",
                     .CF_PARAM_CONTROL_ARGS)))
  }
  if (method == "bcf") {
    if (!requireNamespace("bcf", quietly = TRUE)) return(NULL)
    # `nburn` and `nsim` are upstream's own names and are taken by the CATEomics
    # fitter because bcf has no default for either; `n_chains` has no upstream
    # counterpart because CATEomics runs the chains itself.
    return(setdiff(union(own, names(formals(bcf::bcf))),
                   c("y", "z", "x_control", "x_moderate", "pihat",
                     .CF_PARAM_CONTROL_ARGS)))
  }
  # The meta / DR learners consume their nested blocks only, so a flat entry
  # would have nowhere to go: check it here rather than dropping it.
  if (method %in% c("slearner", "tlearner", "xlearner", "drlearner")) return(own)
  # ganite is implemented in this package, so its hyperparameters are its own
  # formals and are checked here. The other neural methods forward flat entries
  # to an upstream constructor that raises on an unknown keyword, so they stay
  # unchecked on the R side apart from the seed names.
  if (method == "ganite") return(own)
  NULL
}

#' Refuse a random seed supplied through `params`
#'
#' Applies to every method, including the ones with no allow-list: their
#' backend would accept a seed keyword and silently take over from
#' `random_state`.
#'
#' @keywords internal
#' @noRd
.cf_params_reject_rng <- function(plist, where) {
  hit <- intersect(names(plist %||% list()), .CF_PARAM_RNG_ARGS)
  if (!length(hit)) return(invisible(TRUE))
  rlang::abort(
    c(sprintf("`%s` cannot set %s.", where, paste(hit, collapse = " / ")),
      i = "Every draw follows `random_state`; set that instead."),
    class = "cateomics_rng_param"
  )
}

#' Names a grf nuisance forest accepts
#' @keywords internal
#' @noRd
.cf_params_allowed_slot <- function(method, slot) {
  if (method == "grf" && identical(slot, "nuisance")) {
    if (!requireNamespace("grf", quietly = TRUE)) return(NULL)
    return(setdiff(names(formals(grf::regression_forest)), c("X", "Y")))
  }
  NULL  # outcome / final go to sklearn, which validates on set_params()
}

#' Validate one resolved params block against its method
#' @keywords internal
#' @noRd
.cf_params_check <- function(method, plist, arg = "params") {
  plist <- .cf_params_shape(plist, arg)
  if (length(plist) == 0L) return(invisible(TRUE))
  slots_ok <- .CF_PARAM_SLOTS[[method]] %||% character(0)
  reserved <- unique(unlist(.CF_PARAM_SLOTS, use.names = FALSE))
  used_slots <- intersect(names(plist), reserved)
  bad_slots <- setdiff(used_slots, slots_ok)
  if (length(bad_slots)) {
    rlang::abort(
      c(sprintf("Method '%s' has no %s model.", method,
                paste(bad_slots, collapse = " / ")),
        i = if (length(slots_ok))
              sprintf("It accepts: %s.", paste(slots_ok, collapse = ", "))
            else sprintf("`%s$%s` takes flat hyperparameters only.", arg, method)),
      class = "cateomics_unknown_params_slot"
    )
  }
  for (s in used_slots) {
    sub <- .cf_params_shape(plist[[s]], paste0(arg, "$", s))
    # outcome / final go straight to a scikit-learn estimator, which accepts
    # `random_state`, so the seed names are refused here as well.
    .cf_params_reject_rng(sub, sprintf("%s$%s$%s", arg, method, s))
    allowed <- .cf_params_allowed_slot(method, s)
    if (!is.null(allowed)) .cf_params_reject(setdiff(names(sub), allowed),
                                             allowed, paste0(arg, "$", s))
  }
  flat <- plist[setdiff(names(plist), reserved)]
  # Before the generic control-argument message, which would point at
  # cf_fit()'s own arguments and none of them is a seed.
  .cf_params_reject_rng(flat, sprintf("%s (method '%s')", arg, method))
  # CATEomics supplies these itself; letting one through `params` would quietly
  # override the corresponding argument of cf_fit / cf_compare.
  clash <- intersect(names(flat), .CF_PARAM_CONTROL_ARGS)
  if (length(clash)) {
    rlang::abort(
      c(sprintf("`%s` must not set CATEomics-controlled argument(s): %s.",
                arg, paste(clash, collapse = ", ")),
        i = "Pass them directly to cf_fit() / cf_compare() instead."),
      class = "cateomics_params_control_clash"
    )
  }
  is_block <- vapply(flat, is.list, logical(1))
  if (any(is_block)) {
    rlang::abort(
      c(sprintf("`%s` entries must be scalars or vectors, not lists: %s.",
                arg, paste(names(flat)[is_block], collapse = ", ")),
        i = "Nested lists are reserved for per-method blocks and for outcome / final / nuisance."),
      class = "cateomics_invalid_params"
    )
  }
  allowed <- .cf_params_allowed(method)
  if (!is.null(allowed)) {
    .cf_params_reject(setdiff(names(flat), allowed), allowed,
                      sprintf("%s (method '%s')", arg, method),
                      slots = slots_ok)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.cf_params_reject <- function(unknown, allowed, where, slots = character(0)) {
  if (!length(unknown)) return(invisible(TRUE))
  hint <- if (length(allowed)) {
    sprintf("Accepted here: %s.", paste(sort(allowed), collapse = ", "))
  } else if (length(slots)) {
    sprintf("This method takes no flat hyperparameters; put them in %s.",
            paste(sprintf("`%s`", slots), collapse = " / "))
  } else {
    "This method takes no hyperparameters of that kind."
  }
  rlang::abort(
    c(sprintf("Unknown hyperparameter(s) in `%s`: %s.",
              where, paste(unknown, collapse = ", ")),
      i = hint),
    class = "cateomics_unknown_param"
  )
}

#' Split a resolved params block against a fitter's formals and call it
#'
#' Names the fitter declares become named arguments, so the CATEomics defaults in
#' its signature still apply when the user says nothing. The remainder travels
#' as `params`, which the fitter forwards to its backend.
#'
#' @keywords internal
#' @noRd
.cf_apply_params <- function(fun, base_args, plist) {
  plist <- plist %||% list()
  known <- names(formals(fun))
  hit   <- plist[intersect(names(plist), setdiff(known, "params"))]
  rest  <- plist[setdiff(names(plist), names(hit))]
  if (length(rest)) base_args$params <- rest
  do.call(fun, c(base_args, hit))
}

#' Validate hyperparameters aimed at an R-side base learner
#'
#' Used for `eval_params` and for `propensity_params` when the propensity model
#' is fitted in R. ranger warns on unknown arguments instead of failing and
#' glmnet drops them silently, so unknown names are rejected here.
#'
#' @keywords internal
#' @noRd
.cf_learner_params_check <- function(learner, plist, arg, type = NULL) {
  plist <- .cf_params_shape(plist, arg)
  if (length(plist) == 0L) return(invisible(TRUE))
  .cf_params_reject_rng(plist, arg)
  # Derived from the upstream functions rather than listed by hand, so the
  # set is whatever that version of the package accepts. A function that
  # forwards through `...` does not name the forwarded arguments in its own
  # formals, so the entry point is unioned with what it forwards to:
  # cv.glmnet passes to glmnet (nfolds lives on the former, penalty.factor on
  # the latter), glm passes to glm.control, and xgb.train takes the booster
  # hyperparameters inside `params`, enumerated by xgb.params.
  allowed <- switch(learner,
    # `.fit_linear()` runs lm() for a regression and glm() for a
    # classification, and they take different arguments: glm's own and
    # glm.control's are never read by a regression, and lm's `qr` is not read
    # by a classification. So the set follows the resolved outcome type where
    # the caller knows it. `qr` is excluded either way: glm forwards it to
    # glm.control(), which has no `...`, and a regression fitted with
    # qr = FALSE then fails in predict.lm().
    linear = ,
    logistic = {
      glm_names <- union(names(formals(stats::glm)),
                         names(formals(stats::glm.control)))
      pool <- switch(type %||% "unknown",
        regression     = names(formals(stats::lm)),
        classification = glm_names,
        union(names(formals(stats::lm)), glm_names))
      setdiff(pool, c("formula", "data", "family", "y", "x", "qr", "..."))
    },
    lasso = ,
    ridge = ,
    elastic_net = setdiff(union(names(formals(glmnet::cv.glmnet)),
                                names(formals(glmnet::glmnet))),
                          c("x", "y", "alpha", "family", "standardize", "...")),
    rf = setdiff(names(formals(ranger::ranger)),
                 c("x", "y", "formula", "data", "probability", "seed", "...")),
    xgboost = setdiff(union(names(formals(xgboost::xgb.train)),
                            names(formals(xgboost::xgb.params))),
                      c("params", "data", "evals", "objective", "xgb_model",
                        "seed", "...")),
    NULL)
  if (is.null(allowed)) return(invisible(TRUE))
  .cf_params_reject(setdiff(names(plist), allowed), allowed,
                    sprintf("%s (learner '%s')", arg, learner))
  invisible(TRUE)
}

#' Drop NULL entries so a fitter's own defaults survive
#' @keywords internal
#' @noRd
.cf_drop_nulls <- function(x) x[!vapply(x, is.null, logical(1))]

#' Dispatch one fit to its low-level `cf_fit_<method>()`
#'
#' Single entry point for both the in-sample path and the per-fold / full-data
#' path, so the two cannot drift apart. `control` holds the arguments CATEomics
#' supplies itself; each method receives only the ones its signature declares
#' (grf has no `binary_outcome`, the meta-learners have no `device`, and so on).
#'
#' @keywords internal
#' @noRd
.cf_dispatch_fit <- function(method, X, T, Y, cov_names,
                             control = list(), params = list()) {
  fitter <- paste0("cf_fit_", method)
  if (!exists(fitter, mode = "function")) {
    rlang::abort(sprintf("Unknown method: %s", method),
                 class = "cateomics_unknown_method")
  }
  fun   <- get(fitter, mode = "function")
  known <- names(formals(fun))
  # bcf translates a known `propensity_score` into its native `pihat` itself.
  # A NULL seed would leave the Python backends unseeded, so it becomes 0
  # here exactly as the previous per-method dispatch did.
  control$random_state <- as.integer(control$random_state %||% 0L)
  base <- c(list(X = X, T = T, Y = Y, covariate_names = cov_names),
            .cf_drop_nulls(control))
  # Arguments the caller set explicitly but this backend has no formal for
  # would otherwise be filtered out without a word. bcf is the one exception:
  # it takes a known propensity through its native `pihat`.
  .cf_warn_unreachable_control(method, names(base), known)
  base <- base[intersect(names(base), known)]
  .cf_apply_params(fun, base, params)
}

#' Take the eval-nuisance block that applies to one slot
#'
#' `eval_params` is either `list(outcome = , propensity = )` or a single block
#' that applies to both models. Returns that block unvalidated, so the single
#' learner and the Super Learner can each check it against what they fit.
#'
#' @keywords internal
#' @noRd
.cf_eval_params_block <- function(eval_params, slot) {
  eval_params <- .cf_params_shape(eval_params, "eval_params")
  if (length(eval_params) == 0L) return(list())
  keyed <- intersect(names(eval_params), c("outcome", "propensity"))
  if (length(keyed) && length(setdiff(names(eval_params), keyed))) {
    rlang::abort(
      c("`eval_params` mixes slot blocks with flat entries.",
        i = "Use either list(outcome = , propensity = ) or a single flat block."),
      class = "cateomics_invalid_params"
    )
  }
  if (length(keyed)) {
    .cf_params_shape(eval_params[[slot]] %||% list(),
                     paste0("eval_params$", slot))
  } else {
    eval_params
  }
}

#' Pick the eval-nuisance hyperparameters for one slot
#'
#' Single-learner form: the whole slot block goes to that one learner and is
#' validated against it.
#'
#' @keywords internal
#' @noRd
.cf_eval_params_slot <- function(eval_params, slot, learner, type = NULL) {
  out <- .cf_eval_params_block(eval_params, slot)
  if (!length(out)) return(list())
  if (identical(learner, "known")) return(list())
  # The Super Learner form keys by member, so a learner name here is the
  # mirror of that mistake. Caught before the allow-list reports the label as
  # an unknown hyperparameter of the learner the caller did name.
  named_learner <- intersect(names(out),
                             c(names(.CF_SL_MEMBERS), .CF_EVAL_LEARNERS))
  if (length(named_learner)) {
    rlang::abort(
      c(sprintf("`eval_params$%s` contains a learner block: %s.",
                slot, paste(named_learner, collapse = ", ")),
        i = "With a single learner the hyperparameters go in directly, without a learner block.",
        i = sprintf("Set `eval_nuisance = \"%s\"` to key them by library member instead.",
                    .CF_SL_TOKEN)),
      class = "cateomics_eval_params_learner_block"
    )
  }
  .cf_learner_params_check(learner, out, sprintf("eval_params$%s", slot),
                           type = type)
  out
}

#' Split the eval-nuisance hyperparameters across Super Learner members
#'
#' Every key names a library member. A bare hyperparameter is refused rather
#' than applied to the members that happen to accept it: the seven members
#' share almost no argument name (`max_depth` belongs to xgboost alone,
#' `nfolds` to glmnet, `num.trees` to ranger), so a shared entry would abort on
#' whichever member rejects it first and the message would name a member the
#' caller never mentioned. Member labels are looked up in a closed set that
#' shares no name with any learner's arguments, so no key has to be guessed at.
#'
#' @param members Character vector of library member labels.
#' @param type Resolved outcome type for this slot (the propensity slot is
#'   always a classification).
#' @return Named list, one entry per member, each a (possibly empty) parameter
#'   list.
#' @keywords internal
#' @noRd
.cf_eval_params_members <- function(eval_params, slot, members, type = NULL) {
  arg <- sprintf("eval_params$%s", slot)
  empty <- stats::setNames(rep(list(list()), length(members)), members)
  block <- .cf_eval_params_block(eval_params, slot)
  if (!length(block)) return(empty)

  stray <- setdiff(names(block), members)
  if (length(stray)) {
    off_library <- intersect(stray, names(.CF_SL_MEMBERS))
    # Name a member that actually takes settings: `mean` heads the default
    # library and takes none, so an example built from the first entry would
    # itself be refused.
    tunable <- Filter(function(m) !identical(.CF_SL_MEMBERS[[m]]$learner, "mean"),
                      members)
    example <- if (length(tunable)) tunable[[1L]] else members[[1L]]
    rlang::abort(
      c(sprintf("`eval_nuisance` is \"%s\", so `%s` must name a library member: %s.",
                .CF_SL_TOKEN, arg, paste(stray, collapse = ", ")),
        i = sprintf("Members in this library: %s.", paste(members, collapse = ", ")),
        i = if (length(off_library)) {
          sprintf("%s names a member that this `eval_library` leaves out.",
                  paste(off_library, collapse = ", "))
        } else {
          sprintf("For example `%s = list(%s = list(...))`, or name a single learner in `eval_nuisance` to tune it directly.",
                  arg, example)
        }),
      class = "cateomics_sl_params_need_member"
    )
  }

  out <- empty
  for (m in intersect(members, names(block))) {
    p <- .cf_params_shape(block[[m]], paste0(arg, "$", m))
    if (!length(p)) next
    spec  <- .CF_SL_MEMBERS[[m]]
    where <- sprintf("%s$%s", arg, m)
    # The label is the name of the settings it carries, so taking them again
    # here would give one knob two names.
    clash <- intersect(names(p), spec$fixes)
    if (length(clash)) {
      rlang::abort(
        c(sprintf("Member '%s' already fixes %s.", m, paste(clash, collapse = ", ")),
          i = sprintf("Name a member whose label carries the value you want, or set `eval_nuisance = \"%s\"` to tune one learner.",
                      spec$learner)),
        class = "cateomics_duplicate_param"
      )
    }
    if (identical(spec$learner, "mean")) {
      rlang::abort(
        c(sprintf("Member 'mean' takes no hyperparameters (in `%s`): %s.",
                  where, paste(names(p), collapse = ", ")),
          i = "It predicts the mean of the training rows."),
        class = "cateomics_unknown_param"
      )
    }
    .cf_learner_params_check(spec$learner, p, where, type = type)
    out[[m]] <- p
  }
  out
}

#' Match several values, refusing the ones that match nothing.
#'
#' `match.arg(several.ok = TRUE)` drops a value that matches no choice as long
#' as one other value matches, so a misspelling in a vector is silently
#' ignored while the same misspelling alone is an error. This keeps the
#' partial matching and makes an unmatched value an error either way.
#'
#' @keywords internal
#' @noRd
.cf_match_several <- function(x, choices, arg, fn) {
  if (is.null(x) || !length(x)) {
    rlang::abort(sprintf("%s: `%s` must name at least one of: %s.",
                         fn, arg, paste(choices, collapse = ", ")),
                 class = "cateomics_invalid_input")
  }
  hit <- pmatch(as.character(x), choices, duplicates.ok = TRUE)
  bad <- as.character(x)[is.na(hit)]
  if (length(bad) > 0L) {
    rlang::abort(
      c(sprintf("%s: `%s` does not accept %s.", fn, arg,
                paste(sprintf("'%s'", bad), collapse = ", ")),
        i = sprintf("Choose from: %s.", paste(choices, collapse = ", "))),
      class = "cateomics_invalid_input")
  }
  unique(choices[hit])
}


#' Reject named arguments that reached a `...` meant to stay empty
#'
#' Functions that keep `...` for S3 reasons would otherwise swallow a
#' misspelled argument and run with the default instead.
#'
#' @param ... The dots to inspect.
#' @param .fn The calling function's name, for the message.
#' @return `TRUE` invisibly when the dots are empty; otherwise an error.
#' @keywords internal
#' @noRd
.cf_reject_dots <- function(..., .fn) {
  dots <- list(...)
  nm <- names(dots)
  nm <- if (is.null(nm)) character(0) else nm[nzchar(nm)]
  if (!length(nm)) return(invisible(TRUE))
  # Names a reader is likely to try, with what the package calls the same
  # thing. Phrased as the present state, not as a change: nobody has used an
  # earlier spelling of these.
  known <- c(min_fidelity = "not an argument",
             consensus_q = "not an argument",
             benefit_top_q = "not an argument",
             outlier_mad = "not an argument",
             honest = "not an argument: a fit that falls back to in-sample tau-hat says so itself",
             stability_filter = "not an argument: pass `consensus_min_agreement` to switch the filter on",
             stability_metric = "not an argument",
             stability_min_agreement = "spelled `consensus_min_agreement`",
             assay_name = "part of `omics`",
             experiments = "part of `omics`",
             assay_cols = "part of `omics`",
             assay_pattern = "not an argument: name the columns in `omics`",
             trim = "spelled `clip`")
  hits <- intersect(nm, names(known))
  rlang::abort(
    c(sprintf("%s does not accept: %s.", .fn, paste(nm, collapse = ", ")),
      if (length(hits))
        stats::setNames(sprintf("`%s` is %s.", hits, known[hits]),
                        rep("i", length(hits)))),
    class = "cateomics_unknown_argument"
  )
}


#' Warn when a set control argument cannot reach the chosen backend
#'
#' `.cf_dispatch_fit()` keeps only the arguments the per-method fitter
#' declares. For most of the control block that is harmless plumbing, but a
#' user who sets a known propensity or a base learner has asked for something
#' specific; dropping it silently leaves them reading a default fit as if it
#' were theirs.
#'
#' @param method Method name.
#' @param supplied Names present in the assembled call (NULLs already dropped).
#' @param known Formal names of the per-method fitter.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.cf_warn_unreachable_control <- function(method, supplied, known) {
  watched <- .cf_watched_control_args(method)
  lost <- setdiff(intersect(supplied, watched), known)
  if (!length(lost)) return(invisible(NULL))
  # The dispatcher runs once per fold plus once for the full fit. The repeats
  # are collapsed where the fold loop is driven, so this says it plainly.
  rlang::warn(
    c(sprintf("Method '%s' has no %s; the value(s) are not used.",
              method, paste0("`", lost, "`", collapse = " / ")),
      .cf_control_arg_hints()),
    class = "cateomics_control_arg_unreachable"
  )
  invisible(NULL)
}

#' Control arguments worth reporting when a backend cannot take them
#'
#' Every name here has a NULL default at the public entry and is present in
#' the assembled call only when the caller set it. An argument whose default
#' is not NULL, or one the entry reshapes into an empty list, would be present
#' always and would report on every fit -- which is why `propensity_params` is
#' not here: it becomes `list()` at the entry, and it has its own abort for
#' the methods that cannot use it.
#'
#' @param method Method name, or NULL for the whole set.
#' @return Character vector of argument names.
#' @keywords internal
#' @noRd
.cf_watched_control_args <- function(method = NULL) {
  watched <- c("propensity_score", "model_outcome", "model_final",
               "model_propensity", "clip")
  # bcf translates a known propensity into its native `pihat` itself.
  if (!is.null(method) && identical(method, "bcf"))
    watched <- setdiff(watched, "propensity_score")
  watched
}

#' @keywords internal
#' @noRd
.cf_control_arg_hints <- function() {
  c(i = "Methods taking a propensity: grf, bcf, xlearner, drlearner.",
    i = "Methods taking a propensity learner: bcf, xlearner, drlearner.",
    i = "Methods taking an outcome learner: slearner, tlearner, xlearner, drlearner.",
    i = "Methods taking a final learner: xlearner, drlearner.",
    i = "Methods bounding a propensity with `clip`: grf, bcf, xlearner, drlearner.")
}

#' Report a control argument no compared method can take
#'
#' `cf_fit` reports this from the dispatcher. A comparison legitimately mixes
#' methods that differ in what they take, so the honest unit of reporting is
#' the whole set: a value that reaches at least one compared method was used,
#' and only one that reaches none of them was asked for and ignored. That is
#' the same shape `propensity_params` already uses at this entry. The
#' per-candidate reports the dispatcher still makes are collected and kept in
#' `metadata$config_warnings`; the line this raises is not repeated there.
#'
#' @param methods Character vector of the compared method names.
#' @param supplied Names of control arguments the caller actually set.
#' @return The condition class it signalled, invisibly, or an empty character
#'   vector when it had nothing to report.
#' @keywords internal
#' @noRd
.cf_check_unreachable_control <- function(methods, supplied) {
  supplied <- intersect(supplied, .cf_watched_control_args())
  if (!length(supplied)) return(invisible(NULL))
  reaches <- function(nm) {
    any(vapply(methods, function(m) {
      if (!nm %in% .cf_watched_control_args(m)) return(TRUE)  # bcf x propensity
      fun <- .cf_method_fitter(m)
      !is.null(fun) && nm %in% names(formals(fun))
    }, logical(1)))
  }
  lost <- supplied[!vapply(supplied, reaches, logical(1))]
  if (!length(lost)) return(invisible(character(0)))
  rlang::warn(
    c(sprintf("No compared method takes %s; the value(s) are not used.",
              paste0("`", lost, "`", collapse = " / ")),
      i = sprintf("Compared: %s.", paste(methods, collapse = ", ")),
      .cf_control_arg_hints()),
    class = "cateomics_control_arg_unreachable"
  )
  # The class it reported, so a caller that collects the same condition from
  # inside the fits can keep the record without saying the fact twice.
  invisible("cateomics_control_arg_unreachable")
}
