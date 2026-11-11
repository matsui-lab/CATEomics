# utils_bootstrap_ci.R
#
# Where `bootstrap_ci` is accepted, warned about or refused.
#
# This file holds that policy and nothing else. The resampling itself lives
# inside each fitter's own worker (the metalearners and the DR learner are the
# ones whose formals take `bootstrap_ci`), which produces `res$ate_ci` and
# `res$tau_ci` by the percentile method.



#' The interval a method reports when `bootstrap_ci = FALSE`.
#'
#' One table for the whole package: `cf_fit` and `cf_compare` both decide the
#' CI policy before any fit work, and a second copy of this map would go stale
#' the moment a backend gains or loses a native interval.
#'
#' @keywords internal
#' @noRd
.cf_ci_native_label <- function(method, model_final = NULL) {
  labels <- list(
    grf   = "asymptotic Wald",
    bcf   = "Bayesian posterior",
    cevae = "Bayesian posterior",
    # econml's StatsModelsInference needs a linear final stage; with any other
    # model_final the DR-learner has no native CI to prefer.
    drlearner = if (identical(model_final, "linear"))
                  "asymptotic Wald" else "unavailable"
  )
  labels[[method]] %||% "unavailable"
}

#' Whether a method's own fit can resample.
#'
#' Three backends are named outright, because the cost of resampling them is
#' the reason they are refused rather than anything in their signature. Every
#' other method is read from its fitter's formals, so that half stays in step
#' with the fitters.
#'
#' @keywords internal
#' @noRd
.cf_can_bootstrap <- function(method) {
  if (method %in% c("tarnet", "dragonnet", "ganite")) return(FALSE)
  fitter <- .cf_method_fitter(method)
  is.null(fitter) || "bootstrap_ci" %in% names(formals(fitter))
}

#' Refuse a bootstrap request a backend cannot serve.
#'
#' Abort-only; `.cf_check_ci_policy()` adds the native-CI override warning
#' for the methods whose request is served.
#'
#' @keywords internal
#' @noRd
.cf_check_bootstrap_supported <- function(method, native_label) {
  nn_methods <- c("tarnet", "dragonnet", "ganite")
  if (method %in% nn_methods) {
    rlang::abort(
      sprintf(
        "Bootstrap CI is not supported for method '%s' (prohibitive compute cost for this backend).",
        method
      ),
      class = "cateomics_bootstrap_unsupported_nn"
    )
  }
  # Every other method is asked of its own fitter: a method whose fitter has no
  # `bootstrap_ci` formal cannot resample, so the request is refused here
  # instead of being dropped on the way to the backend. Reading the formals
  # keeps this in step with the fitters (a hard-coded method list would go
  # stale the moment one of them gains or loses the argument).
  fitter <- .cf_method_fitter(method)
  if (!is.null(fitter) && !("bootstrap_ci" %in% names(formals(fitter)))) {
    rlang::abort(
      c(
        sprintf("`bootstrap_ci = TRUE` is not supported for method '%s'.",
                method),
        i = sprintf("This backend does not resample; its interval is the %s.",
                    native_label),
        i = "Drop `bootstrap_ci` or choose a method that resamples."
      ),
      class = "cateomics_bootstrap_unsupported"
    )
  }
  invisible(NULL)
}

#' Decide the CI policy at an entry point, before any fit work.
#'
#' Called with one method from `cf_fit` and with the whole candidate set from
#' `cf_compare`, so a request a backend cannot serve is refused at the same
#' point by both entries. Taking the set in one call is what makes the partial
#' case decidable: whether a request can be served "where it can be" is a
#' property of the set, not of one candidate.
#'
#' `model_final` decides whether a method has a native interval to override, and
#' a grid sets it per candidate, so it is accepted either as one value for every
#' method or as a list parallel to `methods`. `labels` names the candidates in
#' the report when they are not simply the method names.
#'
#' @keywords internal
#' @noRd
.cf_check_ci_policy <- function(methods, bootstrap_ci, model_final = NULL,
                                allow_partial = FALSE, labels = NULL) {
  if (!isTRUE(bootstrap_ci)) return(invisible(NULL))
  methods <- as.character(methods)
  if (is.null(labels)) labels <- methods
  mf <- if (is.list(model_final)) model_final else
    rep(list(model_final), length(methods))
  can <- vapply(methods, .cf_can_bootstrap, logical(1))
  # One method is asked, or none of several can resample: the request cannot be
  # served at all, so it is refused. Several methods of which some can: the
  # request is served where it can be, and the rest are named -- once, with
  # every method that the request does not reach.
  keep <- rep(TRUE, length(methods))
  if (allow_partial && any(can) && !all(can)) {
    out <- unique(methods[!can])
    rlang::warn(
      c(sprintf("`bootstrap_ci = TRUE` does not reach %d of %d method(s): %s.",
                length(out), length(unique(methods)), paste(out, collapse = ", ")),
        i = "Those report the interval their own backend provides.",
        i = "The rest resample as asked."),
      class = "cateomics_bootstrap_partial"
    )
    keep <- can
  }
  # Abort for a request that is not served at all.
  for (i in which(keep)) {
    .cf_check_bootstrap_supported(methods[i],
                                  .cf_ci_native_label(methods[i], mf[[i]]))
  }
  # The request is served and replaces an interval the backend would have
  # reported itself. One line for the call, naming each candidate it applies to:
  # a grid asks the same question of every candidate, and the answer differs
  # only by which interval is being replaced.
  ovr <- character(0)
  for (i in which(keep)) {
    label <- .cf_ci_native_label(methods[i], mf[[i]])
    if (!identical(label, "unavailable")) {
      ovr <- c(ovr, sprintf("%s (%s)", labels[i], label))
    }
  }
  if (length(ovr) > 0L) {
    rlang::warn(
      c(sprintf("A nonparametric bootstrap replaces the native CI of %d candidate(s).",
                length(ovr)),
        i = paste(unique(ovr), collapse = "; ")),
      class = "cateomics_bootstrap_overrides_native"
    )
  }
  invisible(NULL)
}

#' Say so when an interval does not contain the estimate it is reported with.
#'
#' A percentile interval is built from where the refitted replicates landed,
#' not around the full-sample estimate, so a biased estimator can put the
#' estimate outside its own interval. That reads like a mistake, and silence
#' would leave a reader to guess, so the fit says which way the replicates
#' went. Checked on the result the caller receives, at each entry, rather than
#' inside the dispatch: with `partition = "oof"` the dispatch also runs once
#' per fold, and those intermediate fits are not what anyone reads.
#'
#' @param res A fitter's return value.
#' @return `res`, invisibly.
#' @keywords internal
#' @noRd
.cf_warn_ci_excludes_estimate <- function(res) {
  if (!is.list(res)) return(invisible(res))
  # A fitter returns either the flat fields or a `res` bucket holding them,
  # depending on which assembler it goes through; both reach this one check.
  src <- if (is.list(res$res)) res$res else res
  ate <- src$ate
  ci  <- src$ate_ci
  if (is.null(ate) || is.null(ci) || length(ate) != 1L || length(ci) != 2L ||
      !all(is.finite(c(ate, ci)))) {
    return(invisible(res))
  }
  if (ate >= ci[1] && ate <= ci[2]) return(invisible(res))
  side <- if (ate > ci[2]) "below" else "above"
  resampled <- is.character(src$ci_method) &&
    any(grepl("bootstrap", src$ci_method, fixed = TRUE))
  why <- if (resampled) {
    sprintf(paste("The refitted replicates sit systematically %s the",
                  "full-sample fit, which happens when a penalised learner",
                  "is refitted on a resample holding fewer distinct rows."),
            side)
  } else {
    sprintf("The interval lies entirely %s the estimate.",
            if (identical(side, "below")) "below" else "above")
  }
  rlang::warn(
    c("The interval does not contain the point estimate reported with it.",
      i = sprintf("ATE %.4g, interval [%.4g, %.4g].", ate, ci[1], ci[2]),
      i = why,
      i = paste("Both are reported as computed; an interval taken from the",
                "spread of the replicates is not built around the estimate.")),
    class = "cateomics_bootstrap_ci_excludes_estimate")
  invisible(res)
}
