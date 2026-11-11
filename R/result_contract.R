#' @title Result Contract: cf_model Object Structure
#'
#' @description
#' This file documents the formal structure (contract) for \code{cf_model}
#' objects returned by \code{\link{cf_fit}}. All causal inference methods in
#' CATEomics MUST return objects conforming to this structure.
#'
#' \strong{Four-pillar layout.} A \code{cf_model} is a flat list with four
#' pillars plus \code{method}:
#'
#' \describe{
#'   \item{method}{Character string identifying the method (e.g. "grf").}
#'   \item{res}{Numeric results (see Res Structure).}
#'   \item{model}{Method-specific fitted backend container (structure varies
#'     by method).}
#'   \item{meta}{Metadata (see Meta Structure).}
#'   \item{data}{List with \code{X_preprocessed}, \code{Y}, \code{T} (the
#'     model-visible training data; used by cf_interpret / cf_compare /
#'     cf_evaluate / cf_rule for nuisance recomputation).}
#'   \item{fold_cache}{Optional: K per-fold fits when
#'     \code{keep_kfold_models = TRUE}.}
#' }
#'
#' \strong{The ensemble is a different shape.} \code{\link{cf_ensemble}}
#' returns class \code{c("cateomics_ensemble", "cf_model")} and departs from the
#' contract above in three ways: \code{res$tau_hat} is \code{NULL} (it combines
#' per-method out-of-fold estimates, so only \code{res$tau_hat_oof} carries
#' values), \code{model} holds the members that carry non-zero weight instead
#' of a single backend, and a fifth element \code{weights} holds the
#' combination weights.
#' \code{meta} is built by the same assembler a \code{cf_fit} result uses, so
#' it carries \code{provenance}, \code{estimator} and \code{named_args}; it
#' has no \code{spec}, since the members hold their own.
#' \code{validate_cateomics_result()} is called by \code{cf_fit} and does not
#' run on it.
#'
#' @section Res Structure:
#'
#' The \code{res} pillar MUST include:
#'
#' \describe{
#'   \item{tau_hat}{Numeric vector of Conditional Average Treatment Effects
#'     (individual \eqn{\hat\tau_i}), length \code{n}.}
#'   \item{ate}{Numeric scalar, Average Treatment Effect.}
#'   \item{mu0_hat, mu1_hat}{Numeric vectors of \emph{modelled} counterfactual
#'     outcomes under \eqn{do(T=0)} / \eqn{do(T=1)}, length \code{n}.}
#'   \item{summary}{List of summary statistics (see Summary Structure).}
#' }
#'
#' Optional \code{res} fields (populated when available):
#' \code{tau_hat_oof} (K-fold OOF \eqn{\hat\tau}), \code{fold_id} (integer
#' vector of length \code{n} giving the fold that held each row out; present
#' with \code{tau_hat_oof} and \code{NULL} under \code{partition = "none"}),
#' \code{ate_ci} (length-2
#' numeric), \code{tau_ci} (data.frame \code{lower}/\code{upper}),
#' \code{ci_method} (character label); both intervals are computed for the
#' full-data fit, so they go with \code{ate} and \code{tau_hat}, and no
#' interval is computed for \code{tau_hat_oof}, \code{nuisance} (list
#' \code{mu0_oof}/\code{mu1_oof}/\code{pi_oof}, plus \code{m_oof}, the
#' pooled regression of Y on X; the nuisance has these four elements under
#' either \code{partition}, out of fold under \code{"oof"} and in sample
#' under \code{"none"} -- the \code{*_oof} names are kept in both cases so
#' that downstream functions read one shape),
#' \code{ps} with
#' \code{ps_type} (the propensity the fit used and where it comes from), and
#' method-specific numeric extras.
#'
#' @section Summary Structure:
#'
#' The \code{res$summary} list MUST contain (IQR form):
#'
#' \describe{
#'   \item{ate}{Numeric, point estimate of ATE.}
#'   \item{ate_ci_lower, ate_ci_upper}{Numeric CI bounds (\code{NA_real_}
#'     when no CI is available).}
#'   \item{tau_mean, tau_sd}{Mean / SD of \eqn{\hat\tau}.}
#'   \item{tau_iqr}{List with \code{q0}, \code{q25}, \code{q50}, \code{q75},
#'     \code{q100} quantiles of \eqn{\hat\tau} (min / Q1 / median / Q3 / max).}
#' }
#'
#' @section Meta Structure:
#'
#' The \code{meta} pillar contains identification fields (\code{n},
#' \code{p}, \code{covariate_names}; the variables the caller named are in
#' \code{named_args}, and the method name is the top-level \code{method})
#' plus:
#'
#' \describe{
#'   \item{provenance}{List with \code{backend}, \code{rng_seed},
#'     \code{software_versions}, \code{device_requested},
#'     \code{device_effective} and \code{sample_names_meaningful}.}
#'   \item{estimator}{List with \code{bootstrap}, \code{preprocess},
#'     \code{outcome_type},
#'     \code{eval_nuisance_resolved} (the learner each nuisance slot ran,
#'     which is \code{"known"} for a propensity that was supplied; the
#'     learner that was asked for is recorded once, on
#'     \code{cf_compare}'s \code{metadata}),
#'     \code{eval_library},
#'     \code{eval_sl_weights}, (when known) \code{ci_method}, and
#'     \code{training_trace} for the backends that report their losses
#'     (\code{ganite}, \code{cevae}).}
#' }
#'
#' Preprocessing slots consumed by \code{\link{predict.cf_model}} on newdata
#' (\code{x_zscore}, \code{y_zscore}, \code{factor_levels},
#' \code{input_type}, \code{named_args}, \code{propensity_X}) are retained at
#' \code{meta} top level, as is \code{spec} (the resolved \code{method},
#' \code{clip}, \code{model_*} and the parameter blocks), which
#' \code{print} reads. The PCA rotation sits in the layer it was fitted on,
#' and \code{partition} with the other estimator settings.
#'
#' The same assembler builds this pillar for a \code{cf_fit} result, for each
#' configuration of \code{\link{cf_compare}} and for
#' \code{\link{cf_ensemble}} (which adds the ensemble method, the loss, the
#' members used and dropped, and the consensus members). On an ensemble the
#' replay slots \code{x_zscore}, \code{y_zscore}, \code{factor_levels},
#' \code{input_type} and \code{propensity_X} are \emph{absent} rather than
#' empty, and so is \code{spec}: each member carries its own preprocessing
#' record and its own resolved settings, and
#' \code{predict(object, newdata = )} delegates to the members, so the
#' ensemble has nothing of its own to record there. \code{named_args} is
#' kept, inherited from the comparison, so the object can still say what it
#' was fitted on.
#'
#' @section S3 Methods:
#'
#' \describe{
#'   \item{\code{predict(object, type, ...)}}{Extract predictions (all, ate,
#'     ate_ci, tau, tau_ci, mu0, mu1, ps, summary).}
#'   \item{\code{summary(object, ...)}}{Display summary with CIs.}
#'   \item{\code{print(object, ...)}}{Print overview of results.}
#' }
#'
#' The package defines no \code{plot} method. \code{\link{cf_plot_cate}()} is
#' a named function that draws the \eqn{\hat\tau} distribution.
#'
#' @section Class Hierarchy:
#'
#' Objects carry class \code{"cf_model"}.
#'
#' @name cateomics_result
#' @aliases cf_model
#' @return This is a documentation-only page describing the structure of the
#'   \code{cateomics_result} (\code{cf_model}) object returned by
#'   \code{\link{cf_fit}}; it does not define a function and has no return value.
#' @seealso \code{\link{cf_fit}}, \code{\link{predict.cf_model}},
#'   \code{\link{summary.cf_model}}
NULL

#' Validate cateomics_result Structure
#'
#' Internal function to validate that a result object conforms to the
#' four-pillar cateomics_result contract.
#'
#' @param x Object to validate
#' @return \code{invisible(TRUE)} if valid; aborts with class
#'   \code{cateomics_invalid_result} otherwise
#' @keywords internal
validate_cateomics_result <- function(x) {
  if (!inherits(x, "cf_model")) {
    rlang::abort(
      "Object must have class 'cf_model'",
      class = "cateomics_invalid_result"
    )
  }

  # Four-pillar top level. `method` + `res` + `meta` are mandatory; `model` /
  # `data` are not validated -- cf_fit always sets them, but the internal
  # K-fold wrappers build lighter objects.
  required_toplevel <- c("method", "res", "meta")
  missing_toplevel <- setdiff(required_toplevel, names(x))
  if (length(missing_toplevel) > 0) {
    rlang::abort(
      sprintf("Missing required components: %s", paste(missing_toplevel, collapse = ", ")),
      class = "cateomics_invalid_result"
    )
  }

  if (!is.character(x$method) || length(x$method) != 1) {
    rlang::abort(
      "'method' must be a single character string",
      class = "cateomics_invalid_result"
    )
  }

  # meta is a single flat metadata pillar with `provenance` / `estimator`
  # buckets (diagnostics are on-demand via cf_diagnostics).
  # `named_args` is in the required set because the newdata path reads it to
  # know which columns were the outcome, the treatment and the baseline
  # covariates; a result without it cannot be scored on new samples.
  required_meta <- c("n", "p", "covariate_names", "named_args",
                     "provenance", "estimator")
  missing_meta <- setdiff(required_meta, names(x$meta))
  if (length(missing_meta) > 0) {
    rlang::abort(
      sprintf("Missing required meta fields: %s", paste(missing_meta, collapse = ", ")),
      class = "cateomics_invalid_result"
    )
  }

  # Results live in the flat `res` pillar.
  res <- x$res
  if (is.null(res) || !is.list(res)) {
    rlang::abort(
      "`res` is NULL or not a list; all CATEomics methods must return a res pillar.",
      class = "cateomics_invalid_result"
    )
  }

  required_res <- c("tau_hat", "ate", "mu0_hat", "mu1_hat", "summary")
  missing_res <- setdiff(required_res, names(res))
  if (length(missing_res) > 0) {
    rlang::abort(
      sprintf("Missing required res fields: %s",
              paste(missing_res, collapse = ", ")),
      class = "cateomics_invalid_result"
    )
  }

  n <- x$meta$n
  .is_num_or_all_na <- function(v, label, expected_len = NULL) {
    if (!(is.numeric(v) || all(is.na(v)))) {
      rlang::abort(
        sprintf("`res$%s` must be numeric or all NA.", label),
        class = "cateomics_invalid_result"
      )
    }
    if (!is.null(expected_len) && length(v) != expected_len) {
      rlang::abort(
        sprintf("`res$%s` must have length n = %d (got %d).",
                label, expected_len, length(v)),
        class = "cateomics_invalid_result"
      )
    }
  }
  .is_num_or_all_na(res$tau_hat, "tau_hat", expected_len = n)
  .is_num_or_all_na(res$mu0_hat, "mu0_hat", expected_len = n)
  .is_num_or_all_na(res$mu1_hat, "mu1_hat", expected_len = n)

  if (!((is.numeric(res$ate) && length(res$ate) == 1L) || is.na(res$ate))) {
    rlang::abort(
      "`res$ate` must be a single numeric value or NA.",
      class = "cateomics_invalid_result"
    )
  }

  # Every method estimates a conditional effect, so a result with no per-unit
  # estimate at all is not one this contract can describe.
  if (all(is.na(res$tau_hat))) {
    rlang::abort(
      sprintf("Method '%s' returned all-NA tau_hat.", x$method),
      class = "cateomics_invalid_result"
    )
  }

  required_summary <- c(
    "ate", "ate_ci_lower", "ate_ci_upper",
    "tau_mean", "tau_sd", "tau_iqr"
  )
  missing_summary <- setdiff(required_summary, names(res$summary))
  if (length(missing_summary) > 0) {
    rlang::abort(
      sprintf("Missing required summary fields: %s",
              paste(missing_summary, collapse = ", ")),
      class = "cateomics_invalid_result"
    )
  }
  if (!is.list(res$summary$tau_iqr) ||
      !all(c("q0", "q25", "q50", "q75", "q100") %in% names(res$summary$tau_iqr))) {
    rlang::abort(
      "`res$summary$tau_iqr` must be a list with q0/q25/q50/q75/q100.",
      class = "cateomics_invalid_result"
    )
  }

  # ate_ci / tau_ci unified location (both optional).
  if (!is.null(res$ate_ci)) {
    if (!is.numeric(res$ate_ci) || length(res$ate_ci) != 2L ||
        anyNA(res$ate_ci) || res$ate_ci[1] > res$ate_ci[2]) {
      rlang::abort(
        "`res$ate_ci` must be numeric length-2 with ate_ci[1] <= ate_ci[2].",
        class = "cateomics_invalid_result"
      )
    }
  }
  if (!is.null(res$tau_ci)) {
    if (!is.data.frame(res$tau_ci) ||
        !all(c("lower", "upper") %in% colnames(res$tau_ci)) ||
        nrow(res$tau_ci) != length(res$tau_hat)) {
      rlang::abort(
        "`res$tau_ci` must be a data.frame with columns 'lower'/'upper' and nrow == length(tau_hat).",
        class = "cateomics_invalid_result"
      )
    }
  }

  invisible(TRUE)
}

#' Build the standard res$summary list (IQR form)
#'
#' Centralised so every method + cf_fit produce an identical summary shape.
#' Called by cf_fit after lifting the per-method res into the flat pillar.
#'
#' @param ate Numeric scalar ATE.
#' @param ate_ci Length-2 numeric \code{c(lower, upper)} or NULL.
#' @param tau_hat Numeric vector of \eqn{\hat\tau}.
#' @return A list with ate / ate_ci_lower / ate_ci_upper / tau_mean / tau_sd /
#'   tau_iqr(q0,q25,q50,q75,q100).
#' @keywords internal
#' @noRd
.cf_build_summary <- function(ate, ate_ci = NULL, tau_hat = NA_real_) {
  ci_lower <- if (!is.null(ate_ci) && length(ate_ci) == 2L) ate_ci[1] else NA_real_
  ci_upper <- if (!is.null(ate_ci) && length(ate_ci) == 2L) ate_ci[2] else NA_real_

  has_tau <- !all(is.na(tau_hat))
  if (has_tau) {
    tv <- as.numeric(tau_hat)
    qs <- stats::quantile(tv, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
    tau_iqr <- list(q0 = qs[1], q25 = qs[2], q50 = qs[3], q75 = qs[4], q100 = qs[5])
    tau_mean <- mean(tv)
    tau_sd   <- stats::sd(tv)
  } else {
    tau_iqr <- list(q0 = NA_real_, q25 = NA_real_, q50 = NA_real_,
                    q75 = NA_real_, q100 = NA_real_)
    tau_mean <- NA_real_
    tau_sd   <- NA_real_
  }


  list(
    ate          = as.numeric(ate),
    ate_ci_lower = as.numeric(ci_lower),
    ate_ci_upper = as.numeric(ci_upper),
    tau_mean     = tau_mean,
    tau_sd       = tau_sd,
    tau_iqr      = tau_iqr
  )
}

#' Put `ci_method` where every reader looks for it
#'
#' A method returns the interval's provenance inside `res$metadata`, the small
#' bag of predict-config it also uses. A reader should not have to know that,
#' so it is lifted to `res$ci_method` and removed from the bag. The rest of the
#' bag -- bcf's tree directory, cevae's hyperparameters, drlearner's final
#' learner -- is read by the per-method predictors and stays. A bag that held
#' nothing else is dropped, which is the case for grf; bcf keeps its chain
#' settings there, so its bag survives with five fields.
#'
#' Both entries call this, which is what makes a comparison's candidates and a
#' single fit answer the question from the same field.
#'
#' @param res The `res` pillar as a method returned it.
#' @return The same list, with `ci_method` at its top level.
#' @keywords internal
#' @noRd
.cf_lift_ci_method <- function(res) {
  if (is.null(res$metadata)) return(res)
  if (is.null(res$ci_method)) res$ci_method <- res$metadata$ci_method
  res$metadata$ci_method <- NULL
  if (length(res$metadata) == 0L) res$metadata <- NULL
  res
}

#' Metadata for one cf_compare configuration.
#'
#' cf_compare does not call `cf_fit`, so its per-config models are assembled
#' here, through the same assembler. That is what lets
#' `predict(newdata = <container>)` rebuild the design matrix from a
#' comparison row: the record of the input and of the preprocessing is the
#' one a `cf_fit` result carries. The
#' same assembler as `cf_fit` is used here instead, with the preprocessing
#' record the shared Step 1-5 produced.
#'
#' @param prep Return value of `.cf_preprocess_input()`.
#' @param cov_names Covariate names of the preprocessed design matrix.
#' @param n Sample count.
#' @param method Method name of this configuration.
#' @param res The `res` pillar of this configuration (read for `y_zscore`).
#' @param random_state Seed recorded in the provenance.
#' @param bootstrap_ci Whether a bootstrap interval was requested.
#' @param propensity_X Propensity covariate subset.
#' @param named_args The names the caller gave.
#' @param spec Resolved hyperparameters of this configuration.
#' @param device,device_resolved Device asked for, and the one resolved for
#'   this method.
#' @param outcome_type Resolved outcome type.
#' @param share Environment holding the large preprocessing members, so the
#'   configurations serialise one copy between them.
#' @param eval_nuisance_resolved,eval_library,eval_sl_weights
#'   The shared evaluation nuisance this comparison scored with.
#' @return A `meta` list with the same contract as a `cf_fit` result, plus the
#'   one field a comparison row carries (`partition`). The learners this
#'   configuration ran with are in `spec`, which records each slot separately --
#'   a single token cannot, since a config may set the outcome and the final
#'   stage to different learners.
#' @keywords internal
#' @noRd
.cf_compare_model_meta <- function(prep, cov_names, n, method, res,
                                   random_state, bootstrap_ci, propensity_X,
                                   named_args, spec,
                                   device = NULL, device_resolved = NULL,
                                   outcome_type = NULL,
                                   eval_nuisance_resolved = NULL,
                                   eval_library = NULL, eval_sl_weights = NULL,
                                   share = NULL) {
  meta <- create_cf_meta(
    n               = n,
    p               = length(cov_names),
    covariate_names = cov_names,
    method          = method,
    bootstrap       = isTRUE(bootstrap_ci),
    random_state    = random_state,
    preprocess      = .cf_preprocess_record(prep, share = share),
    input_type      = prep$input_type,
    factor_levels   = prep$factor_levels,
    x_zscore        = prep$x_zscore,
    # NN backends standardise Y on the data they receive and return the stats.
    y_zscore        = res$y_zscore,
    named_args      = named_args,
    propensity_X    = propensity_X
  )
  # Whether the row names name anything is only knowable where the design was
  # built, so it is carried rather than guessed when values are written back.
  meta$provenance$sample_names_meaningful <-
    isTRUE(prep$sample_names_meaningful)
  meta$provenance <- .cf_record_device(meta$provenance, device, device_resolved,
                                       res$device_effective)
  # The rest of the estimator bucket, which `cf_fit` fills after the assembler
  # returns. Without it a comparison row does not satisfy the contract the
  # result documentation states, and `summary()` reports no interval method.
  meta$estimator$ci_method              <- res$ci_method
  # A backend that logs its losses returns them on `res`; they belong with the
  # other estimator settings, and this entry point has to put them in the same
  # place `cf_fit` does or the same fit would answer differently depending on
  # which function produced it.
  meta$estimator$training_trace         <- res$training_trace
  meta$estimator$outcome_type           <- outcome_type
  meta$estimator$eval_nuisance_resolved <- eval_nuisance_resolved
  meta$estimator$eval_library           <- eval_library
  meta$estimator$eval_sl_weights        <- eval_sl_weights
  # cf_compare always scores out of fold.
  meta$estimator$partition <- "oof"
  meta$spec      <- spec
  meta
}

#' Metadata for a cf_ensemble result.
#'
#' Uses the same assembler as `cf_fit` for the fields every `cf_model` reader
#' expects, then adds the ones only an ensemble has. The slots the newdata path
#' reads on a single fit are absent rather than empty: an ensemble has no
#' preprocessing of its own, and `predict(ens, newdata = )` delegates to the
#' members, each of which carries its own record. Holding them as NULL would
#' read as "nothing was done" instead of "not applicable here".
#'
#' @param compare_result The comparison the members came from.
#' @param cov_names Covariate names.
#' @param n,p Sample and covariate counts.
#' @param extra The ensemble-specific fields (method, loss, weights diagnostics).
#' @return A `meta` list.
#' @keywords internal
#' @noRd
.cf_ensemble_meta <- function(compare_result, cov_names, n, p, extra) {
  member <- compare_result$models[[1L]]$meta
  meta <- create_cf_meta(
    n               = n,
    p               = p,
    covariate_names = cov_names,
    method          = "ensemble",
    bootstrap       = FALSE,
    random_state    = compare_result$metadata$random_state
  )
  # Carried from the comparison that produced the members: writing per-sample
  # values back to a container needs to know whether the row names name
  # anything, and that is only knowable where the design was built.
  meta$provenance$sample_names_meaningful <-
    .cateomics_compare_names_meaningful(compare_result)
  meta$estimator$partition <- "oof"
  # The names the members were fitted under: the same for all of them, since
  # one comparison built one design, so the ensemble can say on its own what it
  # was fitted on.
  meta$named_args <- member$named_args
  # The ensemble reads the comparison's shared nuisance, so the weights that
  # produced it are the ones to report. `cf_evaluate` looks here, and an
  # ensemble is one of its documented inputs.
  meta$estimator$eval_sl_weights <- compare_result$metadata$eval_sl_weights
  meta$estimator$eval_library    <- compare_result$metadata$eval_library
  # The members were fitted on one design, under one nuisance, from one
  # comparison, so the ensemble records the same settings rather than leaving
  # them empty and letting `preprocess = NULL` read as "no preprocessing".
  # `ci_method` is not among them: an ensemble computes no interval.
  meta$estimator$preprocess             <- member$estimator$preprocess
  meta$estimator$outcome_type           <- member$estimator$outcome_type
  meta$estimator$eval_nuisance_resolved <- member$estimator$eval_nuisance_resolved
  meta$estimator$n_folds                <- compare_result$metadata$n_folds
  # An ensemble has no container and no design of its own: it combines the
  # members' out-of-fold estimates, and `predict()` with new samples delegates
  # to them. The slots the newdata path would read are therefore dropped rather
  # than carried empty, which would read as "nothing was done".
  for (nm in c("propensity_X", "input_type", "factor_levels",
               "x_zscore", "y_zscore")) meta[[nm]] <- NULL
  for (nm in names(extra)) meta[[nm]] <- extra[[nm]]
  meta
}

#' Turn a backend's training log into a table.
#'
#' A backend that reports its losses returns one record per logged step. They
#' are stored as a data frame in `meta$estimator$training_trace`, so a run can
#' be inspected afterwards without the backend printing anything. The columns
#' are whatever the backend reported, in the order of its first record.
#'
#' @param trace List of records (each a named list of scalars), or NULL.
#' @return A data frame, or NULL when there is nothing to record.
#' @keywords internal
#' @noRd
.cf_trace_to_frame <- function(trace) {
  if (is.null(trace) || !length(trace)) return(NULL)
  if (is.data.frame(trace)) return(trace)
  nms <- names(trace[[1L]])
  if (is.null(nms)) return(NULL)
  cols <- lapply(nms, function(nm) {
    v <- lapply(trace, function(rec) rec[[nm]])
    v[vapply(v, is.null, logical(1))] <- NA
    unlist(v, use.names = FALSE)
  })
  names(cols) <- nms
  as.data.frame(cols, stringsAsFactors = FALSE)
}

#' Thin a per-step loss series into a trace.
#'
#' A variational fit reports one value per optimisation step, which is far more
#' than is worth carrying. The series is thinned to about `rows` points, always
#' including the first and the step the fit ended on. The result has the same
#' shape as any other training trace: an `iteration` column and the reported
#' loss.
#'
#' Thinning to a count rather than to a fixed cadence is what keeps the two
#' backends that fill this field comparable: the number of optimisation steps
#' depends on the epochs, the batch size and n, so one cadence leaves a
#' different number of rows for each fit.
#'
#' The column is named `loss` because that is what the value is: pyro's
#' `SVI.step` returns the negative of the evidence lower bound, so it falls as
#' the fit converges.
#'
#' @param loss Numeric vector of per-step values, or NULL.
#' @param rows Approximate number of points to keep.
#' @return A data frame with `iteration` and `loss`, or NULL.
#' @keywords internal
#' @noRd
.cf_loss_trace <- function(loss, rows = 20L) {
  if (is.null(loss) || !length(loss)) return(NULL)
  n <- length(loss)
  every <- max(1L, as.integer(ceiling(n / max(1L, as.integer(rows)))))
  idx <- unique(c(seq.int(1L, n, by = every), n))
  .cf_trace_to_frame(lapply(idx, function(i)
    list(iteration = i - 1L, loss = as.numeric(loss[[i]]))))
}

#' Create Standard Meta Object
#'
#' Internal helper to create a standardized \code{meta} pillar for
#' \code{cateomics_result} (a single flat metadata pillar).
#'
#' @param n Number of observations
#' @param p Number of covariates
#' @param covariate_names Character vector of covariate names
#' @param method Character string of method name
#' @param bootstrap Logical, whether bootstrap was used
#' @param random_state Random seed used
#' @param preprocess List describing preprocessing steps
#' @param input_type "data.frame" / "SE" / "MAE"
#' @param factor_levels Named list of factor column train levels (newdata predict)
#' @param x_zscore list(X_mean, X_sd) for X z-score inverse
#' @param y_zscore list(Y_mean, Y_sd) for Y inverse transform (NN methods)
#' @param named_args The names the caller gave (`Y` / `T` / `X_baseline` / `omics`)
#' @param propensity_X PS-specific covariate subset
#' @return A list with the standardized meta structure
#' @keywords internal
create_cf_meta <- function(n,
                           p,
                           covariate_names,
                           method,
                           bootstrap = FALSE,
                           random_state = NULL,
                           preprocess = NULL,
                           input_type = NULL,
                           factor_levels = NULL,
                           x_zscore = NULL,
                           y_zscore = NULL,
                           named_args = NULL,
                           propensity_X = NULL) {

  # Capture software versions
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  pkg_versions <- list(
    CATEomics = tryCatch(as.character(utils::packageVersion("CATEomics")),
                       error = function(e) NA_character_),
    reticulate = tryCatch(as.character(utils::packageVersion("reticulate")),
                          error = function(e) NA_character_)
  )

  method_pkg <- switch(method,
    grf       = "grf",
    bcf       = "bcf",
    NULL
  )
  if (!is.null(method_pkg)) {
    pkg_versions[[method_pkg]] <- tryCatch(
      as.character(utils::packageVersion(method_pkg)),
      error = function(e) NA_character_
    )
  }

  list(
    n = as.integer(n),
    p = as.integer(p),
    covariate_names = covariate_names,
    # `method` is not stored here: the result carries it at the top level
    # (`fit$method`), and a second copy is a second thing to keep in step. The
    # argument is still read, for the backend name and the package version.

    # ----- Bucket 2: provenance (reproducibility metadata) -----
    provenance = list(
      backend           = .infer_backend(method),
      rng_seed          = random_state,
      software_versions = list(
        R         = r_version,
        timestamp = Sys.time(),
        # CATEomics and the R backend of the method that ran. The learners,
        # the preprocessing packages and the Python stack are not here: the
        # Python side is pinned in R/basilisk.R, and the rest are declared
        # dependencies whose versions sessionInfo() reports.
        packages  = pkg_versions
      )
    ),

    # ----- Bucket 3: estimator (method-specific configuration) -----
    estimator = list(
      bootstrap  = bootstrap,
      preprocess = preprocess
    ),

    # ----- named-args facade record (top-level for easy access) -----
    named_args   = named_args,
    propensity_X = propensity_X,

    # ----- preprocessing slots consumed by cf_predict newdata pipeline -----
    # rowData is not retained. The PCA rotation is not here: it belongs to the
    # layer it was fitted on, and sits in that layer's record.
    input_type     = input_type,
    factor_levels  = factor_levels,
    x_zscore       = x_zscore,
    y_zscore       = y_zscore
  )
}

#' Infer the implementation backend used by a method (provenance bucket).
#' @keywords internal
#' @noRd
.infer_backend <- function(method) {
  switch(method,
    grf       = "grf (R)",
    bcf       = "bcf (R)",
    # The S-learner does not go through econml: CATEomics augments the design
    # with the interactions and fits one sklearn model directly.
    slearner  = "CATEomics slearner (sklearn)",
    tlearner  = ,
    xlearner  = "econml.metalearners (Python)",
    drlearner = "econml.dr (Python)",
    ganite    = "CATEomics ganite_impl (PyTorch)",
    cevae     = "pyro-ppl (Python)",
    tarnet    = ,
    dragonnet = "catenets (Python)",
    ensemble  = "CATEomics (R)",
    "unknown"
  )
}
