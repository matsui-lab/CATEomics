#' Attach CATEomics results back to SE / MAE / data.frame
#'
#' @description
#' Writes \code{cf_model} per-sample fields (the estimated effect, the
#' benefit group) into the colData of a \code{SummarizedExperiment} /
#' \code{MultiAssayExperiment} (or adds columns to a data.frame), and the
#' variable importance onto the feature annotation.
#'
#' @param data SE / MAE / data.frame to augment.
#' @param model \code{cf_model}.
#' @param what Character vector subset of \code{c("all", "tau", "group",
#'   "vi")} (default \code{"all"}). \code{"all"} means everything that can
#'   be written from what was passed: \code{tau} always, the group when
#'   \code{groups} is given, the importance when \code{itp} is given and the
#'   object has a feature annotation. \code{"tau"} writes the out-of-fold
#'   \eqn{\hat\tau}, or the in-sample one when the fit has no out-of-fold
#'   estimate. A field named explicitly whose input is missing is an error;
#'   reached through \code{"all"} it is reported and skipped. \code{"vi"}
#'   needs a feature annotation to write to, so naming it for a data frame is
#'   an error as well.
#' @param groups Optional \code{cateomics_groups} for the \code{"group"}
#'   field.
#' @param itp Optional \code{cateomics_interpret} for the \code{"vi"}
#'   field (only when \code{data} is SE/MAE and feature-level write
#'   into rowData is requested via top-k VI markers).
#' @param surrogate Optional name of the surrogate whose variable
#'   importance is written for \code{"vi"} (one of the surrogates in
#'   \code{itp}, e.g. \code{"ridge"} or \code{"xgb"}). \code{NULL}
#'   (default) uses the first surrogate present in \code{itp}.
#' @param overwrite Logical (default \code{FALSE}); overwrite existing
#'   columns of the same name.
#' @return Augmented \code{data} (same class as input). The columns written
#'   are \code{cf_tau} and \code{cf_group} and, for a
#'   \code{SummarizedExperiment} or \code{MultiAssayExperiment},
#'   \code{cf_vi_importance} and \code{cf_vi_rank} on the feature
#'   annotation. The average effect is one number for the whole fit
#'   (\code{model$res$ate}), so it is not written as a column.
#' @examples
#' # Write the estimated individual effects back onto the SummarizedExperiment.
#' se  <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4, seed = 1)
#' fit <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'               X_baseline = c("age", "sex"), method = "grf",
#'               n_features = 40, partition = "oof")
#' se <- cf_attach_results(se, fit, what = "tau")
#' head(SummarizedExperiment::colData(se)$cf_tau)
#' @export
cf_attach_results <- function(data, model,
                              what      = "all",
                              groups    = NULL,
                              itp       = NULL,
                              surrogate = NULL,
                              overwrite = FALSE) {
  # The default is the mode `"all"`, not the list of every item: the guard below
  # treats what it is given as named by the caller, and a default holding
  # "group" and "vi" would read as a request for both.
  what <- .cf_match_several(what, c("all", "tau", "group", "vi"),
                            "what", "cf_attach_results()")
  vi_requested <- "vi" %in% what      # asked for by name, not via "all"
  # A field that was named needs its input: writing nothing and saying nothing
  # is the one outcome the caller cannot act on.
  needs <- c(group = "groups", vi = "itp")
  have  <- c(group = !is.null(groups), vi = !is.null(itp))
  named_missing <- intersect(names(have)[!have], what)
  if (length(named_missing) > 0L) {
    rlang::abort(
      c(sprintf("what = %s was asked for, but %s.",
                paste(sprintf("'%s'", named_missing), collapse = " and "),
                paste(sprintf("`%s` is missing", needs[named_missing]),
                      collapse = " and ")),
        i = sprintf("Pass %s, or drop %s from `what`.",
                    paste(sprintf("`%s =`", needs[named_missing]),
                          collapse = " and "),
                    paste(sprintf("'%s'", named_missing), collapse = " and "))),
      class = "cateomics_attach_missing_input"
    )
  }
  if ("all" %in% what) {
    skipped <- names(have)[!have]
    if (length(skipped) > 0L) {
      rlang::inform(
        sprintf("what = 'all': %s not written (%s not given).",
                paste(sprintf("cf_%s", skipped), collapse = ", "),
                paste(sprintf("`%s`", needs[skipped]), collapse = ", ")),
        class = "cateomics_attach_all_partial")
    }
    what <- c("tau", names(have)[have])
  }
  if (!inherits(model, "cf_model")) {
    rlang::abort("model must be a cf_model.",
                  class = "cateomics_invalid_input")
  }

  add_col <- function(cd, name, value) {
    if (!isTRUE(overwrite) && name %in% colnames(cd)) {
      rlang::warn(
        sprintf("Column '%s' exists; skipping (use overwrite = TRUE).",
                name),
        class = "cateomics_attach_skipped"
      )
      return(cd)
    }
    cd[[name]] <- value
    cd
  }

  # The model's per-sample results are in the order the fit saw the samples,
  # which is not always the order of the object being written to: a
  # MultiAssayExperiment is joined through its sampleMap. Matching on the
  # sample names is what keeps each estimate on its own sample; equal lengths
  # alone would let a permutation through in silence.
  model_ids <- rownames(model$data$X_preprocessed %||% model$data$X)
  # Whether those row names name anything. A data.frame fit carries R's
  # automatic row numbers, which name nothing, and "1" .. "n" is also a
  # legitimate set of integer sample identifiers, so the two cannot be told
  # apart from the strings. `cf_fit` / `cf_compare` record the answer, which
  # is where it is still knowable. An object without that record falls back
  # to reading the shape.
  named <- model$meta$provenance$sample_names_meaningful
  if (is.null(named)) {
    named <- !identical(as.character(model_ids),
                        as.character(seq_along(model_ids)))
  }
  if (!isTRUE(named)) model_ids <- NULL
  align <- function(v, df) {
    if (is.null(v)) return(NULL)
    ids <- rownames(df)
    if (is.null(model_ids) || is.null(ids)) {
      # Nothing names the samples on one side, so position is all there is.
      return(if (length(v) == nrow(df)) v else NULL)
    }
    if (length(v) != length(model_ids)) return(NULL)
    if (all(model_ids %in% ids)) {
      # Covers both the full set and a subset: a MultiAssayExperiment whose
      # sampleMap does not reach every primary is fitted on the complete ones
      # only, so the rest have no estimate. `match` yields NA for those, in
      # whatever type `v` is.
      n_missing <- sum(!ids %in% model_ids)
      if (n_missing > 0L) {
        rlang::inform(
          sprintf(paste("%d of the %d samples were not in the fit and are",
                        "written as NA."), n_missing, length(ids)),
          class = "cateomics_attach_partial"
        )
      }
      return(v[match(ids, model_ids)])
    }
    rlang::abort(
      c("The fitted model names samples the object does not have.",
        i = "cf_attach_results writes each estimate onto its own sample, so every fitted sample has to be present."),
      class = "cateomics_attach_sample_mismatch"
    )
  }

  attach_to_df <- function(df) {
    if ("tau" %in% what) {
      # Inside the branch: aligning can abort, and a caller who did not ask
      # for the estimate should not pay for that.
      # The out-of-fold estimate comes first. The column feeds downstream
      # description, and cf_benefit_groups and cf_interpret read the
      # out-of-fold tau-hat too, so writing the in-sample one here would put a
      # different quantity next to theirs in the same object. A fit with
      # partition = "none" has no out-of-fold estimate and falls back.
      tau_vec <- align(model$res$tau_hat_oof %||% model$res$tau_hat, df)
      if (!is.null(tau_vec)) df <- add_col(df, "cf_tau", tau_vec)
    }
    if ("group" %in% what) {
      g <- align(groups$groups, df)
      if (!is.null(g)) {
        df <- add_col(df, "cf_group",
                       factor(g, levels = seq_along(groups$labels),
                               labels = groups$labels))
      }
    }
    df
  }

  if (is_se(data)) {
    cd <- SummarizedExperiment::colData(data)
    cd <- attach_to_df(as.data.frame(cd))
    SummarizedExperiment::colData(data) <-
      S4Vectors::DataFrame(cd, check.names = FALSE)
    if ("vi" %in% what && length(itp$vi) > 0L) {
      vi_df <- .cf_vi_table(itp, surrogate)
      SummarizedExperiment::rowData(data) <-
        .cf_rowdata_with_vi(SummarizedExperiment::rowData(data), vi_df,
                            overwrite = overwrite)
    }
    return(data)
  }
  if (is_mae(data)) {
    cd <- MultiAssayExperiment::colData(data)
    cd <- attach_to_df(as.data.frame(cd))
    MultiAssayExperiment::colData(data) <-
      S4Vectors::DataFrame(cd, check.names = FALSE)
    if ("vi" %in% what && length(itp$vi) > 0L) {
      vi_df <- .cf_vi_table(itp, surrogate)
      for (ex in names(MultiAssayExperiment::experiments(data))) {
        se_ex <- data[[ex]]
        if (!inherits(se_ex, "SummarizedExperiment")) next
        # A feature can be labelled `<experiment>_<feature>` (the MAE
        # covariate naming) or bare; try the prefixed form first.
        SummarizedExperiment::rowData(se_ex) <-
          .cf_rowdata_with_vi(SummarizedExperiment::rowData(se_ex), vi_df,
                              prefix = ex, overwrite = overwrite)
        data[[ex]] <- se_ex
      }
    }
    return(data)
  }
  if (is.data.frame(data)) {
    if ("vi" %in% what) {
      msg <- c(
        "A data.frame has no per-feature slot to write variable importance into.",
        i = "Read it from the interpretation itself: itp$vi[[surrogate]]$features.",
        i = "The per-sample results (tau, group) are attached as columns.")
      # An explicit `what = "vi"` cannot be honoured; the same field reached
      # via `what = "all"` is a convenience, so the rest is still attached.
      if (vi_requested) {
        rlang::abort(msg, class = "cateomics_attach_vi_unsupported")
      }
      rlang::warn(msg, class = "cateomics_attach_vi_skipped")
    }
    return(attach_to_df(data))
  }
  rlang::abort("Unsupported data class for cf_attach_results.",
                class = "cateomics_invalid_input")
}

#' Pick one surrogate's feature table from an interpretation
#' @keywords internal
#' @noRd
.cf_vi_table <- function(itp, surrogate) {
  sur <- if (is.null(surrogate)) names(itp$vi)[1L] else surrogate
  if (!sur %in% names(itp$vi)) {
    rlang::abort(
      sprintf("surrogate '%s' not in the interpretation (available: %s).",
              sur, paste(names(itp$vi), collapse = ", ")),
      class = "cateomics_invalid_surrogate"
    )
  }
  itp$vi[[sur]]$features
}

#' Add importance / rank columns to a rowData frame
#'
#' Rows with no matching feature get `NA`. With `prefix`, feature names of the
#' form `<prefix>_<feature>` are matched as well, which is how a MAE labels a
#' covariate that came from one of its experiments. Existing columns of the
#' same name follow the same rule as the colData writes: kept unless
#' `overwrite` is TRUE.
#' @keywords internal
#' @noRd
.cf_rowdata_with_vi <- function(rd, vi_df, prefix = NULL, overwrite = FALSE) {
  if (is.null(vi_df) || nrow(vi_df) == 0L) return(rd)
  keys <- rownames(rd)
  ord  <- match(keys, vi_df$feature_name)
  if (!is.null(prefix)) {
    ord2 <- match(paste0(prefix, "_", keys), vi_df$feature_name)
    ord[is.na(ord)] <- ord2[is.na(ord)]
  }
  if (!any(!is.na(ord))) return(rd)
  rd_df <- as.data.frame(rd)
  values <- list(cf_vi_importance = vi_df$importance[ord],
                 cf_vi_rank       = vi_df$rank[ord])
  # The two columns are one table's importance and its rank, so writing one
  # over a stale copy of the other would leave them disagreeing.
  present <- intersect(names(values), colnames(rd_df))
  if (!isTRUE(overwrite) && length(present) > 0L) {
    rlang::warn(
      sprintf("Column(s) %s exist; skipping (use overwrite = TRUE).",
              paste(sprintf("'%s'", present), collapse = " and ")),
      class = "cateomics_attach_skipped"
    )
    return(rd)
  }
  for (nm in names(values)) rd_df[[nm]] <- values[[nm]]
  S4Vectors::DataFrame(rd_df, check.names = FALSE)
}
