#' Attach cfomics results back to SE / MAE / data.frame (E-1)
#'
#' @description
#' Writes \code{cf_model} per-sample fields (\code{tau}, group
#' assignment, ATE constant, etc.) into the colData of a
#' \code{SummarizedExperiment} / \code{MultiAssayExperiment} (or
#' adds columns to a data.frame).
#'
#' @param data SE / MAE / data.frame to augment.
#' @param model \code{cf_model}.
#' @param what Character vector subset of \code{c("all", "tau", "group",
#'   "vi", "ate")} (default \code{"all"}).
#' @param groups Optional \code{cfomics_groups} for the \code{"group"}
#'   field.
#' @param itp Optional \code{cfomics_interpret} for the \code{"vi"}
#'   field (only when \code{data} is SE/MAE and feature-level write
#'   into rowData is requested via top-k VI markers).
#' @param surrogate Optional name of the surrogate whose variable
#'   importance is written for \code{"vi"} (one of the surrogates in
#'   \code{itp}, e.g. \code{"ridge"} or \code{"xgb"}). \code{NULL}
#'   (default) uses the first surrogate present in \code{itp}.
#' @param overwrite Logical (default \code{FALSE}); overwrite existing
#'   columns of the same name.
#' @return Augmented \code{data} (same class as input).
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
                              what      = c("all", "tau", "group",
                                             "vi", "ate"),
                              groups    = NULL,
                              itp       = NULL,
                              surrogate = NULL,
                              overwrite = FALSE) {
  what <- match.arg(what, c("all", "tau", "group", "vi", "ate"),
                     several.ok = TRUE)
  if ("all" %in% what) {
    what <- c("tau", "group", "ate")
    if (!is.null(itp)) what <- c(what, "vi")
  }
  if (!inherits(model, "cf_model")) {
    rlang::abort("model must be a cf_model.",
                  class = "cfomics_invalid_input")
  }

  add_col <- function(cd, name, value) {
    if (!isTRUE(overwrite) && name %in% colnames(cd)) {
      rlang::warn(
        sprintf("Column '%s' exists; skipping (use overwrite = TRUE).",
                name),
        class = "cfomics_attach_skipped"
      )
      return(cd)
    }
    cd[[name]] <- value
    cd
  }

  attach_to_df <- function(df) {
    n <- nrow(df)
    tau_vec <- model$res$tau_hat
    if ("tau" %in% what && length(tau_vec) == n) {
      df <- add_col(df, "cf_tau", tau_vec)
    }
    if ("ate" %in% what) {
      df <- add_col(df, "cf_ate", rep(model$res$ate, n))
    }
    if ("group" %in% what && !is.null(groups)) {
      g <- groups$groups
      if (length(g) == n) {
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
    if ("vi" %in% what && !is.null(itp) && length(itp$vi) > 0L) {
      rd <- SummarizedExperiment::rowData(data)
      # cf_interpret$vi is a named list, one entry per surrogate. Default to the
      # first surrogate present; `surrogate =` selects a specific one.
      sur <- if (is.null(surrogate)) names(itp$vi)[1L] else surrogate
      if (!sur %in% names(itp$vi)) {
        rlang::abort(
          sprintf("surrogate '%s' not in the interpretation (available: %s).",
                  sur, paste(names(itp$vi), collapse = ", ")),
          class = "cfomics_invalid_surrogate"
        )
      }
      vi_df <- itp$vi[[sur]]$features
      ord <- if (is.null(vi_df)) NA_integer_
             else match(rownames(rd), vi_df$feature_name)
      if (any(!is.na(ord))) {
        rd_df <- as.data.frame(rd)
        rd_df$cf_vi_importance <- vi_df$importance[ord]
        rd_df$cf_vi_rank       <- vi_df$rank[ord]
        SummarizedExperiment::rowData(data) <-
          S4Vectors::DataFrame(rd_df, check.names = FALSE)
      }
    }
    return(data)
  }
  if (is_mae(data)) {
    cd <- MultiAssayExperiment::colData(data)
    cd <- attach_to_df(as.data.frame(cd))
    MultiAssayExperiment::colData(data) <-
      S4Vectors::DataFrame(cd, check.names = FALSE)
    return(data)
  }
  if (is.data.frame(data)) {
    return(attach_to_df(data))
  }
  rlang::abort("Unsupported data class for cf_attach_results.",
                class = "cfomics_invalid_input")
}
