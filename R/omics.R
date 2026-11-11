# Resolution of the single `omics` argument into the per-container assay spec.
#
# One user-facing argument describes which measurements are the omics block. `omics` says WHICH layer(s) and WHICH features; `omics_type`
# says how a layer should be treated when it matters (currently: whether the
# layer is RNA-seq and therefore gets library-size normalisation).
#
#   omics = "all"                       every layer, every feature
#   omics = "none"                      no omics block (clinical only)
#   omics = c("BAX", "CASP3")           those features of the single layer
#   omics = list(expr = "all")          named layer, every feature
#   omics = list(methyl = c("cg001"))   named layer, those features
#
# Layer names are SE assay names, MAE experiment names, or free labels for a
# data.frame.

.CF_RNASEQ_LABELS <- c("rnaseq", "rna", "rna_seq", "rna-seq")

#' Is this layer RNA-seq (and therefore internally log-normalised)?
#' @keywords internal
#' @noRd
.cf_layer_is_rnaseq <- function(label, omics_type, n_layers, input_type) {
  tp <- omics_type[[label]]
  if (!is.null(tp) && !is.na(tp)) return(tolower(tp) %in% .CF_RNASEQ_LABELS)
  # A single omics layer is what `preprocess` unambiguously refers to, so it
  # keeps the library-size normalisation regardless of its name. With several
  # layers the name decides, because normalising a methylation or mutation
  # matrix as if it were counts would be wrong.
  if (n_layers <= 1L && !identical(input_type, "MAE")) return(TRUE)
  tolower(label) %in% .CF_RNASEQ_LABELS
}

#' Normalise `omics` / `omics_type` into a per-container assay spec
#'
#' @param omics User value (NULL / "all" / "none" / character / named list).
#' @param omics_type Optional layer-type override; a single string or a vector
#'   named by layer.
#' @param data The input container.
#' @param Y,T,X_baseline Used only to decide what `"all"` means for a
#'   data.frame (every column that is not the outcome, treatment or baseline).
#' @section omics_type: normalised by `.cf_normalise_omics_type()` into a named
#'   character vector, one entry per layer, so the per-layer lookup cannot fall
#'   off the end.
#' @return List with `layers` (named list of "all" or a character vector),
#'   `types` (named logical: is the layer RNA-seq), plus the container-specific
#'   fields the downstream helpers consume.
#' @keywords internal
#' @noRd
.cf_resolve_omics <- function(omics, omics_type = NULL, data = NULL,
                              Y = NULL, T = NULL, X_baseline = NULL) {
  if (is.null(data)) {
    rlang::abort("`data` argument is required.", class = "cateomics_invalid_data")
  }
  input_type <- if (is_mae(data)) "MAE" else if (is_se(data)) "SE" else
    if (is.data.frame(data)) "data.frame" else
      rlang::abort(
        sprintf("Unsupported data class: %s", paste(class(data), collapse = "/")),
        class = "cateomics_unsupported_data_class"
      )

  # ---- shape validation -----------------------------------------------------
  if (!is.null(omics) && !is.character(omics) && !is.list(omics)) {
    rlang::abort(
      sprintf("`omics` must be \"all\", \"none\", a character vector or a named list (got %s).",
              class(omics)[1L]),
      class = "cateomics_invalid_omics"
    )
  }
  if (is.list(omics) &&
      (is.null(names(omics)) || any(!nzchar(names(omics))))) {
    rlang::abort(
      c("Every element of `omics` must be named with its layer.",
        i = "For example omics = list(rnaseq = \"all\", methyl = c(\"cg001\"))."),
      class = "cateomics_invalid_omics"
    )
  }

  # ---- default --------------------------------------------------------------
  if (is.null(omics)) {
    omics <- if (identical(input_type, "data.frame")) "none" else "all"
  }
  if (is.character(omics) && length(omics) == 1L && identical(omics, "none")) {
    return(list(layers = list(), types = logical(0), input_type = input_type,
                assay_name = NULL, experiments = NULL, assay_cols = NULL,
                features = list()))
  }

  avail <- .cf_available_layers(data, input_type)

  # ---- expand the shorthand forms into list(layer = features) ---------------
  layers <- if (is.list(omics)) {
    lapply(names(omics), function(nm) {
      v <- omics[[nm]]
      if (identical(v, "all")) return("all")
      v <- as.character(v)
      if (length(v) == 0L) {
        rlang::abort(
          c(sprintf("`omics$%s` is empty.", nm),
            i = "Use \"all\" for the whole layer, or drop the layer."),
          class = "cateomics_invalid_omics"
        )
      }
      v
    })
  } else if (identical(omics, "all")) {
    if (length(avail) == 0L) {
      # data.frame: one unnamed layer holding every remaining column.
      list(omics = "all")
    } else {
      stats::setNames(rep(list("all"), length(avail)), avail)
    }
  } else if (length(omics) == 1L && omics %in% avail) {
    # A bare layer name means that whole layer.
    stats::setNames(list("all"), omics)
  } else {
    # A bare feature vector: it belongs to the only layer there is.
    if (length(avail) > 1L) {
      rlang::abort(
        c(sprintf("`omics` is a feature vector but %s has several layers: %s.",
                  input_type, paste(avail, collapse = ", ")),
          i = "Name the layer, e.g. omics = list(rnaseq = c(...))."),
        class = "cateomics_omics_layer_ambiguous"
      )
    }
    stats::setNames(list(as.character(omics)),
                    if (length(avail)) avail[[1L]] else "omics")
  }

  if (is.list(omics)) names(layers) <- names(omics)

  # ---- container rules ------------------------------------------------------
  unknown <- setdiff(names(layers), avail)
  if (length(unknown) && !identical(input_type, "data.frame")) {
    rlang::abort(
      c(sprintf("`omics` names a layer absent from the %s: %s.",
                input_type, paste(unknown, collapse = ", ")),
        i = sprintf("Available: %s.", paste(avail, collapse = ", "))),
      class = "cateomics_omics_unknown_layer"
    )
  }
  if (identical(input_type, "SE") && length(avail) > 1L &&
      is.character(omics) && identical(omics, "all")) {
    rlang::abort(
      c(sprintf("This SummarizedExperiment has %d assays, so `omics = \"all\"` is ambiguous.",
                length(avail)),
        i = sprintf("Pick one: %s.",
                    paste(sprintf("omics = \"%s\"", avail), collapse = " / "))),
      class = "cateomics_omics_assay_ambiguous"
    )
  }
  if (identical(input_type, "SE") && length(layers) > 1L) {
    rlang::abort(
      c("A SummarizedExperiment holds one measurement per sample per feature, so its assays cannot be combined.",
        i = sprintf("Pick one: %s.", paste(avail, collapse = ", "))),
      class = "cateomics_omics_multiple_assays"
    )
  }

  omics_type <- .cf_normalise_omics_type(omics_type, names(layers))
  types <- vapply(names(layers), .cf_layer_is_rnaseq, logical(1),
                  omics_type = omics_type, n_layers = length(layers),
                  input_type = input_type)
  # Explicit feature lists (i.e. not "all") are the ones that must be resolved
  # by name before the formula is built.
  features <- layers[vapply(layers, function(v) !identical(v, "all"), logical(1))]

  out <- list(layers = layers, types = types, input_type = input_type,
              # Whether the type came from the layer's name rather than from
              # `omics_type`. The preprocessing step reports what it did for
              # the layers the caller did not name a type for.
              types_from_name = is.na(omics_type[names(layers)]),
              features = features,
              assay_name = NULL, experiments = NULL, assay_cols = NULL)
  if (identical(input_type, "SE")) {
    lab <- names(layers)[[1L]]
    # With unnamed assays the label is a position, and `assay()` needs it as a
    # number: the character method would look for an assay literally called "1".
    real_names <- SummarizedExperiment::assayNames(data)
    out$assay_name <- if (!is.null(real_names) && lab %in% real_names) lab
                      else suppressWarnings(as.integer(lab))
    out$assay_cols <- if (length(features)) features[[1L]] else NULL
  } else if (identical(input_type, "MAE")) {
    out$experiments <- names(layers)
  } else {
    cols <- unlist(lapply(names(layers), function(lab) {
      v <- layers[[lab]]
      if (identical(v, "all")) setdiff(names(data), c(Y, T, X_baseline)) else v
    }), use.names = FALSE)
    missing_cols <- setdiff(cols, names(data))
    if (length(missing_cols)) {
      rlang::abort(
        c("`omics` references column(s) absent from data:",
          sprintf("  %s", paste(utils::head(missing_cols, 5L), collapse = ", "))),
        class = "cateomics_omics_columns_missing"
      )
    }
    cols <- unique(cols)
    if (length(cols) == 0L) {
      return(list(layers = list(), types = logical(0), input_type = input_type,
                  assay_name = NULL, experiments = NULL, assay_cols = NULL,
                  features = list()))
    }
    out$assay_cols <- cols
    # Re-expand "all" so the per-layer feature map is explicit downstream.
    out$layers <- stats::setNames(lapply(names(layers), function(lab) {
      v <- layers[[lab]]
      if (identical(v, "all")) setdiff(names(data), c(Y, T, X_baseline)) else v
    }), names(layers))
    out$features <- out$layers
  }
  out
}

#' Layer names available in a container
#' @keywords internal
#' @noRd
.cf_available_layers <- function(data, input_type) {
  if (identical(input_type, "SE")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(character(0))
    nm <- SummarizedExperiment::assayNames(data)
    if (is.null(nm) || !length(nm)) nm <- as.character(seq_len(length(
      SummarizedExperiment::assays(data))))
    return(nm)
  }
  if (identical(input_type, "MAE")) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) return(character(0))
    return(names(MultiAssayExperiment::experiments(data)))
  }
  character(0)
}

#' Normalise `omics_type` into one entry per layer
#'
#' Accepts a single unnamed value (applies to every layer) or a value named by
#' layer. An unnamed vector longer than one is rejected rather than having its
#' first element quietly applied to everything.
#'
#' @keywords internal
#' @noRd
.cf_normalise_omics_type <- function(omics_type, layer_names) {
  if (is.null(omics_type)) {
    return(stats::setNames(rep(NA_character_, length(layer_names)), layer_names))
  }
  tp <- unlist(omics_type, use.names = TRUE)
  if (is.null(names(tp)) || any(!nzchar(names(tp)))) {
    if (length(tp) != 1L) {
      rlang::abort(
        c("`omics_type` must be a single value or be named by layer.",
          i = sprintf("Layers: %s.", paste(layer_names, collapse = ", "))),
        class = "cateomics_invalid_omics_type"
      )
    }
    return(stats::setNames(rep(as.character(tp), length(layer_names)),
                           layer_names))
  }
  bad <- setdiff(names(tp), layer_names)
  if (length(bad)) {
    rlang::abort(
      c(sprintf("`omics_type` names a layer that is not in `omics`: %s.",
                paste(bad, collapse = ", ")),
        i = sprintf("Layers: %s.", paste(layer_names, collapse = ", "))),
      class = "cateomics_invalid_omics_type"
    )
  }
  out <- stats::setNames(rep(NA_character_, length(layer_names)), layer_names)
  out[names(tp)] <- as.character(tp)
  out
}
