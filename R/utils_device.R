# =============================================================================
# GPU device propagation helpers
#
# device resolution contract:
# - device = "auto" (default): "gpu" if cuda is available, else "cpu"
# - device = "cpu": force CPU
# - device = "gpu": request GPU. The six methods with no GPU path fall back to
#                   CPU, and the entry that was called says so (cf_fit per
#                   method, cf_compare once for the set); the NN methods
#                   attempt GPU in the backend with an automatic fallback.
#
# The caller's argument, before resolution, is recorded in
# meta$provenance$device_requested; the device the backend reports goes to
# meta$provenance$device_effective, and can be a value such as
# "cuda:0".
# =============================================================================

#' Detect CUDA / GPU availability via Torch
#'
#' Treats GPU as available when `torch.cuda.is_available()` is True in the
#' Python environment; treated as CPU when the Python env is unset or the
#' check fails.
#'
#' Only torch.cuda is queried: every NN method (cevae / tarnet / dragonnet /
#' ganite) runs on the PyTorch stack.
#'
#' @return Logical scalar: TRUE if any GPU detected.
#' @keywords internal
#' @noRd
.detect_gpu <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    if (reticulate::py_module_available("torch")) {
      torch <- reticulate::import("torch", convert = TRUE)
      isTRUE(torch$cuda$is_available())
    } else FALSE
  }, error = function(e) FALSE)
}

#' Methods with no GPU path.
#'
#' The R-native backends and the econml meta / DR learners run on the CPU
#' whatever the caller asks for; the four Python neural methods decide
#' backend-side via `torch.cuda.is_available()`. Shared so `cf_fit` and
#' `cf_compare` resolve a device the same way.
#'
#' @keywords internal
#' @noRd
.CF_GPU_UNSUPPORTED <- c("grf", "bcf", "slearner", "tlearner", "xlearner",
                         "drlearner")

#' Resolve `device` argument to canonical "cpu" / "gpu"
#'
#' The methods without GPU support resolve a "gpu" request to "cpu". The NN
#' methods fall back automatically in the backend; here the user's intent is
#' forwarded as-is (actual GPU absence is decided by the backend).
#'
#' @param device User-specified `device` ("auto" / "cpu" / "gpu").
#' @param method cf_fit method name.
#' @param r_native_methods Character vector of method names with no GPU path.
#'   The caller passes the six methods that have none (grf and bcf are
#'   R-native; the four econml ones run on CPU), so the name is narrower than
#'   the contents.
#' @return Resolved device: "cpu" or "gpu" ("auto" resolves to one of them).
#' @keywords internal
#' @noRd
.resolve_device <- function(device, method, r_native_methods) {
  device <- match.arg(device, c("auto", "cpu", "gpu"))
  if (identical(device, "cpu")) {
    # This resolver stays pure: the resolved value is what each Python
    # backend applies itself, passed as an argument on every call because the
    # basilisk worker is reused and inherits the environment only once.
    return("cpu")
  }
  if (identical(device, "auto")) {
    # R-native is pinned to CPU; Python resolves to "gpu" when a GPU is detected.
    if (method %in% r_native_methods) {
      return("cpu")
    }
    return(if (.detect_gpu()) "gpu" else "cpu")
  }
  # device == "gpu"
  if (method %in% r_native_methods) {
    # Both entries say so before calling this: cf_fit for its one method,
    # cf_compare once for the compared set.
    return("cpu")
  }
  "gpu"
}

#' Device a saved neural model should be reloaded and predicted on.
#'
#' The fit records where it ran in `res$device_effective` (e.g. `"cpu"` or
#' `"cuda:0"`). Reusing it keeps `predict(newdata = )` on the same hardware
#' as the fit; without it the reload takes whatever GPU is visible, so one
#' saved model answers differently on different hosts. Falls back to the CPU
#' when nothing was recorded, which is the host-independent choice. Shared by
#' every backend that saves its model to a file and reloads it for newdata.
#'
#' @param object A fitted cf_model with a torch-backed backend.
#'
#' @return `"cpu"` or `"cuda"`.
#'
#' @keywords internal
#' @noRd
.cf_predict_device <- function(object) {
  dev <- object$res$device_effective %||%
         object$meta$provenance$device_effective
  # The backends record NA when the probe itself failed, and `%||%` only
  # replaces NULL, so NA reaches here. `is.na()` is vectorised, so a record
  # that is absent (character(0)) or carries more than one entry would make
  # the condition zero- or two-length and the call would stop on it rather
  # than fall back. Anything that is not a single usable string is cpu.
  if (length(dev) != 1L || is.na(dev) || !nzchar(dev)) return("cpu")
  if (isTRUE(startsWith(as.character(dev), "cuda"))) "cuda" else "cpu"
}

#' Record which device a fit asked for and which one it ran on.
#'
#' `device_requested` is the raw argument and `device_effective` what the
#' backend reports, falling back to the resolved default when a backend reports
#' nothing (a character(0) would otherwise leave both places the predict path
#' reads empty). Shared by `cf_fit` and `cf_compare` so a comparison row
#' carries the same record as a single fit.
#'
#' @param provenance The provenance bucket to add to.
#' @param device Raw `device` argument.
#' @param device_resolved Device after resolution.
#' @param device_effective What the backend reported, or NULL.
#' @return The provenance bucket with the two device fields.
#' @keywords internal
#' @noRd
.cf_record_device <- function(provenance, device, device_resolved,
                              device_effective = NULL) {
  dev_eff <- device_effective
  if (length(dev_eff) != 1L || is.na(dev_eff) || !nzchar(dev_eff)) {
    dev_eff <- device_resolved
  }
  provenance$device_requested <- device
  provenance$device_effective <- dev_eff
  provenance
}
