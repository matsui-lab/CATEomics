# =============================================================================
# GPU device propagation helpers
#
# device resolution contract:
# - device = "auto" (default): "gpu" if cuda is available, else "cpu"
# - device = "cpu": force CPU (sets CUDA_VISIBLE_DEVICES="")
# - device = "gpu": request GPU (R-native methods warn + fall back to CPU;
#                   Python NN methods attempt GPU in the backend + auto fallback)
#
# The resolved value ("cpu" or "gpu") is recorded in meta$provenance$device.
# =============================================================================

#' Detect CUDA / GPU availability via Torch
#'
#' Treats GPU as available when `torch.cuda.is_available()` is True in the
#' Python environment; treated as CPU when the Python env is unset or the
#' check fails.
#'
#' Only torch.cuda is queried. All NN methods (cevae / tarnet / dragonnet /
#' ganite) share the PyTorch stack, and importing TensorFlow unnecessarily can
#' trigger a Blackwell sm_120 PTX JIT warning and a segfault from TF + PyTorch
#' CUDA-context coexistence, so the check is restricted to a single torch.cuda
#' query.
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

#' Resolve `device` argument to canonical "cpu" / "gpu"
#'
#' R-native methods do not support GPU, so a "gpu" request resolves to "cpu".
#' Python methods fall back automatically in the backend; here the user's
#' intent is forwarded as-is (actual GPU absence is decided by the backend).
#'
#' @param device User-specified `device` ("auto" / "cpu" / "gpu").
#' @param method cf_fit method name (R-native vs Python distinction).
#' @param r_native_methods Character vector of method names that are R-only.
#' @return Resolved device: "cpu" or "gpu" (auto → 自動判定後の値).
#' @keywords internal
#' @noRd
.resolve_device <- function(device, method, r_native_methods) {
  device <- match.arg(device, c("auto", "cpu", "gpu"))
  if (identical(device, "cpu")) {
    # This resolver stays pure. The CUDA_VISIBLE_DEVICES forcing
    # (and its restore on exit) is done by the caller via
    # `.with_cpu_cuda_hidden()` so a CPU fit does not permanently pollute
    # the session's CUDA visibility for later GPU fits.
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
    return("cpu")  # the cf_fit caller has already emitted a warning
  }
  "gpu"
}
