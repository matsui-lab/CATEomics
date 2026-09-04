## helper-gpu.R — VRAM budget enforcement for integration tests
##
## GPU VRAM budget = 10-20 GB (20-40% of a 48 GB card). Docker cannot cap VRAM,
## so it is controlled per-process on the application side (PyTorch /
## TensorFlow).
##
## Loads automatically before all test_that blocks (testthat helper convention).
## Idempotent — safe to call multiple times in same R session.
##
## Effects:
##   - PyTorch: torch.cuda.set_per_process_memory_fraction(0.4)
##     → 48 GB * 0.4 ≈ 19.2 GB cap (within 10-20 GB target)
##   - TensorFlow: set_memory_growth(True) for all visible GPUs
##     → grow on demand instead of pre-allocate (TF defaults to all VRAM)
##     → on this env TF may not recognise the CUDA runtime and fall back to
##       CPU, making this a no-op
##
## Note: This is test-only setup. Production users rely on basilisk (automatic)
## for the same behaviour outside tests.

.cfomics_test_setup_gpu <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
  if (!isTRUE(reticulate::py_available(initialize = FALSE))) {
    ## Skip if Python not initialized — avoid forcing init here.
    return(invisible(FALSE))
  }

  ## --- PyTorch VRAM fraction ---
  ok_torch <- tryCatch({
    torch <- reticulate::import("torch", convert = TRUE)
    if (isTRUE(torch$cuda$is_available())) {
      ## 0.4 = 40% of total VRAM ≈ 19.2 GB on RTX PRO 5000 Blackwell (48 GB)
      torch$cuda$set_per_process_memory_fraction(0.4, device = 0L)
      message(sprintf(
        "[helper-gpu] torch.cuda.set_per_process_memory_fraction(0.4) on device %s",
        torch$cuda$get_device_name(0L)
      ))
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    message("[helper-gpu] torch VRAM setup skipped: ", conditionMessage(e))
    FALSE
  })

  ## --- TensorFlow GPU hide (no longer needed) ---
  ## ganite uses the PyPI ganite (PyTorch) backend, so TF is not on cfomics'
  ## critical path; ganite runs on PyTorch CUDA natively. We deliberately do
  ## NOT import tensorflow here — importing TF
  ## initializes its CUDA context which conflicts with PyTorch's, causing
  ## segfaults in subsequent torch cuda ops (verified with cf_has_python). If a
  ## stale TF is sitting in sys.modules from a parallel test (cf_interpret IG),
  ## that is on the caller to manage.
  invisible(c(torch = ok_torch))
}

.cfomics_test_setup_gpu()
