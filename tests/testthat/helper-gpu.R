## helper-gpu.R — VRAM budget enforcement for integration tests
##
## Docker cannot cap VRAM, so a test that wants a bound sets one per process
## on the application side (PyTorch).
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

.cateomics_test_setup_gpu <- function() {
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

  ## No backend uses TensorFlow -- every Python method runs on PyTorch or
  ## sklearn -- and TensorFlow is deliberately not imported here: importing it
  ## initialises a CUDA context that conflicts with PyTorch's and segfaults the
  ## torch cuda calls that follow.
  invisible(c(torch = ok_torch))
}

.cateomics_test_setup_gpu()
