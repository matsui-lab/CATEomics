# Test helper for the basilisk-managed Python backends.
#
# Migrated Python methods run their code through `.cf_py_run()` (a separate
# basilisk process), so tests no longer probe the main-process reticulate
# session. This helper skips a test when the cfomics basilisk environment
# cannot be started (e.g. basilisk not installed, or the conda env has not
# been provisioned on this machine). The liveness probe is memoised so the
# basilisk process is started at most once per test session.

.cf_basilisk_ok <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ok <- tryCatch(
      requireNamespace("basilisk", quietly = TRUE) &&
        isTRUE(.cf_py_run(function() TRUE)),
      error = function(e) FALSE
    )
    cached <<- ok
    ok
  }
})

skip_if_no_cf_python <- function() {
  testthat::skip_if_not(
    .cf_basilisk_ok(),
    "cfomics basilisk Python environment not available"
  )
}
