# Opt-in skips for the heaviest test groups, to speed up dev iteration
# (full R CMD check is dominated by bcf MCMC + Python backends).
#
#   CFOMICS_SKIP_BCF=1     -> skip bcf MCMC tests
#   CFOMICS_SKIP_PYTHON=1  -> skip Python-backend tests
#
# These are dev-only accelerators; CI / release checks run with them unset so
# full coverage still executes.

skip_if_skip_bcf <- function() {
  if (identical(Sys.getenv("CFOMICS_SKIP_BCF"), "1")) {
    testthat::skip("CFOMICS_SKIP_BCF=1 (bcf MCMC tests skipped)")
  }
}

skip_if_skip_python <- function() {
  if (identical(Sys.getenv("CFOMICS_SKIP_PYTHON"), "1")) {
    testthat::skip("CFOMICS_SKIP_PYTHON=1 (Python tests skipped)")
  }
}
