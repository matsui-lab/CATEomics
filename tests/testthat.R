# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(cfomics)

## Sub-commit 6 (2026-05-24): also write the full per-test stream to a
## sibling log file so that R CMD check's `Last 13 lines` summary can be
## augmented when diagnosing fails under the installed-package context.
## Default `summary` reporter still controls the stdout the check
## harness consumes; the extra log is only used during fix iteration.
log_path <- Sys.getenv("CFOMICS_TESTTHAT_LOG", "")
if (nzchar(log_path)) {
  reporter <- testthat::MultiReporter$new(reporters = list(
    testthat::SummaryReporter$new(),
    testthat::TapReporter$new(file = log_path)
  ))
  test_check("cfomics", reporter = reporter)
} else {
  test_check("cfomics")
}
