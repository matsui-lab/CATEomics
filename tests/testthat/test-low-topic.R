# ensure_*_env are not exported.
#
# There is no on-disk config subsystem, so `save_config` does not exist. The
# absence is asserted in test-onload-no-sideeffects.R ("the removed config
# subsystem is not exported").

test_that("ensure_*_env are not exported", {
  ns_exports <- getNamespaceExports("CATEomics")
  expect_false("ensure_conda_env" %in% ns_exports)
  expect_false("ensure_python_env" %in% ns_exports)
})
