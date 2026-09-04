# Fork-safety of the basilisk bridge (.cf_py_run).
#
# When a foreign Python is loaded in the process, basilisk uses socket mode
# (parallel::makePSOCKcluster), whose default port is inherited by every fork,
# so concurrent workers collide on one port and the Python fits are dropped.
# .cf_basilisk_start() gives each process a distinct port + retries. We force
# socket mode here with setBasiliskFork(FALSE) so the test exercises the socket
# path regardless of whether Python is already loaded, then run several
# concurrent .cf_py_run() calls under mclapply and require all to succeed.

test_that(".cf_py_run survives concurrent forks in socket mode", {
  skip_on_os("windows")
  skip_if_no_cf_python()
  skip_if_not_installed("parallel")
  # mclapply forking is a no-op on Windows; and single-core machines cannot
  # exhibit the collision.
  if (parallel::detectCores() < 2L) skip("needs >= 2 cores")

  old_fork <- basilisk::getBasiliskFork()
  basilisk::setBasiliskFork(FALSE)              # force PSOCK socket mode
  on.exit(basilisk::setBasiliskFork(old_fork), add = TRUE)

  # R CMD check sets _R_CHECK_LIMIT_CORES_, which caps mclapply at 2 cores
  # (parallel:::.check_ncores() errors above that). Two concurrent forks still
  # exercise the shared-port collision; dev runs (devtools::test) use 6.
  n <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else 6L
  res <- parallel::mclapply(seq_len(n), function(i) {
    tryCatch(
      cfomics:::.cf_py_run(function() 2L + 2L),
      error = function(e) conditionMessage(e)
    )
  }, mc.cores = n)

  ok <- vapply(res, function(x) identical(x, 4L), logical(1))
  # Without the per-process port fix, some workers fail with
  # "creation of server socket failed: port ... cannot be opened".
  failed <- res[!ok]
  expect_true(
    all(ok),
    info = paste0(sum(!ok), "/", n, " concurrent basilisk starts failed: ",
                  paste(unique(unlist(failed)), collapse = " | "))
  )
})
