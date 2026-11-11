# =============================================================================
# basilisk-managed Python environment for all Python backends.
#
# Every Python-backed method (s/t/x-learner, drlearner, cevae, ganite,
# tarnet, dragonnet) runs its Python code inside a single reproducible conda
# environment provisioned by basilisk. basilisk isolates the
# environment (its own miniconda + pinned packages), so the package works out of
# the box on a clean Bioconductor machine and coexists with any RETICULATE_PYTHON
# the user may have set (basilisk overrides it before reticulate init).
#
# basilisk may or may not run Python in a separate process: under its shared
# mode it binds reticulate inside the calling R process, and it falls back to a
# socket worker when a foreign Python is already initialised there. A fitted
# Python model therefore cannot be relied on to survive as a live R object.
# Each method serialises its fitted model to a file during fit (joblib for
# sklearn/econml, torch.save for the NN backends) and reloads it at predict
# time. Only R-serialisable values (numeric result vectors + the model path +
# metadata) cross back.
#
# Because the interpreter can be the caller's own, the Python-side global
# generators (numpy / random / torch) are snapshotted and restored around
# every call by `.cf_py_rng_guarded()`.
#
# The environment is a conda `python` base plus pip-pinned packages, which
# pins the versions the package was developed and tested against.
# =============================================================================

#' Pinned pip packages for the CATEomics unified environment.
#' @keywords internal
#' @noRd
.cateomics_basilisk_pip <- c(
  "numpy==2.4.6", "scipy==1.15.3", "scikit-learn==1.6.1", "pandas==3.0.3",
  "statsmodels==0.14.6", "Cython==3.2.4",
  "torch==2.12.0",
  "econml==0.16.0", "catenets==0.2.4", "pyro-ppl==1.9.1",
  "shap==0.48.0", "xgboost==2.1.4"
)

#' The CATEomics unified basilisk environment.
#'
#' Defined at load time (a cheap S4 object); the heavy conda environment is
#' built lazily by basilisk on the first \code{basiliskStart()}.
#' @keywords internal
#' @noRd
.cateomics_env <- basilisk::BasiliskEnvironment(
  envname  = "cateomics_unified",
  pkgname  = "CATEomics",
  packages = "python=3.12.3",
  pip      = .cateomics_basilisk_pip
)

#' Run a function inside the CATEomics basilisk environment.
#'
#' Starts (and stops) a basilisk worker bound to the CATEomics unified conda
#' environment and evaluates \code{fun(...)} inside it; under basilisk's shared
#' mode that worker is the calling R process itself. \code{fun} may call
#' \code{reticulate::import()} freely; only its return value (which must be
#' R-serialisable) crosses back. This is the single entry point every Python
#' backend uses instead of importing modules in the main R process, and it
#' restores the Python-side global generators afterwards.
#'
#' @param fun A function to evaluate in the basilisk worker.
#' @param ... Arguments forwarded to \code{fun} (must be R-serialisable).
#' @return The value returned by \code{fun}.
#' @keywords internal
#' @noRd
.cf_py_run <- function(fun, ...) {
  # Keep the torch/inductor compile cache inside the session tempdir (inherited
  # by the basilisk worker) so it is cleaned up with the session instead of
  # leaving a torchinductor_<user> directory in the system temp.
  if (!nzchar(Sys.getenv("TORCHINDUCTOR_CACHE_DIR"))) {
    ti_dir <- file.path(tempdir(), "cateomics_torchinductor")
    dir.create(ti_dir, showWarnings = FALSE, recursive = TRUE)
    Sys.setenv(TORCHINDUCTOR_CACHE_DIR = ti_dir)
    on.exit(Sys.unsetenv("TORCHINDUCTOR_CACHE_DIR"), add = TRUE)
  }
  proc <- .cf_basilisk_start(.cateomics_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  tryCatch(
    basilisk::basiliskRun(proc, .cf_py_rng_guarded, .cf_fun = fun,
                          .cf_pypath = system.file("python", package = "CATEomics"),
                          ...),
    error = function(e) {
      doubles <- .cf_whole_doubles(list(...))
      if (!length(doubles)) stop(e)
      rlang::abort(
        c(conditionMessage(e),
          i = sprintf("Passed as a double: %s.", paste(doubles, collapse = "; ")),
          i = paste("A hyperparameter the Python side requires to be whole",
                    "has to be an R integer: write 10L, not 10.")),
        class = "cateomics_python_param_type", parent = e)
    })
}

#' Whole-number doubles among the hyperparameter blocks of a Python call
#'
#' R writes `10` as a double, and the Python libraries check the type of an
#' integer hyperparameter, so the natural R literal is refused. The refusal
#' reads as an upstream type error (scikit-learn) or, for a keyword that
#' reaches torch, as an error about a tensor, neither of which names the fix.
#' This collects the candidates so the failure can.
#'
#' Only lists are walked. The design matrix, the outcome and the treatment
#' arrive as a matrix and as vectors, so they are never scanned; the blocks a
#' caller supplies (`params`, the per-slot blocks, and the leftover keywords
#' the neural methods forward) are lists.
#'
#' @param args The arguments forwarded to the worker.
#' @return Character vector of `name = value` strings, possibly empty.
#' @keywords internal
#' @noRd
.cf_whole_doubles <- function(args) {
  out <- character(0)
  walk <- function(x, path) {
    if (is.list(x)) {
      nms <- names(x) %||% rep("", length(x))
      for (i in seq_along(x)) {
        p <- if (nzchar(nms[i])) paste(c(path, nms[i]), collapse = "$") else path
        walk(x[[i]], p)
      }
      return(invisible(NULL))
    }
    # A vector mixing whole and fractional values is a float argument whatever
    # its whole entries look like, so it is left alone.
    if (nzchar(path) && is.double(x) && length(x) >= 1L && length(x) <= 10L &&
        all(is.finite(x)) && all(x == trunc(x))) {
      out <<- c(out, sprintf("%s = %s", path,
                             paste(format(x, trim = TRUE), collapse = ", ")))
    }
    invisible(NULL)
  }
  for (a in args) if (is.list(a)) walk(a, character(0))
  unique(out)
}

#' Evaluate a worker body with the Python RNG state restored afterwards.
#'
#' basilisk runs the body in the calling R process under its shared mode, so
#' the seeding each fitter does from `random_state` would otherwise persist in
#' the user's own Python session. This wrapper takes a snapshot of the numpy /
#' random / torch generators, runs the body, and puts the snapshot back,
#' mirroring the `on.exit()` save-restore the R-side methods already do.
#' Placed at the single entry point rather than in each backend so a new
#' backend cannot forget it.
#'
#' @param .cf_fun The worker body.
#' @param .cf_pypath Path to the bundled Python modules.
#' @param ... Arguments forwarded to `.cf_fun`.
#' @return The value returned by `.cf_fun`.
#' @keywords internal
#' @noRd
.cf_py_rng_guarded <- function(.cf_fun, .cf_pypath, ...) {
  guard <- NULL
  if (is.character(.cf_pypath) && length(.cf_pypath) == 1L && nzchar(.cf_pypath)) {
    guard <- tryCatch(
      reticulate::import_from_path("rng_guard", path = .cf_pypath,
                                   convert = FALSE),
      error = function(e) NULL
    )
  }
  if (!is.null(guard)) {
    state <- tryCatch(guard$snapshot(), error = function(e) NULL)
    if (!is.null(state)) {
      on.exit(try(guard$restore(state), silent = TRUE), add = TRUE)
    }
  }
  # One deprecation of a pinned pair of dependencies is filtered here, at the
  # same single entry point; `inst/python/warning_filters.py` says which and
  # why. Every other Python-side warning still reaches the console.
  if (is.character(.cf_pypath) && length(.cf_pypath) == 1L && nzchar(.cf_pypath)) {
    try(reticulate::import_from_path("warning_filters", path = .cf_pypath,
                                     convert = FALSE)$apply(), silent = TRUE)
  }
  .cf_fun(...)
}

#' Start a basilisk process, robust to concurrent PSOCK-port collisions.
#'
#' basilisk falls back to \emph{socket} mode (\code{parallel::makePSOCKcluster})
#' whenever a foreign Python is already initialised in the calling process
#' (\code{reticulate::py_available()} is TRUE and the loaded interpreter is not
#' the CATEomics basilisk environment). The PSOCK port is \code{parallel}'s
#' default, which is computed once when \code{parallel} loads and is then
#' \emph{inherited unchanged by every fork}. So when \code{cf_fit} / cf_compare
#' Python fits run under fork-based parallelism
#' (\code{BiocParallel::MulticoreParam}, \code{parallel::mclapply},
#' \code{future::multicore}), all workers try to open the same port, all but one
#' fail with "server socket ... cannot be opened", and the affected per-method
#' fits are silently dropped. \code{grf} / \code{bcf} are R-native and use no
#' socket, so they are unaffected.
#'
#' This helper gives each process a distinct default PSOCK port (restoring the
#' caller's default afterwards) and retries on the residual collision. In fork
#' mode no socket is opened, so it is a harmless no-op there.
#'
#' @param env A \code{BasiliskEnvironment}.
#' @param attempts Maximum port attempts before giving up.
#' @return The basilisk process handle from \code{basiliskStart}.
#' @keywords internal
#' @noRd
.cf_basilisk_start <- function(env, attempts = 10L) {
  set_opt <- utils::getFromNamespace("setDefaultClusterOptions", "parallel")
  get_opt <- utils::getFromNamespace("getClusterOption", "parallel")
  old_port <- tryCatch(get_opt("port"), error = function(e) NULL)
  restore  <- function() {
    if (!is.null(old_port)) try(set_opt(port = old_port), silent = TRUE)
  }
  # Process-unique base port; each retry shifts it so two processes that collide
  # (same pid modulus) diverge deterministically without touching the RNG.
  base_port <- 11000L + (Sys.getpid() %% 40000L)
  last_err  <- NULL
  for (i in seq_len(attempts)) {
    try(set_opt(port = base_port + (i - 1L) * 131L), silent = TRUE)
    proc <- tryCatch(basilisk::basiliskStart(env, testload = NULL),
                     error = function(e) e)
    if (!inherits(proc, "error")) {
      restore()
      return(proc)
    }
    last_err <- proc
    # Only a socket/port contention error is worth retrying; anything else
    # (env build failure, etc.) should surface immediately.
    if (!grepl("socket|port", conditionMessage(proc), ignore.case = TRUE)) {
      restore()
      stop(proc)
    }
  }
  restore()
  stop(last_err)
}

#' Directory where fitted Python models are cached.
#'
#' Uses \code{tools::R_user_dir()} (never the working directory or HOME root) so
#' serialised models persist across sessions without polluting the user's
#' space.
#'
#' @keywords internal
#' @noRd
.cf_model_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("CATEomics", which = "cache"), "py_models")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Allocate a fresh path for a serialised Python model.
#'
#' Uses \code{tempfile()} (an internal counter, not the R RNG) so allocating a
#' path does not perturb the caller's \code{.Random.seed}.
#' @keywords internal
#' @noRd
.cf_new_model_path <- function(method, ext) {
  tempfile(pattern = paste0(method, "_"),
           tmpdir  = .cf_model_cache_dir(),
           fileext = paste0(".", ext))
}
