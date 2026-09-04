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
# Because basilisk runs Python in a separate process, a fitted Python model
# cannot be handed back to R as a live object. Each method therefore serialises
# its fitted model to a file during fit (joblib for sklearn/econml,
# torch.save for the NN backends) and reloads it at predict time. Only
# R-serialisable values (numeric result vectors + the model path + metadata)
# cross the process boundary.
#
# Bleeding-edge versions (numpy 2.x / pandas 3.x / torch 2.12) are not on
# conda-forge, so the environment is a conda `python` base plus pip-pinned
# packages; this pins the versions the package was developed and tested
# against.
# =============================================================================

#' Pinned pip packages for the cfomics unified environment.
#' @keywords internal
#' @noRd
.cfomics_basilisk_pip <- c(
  "numpy==2.4.6", "scipy==1.15.3", "scikit-learn==1.6.1", "pandas==3.0.3",
  "statsmodels==0.14.6", "Cython==3.2.4",
  "torch==2.12.0",
  "econml==0.16.0", "catenets==0.2.4", "pyro-ppl==1.9.1",
  "ganite==0.1.2",
  "shap==0.48.0", "xgboost==2.1.4"
)

#' The cfomics unified basilisk environment.
#'
#' Defined at load time (a cheap S4 object); the heavy conda environment is
#' built lazily by basilisk on the first \code{basiliskStart()}.
#' @keywords internal
#' @noRd
.cfomics_env <- basilisk::BasiliskEnvironment(
  envname  = "cfomics_unified",
  pkgname  = "cfomics",
  packages = "python=3.12.3",
  pip      = .cfomics_basilisk_pip
)

#' Run a function inside the cfomics basilisk environment.
#'
#' Starts (and stops) a basilisk process bound to the cfomics unified conda
#' environment and evaluates \code{fun(...)} inside it. \code{fun} may call
#' \code{reticulate::import()} freely; only its return value (which must be
#' R-serialisable) crosses back. This is the single entry point every Python
#' backend uses instead of importing modules in the main R process.
#'
#' @param fun A function to evaluate in the basilisk process.
#' @param ... Arguments forwarded to \code{fun} (must be R-serialisable).
#' @return The value returned by \code{fun}.
#' @keywords internal
#' @noRd
.cf_py_run <- function(fun, ...) {
  # Keep the torch/inductor compile cache inside the session tempdir (inherited
  # by the basilisk worker process) so it is cleaned up with the session. This
  # avoids the R CMD check "detritus in the temp directory" NOTE from a
  # persistent torchinductor_<user> directory in the system temp.
  if (!nzchar(Sys.getenv("TORCHINDUCTOR_CACHE_DIR"))) {
    ti_dir <- file.path(tempdir(), "cfomics_torchinductor")
    dir.create(ti_dir, showWarnings = FALSE, recursive = TRUE)
    Sys.setenv(TORCHINDUCTOR_CACHE_DIR = ti_dir)
    on.exit(Sys.unsetenv("TORCHINDUCTOR_CACHE_DIR"), add = TRUE)
  }
  proc <- .cf_basilisk_start(.cfomics_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  basilisk::basiliskRun(proc, fun, ...)
}

#' Start a basilisk process, robust to concurrent PSOCK-port collisions.
#'
#' basilisk falls back to \emph{socket} mode (\code{parallel::makePSOCKcluster})
#' whenever a foreign Python is already initialised in the calling process
#' (\code{reticulate::py_available()} is TRUE and the loaded interpreter is not
#' the cfomics basilisk environment). The PSOCK port is \code{parallel}'s
#' default, which is computed once when \code{parallel} loads and is then
#' \emph{inherited unchanged by every fork}. So when \code{cf_fit} / cf_compare
#' Python fits run under fork-based parallelism
#' (\code{BiocParallel::MulticoreParam}, \code{parallel::mclapply},
#' \code{future::multicore}), all workers try to open the same port, all but one
#' fail with "server socket ... cannot be opened", and the affected per-method
#' fits are silently dropped. (\code{grf} / \code{bcf} are R-native and use no
#' socket, so they are immune -- which is why the failure looks method-specific.)
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
#' serialised models persist across sessions without polluting the user's space,
#' matching the Bioconductor convention for on-disk caches.
#'
#' @keywords internal
#' @noRd
.cf_model_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("cfomics", which = "cache"), "py_models")
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
