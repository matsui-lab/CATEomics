# utils_bootstrap_ci.R
#
# `bootstrap_ci` opt-in common machinery for CATE-capable methods.
# Shared by s/t/x-learner and drlearner (non-linear final): bootstrap
# percentile CIs for tau when native CI is unavailable.
#
# The wrapper accepts a `fit_fn(X_b, T_b, Y_b, random_state = ..., ...)`
# closure that returns a cfomics_result-like list with `$res$ite` (the backend
# tau field) and `$res$ate`, plus a `predict_fn(model, newdata)` closure for
# per-row tau evaluation on the original sample. The percentile method yields
# `res$ate_ci` and `res$tau_ci_lower` / `res$tau_ci_upper`.


#' Validate `bootstrap_ci` request for a method
#'
#' Centralised policy enforcement (hybrid policy: accept / warn / error).
#' Returns the canonical `ci_method` label string the caller should record
#' under `meta$estimator$ci_method`.
#'
#' @param method Character. cfomics method name.
#' @param bootstrap_ci Logical. User's request.
#' @param native_label Character. The `ci_method` label that would apply
#'   when `bootstrap_ci = FALSE` (e.g., `"asymptotic Wald"`,
#'   `"Bayesian posterior"`, `"EIF asymptotic"`, `"sandwich SE"`,
#'   `"unavailable"`).
#' @return Character. Either `native_label` (when `bootstrap_ci = FALSE`)
#'   or `"nonparametric bootstrap"`. Errors when the request is structurally
#'   unsupported.
#' @keywords internal
#' @noRd
.resolve_ci_method <- function(method, bootstrap_ci, native_label) {
  if (!isTRUE(bootstrap_ci)) {
    return(native_label)
  }
  nn_methods <- c("tarnet", "dragonnet", "ganite")
  if (method %in% nn_methods) {
    rlang::abort(
      sprintf(
        "Bootstrap CI is not supported for method '%s' (prohibitive compute cost for this backend).",
        method
      ),
      class = "cfomics_bootstrap_unsupported_nn"
    )
  }
  has_native <- !identical(native_label, "unavailable") &&
                !identical(native_label, "nonparametric bootstrap")
  if (has_native) {
    rlang::warn(
      sprintf(
        "Native CI (%s) is faster; bootstrap is for research comparison.",
        native_label
      ),
      class = "cfomics_bootstrap_overrides_native"
    )
  }
  "nonparametric bootstrap"
}
