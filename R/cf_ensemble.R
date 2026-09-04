#' Combine multiple CATE estimators into a tau-ensemble.
#'
#' \code{cf_ensemble()} is the third leg of the comparison / selection /
#' integration triplet (alongside \code{\link{cf_compare}} and
#' \code{\link{cf_select}}). It consumes a \code{cfomics_compare} object and
#' returns a CATE estimator that combines the per-method OOF tau hats
#' without refitting any nuisance models.
#'
#' Five combination rules are available, all powered by the cached
#' \code{oof_tau_matrix} and \code{nuisance_oof} on the input object:
#'
#' \describe{
#'   \item{\code{"mean"}}{Equal weights \eqn{w_k = 1/K}. The simplest
#'     baseline ensemble.}
#'   \item{\code{"softmax"}}{Score-based softmax weights
#'     \eqn{w_k = \exp(\eta s_k) / \sum_j \exp(\eta s_j)} with
#'     \eqn{s_k = -\text{loss}_k}. Reduces to mean as \code{eta -> 0}
#'     and to a single-best selection as \code{eta -> Inf}.}
#'   \item{\code{"r_stacking"}}{Non-negative least-squares (NNLS) stacking
#'     against the R-loss residuals, solved as a small QP under
#'     \eqn{w \ge 0} only (no sum-to-one constraint). This matches the
#'     R-loss "R-Stacking" formulation (Nie & Wager R-loss objective): the
#'     residual-on-residual regression \eqn{Yres \approx \sum_k w_k\, Tres\,
#'     \tau_k} sets the overall scale itself, so forcing \eqn{\sum w = 1}
#'     would over-constrain it (and block legitimate shrinkage). The convex
#'     \eqn{\sum w = 1} constraint belongs instead to the DR pseudo-outcome
#'     objective (Han & Wu causal stacking); see \code{dr_stacking}.}
#'   \item{\code{"dr_stacking"}}{Convex stacking (\eqn{w \ge 0, \sum w = 1})
#'     targeting the DR pseudo-outcome \eqn{\hat\tau^{PO}}, i.e. Han & Wu
#'     (2022) causal stacking \eqn{\arg\min_w \sum_i (\hat\tau^{PO}_i -
#'     \sum_k w_k \tau_k(x_i))^2}. The convex constraint keeps the ensemble
#'     calibrated to the \eqn{\tau} scale.}
#'   \item{\code{"cba"}}{Consensus-Based Averaging (Machluf et al. 2026
#'     \emph{arXiv:2407.03690}). Ranks methods by mean off-diagonal
#'     correlation (the correlation-knee helper) and averages the consensus subset above the
#'     kneedle knee.}
#' }
#'
#' Refit cost is zero in all cases: every method's full-data and OOF
#' artifacts already live on the input \code{cfomics_compare} object.
#'
#' @param compare_result A \code{cfomics_compare} object from
#'   \code{cf_compare()}. Single \code{cf_model} inputs are rejected for
#'   consistency with \code{cf_select()}.
#' @param method One of \code{"mean"}, \code{"softmax"}, \code{"r_stacking"},
#'   \code{"dr_stacking"}, \code{"cba"}.
#' @param loss Loss column used by \code{method = "softmax"} only; one of
#'   \code{"r_loss"}, \code{"dr_loss"}, \code{"factual_loss"} (pulled from the
#'   cf_compare scoreboard to weight methods). \code{r_stacking} and
#'   \code{dr_stacking} ignore \code{loss} — they hard-code their objective
#'   (R-loss / DR pseudo-outcome respectively), and \code{mean} / \code{cba}
#'   do not use a loss at all.
#' @param eta Softmax temperature; only used when \code{method = "softmax"}.
#'   \code{eta = 0} reduces to equal weights, \code{eta -> Inf} concentrates
#'   on the lowest-loss method.
#' @param cor_method Correlation method for \code{method = "cba"}; one of
#'   \code{"kendall"} (default), \code{"spearman"}, \code{"pearson"}.
#' @param min_cba Minimum number of methods to retain in the cba consensus
#'   subset (clamps the kneedle output from below).
#' @param ... Currently unused; reserved for future extensions.
#'
#' @return An object of class
#'   \code{c("cfomics_ensemble", "cf_model")} with the
#'   minimum contract required by downstream functions:
#'   \itemize{
#'     \item \code{method} = \code{"ensemble"} (so \code{cf_evaluate} /
#'           \code{cf_benefit_groups} treat it as a CATE-capable result);
#'     \item \code{tau_hat_oof} = the combined OOF tau hat;
#'     \item \code{nuisance} = the cached \code{nuisance_oof};
#'     \item \code{data} = list with \code{Y}, \code{T}, \code{X};
#'     \item \code{meta} = list with \code{partition = "oof"},
#'           \code{ensemble_method}, \code{loss}, \code{eta},
#'           \code{covariate_names}, \code{n}, \code{p}, \code{methods_used};
#'     \item \code{fit$res} = list with \code{ite}, \code{ate},
#'           \code{summary}, \code{mu0_hat}, \code{mu1_hat} (Robinson
#'           reconstruction from nuisance + ensemble tau);
#'     \item \code{weights} = named numeric vector of length \eqn{K}
#'           giving the combination weights (zeros for methods not
#'           included by cba).
#'   }
#'
#' @section Downstream:
#' The minimum contract above is enough for
#' \code{cf_evaluate(model, ...)}, \code{cf_benefit_groups(model, ...)},
#' \code{cf_rule(model, ...)}, and \code{cf_interpret(model, ...)} to run
#' against the ensemble. The combination weights are stored on the result.
#'
#' @section References:
#' Machluf, Frostig, Milo, Carlin, Berkman & Pryluk (2026) Robust CATE
#' estimation using novel ensemble methods (CBA), \emph{J. Biopharmaceutical
#' Statistics}, doi:10.1080/10543406.2026.2668583.
#' Han & Wu (2022) Causal stacking, \emph{arXiv:2202.12445}.
#' Nie & Wager (2021) Quasi-Oracle Estimation of HTE, \emph{Biometrika}.
#'
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf), then combine their
#' # per-patient effects. Runs without Python.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   n_burn = 10L, n_iter = 10L, n_chains = 1L,  # 1 chain keeps bcf fast
#'                   propensity_score = 0.5, random_state = 1)
#' ens <- cf_ensemble(cmp, method = "mean")
#' ens
#'
#' @export
cf_ensemble <- function(compare_result,
                        method     = c("mean", "softmax",
                                       "r_stacking", "dr_stacking", "cba"),
                        loss       = c("r_loss", "dr_loss", "factual_loss"),
                        eta        = 1.0,
                        cor_method = c("kendall", "spearman", "pearson"),
                        min_cba    = 2L,
                        ...) {
  method     <- match.arg(method)
  loss       <- match.arg(loss)
  cor_method <- match.arg(cor_method)
  eta        <- as.numeric(eta)[1L]
  min_cba    <- as.integer(min_cba)[1L]

  if (!inherits(compare_result, "cfomics_compare")) {
    if (inherits(compare_result, "cf_model")) {
      rlang::abort(
        paste0("cf_ensemble() requires a cfomics_compare object ",
               "(from cf_compare), not a single cf_model."),
        class = "cfomics_ensemble_single_model"
      )
    }
    rlang::abort(
      "compare_result must be a cfomics_compare object.",
      class = "cfomics_invalid_compare_result"
    )
  }

  # ----- Step 1: clean oof_tau_matrix via the OOF-matrix cleaner, abort on degenerate input.
  cleaned <- .cfomics_clean_oof_tau_matrix(
    compare_result$metadata$oof_tau_matrix
  )
  if (cleaned$status %in% c("no_input", "no_columns", "no_rows")) {
    rlang::abort(
      sprintf(paste0("cf_ensemble: oof_tau_matrix is unusable (status = '%s'). ",
                     "Need at least one surviving method and one complete row."),
              cleaned$status),
      class = "cfomics_ensemble_empty_input"
    )
  }
  if (cleaned$status == "single_column") {
    rlang::warn(
      paste0("cf_ensemble: only one method survived cleaning ('",
             cleaned$kept_methods, "'); ensemble reduces to that method."),
      class = "cfomics_ensemble_single_method"
    )
  }

  if (length(cleaned$dropped_methods) > 0L) {
    rlang::inform(
      sprintf("cf_ensemble: dropped %d failed method(s): %s",
              length(cleaned$dropped_methods),
              paste(cleaned$dropped_methods, collapse = ", ")),
      class = "cfomics_ensemble_dropped_methods"
    )
  }

  # We work on the *original* oof_tau_matrix (not the cleaning matrix,
  # which dropped rows). Stacking / softmax need per-obs alignment with Y,
  # T, and nuisance, all of which are length n. Reasoning: rows dropped by
  # complete.cases are only those where SOME surviving method has NA at
  # that row; in cf_compare's failed-method semantics this should be rare
  # (failed methods are whole columns of NA, not isolated rows). We still
  # honour the row count for cba (which uses the cleaned matrix directly).
  oof_tau_full <- compare_result$metadata$oof_tau_matrix
  kept_methods <- cleaned$kept_methods
  oof_tau_kept <- oof_tau_full[, kept_methods, drop = FALSE]

  # ----- Step 2: compute weights per method dispatch.
  weights_kept <- switch(method,
    mean       = .cfomics_ensemble_weights_mean(kept_methods),
    softmax    = .cfomics_ensemble_weights_softmax(compare_result, kept_methods,
                                                    loss, eta),
    r_stacking = .cfomics_ensemble_weights_stacking(compare_result, kept_methods,
                                                    target = "r_loss"),
    dr_stacking = .cfomics_ensemble_weights_stacking(compare_result, kept_methods,
                                                     target = "dr_loss"),
    cba        = .cfomics_ensemble_weights_cba(cleaned,
                                                cor_method = cor_method,
                                                min_cba    = min_cba)
  )

  # ----- Step 3: assemble combined OOF tau hat. NA-safe weighted sum.
  oof_ite <- .cfomics_weighted_row_sum(oof_tau_kept, weights_kept$weights)

  # Re-expand weights to the original method set so consumers see a stable
  # length-K vector (failed methods get weight 0).
  all_methods <- colnames(oof_tau_full)
  if (is.null(all_methods)) all_methods <- as.character(seq_len(ncol(oof_tau_full)))
  weights_full <- setNames(rep(0, length(all_methods)), all_methods)
  weights_full[kept_methods] <- weights_kept$weights

  # ----- Step 4: extract cached cache for downstream fit$res + nuisance.
  cache <- .cfomics_extract_compare_cache(compare_result)

  # Reconstruct the ensemble's per-arm counterfactual means from the marginal
  # outcome regression E[Y|X] and the ensemble tau.
  mu_bar <- cache$m_oof
  mu0_hat <- mu_bar - cache$pi_oof * oof_ite
  mu1_hat <- mu_bar + (1 - cache$pi_oof) * oof_ite

  ate_value <- mean(oof_ite, na.rm = TRUE)

  # ----- Step 5: build output object (cf_model + cfomics_ensemble class).
  cov_names <- compare_result$metadata$parsed$covariate_names
  if (is.null(cov_names)) {
    cov_names <- colnames(cache$oof_tau_matrix)  # fallback (rare)
  }
  X_data <- NULL
  for (m in names(compare_result$models)) {
    mdl <- compare_result$models[[m]]
    if (!is.null(mdl$data$X)) {
      X_data <- mdl$data$X
      break
    }
  }

  out <- structure(
    list(
      method   = "ensemble",
      # Four-pillar: ensembled tau + OOF tau + shared nuisance live in res.
      res = list(
        tau_hat     = oof_ite,
        tau_hat_oof = oof_ite,
        ate         = ate_value,
        mu0_hat     = mu0_hat,
        mu1_hat     = mu1_hat,
        nuisance    = compare_result$metadata$nuisance_oof,
        summary     = .cf_build_summary(ate_value, ate_ci = NULL, tau_hat = oof_ite)
      ),
      model    = NULL,   # no single backend for an ensemble
      data     = list(Y = cache$Y, T = cache$T, X = X_data),
      meta     = list(
        partition         = "oof",
        ensemble_method   = method,
        loss              = loss,
        eta               = if (method == "softmax") eta else NA_real_,
        cor_method        = if (method == "cba") cor_method else NA_character_,
        methods_used      = kept_methods,
        methods_dropped   = cleaned$dropped_methods,
        consensus_members = weights_kept$consensus_members,
        covariate_names   = cov_names,
        n                 = length(cache$Y),
        p                 = if (!is.null(X_data)) ncol(X_data) else NA_integer_,
        n_folds           = compare_result$metadata$n_folds,
        formula           = compare_result$metadata$formula,
        random_state      = compare_result$metadata$random_state
      ),
      weights = weights_full
    ),
    class = c("cfomics_ensemble", "cf_model")
  )

  out
}


# -------------------------------------------------------------------------
# Weight solvers
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_ensemble_weights_mean <- function(kept_methods) {
  K <- length(kept_methods)
  w <- rep(1 / K, K)
  names(w) <- kept_methods
  list(weights = w, consensus_members = kept_methods)
}


#' @keywords internal
#' @noRd
.cfomics_ensemble_weights_softmax <- function(compare_result, kept_methods,
                                              loss, eta) {
  # Pull scalar per-method loss; prefer scoreboard if the column exists, else
  # recompute via the per-observation loss helper.
  losses <- .cfomics_get_method_losses(compare_result, kept_methods, loss)
  scores <- -losses
  if (!is.finite(eta) || eta <= 0) {
    w <- rep(1 / length(kept_methods), length(kept_methods))
  } else {
    # Numerically stable softmax.
    z <- eta * (scores - max(scores, na.rm = TRUE))
    e <- exp(z)
    w <- e / sum(e)
  }
  names(w) <- kept_methods
  list(weights = w, consensus_members = kept_methods)
}


#' @keywords internal
#' @noRd
.cfomics_ensemble_weights_stacking <- function(compare_result, kept_methods,
                                                target = c("r_loss", "dr_loss")) {
  target <- match.arg(target)
  mats   <- .cfomics_compute_stacking_materials(compare_result)
  tau_mat <- mats$oof_tau_matrix[, kept_methods, drop = FALSE]
  K       <- ncol(tau_mat)

  # Build (design, target) per the loss type.
  if (target == "r_loss") {
    # Linear-in-w: (Yres - sum_k w_k Tres * tau_k)^2
    design <- tau_mat * mats$Tres   # n x K, columnwise multiply by Tres
    target_vec <- mats$Yres
  } else {
    # Linear-in-w: (sum_k w_k tau_k - Y_DR)^2
    design     <- tau_mat
    target_vec <- mats$Y_DR
  }

  good <- stats::complete.cases(design, target_vec) &
          is.finite(target_vec) &
          apply(is.finite(design), 1, all)
  design     <- design[good, , drop = FALSE]
  target_vec <- target_vec[good]

  if (nrow(design) < K) {
    rlang::warn(
      paste0("cf_ensemble (", target, "_stacking): only ", nrow(design),
             " usable rows for K = ", K,
             " methods; falling back to equal weights."),
      class = "cfomics_ensemble_stacking_insufficient_rows"
    )
    w <- rep(1 / K, K)
    names(w) <- kept_methods
    return(list(weights = w, consensus_members = kept_methods))
  }

  # Set up QP.
  # min  1/2 w' D w - d' w
  # s.t. A' w = 1 (equality), w_k >= 0 (inequalities)
  D <- crossprod(design)
  # Regularise diagonal for numerical PD (small jitter).
  reg <- 1e-8 * max(diag(D), 1)
  diag(D) <- diag(D) + reg
  d <- as.numeric(crossprod(design, target_vec))

  # Constraint depends on the stacking objective:
  #   r_loss  -> NNLS:   w >= 0 only (no sum-to-one). The R-loss residual-on-
  #              residual regression sets its own scale; sum-to-one would
  #              over-constrain it. = Nie-Wager R-loss "R-Stacking".
  #   dr_loss -> convex: sum(w) = 1 AND w >= 0. = Han & Wu (2022) causal
  #              stacking on the DR pseudo-outcome (keeps tau-scale calibrated).
  # solve.QP expects Amat in K-by-(meq + nineq) form; meq leading columns are
  # equality constraints.
  convex <- (target == "dr_loss")
  if (convex) {
    Amat <- cbind(rep(1, K), diag(K))
    bvec <- c(1, rep(0, K))
    meq  <- 1L
  } else {
    Amat <- diag(K)        # w_k >= 0 only (NNLS)
    bvec <- rep(0, K)
    meq  <- 0L
  }

  qp <- tryCatch(
    quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(qp)) {
    rlang::warn(
      paste0("cf_ensemble (", target,
             "_stacking): QP solver failed; falling back to equal weights."),
      class = "cfomics_ensemble_stacking_qp_failed"
    )
    w <- rep(1 / K, K)
  } else {
    w <- qp$solution
    w[w < 0] <- 0           # clean tiny negative numerical dust
    if (convex) {
      # convex: renormalise to the simplex (cleanup after clipping).
      if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / K, K)
    } else {
      # NNLS: keep the free-scale solution; only guard the all-zero case.
      if (sum(w) <= 0) w <- rep(1 / K, K)
    }
  }
  names(w) <- kept_methods
  list(weights = w, consensus_members = kept_methods)
}


#' @keywords internal
#' @noRd
.cfomics_ensemble_weights_cba <- function(cleaned,
                                          cor_method = "kendall",
                                          min_cba    = 2L) {
  ck <- .cfomics_tau_correlation_knee(
    mat        = cleaned$mat,
    cor_method = cor_method,
    min_knee   = min_cba
  )
  K <- length(cleaned$kept_methods)
  if (is.na(ck$knee_k) || K < 2L) {
    w <- rep(1 / K, K)
    names(w) <- cleaned$kept_methods
    return(list(weights = w, consensus_members = cleaned$kept_methods))
  }
  ordering         <- ck$ordering
  knee_k           <- ck$knee_k
  consensus_idx    <- ordering[seq_len(knee_k)]
  consensus_names  <- cleaned$kept_methods[consensus_idx]

  w <- rep(0, K)
  names(w) <- cleaned$kept_methods
  w[consensus_names] <- 1 / knee_k
  list(weights = w, consensus_members = consensus_names)
}


# -------------------------------------------------------------------------
# Loss column dispatch (softmax: prefer scoreboard, else per-observation-loss recompute)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_get_method_losses <- function(compare_result, kept_methods, loss) {
  sb <- compare_result$scoreboard
  loss_col <- loss
  if (is.data.frame(sb) && loss_col %in% names(sb)) {
    # Identity column: `label` (method_base) in grid mode, else `method`.
    method_col <- if ("label" %in% names(sb)) sb$label else
                  if ("method" %in% names(sb)) sb$method else rownames(sb)
    losses <- vapply(kept_methods, function(m) {
      idx <- which(method_col == m)
      if (length(idx) == 0L) return(NA_real_)
      val <- sb[[loss_col]][idx[1L]]
      if (is.null(val) || !is.finite(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
    if (all(is.finite(losses))) return(losses)
  }
  # Fallback: recompute scalar loss per method via the per-observation loss helper.
  vapply(kept_methods, function(m) {
    pof <- .cfomics_compute_per_obs_loss_from_compare(compare_result, m, loss)
    mean(pof, na.rm = TRUE)
  }, numeric(1))
}


# -------------------------------------------------------------------------
# NA-safe weighted row sum
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cfomics_weighted_row_sum <- function(mat, weights) {
  # apply the weights on their native scale. r_stacking uses NNLS
  # weights whose sum is NOT constrained to 1 (R-loss "R-Stacking",
  # Nie & Wager 2021; the free scale is intrinsic to the R-loss residual
  # regression); the previous
  # per-row renormalisation `w / sum(w)` forced sum(w) = 1 and silently
  # cancelled that free scale. Complete rows now use the raw weighted sum;
  # only rows with NA/dropped columns are rescaled, and to the intended total
  # weight (preserving the NNLS scale) rather than to 1.
  if (is.null(mat) || ncol(mat) == 0L) return(numeric(0))
  weights <- as.numeric(weights)
  if (length(weights) != ncol(mat)) {
    rlang::abort("weights length must equal ncol(mat).",
                 class = "cfomics_weight_length_mismatch")
  }
  n <- nrow(mat)
  out <- numeric(n)
  pos <- weights > 0
  total_w <- sum(weights[pos])
  for (i in seq_len(n)) {
    row_i <- mat[i, ]
    ok <- is.finite(row_i) & pos
    if (!any(ok)) {
      out[i] <- NA_real_
    } else if (all(ok == pos)) {
      # All positive-weight columns are finite: native-scale weighted sum.
      out[i] <- sum(row_i[ok] * weights[ok])
    } else {
      # Some positive-weight columns are NA: rescale the available weights so
      # the effective weight-sum stays at `total_w` (keeps the NNLS scale).
      out[i] <- sum(row_i[ok] * weights[ok]) * (total_w / sum(weights[ok]))
    }
  }
  out
}


# -------------------------------------------------------------------------
# S3 print / summary
# -------------------------------------------------------------------------

#' @export
print.cfomics_ensemble <- function(x, ...) {
  cat("<cfomics_ensemble>: method = ", x$meta$ensemble_method,
      " (n = ", x$meta$n, ", K = ", length(x$meta$methods_used), ")\n",
      sep = "")
  cat("  ATE = ", format(x$res$ate, digits = 4),
      ", tau summary: mean = ", format(x$res$summary$tau_mean, digits = 4),
      ", sd = ", format(x$res$summary$tau_sd, digits = 4), "\n", sep = "")
  cat("  Methods used (", length(x$meta$methods_used), "): ",
      paste(x$meta$methods_used, collapse = ", "), "\n", sep = "")
  if (length(x$meta$methods_dropped) > 0L) {
    cat("  Methods dropped (", length(x$meta$methods_dropped), "): ",
        paste(x$meta$methods_dropped, collapse = ", "), "\n", sep = "")
  }
  if (x$meta$ensemble_method == "softmax") {
    cat("  loss = ", x$meta$loss, ", eta = ", x$meta$eta, "\n", sep = "")
  }
  if (x$meta$ensemble_method == "cba") {
    cat("  cor_method = ", x$meta$cor_method, ", consensus_members = ",
        length(x$meta$consensus_members), "\n", sep = "")
  }
  w <- x$weights[x$weights > 0]
  if (length(w) > 0L) {
    cat("  Top weights:\n")
    ord <- order(w, decreasing = TRUE)
    for (i in seq_len(min(length(ord), 5L))) {
      idx <- ord[i]
      cat("    ", names(w)[idx], " = ", format(w[idx], digits = 3), "\n",
          sep = "")
    }
  }
  invisible(x)
}


#' @export
summary.cfomics_ensemble <- function(object, ...) {
  res <- list(
    method            = object$meta$ensemble_method,
    n                 = object$meta$n,
    methods_used      = object$meta$methods_used,
    methods_dropped   = object$meta$methods_dropped,
    weights           = object$weights,
    ate               = object$res$ate,
    tau_summary       = object$res$summary,
    consensus_members = object$meta$consensus_members,
    loss              = object$meta$loss,
    eta               = object$meta$eta,
    cor_method        = object$meta$cor_method
  )
  class(res) <- "summary.cfomics_ensemble"
  res
}


#' @export
print.summary.cfomics_ensemble <- function(x, ...) {
  cat("cf_ensemble summary\n")
  cat("  ensemble_method: ", x$method, "\n", sep = "")
  cat("  n: ", x$n, ", K: ", length(x$methods_used), "\n", sep = "")
  cat("  methods_used: ", paste(x$methods_used, collapse = ", "), "\n",
      sep = "")
  if (length(x$methods_dropped) > 0L) {
    cat("  methods_dropped: ", paste(x$methods_dropped, collapse = ", "),
        "\n", sep = "")
  }
  cat("  ATE: ", format(x$ate, digits = 4), "\n", sep = "")
  cat("  tau mean / sd: ", format(x$tau_summary$tau_mean, digits = 4),
      " / ", format(x$tau_summary$tau_sd, digits = 4), "\n", sep = "")
  cat("  weights:\n")
  for (nm in names(x$weights)) {
    cat("    ", nm, " = ", format(x$weights[nm], digits = 3), "\n", sep = "")
  }
  invisible(x)
}
