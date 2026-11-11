#' Combine multiple CATE estimators into a tau-ensemble.
#'
#' \code{cf_ensemble()} is the combination step of the comparison /
#' integration triplet (alongside \code{\link{cf_compare}} and
#' \code{\link{cf_select}}). It consumes a \code{cateomics_compare} object and
#' returns a CATE estimator that combines the per-method OOF tau hats
#' without refitting any nuisance models.
#'
#' Five combination rules are available, all powered by the cached
#' \code{oof_tau_matrix} and \code{nuisance_oof} on the input object:
#'
#' \describe{
#'   \item{\code{"mean"}}{Equal weights \eqn{w_k = 1/K}.}
#'   \item{\code{"softmax"}}{Score-based softmax weights
#'     \eqn{w_k = \exp(\eta s_k) / \sum_j \exp(\eta s_j)} with
#'     \eqn{s_k = -\text{loss}_k}. Reduces to mean as \code{eta -> 0}
#'     and to a single-best selection as \code{eta -> Inf}.}
#'   \item{\code{"r_stacking"}}{Non-negative least squares against the R-loss
#'     residuals, \eqn{Yres \approx \sum_k w_k\, Tres\, \tau_k} under
#'     \eqn{w \ge 0}, solved as a small QP. The loss is the R-loss of
#'     Nie & Wager (2021); minimising it over non-negative weights, with no
#'     sum-to-one constraint, is what Han & Wu (2022) write as their
#'     eq. (4). Nie & Wager's own stacking step also fits an intercept and a
#'     constant treatment-effect term; this does not.}
#'   \item{\code{"dr_stacking"}}{Convex stacking (\eqn{w \ge 0, \sum w = 1})
#'     targeting the DR pseudo-outcome \eqn{\hat\tau^{PO}}, i.e. the causal
#'     stacking Han & Wu (2022) write as their eq. (2),
#'     \eqn{\arg\min_w \sum_i (\hat\tau^{PO}_i -
#'     \sum_k w_k \tau_k(x_i))^2}.}
#'   \item{\code{"cba"}}{Consensus-Based Averaging (Machluf et al. 2024).
#'     Ranks methods by mean off-diagonal correlation and averages the
#'     consensus subset above the knee of that curve.}
#' }
#'
#' No method is refitted: the weights are computed from the OOF artifacts
#' already on the input \code{cateomics_compare} object.
#'
#' @param compare_result A \code{cateomics_compare} object from
#'   \code{cf_compare()}. A single \code{cf_model} is rejected.
#' @param method One of \code{"mean"} (the default), \code{"softmax"},
#'   \code{"r_stacking"}, \code{"dr_stacking"}, \code{"cba"}. The two
#'   stacking rules require more scored rows than they have free weights (one
#'   per kept method, less one for \code{dr_stacking}, whose weights sum to
#'   one) and abort otherwise; they fall back to equal weights, with a
#'   warning, when the quadratic program fails.
#' @param loss Loss column used by \code{method = "softmax"} only; one of
#'   \code{"r_loss"} (the default), \code{"dr_loss"},
#'   \code{"factual_loss"} (pulled from the
#'   cf_compare scoreboard to weight methods). \code{r_stacking} and
#'   \code{dr_stacking} ignore \code{loss} — they hard-code their objective
#'   (R-loss / DR pseudo-outcome respectively), and \code{mean} / \code{cba}
#'   do not use a loss at all.
#' @param eta Softmax temperature (default \code{1}); only used when
#'   \code{method = "softmax"}. Any value that is not positive and finite,
#'   including \code{0} and negatives, gives equal weights;
#'   \code{eta -> Inf} concentrates on the lowest-loss method.
#' @param cor_method Correlation method for \code{method = "cba"}; one of
#'   \code{"kendall"} (default), \code{"spearman"}, \code{"pearson"}.
#' @param min_cba Minimum number of methods to retain in the cba consensus
#'   subset (default \code{2}; clamps the kneedle output from below).
#' @param ... Not used; naming an unknown argument is an error
#'   (\code{cateomics_unknown_argument}).
#'
#' @return An object of class
#'   \code{c("cateomics_ensemble", "cf_model")} with the
#'   minimum contract required by downstream functions:
#'   \itemize{
#'     \item \code{method} = \code{"ensemble"} (so \code{cf_evaluate} /
#'           \code{cf_benefit_groups} treat it as a CATE-capable result);
#'     \item \code{res} = list with \code{tau_hat_oof} (the combined OOF
#'           tau hat), \code{tau_hat = NULL} (an ensemble has no backend of
#'           its own), \code{fold_id} (the shared split that produced the
#'           out-of-fold estimates), \code{ate}, \code{mu0_hat},
#'           \code{mu1_hat}
#'           (Robinson reconstruction from the nuisance and the ensemble tau),
#'           \code{nuisance} (the cached \code{nuisance_oof}) and
#'           \code{summary};
#'     \item \code{model} = the members carrying non-zero weight, which is
#'           what \code{predict(object, newdata = )} scores;
#'     \item \code{data} = list with \code{Y}, \code{T},
#'           \code{X_preprocessed};
#'     \item \code{meta} = the same pillar a \code{\link{cf_fit}} result
#'           carries (see \code{cateomics_result}), with
#'           \code{partition = "oof"} and the ensemble's own fields:
#'           \code{ensemble_method}, \code{loss}, \code{eta},
#'           \code{cor_method}, \code{methods_used},
#'           \code{methods_dropped}, \code{consensus_members},
#'           and \code{n_rows_rescaled} (rows combined from fewer members than
#'           the weights name, because a member's out-of-fold column was
#'           missing there). The fold count and the split are in
#'           \code{estimator}, the seed in \code{provenance$rng_seed}. The
#'           slots the newdata path reads on a single fit are \emph{absent}
#'           here: each member carries its own preprocessing record and
#'           \code{predict(object, newdata = )} delegates to the members, so
#'           the ensemble has nothing of its own to record there;
#'     \item \code{weights} = named numeric vector of length \eqn{K}
#'           giving the combination weights (zero for configs dropped during
#'           cleaning, and for the ones cba leaves out of the consensus
#'           subset).
#'   }
#'
#'   With a binary outcome \code{mu0_hat} / \code{mu1_hat} can fall outside
#'   \eqn{[0, 1]}, with every combination method, since they are reconstructed
#'   from the nuisance and the combined effect rather than predicted directly.
#'
#' @section Downstream:
#' The minimum contract above is enough for
#' \code{cf_evaluate(model, ...)}, \code{cf_benefit_groups(model, ...)},
#' \code{cf_rule(model, ...)}, and \code{cf_interpret(model, ...)} to run
#' against the ensemble. The combination weights are stored on the result.
#'
#' @references
#' Machluf, Frostig, Shoham, Milo, Berkman & Pryluk (2024) Robust CATE
#' Estimation Using Novel Ensemble Methods (CBA), \emph{arXiv:2407.03690}.
#'
#' Han & Wu (2022) Ensemble Method for Estimating Individualized Treatment
#' Effects, \emph{arXiv:2202.12445}.
#'
#' Nie & Wager (2021) Quasi-Oracle Estimation of Heterogeneous Treatment
#' Effects, \emph{Biometrika} 108(2):299-319.
#'
#' @examples
#' # Two R-native CATE methods (grf + a short-MCMC bcf), then combine their
#' # per-patient effects. Runs without Python.
#' se  <- cf_example_data(n = 50, p = 40, n_prog = 8, n_pred = 4, seed = 1)
#' cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
#'                   X_baseline = c("age", "sex"),
#'                   methods = c("grf", "bcf"), n_features = 25,
#'                   # 1 chain keeps bcf fast
#'                   params = list(bcf = list(nburn = 10L, nsim = 10L,
#'                                            n_chains = 1L)),
#'                   propensity_score = 0.5, eval_propensity_score = 0.5,
#'                   random_state = 1)
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
  .cf_reject_dots(..., .fn = "cf_ensemble()")
  method     <- match.arg(method)
  loss       <- match.arg(loss)
  cor_method <- match.arg(cor_method)
  eta        <- as.numeric(eta)[1L]
  min_cba    <- as.integer(min_cba)[1L]

  if (!inherits(compare_result, "cateomics_compare")) {
    if (inherits(compare_result, "cf_model")) {
      rlang::abort(
        paste0("cf_ensemble() requires a cateomics_compare object ",
               "(from cf_compare), not a single cf_model."),
        class = "cateomics_ensemble_single_model"
      )
    }
    rlang::abort(
      "compare_result must be a cateomics_compare object.",
      class = "cateomics_invalid_compare_result"
    )
  }

  # ----- Step 1: clean oof_tau_matrix via the OOF-matrix cleaner, abort on degenerate input.
  cleaned <- .cateomics_clean_oof_tau_matrix(
    compare_result$metadata$oof_tau_matrix
  )
  if (cleaned$status %in% c("no_input", "no_columns", "no_rows")) {
    rlang::abort(
      sprintf(paste0("cf_ensemble: oof_tau_matrix is unusable (status = '%s'). ",
                     "Need at least one surviving method and one complete row."),
              cleaned$status),
      class = "cateomics_ensemble_empty_input"
    )
  }
  if (cleaned$status == "single_column") {
    rlang::warn(
      paste0("cf_ensemble: only one method survived cleaning ('",
             cleaned$kept_methods, "'); ensemble reduces to that method."),
      class = "cateomics_ensemble_single_method"
    )
  }

  if (length(cleaned$dropped_methods) > 0L) {
    rlang::inform(
      sprintf("cf_ensemble: dropped %d failed method(s): %s",
              length(cleaned$dropped_methods),
              paste(cleaned$dropped_methods, collapse = ", ")),
      class = "cateomics_ensemble_dropped_methods"
    )
  }

  # We work on the *original* oof_tau_matrix (not the cleaning matrix,
  # which dropped rows). Stacking / softmax need per-obs alignment with Y,
  # T, and nuisance, all of which are length n. cba uses the cleaned matrix
  # directly, so it honours the reduced row count.
  oof_tau_full <- compare_result$metadata$oof_tau_matrix
  kept_methods <- cleaned$kept_methods
  oof_tau_kept <- oof_tau_full[, kept_methods, drop = FALSE]

  # ----- Step 2: compute weights per method dispatch.
  weights_kept <- switch(method,
    mean       = .cateomics_ensemble_weights_mean(kept_methods),
    softmax    = .cateomics_ensemble_weights_softmax(compare_result, kept_methods,
                                                    loss, eta),
    r_stacking = .cateomics_ensemble_weights_stacking(compare_result, kept_methods,
                                                    target = "r_loss"),
    dr_stacking = .cateomics_ensemble_weights_stacking(compare_result, kept_methods,
                                                     target = "dr_loss"),
    cba        = .cateomics_ensemble_weights_cba(cleaned,
                                                cor_method = cor_method,
                                                min_cba    = min_cba)
  )

  # ----- Step 3: assemble combined OOF tau hat. NA-safe weighted sum.
  oof_ite <- .cateomics_weighted_row_sum(oof_tau_kept, weights_kept$weights)
  # Record how many rows stood in for a smaller member set, then return a bare
  # vector so the stored estimate carries no attributes.
  n_rows_rescaled <- attr(oof_ite, "n_rows_rescaled") %||% 0L
  oof_ite <- as.numeric(oof_ite)

  # Re-expand weights to the original method set so consumers see a stable
  # length-K vector (failed methods get weight 0).
  all_methods <- colnames(oof_tau_full)
  if (is.null(all_methods)) all_methods <- as.character(seq_len(ncol(oof_tau_full)))
  weights_full <- setNames(rep(0, length(all_methods)), all_methods)
  weights_full[kept_methods] <- weights_kept$weights
  if (length(.cf_ensemble_weighted_members(weights_full)) == 0L) {
    # r_stacking solves a non-negative regression with no constraint on the
    # sum, so "no candidate helps" is an available answer. It arrives either as
    # exact zeros or, more often, as the solver's own zero: every weight pinned
    # at its lower bound, of the order 1e-16, which is below the floor that
    # counts a weight as a share. Either way the combined estimate is zero for
    # every sample, which reads as a finding unless it is said out loud.
    rlang::warn(
      c("No candidate carries weight: the combined effect is zero for every sample.",
        i = "This is what the weighting found, not an estimate of no effect."),
      class = "cateomics_ensemble_zero_weights"
    )
  }

  # ----- Step 4: extract cached cache for downstream fit$res + nuisance.
  cache <- .cateomics_extract_compare_cache(compare_result)

  # Reconstruct the ensemble's per-arm counterfactual means from the marginal
  # outcome regression E[Y|X] and the ensemble tau.
  mu_bar <- cache$m_oof
  mu0_hat <- mu_bar - cache$pi_oof * oof_ite
  mu1_hat <- mu_bar + (1 - cache$pi_oof) * oof_ite

  ate_value <- mean(oof_ite, na.rm = TRUE)

  # ----- Step 5: build output object (cf_model + cateomics_ensemble class).
  # The columns of the OOF tau matrix are method labels, not covariates, so
  # there is no fallback here: if the comparison recorded no covariate names
  # the field stays NULL.
  cov_names <- compare_result$metadata$covariate_names
  X_data <- NULL
  for (m in names(compare_result$models)) {
    mdl <- compare_result$models[[m]]
    X_candidate <- mdl$data$X_preprocessed %||% mdl$data$X
    if (!is.null(X_candidate)) {
      X_data <- X_candidate
      break
    }
  }

  out <- structure(
    list(
      method   = "ensemble",
      # Four-pillar: ensembled tau + OOF tau + shared nuisance live in res.
      res = list(
        # An ensemble has no backend to re-predict with, so there is no
        # in-sample tau-hat: the only estimate it owns is the combined OOF one.
        tau_hat     = NULL,
        tau_hat_oof = oof_ite,
        # The combined estimate is a weighted sum over the members, which all
        # share one split, so the split that produced it is that same one.
        fold_id     = compare_result$metadata$fold_id,
        ate         = ate_value,
        mu0_hat     = mu0_hat,
        mu1_hat     = mu1_hat,
        nuisance    = compare_result$metadata$nuisance_oof,
        summary     = .cf_build_summary(ate_value, ate_ci = NULL, tau_hat = oof_ite)
      ),
      # An ensemble has no backend of its own. What it needs to score new
      # observations is the members it actually weights: each one carries its
      # own preprocessing record and backend, so newdata is reconstructed by
      # the member and this object holds no second copy of that record.
      model    = list(
        kind    = "ensemble",
        members = compare_result$models[
          intersect(.cf_ensemble_weighted_members(weights_full),
                    names(compare_result$models))]
      ),
      data     = list(Y = cache$Y, T = cache$T, X_preprocessed = X_data),
      meta     = .cf_ensemble_meta(
        compare_result = compare_result,
        cov_names = cov_names, n = length(cache$Y),
        p = if (!is.null(X_data)) ncol(X_data) else NA_integer_,
        extra = list(
          ensemble_method   = method,
          loss              = loss,
          eta               = if (method == "softmax") eta else NA_real_,
          cor_method        = if (method == "cba") cor_method else NA_character_,
          methods_used      = kept_methods,
          methods_dropped   = cleaned$dropped_methods,
          consensus_members = weights_kept$consensus_members,
          n_rows_rescaled   = n_rows_rescaled
        )
      ),
      weights = weights_full
    ),
    class = c("cateomics_ensemble", "cf_model")
  )

  out
}


# -------------------------------------------------------------------------
# Weight solvers
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.cateomics_ensemble_weights_mean <- function(kept_methods) {
  K <- length(kept_methods)
  w <- rep(1 / K, K)
  names(w) <- kept_methods
  list(weights = w, consensus_members = kept_methods)
}


#' @keywords internal
#' @noRd
.cateomics_ensemble_weights_softmax <- function(compare_result, kept_methods,
                                              loss, eta) {
  # Pull scalar per-method loss; prefer scoreboard if the column exists, else
  # recompute via the per-observation loss helper.
  losses <- .cateomics_get_method_losses(compare_result, kept_methods, loss)
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
.cateomics_ensemble_weights_stacking <- function(compare_result, kept_methods,
                                                target = c("r_loss", "dr_loss")) {
  target <- match.arg(target)
  mats   <- .cateomics_compute_stacking_materials(compare_result)
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

  if (nrow(design) == 0L) {
    rlang::abort(
      c(sprintf("cf_ensemble (%s_stacking): no row has a usable value for every method.",
                target),
        i = "The weights are solved on the rows where all of them predict and the target is finite.",
        i = "cf_compare() records per-method failures in metadata$failed_methods."),
      class = "cateomics_ensemble_stacking_no_rows"
    )
  }
  # The weights are a regression of a target on the K candidate columns, so
  # they need more rows than free parameters. `r_stacking` solves K free
  # weights; `dr_stacking` constrains them to sum to one, which leaves K - 1.
  # This is the identification condition the Super Learner applies to its own
  # library (`cateomics_sl_library_too_large`); the fold structure stays honest
  # below it, so nothing else catches it.
  n_free <- if (identical(target, "dr_loss")) K - 1L else K
  if (nrow(design) <= n_free) {
    rlang::abort(
      c(sprintf("cf_ensemble (%s_stacking): the weights are solved on %d row(s) over %d free weight(s).",
                target, nrow(design), n_free),
        i = "There is no residual degree of freedom unless the rows outnumber the free weights.",
        i = sprintf("Compare fewer methods, or use method = \"mean\", which needs no weights.")),
      class = "cateomics_ensemble_library_too_large"
    )
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
  #              over-constrain it. The loss is the R-loss of Nie & Wager; the
  #              stacking form fitted here, with no intercept and no constant
  #              treatment-effect term, is Han & Wu (2022) eq. 4.
  #   dr_loss -> convex: sum(w) = 1 AND w >= 0. = Han & Wu (2022) eq. 2, causal
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
      class = "cateomics_ensemble_stacking_qp_failed"
    )
    w <- rep(1 / K, K)
  } else {
    w <- qp$solution
    w[w < 0] <- 0           # clean tiny negative numerical dust
    if (convex) {
      # convex: renormalise to the simplex (cleanup after clipping).
      if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / K, K)
    }
    # NNLS keeps the free-scale solution as it came back. "No candidate helps"
    # is an available answer when the weights are only required to be
    # non-negative, and replacing it with equal weights would put an answer in
    # its place that the data did not give -- and would hide it, because the
    # report below only speaks when no weight clears the floor.
  }
  names(w) <- kept_methods
  list(weights = w, consensus_members = kept_methods)
}


#' @keywords internal
#' @noRd
.cateomics_ensemble_weights_cba <- function(cleaned,
                                          cor_method = "kendall",
                                          min_cba    = 2L) {
  ck <- .cateomics_tau_correlation_knee(
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
.cateomics_get_method_losses <- function(compare_result, kept_methods, loss) {
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
    pof <- .cateomics_compute_per_obs_loss_from_compare(compare_result, m, loss)
    mean(pof, na.rm = TRUE)
  }, numeric(1))
}


# -------------------------------------------------------------------------
# NA-safe weighted row sum
# -------------------------------------------------------------------------

#' The members a combination actually leans on
#'
#' A solver leaves residue: a non-negative least squares can return a weight of
#' order 1e-18 for a method it did not use. Carrying it would score that
#' backend on new data for nothing, and if it were the only member that managed
#' to score, the per-row renormalisation would hand it the whole prediction. So
#' the cut is relative to the largest weight rather than at zero.
#'
#' @param w The full weight vector.
#' @return The names whose weight is more than dust.
#' @keywords internal
#' @noRd
.cf_ensemble_weighted_members <- function(w) {
  w <- w[is.finite(w)]
  if (!length(w)) return(character(0))
  # An absolute floor as well as a relative one. A solver that pins every
  # weight at its lower bound leaves values of the order 1e-12, and a tolerance
  # taken from the largest of those would call one of them a member.
  if (max(w, 0) <= .CF_ENSEMBLE_WEIGHT_FLOOR) return(character(0))
  tol <- max(w, 0) * sqrt(.Machine$double.eps)
  names(w)[w > tol]
}


# Below this a weight is the solver's zero, not a share of the combination.
.CF_ENSEMBLE_WEIGHT_FLOOR <- sqrt(.Machine$double.eps)


#' Score new observations with a weighted combination of the ensemble members
#'
#' Each member reconstructs its own design matrix from its own preprocessing
#' record, so `newdata` is passed through untouched and this function never
#' preprocesses it a second time.
#'
#' The potential outcomes are not returned. In-sample they are recovered from
#' the shared evaluation nuisance by solving `m = pi * mu1 + (1 - pi) * mu0`
#' together with `tau = mu1 - mu0`, which keeps `mu1 - mu0 == tau` exact. That
#' nuisance is stored as out-of-fold prediction vectors, not as fitted
#' learners, so there is no `m` for rows it was not fitted on. Weighting the
#' members' own `mu` instead would break the identity, because three of the
#' methods (grf, xlearner, drlearner) do not satisfy `mu1 - mu0 == tau`
#' individually.
#'
#' @keywords internal
#' @noRd
.cf_ensemble_predict_newdata <- function(object, newdata, type) {
  members <- object$model$members
  w       <- object$weights
  if (is.null(members) || !length(members)) {
    # Two different reasons land here, and only one of them is about how the
    # object was built: an ensemble whose weighting found nothing has no member
    # to delegate to either, and telling that caller to rebuild it would point
    # at the wrong thing.
    if (!is.null(w) && length(w) &&
        length(.cf_ensemble_weighted_members(w)) == 0L) {
      rlang::abort(
        c("predict(): no candidate carried weight, so there is no member to score newdata with.",
          i = "This is what the weighting found, not an estimate of no effect: the combined effect is zero for every sample.",
          i = "Read the candidates in the cf_compare() scoreboard, or combine them another way (see `method` in cf_ensemble())."),
        class = "cateomics_ensemble_zero_weights_newdata"
      )
    }
    rlang::abort(
      c("predict(): this ensemble carries no member models.",
        i = "Rebuild it with cf_ensemble() on a cf_compare() result."),
      class = "cateomics_ensemble_no_members"
    )
  }
  use <- names(members)
  tau_list <- vector("list", length(use))
  why      <- list()
  for (j in seq_along(use)) {
    p <- tryCatch(
      as.numeric(stats::predict(members[[use[j]]], newdata = newdata,
                                type = "tau")),
      error = function(e) {why[[use[j]]] <<- conditionMessage(e); NULL}
    )
    if (length(p)) tau_list[[j]] <- p
  }
  len <- vapply(tau_list, function(z) if (is.null(z)) 0L else length(z),
                integer(1))
  got    <- len > 0L
  failed <- use[!got]
  if (!any(got)) {
    rlang::abort(
      c("predict(): no ensemble member could score newdata.",
        i = paste(vapply(use, function(k)
          sprintf("%s: %s", k, why[[k]] %||% "returned no values"),
          character(1)), collapse = " | "),
        i = "Check that newdata carries the columns the members were fitted on."),
      class = "cateomics_ensemble_newdata_failed"
    )
  }
  # The members score the same rows, so a disagreement on how many there are
  # means at least one is wrong and nothing here says which. Picking the
  # commonest or the longest would make that choice silently.
  n_seen <- unique(len[got])
  if (length(n_seen) > 1L) {
    rlang::abort(
      c("predict(): the ensemble members disagree on how many rows newdata has.",
        i = paste(sprintf("%s: %d", use[got], len[got]), collapse = ", "),
        i = "Score the members individually to find which one is wrong."),
      class = "cateomics_ensemble_newdata_length_mismatch"
    )
  }
  n_new <- n_seen
  mat <- matrix(NA_real_, nrow = n_new, ncol = length(use),
                dimnames = list(NULL, use))
  for (j in which(got)) mat[, j] <- tau_list[[j]]
  if (length(failed)) {
    rlang::warn(
      c(sprintf("predict(): %d of %d ensemble member(s) could not score newdata.",
                length(failed), length(use)),
        i = paste(vapply(failed, function(k)
          sprintf("%s: %s", k, why[[k]] %||% "returned no values"),
          character(1)), collapse = " | "),
        i = "Their weight is spread over the members that did predict."),
      class = "cateomics_ensemble_member_newdata_failed"
    )
  }
  # One rule for both entry points: the in-sample assembler renormalises the
  # available weights per row, which here rescales away the dropped members.
  tau <- as.numeric(.cateomics_weighted_row_sum(mat, w[use]))
  ate <- mean(tau, na.rm = TRUE)
  na_n <- rep(NA_real_, n_new)
  switch(
    type,
    "all" = list(
      ate     = ate,
      # NULL for an interval that was not computed; NA for a per-row quantity
      # that exists but has no value here.
      ate_ci  = NULL,
      tau_hat = tau,
      tau_ci  = NULL,
      mu0_hat = na_n,
      mu1_hat = na_n,
      ps      = na_n,
      ps_type = NA_character_
    ),
    "ate"     = ate,
    "ate_ci"  = NA_real_,
    "tau"     = tau,
    "tau_ci"  = NA_real_,
    "mu0"     = na_n,
    "mu1"     = na_n,
    "ps"      = na_n,
    "summary" = rlang::abort(
      c("predict(): `type = \"summary\"` is not available with newdata.",
        i = "Call predict(object, type = \"summary\") for the stored summary, or ask for `tau`."),
      class = "cateomics_unsupported_predict_type"),
    stop(sprintf("Unsupported type '%s'.", type))
  )
}


#' @keywords internal
#' @noRd
.cateomics_weighted_row_sum <- function(mat, weights) {
  # Apply the weights on their native scale: r_stacking does not constrain
  # sum(w) to 1, so renormalising a complete row would cancel the scale the
  # R-loss residual regression chose. Only rows with NA / dropped columns are
  # rescaled, and to the intended total weight rather than to 1.
  if (is.null(mat) || ncol(mat) == 0L) return(numeric(0))
  weights <- as.numeric(weights)
  if (length(weights) != ncol(mat)) {
    rlang::abort("weights length must equal ncol(mat).",
                 class = "cateomics_weight_length_mismatch")
  }
  n <- nrow(mat)
  out <- numeric(n)
  pos <- weights > 0
  total_w <- sum(weights[pos])
  n_rescaled <- 0L
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
      n_rescaled <- n_rescaled + 1L
    }
  }
  # How many rows were combined from fewer members than the weights name. The
  # value is not recoverable from the returned vector, and a row standing in
  # for a different subset of members is not the same estimate as its
  # neighbours.
  attr(out, "n_rows_rescaled") <- n_rescaled
  out
}


# -------------------------------------------------------------------------
# S3 print / summary
# -------------------------------------------------------------------------

#' @export
print.cateomics_ensemble <- function(x, ...) {
  # `%||%` on each piece rather than one guard: a print says what it cannot
  # find rather than stopping part way, and an absent name is the thing a
  # reader of a malformed object most needs to see.
  cat("<cateomics_ensemble>: method = ",
      x$meta$ensemble_method %||% "unrecorded",
      " (n = ", x$meta$n %||% NA_integer_,
      ", K = ", length(x$meta$methods_used), ")\n",
      sep = "")
  # An ensemble owns no in-sample estimate, so what is summarised here is the
  # combined out-of-fold one; the other two screens of the family say which
  # estimate their numbers come from, and this one did not.
  cat("  ATE = ", format(x$res$ate, digits = 4),
      ", out-of-fold tau: mean = ", format(x$res$summary$tau_mean, digits = 4),
      ", sd = ", format(x$res$summary$tau_sd, digits = 4), "\n", sep = "")
  cat("  Methods used (", length(x$meta$methods_used), "): ",
      paste(x$meta$methods_used, collapse = ", "), "\n", sep = "")
  if (length(x$meta$methods_dropped) > 0L) {
    cat("  Methods dropped (", length(x$meta$methods_dropped), "): ",
        paste(x$meta$methods_dropped, collapse = ", "), "\n", sep = "")
  }
  if (identical(x$meta$ensemble_method, "softmax")) {
    cat("  loss = ", x$meta$loss, ", eta = ", x$meta$eta, "\n", sep = "")
  }
  if (identical(x$meta$ensemble_method, "cba")) {
    cat("  cor_method = ", x$meta$cor_method, ", consensus_members = ",
        length(x$meta$consensus_members), "\n", sep = "")
  }
  # How many of the candidates actually ended up in the combination: the line
  # above counts the ones it was solved over, and a member at zero is in that
  # count while contributing nothing and not being kept in `model`.
  carry <- .cf_ensemble_weighted_members(x$weights)
  if (length(carry) == 0L) {
    cat("  no member carries weight: the combination is zero for every sample\n")
  } else {
    w <- x$weights[carry]
    cat(sprintf("  weights (%d of %d members carry weight):\n",
                length(carry), length(x$weights)))
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
summary.cateomics_ensemble <- function(object, ...) {
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
  class(res) <- "summary.cateomics_ensemble"
  res
}


#' @export
print.summary.cateomics_ensemble <- function(x, ...) {
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
  cat("  out-of-fold tau mean / sd: ", format(x$tau_summary$tau_mean, digits = 4),
      " / ", format(x$tau_summary$tau_sd, digits = 4), "\n", sep = "")
  # Same rule as print: a weight the solver left at its lower bound is not a
  # share of the combination, and listing it as 7e-15 invites reading it as one.
  carry <- .cf_ensemble_weighted_members(x$weights)
  if (length(carry) == 0L) {
    cat("  weights: none; the combination is zero for every sample\n")
  } else {
    cat(sprintf("  weights (%d of %d members carry weight):\n",
                length(carry), length(x$weights)))
    for (nm in carry[order(x$weights[carry], decreasing = TRUE)]) {
      cat("    ", nm, " = ", format(x$weights[nm], digits = 3), "\n", sep = "")
    }
  }
  invisible(x)
}

#' Whether the samples behind a comparison carry meaningful names
#'
#' The member models all come from one design, so any of them answers for the
#' comparison. Returns `NA` when no member carries the record (a comparison
#' saved before it existed), which leaves the reader on its fallback.
#'
#' @param compare_result A `cateomics_compare` object.
#' @return `TRUE` / `FALSE`, or `NA` when no member records it.
#' @keywords internal
#' @noRd
.cateomics_compare_names_meaningful <- function(compare_result) {
  models <- compare_result$models
  if (!length(models)) return(NA)
  for (m in models) {
    v <- m$meta$provenance$sample_names_meaningful
    if (!is.null(v)) return(isTRUE(v))
  }
  NA
}
