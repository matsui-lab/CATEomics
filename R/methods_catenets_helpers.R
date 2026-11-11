# =============================================================================
# methods_catenets_helpers.R — shared Y z-score wrapper for tarnet / dragonnet
# =============================================================================
#
# catenets v0.2.x standardises neither X nor Y internally. X is z-scored
# centrally in preprocessing, so this shared helper (tarnet / dragonnet) only
# z-scores Y R-side and inverse-transforms the predictions.
#
# Helpers:
#   .catenets_scale_y(X, Y): compute the Y z-score, return Y-scaled data
#     (X is coerced to a double matrix and passed through unchanged)
#   .catenets_inverse_predictions(...): map catenets output back to the original
#     scale


#' Scale Y for catenets methods (X coerced and passed through)
#'
#' X is standardised centrally in `.cf_preprocess_input` (every column,
#' including binary indicators), so it is passed through unchanged here and
#' `x_zscore` is `NULL` --
#' cf_predict.R `.apply_newdata_preprocessing` re-applies the central train
#' stats to newdata, and the OOF held-out rows are already standardised, so
#' returning `NULL` avoids double-standardising. Only Y is standardised NN-side
#' (catenets v0.2.x standardises neither), and its predictions are
#' inverse-transformed via `y_zscore`. A constant Y falls back to sigma = 1.
#'
#' @param X Numeric matrix (n × p), already centrally z-scored.
#' @param Y Numeric vector (length n).
#' @return List with components:
#'   * X_scaled: X passed through unchanged (already standardised)
#'   * Y_scaled: standardised Y (numeric vector)
#'   * x_zscore: NULL (X standardisation handled centrally)
#'   * y_zscore: list(Y_mean = numeric(1), Y_sd = numeric(1))
#'
#' @keywords internal
#' @noRd
.catenets_scale_y <- function(X, Y) {
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "double"

  Y_vec  <- as.numeric(Y)
  Y_mean <- mean(Y_vec)
  Y_sd   <- stats::sd(Y_vec)
  if (!is.finite(Y_sd) || Y_sd == 0) Y_sd <- 1

  Y_scaled <- (Y_vec - Y_mean) / Y_sd

  list(
    X_scaled = X_mat,
    Y_scaled = Y_scaled,
    x_zscore = NULL,
    y_zscore = list(Y_mean = Y_mean, Y_sd = Y_sd)
  )
}


#' Inverse-transform catenets predictions
#'
#' catenets fit / predict on scaled (X, Y) returns mu0_hat / mu1_hat / ite on
#' the scaled Y space. Convert to the original Y scale: μ_t = scaled * sd + mean,
#' τ = scaled * sd (the mean term cancels in the difference).
#'
#' @param mu0_scaled,mu1_scaled,ite_scaled Numeric vectors from catenets output.
#' @param y_zscore list(Y_mean, Y_sd).
#'
#' @return list(mu0_hat, mu1_hat, ite, ate)
#'
#' @keywords internal
#' @noRd
.catenets_inverse_predictions <- function(mu0_scaled, mu1_scaled, ite_scaled,
                                          y_zscore) {
  y_mean <- y_zscore$Y_mean
  y_sd   <- y_zscore$Y_sd
  mu0_hat <- mu0_scaled * y_sd + y_mean
  mu1_hat <- mu1_scaled * y_sd + y_mean
  ite     <- ite_scaled * y_sd
  ate     <- mean(ite)
  list(mu0_hat = mu0_hat, mu1_hat = mu1_hat, ite = ite, ate = ate)
}
