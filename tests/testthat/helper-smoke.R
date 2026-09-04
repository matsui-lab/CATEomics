# Lightweight simulation helpers shared across smoke tests.

sim_simple_linear <- function(n = 200, p = 5, seed = 42L,
                              tau = 1.5) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  # Confounded treatment: T | X ~ Bern(plogis(0.5 X1)).
  T <- stats::rbinom(n, 1, stats::plogis(0.5 * X[, 1]))
  # Linear outcome with constant treatment effect tau.
  Y <- tau * T + 0.8 * X[, 1] + 0.4 * X[, 2] + stats::rnorm(n)

  data.frame(Y = Y, T = T, X, check.names = FALSE)
}
