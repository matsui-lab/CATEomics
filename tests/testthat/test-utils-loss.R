# Smoke tests for the per-observation loss wrappers.
# Verify cache extraction + stacking materials + per-obs loss dispatch on a
# mock cfomics_compare object. No nuisance refit involved.

set.seed(456)
.make_mock_compare <- function(n = 30L, K = 3L,
                                include_models = TRUE,
                                include_nuisance = TRUE,
                                include_oof = TRUE) {
  Y <- rnorm(n, sd = 1)
  T_vec <- rbinom(n, 1L, prob = 0.5)
  mu0 <- rnorm(n, sd = 0.5)
  mu1 <- mu0 + 0.5 + rnorm(n, sd = 0.2)
  pi  <- pmin(pmax(plogis(rnorm(n, sd = 0.3)), 1e-4), 1 - 1e-4)

  oof_tau_matrix <- matrix(rnorm(n * K), nrow = n, ncol = K)
  colnames(oof_tau_matrix) <- paste0("m", seq_len(K))

  models <- if (include_models) {
    out <- list()
    for (k in seq_len(K)) {
      mdl <- list(data = list(Y = Y, T = T_vec, X = matrix(0, nrow = n, ncol = 2)))
      out[[paste0("m", k)]] <- mdl
    }
    out
  } else {
    list()
  }

  metadata <- list(
    nuisance_oof   = if (include_nuisance) list(mu0_oof = mu0, mu1_oof = mu1,
                                                 pi_oof = pi,
                                                 m_oof = rnorm(n)) else NULL,
    oof_tau_matrix = if (include_oof) oof_tau_matrix else NULL
  )

  structure(
    list(scoreboard = data.frame(method = colnames(oof_tau_matrix)),
         models     = models,
         metadata   = metadata),
    class = c("cfomics_compare", "list")
  )
}


# .cfomics_extract_compare_cache --------------------------------------------

test_that(".cfomics_extract_compare_cache returns Y, T, nuisance, oof", {
  cmp <- .make_mock_compare(n = 25L, K = 3L)
  cache <- .cfomics_extract_compare_cache(cmp)

  expect_named(cache, c("Y", "T", "mu0_oof", "mu1_oof", "pi_oof", "m_oof",
                        "oof_tau_matrix", "methods"))
  expect_equal(length(cache$Y), 25L)
  expect_equal(length(cache$T), 25L)
  expect_equal(length(cache$mu0_oof), 25L)
  expect_equal(dim(cache$oof_tau_matrix), c(25L, 3L))
  expect_setequal(cache$methods, c("m1", "m2", "m3"))
})

test_that(".cfomics_extract_compare_cache rejects non-cfomics_compare", {
  expect_error(
    .cfomics_extract_compare_cache(list(scoreboard = NULL)),
    class = "cfomics_invalid_compare_result"
  )
})

test_that(".cfomics_extract_compare_cache errors when nuisance_oof missing", {
  cmp <- .make_mock_compare(include_nuisance = FALSE)
  expect_error(
    .cfomics_extract_compare_cache(cmp),
    class = "cfomics_missing_nuisance_oof"
  )
})

test_that(".cfomics_extract_compare_cache errors when oof_tau_matrix missing", {
  cmp <- .make_mock_compare(include_oof = FALSE)
  expect_error(
    .cfomics_extract_compare_cache(cmp),
    class = "cfomics_missing_oof_tau_matrix"
  )
})

test_that(".cfomics_extract_compare_cache errors when no surviving model has Y/T", {
  cmp <- .make_mock_compare(include_models = FALSE)
  expect_error(
    .cfomics_extract_compare_cache(cmp),
    class = "cfomics_compare_cache_empty"
  )
})


# .cfomics_compute_stacking_materials ----------------------------------------

test_that(".cfomics_compute_stacking_materials returns Yres, Tres, Y_DR", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  mats <- .cfomics_compute_stacking_materials(cmp)

  expect_true(all(c("Yres", "Tres", "Y_DR") %in% names(mats)))
  expect_equal(length(mats$Yres), 30L)
  expect_equal(length(mats$Tres), 30L)
  expect_equal(length(mats$Y_DR), 30L)
  # Yres = Y - m_oof (marginal outcome regression) by construction.
  expect_true(all(abs((mats$Y - mats$m_oof) - mats$Yres) < 1e-12))
  # Tres = T - pi by construction.
  expect_true(all(abs(as.numeric(mats$T) - mats$pi_oof - mats$Tres) < 1e-12))
})


# .cfomics_compute_per_obs_loss_from_compare ---------------------------------

test_that("per-obs loss wrapper dispatches r_loss / dr_loss / factual_loss", {
  cmp <- .make_mock_compare(n = 30L, K = 2L)

  r  <- .cfomics_compute_per_obs_loss_from_compare(cmp, "m1", "r_loss")
  dr <- .cfomics_compute_per_obs_loss_from_compare(cmp, "m1", "dr_loss")
  fl <- .cfomics_compute_per_obs_loss_from_compare(cmp, "m1", "factual_loss")

  expect_equal(length(r), 30L)
  expect_equal(length(dr), 30L)
  expect_equal(length(fl), 30L)
  expect_true(all(is.finite(r)))
  expect_true(all(is.finite(dr)))
  expect_true(all(is.finite(fl)))
  # All per-obs losses are squared residuals -> non-negative.
  expect_true(all(r >= 0))
  expect_true(all(dr >= 0))
  expect_true(all(fl >= 0))
})

test_that("per-obs loss wrapper errors on unknown method", {
  cmp <- .make_mock_compare(n = 30L, K = 2L)
  expect_error(
    .cfomics_compute_per_obs_loss_from_compare(cmp, "missing", "r_loss"),
    class = "cfomics_method_not_in_compare"
  )
})


# .cfomics_compute_per_obs_loss_for_tau --------------------------------------

test_that("per-obs loss for tau works with raw vector input", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  tau_ens <- rowMeans(cmp$metadata$oof_tau_matrix)

  r_ens <- .cfomics_compute_per_obs_loss_for_tau(cmp, tau_ens, "r_loss")
  expect_equal(length(r_ens), 30L)
  expect_true(all(is.finite(r_ens)))
})

test_that("per-obs loss for tau errors on length mismatch", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  expect_error(
    .cfomics_compute_per_obs_loss_for_tau(cmp, rep(0, 10L), "r_loss"),
    class = "cfomics_tau_length_mismatch"
  )
})


# Sanity equivalence: r_loss(tau_m) from method dispatch == from tau vector
test_that("dispatch-by-method and dispatch-by-vector agree", {
  cmp <- .make_mock_compare(n = 30L, K = 3L)
  tau_m1 <- cmp$metadata$oof_tau_matrix[, "m1"]
  r_by_method <- .cfomics_compute_per_obs_loss_from_compare(cmp, "m1", "r_loss")
  r_by_tau    <- .cfomics_compute_per_obs_loss_for_tau(cmp, tau_m1, "r_loss")
  expect_true(all(abs(r_by_method - r_by_tau) < 1e-12))
})


# theta* R-loss baseline ------------------------------------------------------

test_that("theta* is the R-loss-minimising constant-effect model", {
  set.seed(5)
  n <- 200L
  m_oof  <- stats::rnorm(n)
  pi_oof <- rep(0.5, n)
  T      <- stats::rbinom(n, 1, 0.5)
  Y      <- m_oof + 0.4 * (T - 0.5) + stats::rnorm(n, sd = 0.5)
  Yres <- Y - m_oof
  Tres <- T - pi_oof
  theta_star <- sum(Yres * Tres) / sum(Tres^2)

  loss <- function(th) mean(.compute_r_loss_per_obs(Y, T, m_oof, pi_oof, rep(th, n)))
  # no other constant effect beats theta* under the R-loss
  expect_lte(loss(theta_star), loss(theta_star + 0.5))
  expect_lte(loss(theta_star), loss(theta_star - 0.5))
  expect_lte(loss(theta_star), loss(0))
})


# .compute_r_squared_score: numerator/denominator over the same finite mask ----

test_that(".compute_r_squared_score is a no-op when per_obs is fully finite", {
  set.seed(3)
  per  <- runif(50, 0.1, 0.5)
  base <- runif(50, 0.3, 0.9)
  got  <- cfomics:::.compute_r_squared_score(per, base)
  expect_equal(got, 1 - mean(per) / mean(base))
})

test_that(".compute_r_squared_score masks base to the finite per_obs subset", {
  per  <- c(0.2, NA, 0.4, NA, 0.6)
  base <- c(0.5, 0.5, 0.5, 0.5, 0.5)
  ok   <- is.finite(per)
  # denominator must average base over the SAME finite subset, not all 5 obs.
  expect_equal(cfomics:::.compute_r_squared_score(per, base),
               1 - mean(per[ok]) / mean(base[ok]))
})

test_that(".compute_r_squared_score returns NA when no finite overlap", {
  expect_true(is.na(cfomics:::.compute_r_squared_score(rep(NA_real_, 5), 1:5)))
})
