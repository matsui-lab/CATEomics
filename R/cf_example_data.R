#' Simulate an example omics dataset for cfomics
#'
#' Generates a small semi-synthetic \code{SummarizedExperiment} of bulk RNA-seq
#' counts with a known heterogeneous treatment effect, for use in examples,
#' vignettes and quick experimentation. The design mimics the \eqn{n \ll p}
#' regime that cfomics targets: many genes, few samples, a randomised treatment,
#' and a continuous outcome whose conditional average treatment effect (CATE)
#' depends on a handful of genes and on a clinical effect-modifier.
#'
#' @details
#' The data are **fully simulated** and carry no biological meaning; they exist
#' only to make the workflow runnable end-to-end. The generating process is:
#' \itemize{
#'   \item **Expression.** Per-gene mean counts are drawn from a log-normal
#'     distribution (to give a realistic dynamic range) and counts are sampled
#'     from a negative-binomial distribution (\code{\link[stats]{rnbinom}}), the
#'     standard model for bulk RNA-seq (as in edgeR / DESeq2). Prognostic and
#'     predictive genes are given higher dispersion so they carry across-sample
#'     variation and survive default variance-based feature selection.
#'   \item **Signal scale.** The outcome model is defined on the *normalised*
#'     expression scale: per-gene log2-CPM values are z-scored internally
#'     (\code{z(logCPM)}), matching the log-CPM transform cfomics applies. The
#'     treatment effect is therefore a smooth function of normalised expression,
#'     not of raw counts.
#'   \item **Clinical covariates.** \code{age} (\eqn{\sim N(60, 10)}) is
#'     *prognostic* (enters the baseline \eqn{\mu_0}); \code{sex}
#'     (\code{"F"}/\code{"M"}) is an *effect-modifier* (shifts the treatment
#'     effect \eqn{\tau}), illustrating both a numeric covariate and a factor
#'     that is one-hot expanded downstream.
#'   \item **Treatment.** Randomised (\eqn{\Pr(T=1)=0.5}, no confounding), so
#'     the setting is an RCT.
#'   \item **Outcome.**
#'     \deqn{\mu_0(X) = v_{age}\, z(age) + \sum_{g \in \mathrm{prog}} v_g\, z(logCPM_g)}
#'     \deqn{\tau(X)  = ATE + \delta\, [sex = M] + \sum_{g \in \mathrm{pred}} w_g\, z(logCPM_g)}
#'     \deqn{Y = \mu_0(X) + T\, \tau(X) + \varepsilon, \quad \varepsilon \sim N(0, \sigma^2)}
#'     The prognostic (\eqn{\mu_0}) and predictive (\eqn{\tau}) signals are of
#'     comparable strength (\eqn{\mathrm{Var}(\mu_0) \approx \mathrm{Var}(\tau)}),
#'     and \eqn{\sigma} is set so that
#'     \eqn{\mathrm{Var}(\tau)/\mathrm{Var}(\varepsilon)} equals \code{snr}. A
#'     minority of individuals have \eqn{\tau < 0}, so treatment-rule and
#'     benefit-group analyses are non-trivial.
#' }
#' At the default size (\eqn{n=200 \ll p=500}) the heterogeneous effect is only
#' modestly recoverable by any single estimator -- illustrating exactly why
#' cfomics compares, selects and ensembles methods rather than trusting one.
#' Recovery quality is a property of this data-generating process, not of
#' cfomics, and is deliberately not emphasised in the vignette.
#' The predictive genes are members of a real MSigDB Hallmark pathway
#' (\code{HALLMARK_APOPTOSIS}); this only makes the biology-bridge recipe in the
#' vignette coherent (enrichment against that gene set) and is **not** a claim
#' of biological signal recovery — recovery accuracy is a property of the
#' data-generating process, not of cfomics.
#'
#' @param n Integer, number of samples (columns). Default 200.
#' @param p Integer, number of genes (rows). Default 500. Must exceed
#'   \code{n_prog + n_pred}.
#' @param n_prog Integer, number of prognostic genes (drive \eqn{\mu_0}).
#'   Default 20.
#' @param n_pred Integer, number of predictive genes (drive \eqn{\tau}).
#'   Default 10.
#' @param snr Numeric, target treatment-effect signal-to-noise ratio
#'   \eqn{\mathrm{Var}(\tau)/\mathrm{Var}(\varepsilon)}. Default 1.5.
#' @param seed Integer, random seed. The caller's RNG state is restored on exit.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with one
#'   \code{counts} assay (\code{p} genes x \code{n} samples, real gene symbols as
#'   row names) and:
#'   \describe{
#'     \item{\code{colData}}{\code{treatment} (0/1), \code{outcome} (numeric),
#'       \code{age} (numeric), \code{sex} (factor), and \code{cf_tau_true} (the
#'       true individual treatment effect).}
#'     \item{\code{rowData}}{\code{is_prognostic} and \code{is_predictive}
#'       logical flags marking the signal genes.}
#'     \item{\code{metadata}}{\code{cf_pathway} (the embedded pathway name),
#'       \code{cf_predictive_genes}, \code{cf_prognostic_genes},
#'       \code{cf_true_ate}, \code{cf_snr} (realised), and the call parameters.}
#'   }
#'
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- cf_example_data(n = 60, p = 200, seed = 1)
#'   se
#'   # ground truth is carried alongside the data
#'   mean(SummarizedExperiment::colData(se)$cf_tau_true)   # ~ true ATE
#'   table(SummarizedExperiment::rowData(se)$is_predictive)
#' }
#'
#' @seealso \code{\link{cf_fit}}, \code{\link{cf_compare}}, \code{\link{as_cf_data}}
#' @export
cf_example_data <- function(n = 200L, p = 500L,
                            n_prog = 20L, n_pred = 10L,
                            snr = 1.5, seed = 1L) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    rlang::abort(
      "cf_example_data() requires the 'SummarizedExperiment' and 'S4Vectors' packages.",
      class = "cfomics_missing_package_error"
    )
  }
  n      <- as.integer(n)
  p      <- as.integer(p)
  n_prog <- as.integer(n_prog)
  n_pred <- as.integer(n_pred)
  if (anyNA(c(n, p, n_prog, n_pred)) || n < 4L || p < 10L ||
      n_prog < 1L || n_pred < 1L) {
    rlang::abort("n, p, n_prog and n_pred must be positive integers (n >= 4, p >= 10).",
                 class = "cfomics_arg_error")
  }
  if (n_prog + n_pred >= p) {
    rlang::abort("p must exceed n_prog + n_pred.", class = "cfomics_arg_error")
  }
  if (!is.numeric(snr) || length(snr) != 1L || snr <= 0) {
    rlang::abort("snr must be a positive scalar.", class = "cfomics_arg_error")
  }

  # Restore the caller's RNG state on exit (see CLAUDE.md RNG convention).
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(as.integer(seed))

  genes <- .cfomics_example_genes
  predictive <- genes$predictive[seq_len(min(n_pred, length(genes$predictive)))]
  prognostic <- genes$prognostic[seq_len(min(n_prog, length(genes$prognostic)))]
  # If the caller asks for more signal genes than the curated pools hold, top up
  # from the real noise pool so the counts still line up.
  if (length(predictive) < n_pred) {
    predictive <- c(predictive,
                    sample(genes$noise_pool, n_pred - length(predictive)))
  }
  if (length(prognostic) < n_prog) {
    prognostic <- c(prognostic,
                    sample(setdiff(genes$noise_pool, predictive),
                           n_prog - length(prognostic)))
  }
  signal    <- c(prognostic, predictive)
  noise_n   <- p - length(signal)
  noise_ids <- sample(setdiff(genes$noise_pool, signal), noise_n)
  gene_ids  <- c(signal, noise_ids)                       # rows

  is_prog <- gene_ids %in% prognostic
  is_pred <- gene_ids %in% predictive
  is_sig  <- is_prog | is_pred

  # --- Negative-binomial counts ---------------------------------------------
  mu_gene <- pmin(pmax(exp(stats::rnorm(p, mean = log(200), sd = 0.8)), 10), 5000)
  # Prognostic/predictive genes get higher dispersion -> more across-sample
  # variance -> they survive default variance feature selection.
  disp <- ifelse(is_sig, 0.5, 0.15)
  counts <- matrix(0L, nrow = p, ncol = n,
                   dimnames = list(gene_ids, paste0("sample", seq_len(n))))
  for (g in seq_len(p)) {
    counts[g, ] <- stats::rnbinom(n, size = 1 / disp[g], mu = mu_gene[g])
  }
  storage.mode(counts) <- "integer"

  # --- log2-CPM + per-gene z-score (normalised signal scale) ----------------
  libsize <- pmax(colSums(counts), 1)
  cpm     <- sweep(counts, 2, libsize, "/") * 1e6
  logcpm  <- log2(cpm + 1)
  z_gene  <- t(scale(t(logcpm)))                          # per-gene z, p x n
  z_gene[!is.finite(z_gene)] <- 0
  Zt      <- t(z_gene)                                    # n x p (samples x genes)

  # --- Clinical covariates ---------------------------------------------------
  age  <- stats::rnorm(n, mean = 60, sd = 10)
  z_age <- as.numeric(scale(age))
  sexM <- stats::rbinom(n, 1, 0.5)
  sex  <- factor(ifelse(sexM == 1, "M", "F"), levels = c("F", "M"))

  # --- Outcome model -------------------------------------------------------
  # Prognostic (mu0) and predictive (tau) signals are given comparable strength
  # (Var(mu0) ~ Var(tau)); the treatment-effect heterogeneity is mostly
  # gene-driven with a smaller clinical (sex) modifier. Each additive block is
  # scaled empirically to its target variance so the moments hold for any
  # realisation, then the noise is set from the requested SNR.
  var_tau  <- 1.0     # target Var(tau)
  var_mu0  <- 1.0     # 1:1 prognostic:predictive strength
  het_frac <- 0.7     # gene share of Var(tau) (rest = sex modifier)
  ATE      <- 1.0     # baseline effect; a minority have tau < 0

  w_pred <- sample(c(-1, 1), n_pred, replace = TRUE)      # +/- effect direction
  v_prog <- sample(c(-1, 1), n_prog, replace = TRUE)
  raw_het  <- as.numeric(Zt[, which(is_pred), drop = FALSE] %*% w_pred)
  raw_prog <- as.numeric(Zt[, which(is_prog), drop = FALSE] %*% v_prog)

  het_tau <- sqrt(het_frac * var_tau / stats::var(raw_het)) * raw_het
  delta   <- sqrt((1 - het_frac) * var_tau / 0.25)         # Var(sex term) target
  tau <- ATE + delta * (sexM - 0.5) + het_tau
  mu0 <- sqrt(0.3 * var_mu0) * z_age +
    sqrt(0.7 * var_mu0 / stats::var(raw_prog)) * raw_prog

  # Noise scaled to the requested SNR against the realised Var(tau).
  sigma <- sqrt(stats::var(tau) / snr)
  epsilon <- stats::rnorm(n, 0, sigma)
  Tr <- stats::rbinom(n, 1, 0.5)
  Y  <- mu0 + Tr * tau + epsilon

  # --- Assemble SummarizedExperiment ----------------------------------------
  col_data <- S4Vectors::DataFrame(
    treatment   = as.integer(Tr),
    outcome     = as.numeric(Y),
    age         = as.numeric(age),
    sex         = sex,
    cf_tau_true = as.numeric(tau),
    row.names   = colnames(counts)
  )
  row_data <- S4Vectors::DataFrame(
    is_prognostic = is_prog,
    is_predictive = is_pred,
    row.names     = gene_ids
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = counts),
    colData = col_data,
    rowData = row_data
  )
  S4Vectors::metadata(se) <- list(
    cf_pathway           = genes$pathway_name,
    cf_predictive_genes  = predictive,
    cf_prognostic_genes  = prognostic,
    cf_true_ate          = mean(tau),
    cf_snr               = as.numeric(stats::var(tau) / stats::var(epsilon)),
    cf_params            = list(n = n, p = p, n_prog = n_prog,
                                n_pred = n_pred, snr = snr, seed = seed)
  )
  se
}
