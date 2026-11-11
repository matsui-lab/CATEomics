# The largest `p` cf_example_data() accepts. The pools can reach a little
# further, by an amount that depends on n_prog and n_pred, so one stated
# maximum is published instead (see the comment at the cap check).
.CF_EXAMPLE_MAX_P <- 4000L


#' Simulate an example omics dataset for CATEomics
#'
#' Generates a small semi-synthetic \code{SummarizedExperiment} of bulk RNA-seq
#' counts with a known heterogeneous treatment effect, for use in examples,
#' vignettes and quick experimentation. The design is \eqn{n \ll p}:
#' many genes, few samples, a randomised treatment,
#' and a continuous outcome whose conditional average treatment effect (CATE)
#' depends on a handful of genes and on a clinical effect-modifier.
#'
#' @details
#' The data are \strong{fully simulated} and carry no biological meaning; they exist
#' only to make the workflow runnable end-to-end. Every argument is a constant
#' of the generating process, so the moments below follow from the arguments
#' without any hidden rescaling.
#' \itemize{
#'   \item \strong{Expression.} Per-gene mean counts are drawn from a log-normal
#'     distribution (to give a realistic dynamic range) and counts are sampled
#'     from a negative-binomial distribution (\code{\link[stats]{rnbinom}}).
#'     Prognostic and
#'     predictive genes are given higher dispersion so they carry across-sample
#'     variation and survive default variance-based feature selection.
#'   \item \strong{Signal scale.} The outcome model is defined on the \emph{normalised}
#'     expression scale: per-gene log2-CPM values are z-scored internally
#'     (\eqn{z_g}), on the same scale as the log-CPM transform CATEomics
#'     applies (this generator uses no TMM step and a different offset). The
#'     treatment effect is therefore a smooth function of normalised expression,
#'     not of raw counts.
#'   \item \strong{Clinical covariates.} \code{age} (\eqn{\sim N(60, 10)}) enters
#'     the baseline \eqn{\mu_0}; \code{sex} (\code{"F"}/\code{"M"}) shifts the
#'     treatment effect \eqn{\tau}; \code{sex} is one-hot expanded downstream.
#'   \item \strong{Treatment.} Randomised (\eqn{\Pr(T=1)=0.5}, no confounding), so
#'     the setting is an RCT.
#'   \item \strong{Outcome.} With \eqn{e} = \code{effect_size} and
#'     \eqn{z(age)} the standardised age,
#'     \deqn{\mu_0(X) = e \left( \sum_{g \in \mathrm{prog}} z_g + z(age) \right)}
#'     \deqn{\tau(X)  = \mathrm{ate} + e \left( \sum_{g \in \mathrm{pred}} z_g + 1[sex = M] - 1/2 \right)}
#'     \deqn{Y = \mu_0(X) + T\, \tau(X) + \varepsilon, \quad \varepsilon \sim N(0, \mathrm{noise\_sd}^2)}
#'     so \eqn{e} is the coefficient one signal gene carries. Every gene of a
#'     pathway pushes the same way, so a gene-set score over that pathway
#'     tracks the effect it drives. The spreads follow from the counts:
#'     \eqn{\mathrm{sd}(\tau) \approx e\sqrt{n_{pred} + 0.25}},
#'     \eqn{\mathrm{sd}(\mu_0) \approx e\sqrt{n_{prog} + 1}} and
#'     \eqn{\Pr(\tau < 0) \approx \Phi(-\mathrm{ate} / \mathrm{sd}(\tau))}.
#'     Raising \code{n_pred} therefore raises the amount of heterogeneity rather
#'     than dividing a fixed amount over more genes.
#' }
#' The predictive genes are the first \code{n_pred} members of the MSigDB
#' Hallmark pathway \code{HALLMARK_APOPTOSIS}, and the prognostic genes the
#' first \code{n_prog} members of \code{HALLMARK_OXIDATIVE_PHOSPHORYLATION}
#' (the two sets are disjoint), so the enrichment recipe in the vignette has a
#' gene set to match against. The other \code{p - n_prog - n_pred} rows are
#' real symbols drawn from the rest of the Hallmark collection and carry no
#' signal.
#'
#' @param n Integer, number of samples (columns). Default 150.
#' @param p Integer, number of genes (rows). Default 600. Must exceed
#'   \code{n_prog + n_pred}, and at most 4000. The pools hold more symbols than
#'   that, but how many a call can reach depends on \code{n_prog} and
#'   \code{n_pred}, so one stated maximum below the smallest of them is easier
#'   to rely on.
#' @param n_prog Integer, number of prognostic genes (drive \eqn{\mu_0}).
#'   Default 5. At most 194: the members of the prognostic pathway that are not
#'   also in the predictive one.
#' @param n_pred Integer, number of predictive genes (drive \eqn{\tau}).
#'   Default 5. At most 161, every member of the predictive pathway.
#' @param effect_size Numeric, the coefficient one signal gene carries on the
#'   z-scored expression scale (default 1.5); it scales both \eqn{\mu_0} and
#'   \eqn{\tau}.
#' @param ate Numeric, the constant term of \eqn{\tau}, i.e. the average
#'   treatment effect (default 0, so about half the samples are harmed).
#' @param noise_sd Numeric, the standard deviation of the outcome noise
#'   \eqn{\varepsilon} (default 1).
#' @param seed Integer, random seed (default \code{1}). The caller's RNG state
#'   is restored on exit.
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
#'     \item{\code{metadata}}{\code{cf_pathway} (the pathway the predictive
#'       genes come from), \code{cf_predictive_genes},
#'       \code{cf_prognostic_genes}, \code{cf_true_ate} (the realised mean of
#'       \eqn{\tau}), \code{cf_snr} (the realised
#'       \eqn{\mathrm{Var}(\tau)/\mathrm{Var}(\varepsilon)}), and
#'       \code{cf_params} (the call).}
#'   }
#'
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- cf_example_data(n = 60, p = 200, seed = 1)
#'   se
#'   # ground truth is carried alongside the data
#'   mean(SummarizedExperiment::colData(se)$cf_tau_true)   # ~ ate
#'   table(SummarizedExperiment::rowData(se)$is_predictive)
#' }
#'
#' @section Gene names and their licence:
#' The gene symbols are drawn from the MSigDB hallmark collection (Liberzon et
#' al. 2015). Set names and gene symbols from that collection are redistributed
#' with the package under the Creative Commons Attribution 4.0 International
#' licence (CC BY 4.0); nothing else of MSigDB is included.
#'
#' @seealso \code{\link{cf_fit}}, \code{\link{cf_compare}}
#' @references
#' Liberzon, Birger, Thorvaldsdóttir, Ghandi, Mesirov & Tamayo (2015) The
#' Molecular Signatures Database hallmark gene set collection,
#' \emph{Cell Syst.} 1(6):417-425.
#'
#' @export
cf_example_data <- function(n = 150L, p = 600L,
                            n_prog = 5L, n_pred = 5L,
                            effect_size = 1.5, ate = 0, noise_sd = 1,
                            seed = 1L) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    rlang::abort(
      "cf_example_data() requires the 'SummarizedExperiment' and 'S4Vectors' packages.",
      class = "cateomics_missing_package_error"
    )
  }
  # Every argument is checked before `set.seed()` runs. A seed that cannot be
  # coerced would otherwise raise base R's own error once the RNG-restoring
  # on.exit() had been registered, which adds a second, unrelated warning. The
  # counts are checked for being whole numbers because as.integer() truncates,
  # so n_pred = 5.7 would otherwise plant five genes without saying so.
  .cf_check_scalar <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      rlang::abort(sprintf("%s must be a finite numeric scalar of length 1.", nm),
                   class = "cateomics_arg_error")
    }
  }
  .cf_check_count <- function(x, nm) {
    .cf_check_scalar(x, nm)
    if (x != trunc(x)) {
      rlang::abort(sprintf("%s must be a whole number (got %s).", nm, format(x)),
                   class = "cateomics_arg_error")
    }
  }
  .cf_check_count(n, "n")
  .cf_check_count(p, "p")
  .cf_check_count(n_prog, "n_prog")
  .cf_check_count(n_pred, "n_pred")
  .cf_check_count(seed, "seed")
  .cf_check_scalar(effect_size, "effect_size")
  .cf_check_scalar(ate, "ate")
  .cf_check_scalar(noise_sd, "noise_sd")
  n      <- as.integer(n)
  p      <- as.integer(p)
  n_prog <- as.integer(n_prog)
  n_pred <- as.integer(n_pred)
  seed   <- as.integer(seed)
  if (n < 4L || p < 10L || n_prog < 1L || n_pred < 1L) {
    rlang::abort("n, p, n_prog and n_pred must be positive integers (n >= 4, p >= 10).",
                 class = "cateomics_arg_error")
  }
  if (n_prog + n_pred >= p) {
    rlang::abort("p must exceed n_prog + n_pred.", class = "cateomics_arg_error")
  }
  if (effect_size <= 0) {
    rlang::abort("effect_size must be positive.", class = "cateomics_arg_error")
  }
  if (noise_sd <= 0) {
    rlang::abort("noise_sd must be positive.", class = "cateomics_arg_error")
  }

  genes <- .cateomics_example_genes
  if (n_pred > length(genes$predictive)) {
    rlang::abort(
      c(sprintf("`n_pred` = %d is more than the %d members of %s that are shipped.",
                n_pred, length(genes$predictive), genes$pathway_name),
        i = sprintf("The predictive signal is planted in members of that pathway, so n_pred <= %d.",
                    length(genes$predictive))),
      class = "cateomics_arg_error"
    )
  }
  if (n_prog > length(genes$prognostic)) {
    rlang::abort(
      c(sprintf("`n_prog` = %d is more than the %d members of %s that are shipped.",
                n_prog, length(genes$prognostic), genes$prognostic_pathway_name),
        i = sprintf("The prognostic signal is planted in members of that pathway, so n_prog <= %d.",
                    length(genes$prognostic))),
      class = "cateomics_arg_error"
    )
  }

  # `p` is capped at a stated number rather than at everything the pools can
  # reach. The background pool holds every Hallmark symbol except the curated
  # head of the two signal pools, so how many symbols a call can actually use
  # depends on n_prog and n_pred, and one documented maximum is easier to rely
  # on than that figure. The cap is below the smallest value it can take, so
  # the draw below can never run short; nothing needs the rest, since the
  # default is 600 rows. Checking here also beats `sample()`'s "larger than the
  # population" further down.
  if (p > .CF_EXAMPLE_MAX_P) {
    rlang::abort(
      c(sprintf("`p` = %d is more genes than this generator offers.", p),
        i = sprintf("Use p <= %d.", .CF_EXAMPLE_MAX_P)),
      class = "cateomics_arg_error"
    )
  }
  predictive <- genes$predictive[seq_len(n_pred)]
  prognostic <- genes$prognostic[seq_len(n_prog)]
  signal     <- c(prognostic, predictive)

  # Restore the caller's RNG state on exit (see CLAUDE.md RNG convention).
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

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
  # Every coefficient is `effect_size` and nothing is rescaled afterwards, so
  # the spreads are the ones `@details` states and adding signal genes adds
  # signal instead of thinning it.
  #
  # All the genes of one pathway push the same way. A sign drawn per gene would
  # set members of one pathway against each other, which is odd on its face
  # and has a consequence: a gene-set score sums its members, so a signal whose
  # genes point in opposite directions cancels in it and a pathway-level
  # recipe could not see the planted effect at all. The spreads are unchanged,
  # because the variance of a sum of independent terms does not depend on
  # their signs.
  tau <- ate + effect_size *
    (rowSums(Zt[, which(is_pred), drop = FALSE]) + (sexM - 0.5))
  mu0 <- effect_size *
    (rowSums(Zt[, which(is_prog), drop = FALSE]) + z_age)

  epsilon <- stats::rnorm(n, 0, noise_sd)
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
    cf_params            = list(n = n, p = p, n_prog = n_prog, n_pred = n_pred,
                                effect_size = effect_size, ate = ate,
                                noise_sd = noise_sd, seed = seed)
  )
  se
}
