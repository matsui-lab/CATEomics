# data-raw/cf_example_genes.R
# Build the real-gene-symbol pools embedded in CATEomics (R/sysdata.rda) that
# cf_example_data() draws on. Run manually from the package directory; requires
# msigdbr (data generation only) so that the runtime function has NO annotation
# dependency. The pools are real human gene symbols; the signal is deliberately
# embedded in members of real MSigDB Hallmark pathways so the vignette
# biology-bridge recipe (fgsea / GSVA on those gene sets) is coherent. This does
# NOT claim any real biological discovery -- the data are simulated.
#
# Three pools:
#
#   predictive   the whole HALLMARK_APOPTOSIS pathway
#   prognostic   HALLMARK_OXIDATIVE_PHOSPHORYLATION minus the genes it shares
#                with the predictive pathway, so the two signal sets are disjoint
#   noise_pool   every Hallmark symbol except the curated head of the two signal
#                pools, which is where a default-sized call plants its signal
#
# The background pool is not disjoint from the signal pathways: a pathway member
# beyond the curated head stays in it and is removed from the draw at run time,
# which is what keeps the row names unique. The generator therefore has slightly
# fewer symbols to offer than it holds, by an amount that depends on n_prog and
# n_pred, and cf_example_data() states one fixed maximum for `p` instead of that
# figure.
#
# Each signal pool is ordered so that a hand-picked, well-known head comes
# first, and cf_example_data() takes the first `n_pred` / `n_prog` of each. The
# head symbols are written out here rather than read back from the shipped data:
# a script that rebuilds the data must not depend on the data it rebuilds.

suppressMessages({
  library(msigdbr)
})

# The curated head of each pool, in the order the generator plants them.
HEAD_PRED <- c("BAX", "BID", "BCL2L11", "CASP3", "CASP8",
               "CASP9", "FAS", "CFLAR", "GADD45A", "BIK")
HEAD_PROG <- c("ACAT1", "ATP5F1A", "ATP5MC3", "ATP6V1G1", "COX15",
               "COX7A2", "ECI1", "FH", "FXN", "HADHA",
               "IDH1", "MRPL35", "MRPS11", "MRPS30", "NDUFB6",
               "SDHB", "SLC25A11", "SLC25A20", "TIMM9", "UQCRB")

# --- Hallmark gene sets (real membership) -----------------------------------
h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
sym <- h$gene_symbol
gs  <- h$gs_name

set_members <- function(name) sort(unique(sym[gs == name]))

pathway_name            <- "HALLMARK_APOPTOSIS"
prognostic_pathway_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
apoptosis <- set_members(pathway_name)
oxphos    <- setdiff(set_members(prognostic_pathway_name), apoptosis)

predictive <- c(HEAD_PRED, setdiff(apoptosis, HEAD_PRED))
prognostic <- c(HEAD_PROG, setdiff(oxphos, HEAD_PROG))
noise_pool <- sort(setdiff(unique(sym), c(HEAD_PRED, HEAD_PROG)))

stopifnot(
  # the curated head is intact and still inside its pathway
  all(HEAD_PRED %in% apoptosis), all(HEAD_PROG %in% oxphos),
  identical(predictive[seq_along(HEAD_PRED)], HEAD_PRED),
  identical(prognostic[seq_along(HEAD_PROG)], HEAD_PROG),
  # the two signal sets are disjoint, and the background pool holds neither head
  length(intersect(predictive, prognostic)) == 0L,
  length(intersect(noise_pool, c(HEAD_PRED, HEAD_PROG))) == 0L,
  !anyDuplicated(noise_pool),
  # enough background for a realistic feature set
  length(noise_pool) >= 1000L
)

.cateomics_example_genes <- list(
  pathway_name            = pathway_name,
  prognostic_pathway_name = prognostic_pathway_name,
  predictive              = predictive,   # whole pathway, curated head first
  prognostic              = prognostic,   # whole pathway, curated head first
  noise_pool              = noise_pool    # every Hallmark symbol but the head
)

usethis::use_data(.cateomics_example_genes, internal = TRUE, overwrite = TRUE)

cat("predictive pathway:", pathway_name, "-", length(predictive), "genes",
    "(head:", paste(head(predictive, 3), collapse = ", "), "...)\n")
cat("prognostic pathway:", prognostic_pathway_name, "-", length(prognostic), "genes",
    "(head:", paste(head(prognostic, 3), collapse = ", "), "...)\n")
cat("noise pool size:", length(noise_pool), "\n")
cat("distinct symbols in total:",
    length(unique(c(predictive, prognostic, noise_pool))), "\n")
