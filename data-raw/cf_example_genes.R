# packages/cfomics/data-raw/cf_example_genes.R
# Build the real-gene-symbol pools embedded in cfomics (R/sysdata.rda) that
# cf_example_data() draws on. Run manually; requires msigdbr (Suggests-level,
# data-generation only) so that the runtime function has NO annotation
# dependency. The pools are real human gene symbols; the *predictive* signal is
# deliberately embedded in members of a real MSigDB Hallmark pathway so the
# vignette biology-bridge recipe (fgsea / GSVA on that gene set) is coherent.
# This does NOT claim any real biological discovery -- the data are simulated.

suppressMessages({
  library(msigdbr)
})

set.seed(1L)

# --- Hallmark gene sets (real membership) -----------------------------------
h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
sym <- h$gene_symbol
gs  <- h$gs_name

set_members <- function(name) sort(unique(sym[gs == name]))

# Predictive genes = members of one real pathway (effect-modifying signal).
pathway_name <- "HALLMARK_APOPTOSIS"
apoptosis <- set_members(pathway_name)

# Prefer well-known apoptosis genes for interpretability, keep only those that
# are genuinely in the set (no fabricated membership).
predictive_pref <- c("BAX", "BID", "BCL2L11", "CASP3", "CASP8", "CASP9",
                     "FAS", "TNFRSF10B", "CFLAR", "GADD45A", "BIK", "BMF",
                     "PMAIP1", "DFFA", "CASP7")
predictive <- intersect(predictive_pref, apoptosis)
stopifnot(length(predictive) >= 10)
predictive <- predictive[seq_len(10)]

# Prognostic genes = real symbols from a biologically distinct pathway
# (baseline-outcome signal), disjoint from the predictive/apoptosis set.
prognostic_src <- set_members("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
prognostic_src <- setdiff(prognostic_src, apoptosis)
prognostic <- sort(sample(prognostic_src, 20L))

# Noise pool = union of all Hallmark genes minus every signal gene. Real
# symbols so an enrichment background over the full feature set is realistic.
signal <- union(predictive, prognostic)
noise_pool <- sort(setdiff(unique(sym), signal))
stopifnot(length(noise_pool) >= 1000L)

.cfomics_example_genes <- list(
  pathway_name = pathway_name,
  predictive   = predictive,   # length 10, real pathway members
  prognostic   = prognostic,   # length 20, real, distinct function
  noise_pool   = noise_pool    # real symbols, disjoint from signal
)

usethis::use_data(.cfomics_example_genes, internal = TRUE, overwrite = TRUE)

cat("pathway:", pathway_name, "\n")
cat("predictive (10):", paste(predictive, collapse = ", "), "\n")
cat("prognostic (20):", paste(prognostic, collapse = ", "), "\n")
cat("noise pool size:", length(noise_pool), "\n")
