# CATEomics

## Background

High-throughput omics assays measure thousands of molecular features in a
sample from a single patient. Such data are widely used to predict prognosis,
but choosing a treatment calls for the treatment effect rather than the
prognosis: how far apart the outcomes under the treatments are for patients
with a given molecular profile. Under the assumptions required for causal
identification, that difference is the conditional average treatment effect
(CATE), and estimating it from omics data can identify molecular profiles
associated with differential treatment benefit.

## Why it is needed

Many CATE estimators have been proposed, and which of them is best depends on
the structure of the data, so each analysis has to choose one of them or
combine several. Each patient, however, reveals an outcome under only one
treatment, so the estimation error that would settle that choice cannot be
evaluated directly from the data. Candidates can instead be compared on
*proxy losses*: stand-ins for the true error, computed on patients a model was
not fitted on, and it is this comparison that guides the selection or the
ensembling. Such a comparison is possible, but for the comparison to be valid,
every candidate needs the same preprocessing, the same data partitions and the
same scoring.

The estimates are also not the end point of an omics study. They do their work
once they are back on the original data container and feed the
differential-expression and gene-set analyses in which one can explore which
molecular profiles are associated with differential treatment benefit.

Existing software covers only parts of this workflow. `causalToolbox` provides
an R interface to several CATE estimators, and `EconML` supports estimation,
R-loss-based model selection and ensembling in Python. Neither provides an
end-to-end Bioconductor workflow that maintains consistent preprocessing and
evaluation and returns the CATE estimates to standard omics data containers
for downstream analysis. Comparing estimators therefore depends on custom
cross-language code and is difficult to reproduce.

## What it does

`CATEomics` provides that workflow. It takes a `SummarizedExperiment`, fits
the candidates under the same conditions, compares them, selects one of them
or combines several, and writes the estimates and the benefit groups back
onto the object.

- **Common preprocessing and data partitions.** Ten estimators are fitted
  under the same conditions: the causal forest and the Bayesian causal forest
  in R, and the S-, T- and X-learners, the DR-learner, TARNet, DragonNet,
  CEVAE and GANITE in Python. The four learners can be crossed with their base
  learners as well. RNA-seq counts can be TMM-normalised and log-CPM
  transformed internally; features can be reduced by variance, by principal
  components or by gene-set scoring, or not at all.
- **Comparison and selection on proxy losses.** The true estimation error
  cannot be computed, so candidates are scored on the proxy losses that stand
  in for it -- the R-loss and the DR-loss, with the factual loss reported
  alongside. Because every step up to the scoring is shared, a difference in
  loss is a difference between the estimates themselves.
- **Ensembling.** Five ways to combine the estimates: a simple average,
  softmax weights, stacking against the R-loss or the DR-loss, and an average
  of the models that agree most closely.
- **Evaluation of the estimates.** Individual effects are never observed, but
  a per-patient pseudo-outcome that equals the true CATE on average can be
  built from the observed outcomes. Scored against it, the estimates show
  whether they really rank patients by benefit (TOC, AUTOC, GATES, BLP).
- **Treatment rules and policy value.** A rule for whom to treat, and the mean
  outcome under that rule against treating everyone, treating no one, and
  treating a randomly chosen group of the same size.
- **Interpretation.** SHAP values of a surrogate model rank the molecular
  features that move the estimates.
- **Bioconductor-native.** `SummarizedExperiment` and `MultiAssayExperiment`
  go in as they are (a `data.frame` is accepted too), and `cf_attach_results()`
  writes the estimates and the benefit groups back onto them, so that the
  analysis connects to the established downstream tools: differential
  expression between the benefit groups with `limma`, pathway activity against
  the estimates with `GSVA`.
- **Python provisioned for you.** `basilisk` builds the environment on first
  use; there is nothing to set up.

## Installation

Install the current version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("matsui-lab/CATEomics")
```

Bioconductor dependencies are resolved automatically when `BiocManager` is
available.

## Getting started

```r
library(CATEomics)

se  <- cf_example_data()                       # example SummarizedExperiment
cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
                  X_baseline = c("age", "sex"),
                  methods = c("grf", "tlearner", "xlearner"))
cmp                                            # loss scoreboard
ens <- cf_ensemble(cmp, method = "r_stacking") # combine the estimators
```

The first call to a method with a Python backend builds its environment with
`basilisk`. This takes a few minutes and happens once.

See the vignette for the full workflow — comparing, selecting and combining
estimators, interpreting the effect, and connecting the results to
differential-expression and pathway analyses:

```r
vignette("CATEomics")
```

## Getting help

Please report bugs and feature requests on the
[issue tracker](https://github.com/matsui-lab/CATEomics/issues).

## How this package was written

Much of the code, tests, documentation and vignettes in CATEomics was
written by an AI coding agent (Claude Code, Anthropic), working from the
authors' design and under their direction and review. Provenance is recorded
per commit with a `Co-Authored-By: Claude` trailer. The authors are responsible
for the correctness of the code and for maintaining it.

The contributed code is redistributable under this package's MIT terms, and
`Assisted-by: Claude Code (Anthropic)` records the same provenance as the
trailer above.

No third-party source file is vendored into the package; the only third-party
material it redistributes is the MSigDB hallmark gene symbols carried in
`R/sysdata.rda`, under CC BY 4.0 (see `?cf_example_data`).

## License

MIT + file LICENSE.
