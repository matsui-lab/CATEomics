# cfomics

`cfomics` estimates conditional average treatment effects (CATE) for
high-dimensional omics data (many features, few samples). Many CATE estimators
exist and none is reliably best in advance — the winner depends on the data.
`cfomics` therefore provides one interface to **compare** estimators on a
data-only loss, **select** one, or **combine** them into an ensemble, and to
carry the resulting individual effects into standard Bioconductor analyses. It
ingests `SummarizedExperiment` and `MultiAssayExperiment` objects and wraps 10
R-native and Python estimators behind a single API; the Python backends are
provisioned automatically through `basilisk`.

## Installation

Install the current version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("matsui-lab/cfomics")
```

Bioconductor dependencies are resolved automatically when `BiocManager` is
available. `cfomics` is being submitted to Bioconductor; once accepted it will
also be installable with `BiocManager::install("cfomics")`.

## Getting started

```r
library(cfomics)

se  <- cf_example_data()                       # example SummarizedExperiment
cmp <- cf_compare(data = se, Y = "outcome", T = "treatment",
                  X_baseline = c("age", "sex"),
                  methods = c("grf", "tlearner", "xlearner"))
cmp                                            # loss scoreboard
ens <- cf_ensemble(cmp, method = "r_stacking") # combine the estimators
```

See the vignette for the full workflow — comparing, selecting and combining
estimators, interpreting the effect, and connecting the results to
differential-expression and pathway analyses:

```r
vignette("cfomics")
```

## Getting help

Please report bugs and feature requests on the
[issue tracker](https://github.com/matsui-lab/cfomics/issues).

## License

MIT + file LICENSE.
