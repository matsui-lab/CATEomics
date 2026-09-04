# CHANGES IN VERSION 0.99.0

First Bioconductor release. `cfomics` provides a unified workflow to compare,
select and combine conditional average treatment effect (CATE) estimators for
high-dimensional omics data (many features, few samples).

## NEW FEATURES

- Unified fitting interface `cf_fit()` wrapping 10 CATE estimators behind one
  API: `grf`, `bcf` (R-native) and `slearner`, `tlearner`, `xlearner`,
  `drlearner`, `cevae`, `ganite`, `tarnet`, `dragonnet` (Python, provisioned
  automatically through `basilisk`). Both continuous and binary outcomes are
  supported.

- Bioconductor-native input: `SummarizedExperiment` and `MultiAssayExperiment`
  are accepted directly, with clinical covariates and per-assay omics selection.
  `as_cf_data()` converts these containers, and `cf_attach_results()` writes
  estimated effects and benefit groups back onto a `SummarizedExperiment`.

- Compare, select and combine estimators: `cf_compare()` fits several methods
  with cross-fitting -- optionally across their base learners (e.g. elastic-net
  vs random forest) via `base_learners` -- and ranks them by a data-only CATE
  loss (R-loss / DR-loss); `cf_select()` picks the lowest-loss method;
  `cf_ensemble()` combines the per-sample effects (mean, softmax, R-stacking,
  DR-stacking, CBA).

- Stability diagnostics: `cf_model_stability()` summarises how consistently the
  compared methods rank samples by benefit.

- Evaluation of treatment prioritisation: `cf_evaluate()` reports the
  rank-weighted average treatment effect (RATE), a best linear predictor (BLP)
  and grouped average effects (GATES), using an independent nuisance model.

- Feature interpretation: `cf_interpret()` distils a fitted model or an ensemble
  into a surrogate SHAP explanation (an importance and a direction per feature),
  and `cf_aggregate_vi()` aggregates importances across methods by robust rank
  aggregation.

- Benefit groups and treatment rules: `cf_benefit_groups()` forms high- and
  low-benefit groups from out-of-fold effects, and `cf_rule()` builds and
  evaluates treat / don't-treat rules (optionally across a grid of operating
  points) with a doubly-robust policy value.

- Assumption diagnostics: `cf_diagnostics()` reports propensity overlap,
  standardised mean differences and an E-value sensitivity bound.

- Visualisation helpers built on `ggplot2`: `cf_plot_cate()`,
  `cf_plot_compare()`, `cf_plot_evaluate()` (RATE / GATES / BLP),
  `cf_plot_interpret()` (SHAP importance, beeswarm, dependence),
  `cf_plot_stability()`, `cf_plot_rule()`, `cf_plot_policy_curve()` and
  `cf_plot_diagnostics()`.

- `cf_example_data()` generates a small semi-synthetic `SummarizedExperiment`
  with a known heterogeneous treatment effect for examples and the vignette.

- A workflow vignette walks through fitting, comparing, selecting, combining and
  interpreting estimators, and connects the estimated effects to standard
  Bioconductor differential-expression and pathway-activity analyses.
