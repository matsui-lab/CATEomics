# CHANGES IN VERSION 0.99.0

First public release. `CATEomics` provides a unified workflow to compare,
select and combine conditional average treatment effect (CATE) estimators for
high-dimensional omics data (many features, few samples).

## NEW FEATURES

- Bioconductor-native input: `SummarizedExperiment` and `MultiAssayExperiment`
  are accepted directly (a `data.frame` too), with clinical covariates and
  per-assay omics selection taken from the same object.

- Ten CATE estimators behind one interface: `grf` and `bcf` (R-native), and
  `slearner`, `tlearner`, `xlearner`, `drlearner`, `cevae`, `ganite`,
  `tarnet` and `dragonnet` (Python, provisioned automatically through
  `basilisk`), optionally crossed with their base learners (e.g. elastic-net
  vs random forest) via `base_learners`. `cf_compare()` fits the candidates
  under shared preprocessing and data partitions and ranks them by a proxy
  CATE loss (R-loss / DR-loss); `cf_select()` picks the lowest-loss candidate
  and `cf_ensemble()` combines the per-sample effects (mean, softmax,
  R-stacking, DR-stacking, and consensus-based averaging, an average of the
  models that agree most closely). `cf_fit()` is the same interface for a
  single estimator. Continuous and binary outcomes are supported, with
  effects on binary outcomes estimated as risk differences.

- Stability diagnostics: `cf_model_stability()` measures how similar the
  candidates' per-sample estimates are to one another.

- Evaluation of the estimates: `cf_evaluate()` scores them against a
  per-patient pseudo-outcome built from the observed outcomes, using an
  independent nuisance model. It reports the targeting operator characteristic
  at a given treated fraction (TOC) and its average over the grid (AUTOC),
  grouped average effects (GATES) and the best linear predictor (BLP).

- Feature interpretation: `cf_interpret()` distils a fitted model or an ensemble
  into a surrogate SHAP explanation (an importance and a direction per feature),
  and `cf_aggregate_vi()` aggregates importances across methods by robust rank
  aggregation.

- Benefit groups and treatment rules: `cf_benefit_groups()` forms high- and
  low-benefit groups from out-of-fold effects, and `cf_rule()` builds and
  evaluates treat / don't-treat rules (optionally across a grid of operating
  points) with a doubly-robust policy value.

- Assumption diagnostics: `cf_diagnostics()` reports propensity overlap,
  standardised mean differences, the per-arm distribution of any covariate
  named in `variables`, and an E-value sensitivity bound.

- Back to the container, and on to the usual analyses: `cf_attach_results()`
  writes the estimated effects and the benefit groups onto the
  `SummarizedExperiment` they came from, so the analysis connects to the
  established downstream tools -- differential expression between the benefit
  groups with `limma`, pathway activity against the estimates with `GSVA`.

- Visualisation helpers built on `ggplot2`: `cf_plot_cate()`,
  `cf_plot_compare()`, `cf_plot_evaluate()` (RATE / GATES / BLP),
  `cf_plot_interpret()` (SHAP importance, beeswarm, dependence),
  `cf_plot_stability()`, `cf_plot_policy_curve()` and
  `cf_plot_diagnostics()`. `cf_plot_rule()` draws the policy tree and returns
  a `DiagrammeR` widget.

- `cf_example_data()` generates a small semi-synthetic `SummarizedExperiment`
  with a known heterogeneous treatment effect for examples and the vignette.

- A workflow vignette walks through fitting, comparing, selecting, combining and
  interpreting estimators, and connects the estimated effects to standard
  Bioconductor differential-expression and pathway-activity analyses.
