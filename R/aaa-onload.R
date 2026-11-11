#' @import methods
#' @importFrom stats model.frame as.formula var glm binomial predict rnorm rbinom setNames median reorder prcomp complete.cases plogis sd quantile lm coef confint
#' @importFrom Formula Formula
#' @importFrom rlang abort warn inform .data
#' @importFrom utils packageVersion head
#' @importFrom graphics plot
#' @importFrom SummarizedExperiment colData assay
#' @section Backends:
#' Each method is computed by one backend. The R backends are declared in
#' \code{DESCRIPTION}; the Python ones are pinned in \code{R/basilisk.R} and
#' are provisioned by \pkg{basilisk} on first use.
#'
#' \tabular{lll}{
#'   \strong{Method} \tab \strong{Backend} \tab \strong{Pinned version} \cr
#'   \code{grf} \tab \pkg{grf} (R) \tab declared in \code{DESCRIPTION}, not pinned \cr
#'   \code{bcf} \tab \pkg{bcf} (R) \tab declared in \code{DESCRIPTION}, not pinned \cr
#'   \code{slearner} \tab scikit-learn (Python) \tab 1.6.1 \cr
#'   \code{tlearner}, \code{xlearner}, \code{drlearner} \tab econml (Python) \tab 0.16.0 \cr
#'   \code{tarnet}, \code{dragonnet} \tab catenets (Python) \tab 0.2.4 \cr
#'   \code{cevae} \tab pyro-ppl (Python) \tab 1.9.1 \cr
#'   \code{ganite} \tab torch (Python) \tab 2.12.0 \cr
#'   base learners \tab \pkg{glmnet}, \pkg{ranger}, \pkg{xgboost} (R) \tab declared in \code{DESCRIPTION}, not pinned \cr
#'   base learners \tab scikit-learn, xgboost (Python) \tab 1.6.1, 2.1.4 \cr
#' }
#'
#' \code{meta$provenance$software_versions} records R, the time of the fit and
#' the version of CATEomics with that of the R backend the method used. The
#' Python side is pinned rather than recorded: the versions above are the ones
#' \pkg{basilisk} installs, and neither that record nor \code{sessionInfo()}
#' reports them.
#'
#' @references
#' Huber et al. (2015) Orchestrating high-throughput genomic analysis with
#' Bioconductor, \emph{Nat. Methods} 12(2):115-121 (SummarizedExperiment).
#'
#' Ramos et al. (2017) Software for the integration of multi-omics experiments
#' in Bioconductor, \emph{Cancer Res.} 77(21):e39-42 (MultiAssayExperiment).
#'
#' Lun (2022) basilisk: a Bioconductor package for managing Python
#' environments, \emph{JOSS} 7(79):4742.
#'
#' Pedregosa et al. (2011) Scikit-learn: Machine Learning in Python,
#' \emph{JMLR} 12:2825-2830.
#'
#' Paszke et al. (2019) PyTorch: An Imperative Style, High-Performance Deep
#' Learning Library, \emph{NeurIPS} 32:8024-8035.
#'
#' Friedman, Hastie & Tibshirani (2010) Regularization Paths for Generalized
#' Linear Models via Coordinate Descent, \emph{JSS} 33(1):1-22 (glmnet).
#'
#' Wright & Ziegler (2017) ranger: A Fast Implementation of Random Forests for
#' High Dimensional Data in C++ and R, \emph{JSS} 77(1):1-17.
#'
#' Chen & Guestrin (2016) XGBoost: A Scalable Tree Boosting System,
#' \emph{KDD} 785-794.
#'
#' @section Development:
#' Much of the code, tests and documentation was written by an AI coding
#' agent (Claude Code, Anthropic) from the authors' design and under their
#' direction and review. Provenance is recorded per commit with a
#' \code{Co-Authored-By: Claude} trailer, which records the same provenance as
#' \code{Assisted-by: Claude Code (Anthropic)}. The authors are responsible for
#' the correctness of the code and for maintaining it, and the contributed code
#' is redistributable under this package's MIT terms. No third-party source file
#' is vendored into the package; the only third-party material it redistributes
#' is the MSigDB hallmark gene symbols in \code{R/sysdata.rda} (CC BY 4.0, see
#' \code{?cf_example_data}). The README has the longer form.
#'
#' @keywords internal
"_PACKAGE"

# The import directive for SummarizedExperiment is what loads that namespace
# when CATEomics loads. The package reaches the container through
# `SummarizedExperiment::` calls, which would load it only when one of them
# first runs; until then the S4 methods of the class a user hands in are not
# in R's method tables, and `se$column` fails in a session that attached
# nothing else.

# NOTE: S3 methods are registered in NAMESPACE via roxygen2 @export tags.
# No .onLoad needed for S3 method registration.
# Python initialization is deliberately NOT performed at package load.
# Python-based methods provision their environment automatically via basilisk
# on first use; no user setup is required.
