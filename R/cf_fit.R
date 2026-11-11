#' Fit a Causal Inference Model
#'
#' @description
#' Fits a causal inference model using various backends with a unified API.
#' Supports data.frame, SummarizedExperiment, and MultiAssayExperiment objects.
#'
#' @details
#' \strong{Estimand:}
#'
#' All methods estimate the \strong{Average Treatment Effect (ATE)}:
#' \deqn{ATE = E[Y(1) - Y(0)]}
#' where \eqn{Y(1)} and \eqn{Y(0)} are potential outcomes under treatment
#' and control, respectively. Every method also estimates the
#' \strong{Conditional Average Treatment Effect}:
#' \deqn{\tau(x) = E[Y_i(1) - Y_i(0) \mid X_i = x]}
#'
#' All methods estimate the individual/conditional treatment effect:
#' \code{grf}, \code{bcf}, \code{slearner}, \code{tlearner}, \code{xlearner},
#' \code{drlearner}, \code{ganite}, \code{cevae}, \code{tarnet},
#' \code{dragonnet}. \code{predict()} returns per-sample \eqn{\hat\tau} and
#' potential outcomes for every method.
#'
#' \strong{Counterfactual outcomes \eqn{\hat\mu_0}, \eqn{\hat\mu_1}:}
#' \eqn{\hat\mu_t} is always the \emph{modelled} counterfactual under
#' \eqn{\text{do}(T=t)}, not the observed value for treated / control
#' subsets.
#'
#' \strong{Methods:}
#' \itemize{
#'   \item \strong{grf}: R-native generalized random forest.
#'   \item \strong{bcf}: Bayesian Causal Forests (\code{bcf::bcf}); the
#'     outcome model is a BART forest and the propensity is fitted
#'     separately.
#'   \item \strong{slearner}: one sklearn model on a design augmented with
#'     the treatment interactions, fitted by CATEomics; see
#'     \code{\link{metalearners}}.
#'   \item \strong{tlearner / xlearner}: sklearn-backed metalearners via
#'     econml; see \code{\link{metalearners}}.
#'   \item \strong{drlearner}: Neyman-orthogonal DR-Learner (econml).
#'   \item \strong{ganite}: GAN-based ITE (PyTorch).
#'   \item \strong{cevae}: latent-confounding VAE, requires Pyro.
#'   \item \strong{tarnet / dragonnet}: neural CATE via catenets.
#' }
#'
#' @param data Data source: data.frame, SummarizedExperiment, or
#'   MultiAssayExperiment object.
#' @param Y Outcome column name (character scalar). Numeric; a 0/1 outcome is
#'   a supported mode (risk difference). For a SummarizedExperiment the
#'   outcome must be a column of \code{colData}; one found among the assay
#'   rownames is an error (\code{cateomics_y_in_assay}).
#'
#'   CATEomics treats \strong{larger values of \code{Y} as the better outcome}.
#'   A positive \code{tau_hat} therefore means the treatment moves the outcome
#'   in the better direction, \code{\link{cf_benefit_groups}} calls the units
#'   with the largest \code{tau_hat} the high-benefit group, and
#'   \code{\link{cf_rule}} treats the units with the largest \code{tau_hat}.
#'   If smaller values of your outcome are better -- a severity score, a
#'   symptom count -- analyse the negated outcome (\code{-Y}, or
#'   \code{1 - Y} for a binary outcome) so that this holds.
#' @param T Treatment column name (character scalar). Binary: numeric 0/1, a
#'   logical, or a factor / character whose values are exactly "0" and "1".
#'   Any other labelling is an error, because which arm is the treated one is
#'   not inferred (\code{cateomics_treatment_not_binary}).
#' @param X_baseline Character vector of non-omics (clinical / demographic)
#'   covariate column names, for the named-args API.
#' @param omics Which measurements form the omics block. \code{"all"} (the
#'   default for SummarizedExperiment / MultiAssayExperiment) takes every
#'   feature of every layer, \code{"none"} (the default for a data.frame) takes
#'   none, a character vector names the features of the single layer (or names
#'   a whole layer), and a named list selects per layer, e.g.
#'   \code{list(rnaseq = "all", methyl = c("cg001", "cg002"))}. Layer names are
#'   assay names (SE), experiment names (MAE) or free labels (data.frame). A
#'   SummarizedExperiment with several assays must name the one to use; they
#'   are not combined. For a MultiAssayExperiment each column carries its
#'   experiment's name (\code{<experiment>_<feature>}). A bare feature vector
#'   is refused when the container holds several layers, since it does not say
#'   which layer they belong to; name the layer or use a named list. A
#'   data.frame has no layers to name, so a single string is read as a column
#'   name and the names of a named list are free labels. Naming the same
#'   column in both \code{X_baseline} and \code{omics} is an error. The omics
#'   block can carry the covariates on its own: \code{X_baseline} may be left
#'   out. A missing value in an omics feature is dropped with a warning naming
#'   it when the layer is taken whole (\code{"all"}, a layer name, or a
#'   data.frame's omics columns). Naming the features one by one makes them
#'   covariates of the model, and a missing value in a named covariate --
#'   whether it came from \code{omics} or from \code{X_baseline} -- is an
#'   error.
#' @param omics_type Optional layer type, as a single string or a vector named
#'   by layer. Only the RNA-seq layers are library-size normalised. With a
#'   single omics layer that is the one \code{preprocess} refers to, whatever it
#'   is called. With several layers -- and always for a
#'   \code{MultiAssayExperiment} -- the layer's name decides, and the names that
#'   count as RNA-seq are \code{"rnaseq"}, \code{"rna"}, \code{"rna_seq"} and
#'   \code{"rna-seq"} (case is ignored). Any other name is treated as not
#'   RNA-seq, so use this argument when the name does not say what the layer is
#'   (e.g. \code{omics_type = c(counts = "rnaseq")}). With more than one layer,
#'   which layers were normalised is reported once per call.
#' @param method Character, the causal inference method to use
#'   (default \code{"grf"}):
#'   \describe{
#'     \item{"grf"}{Generalized Random Forest (R-native)}
#'     \item{"bcf"}{Bayesian Causal Forests (R-native, requires bcf)}
#'     \item{"slearner"}{one sklearn model on an interaction-augmented
#'       design, fitted by CATEomics (Python)}
#'     \item{"tlearner" / "xlearner"}{sklearn-backed metalearners via econml
#'       (Python)}
#'     \item{"drlearner"}{Neyman-orthogonal DR-Learner (econml, Python)}
#'     \item{"cevae"}{Causal Effect VAE (Python/Pyro)}
#'     \item{"ganite"}{GAN-based ITE estimation (Python/PyTorch)}
#'     \item{"tarnet" / "dragonnet"}{representation-learning neural CATE
#'       (catenets, Python)}
#'   }
#' @param n_bootstrap Integer, bootstrap replicates used when
#'   \code{bootstrap_ci = TRUE}, for the methods that accept it
#'   (s/t/x-learner, drlearner). Default 100.
#' @param random_state Integer random seed for reproducibility (default
#'   \code{0}). With the same inputs and the same thread count it yields
#'   bit-identical results for every method. Three methods also depend on
#'   the thread count, which the seed does not reach: \code{bcf::bcf()}
#'   takes its thread count from the machine and combines each iteration's
#'   partial results in whatever order the threads finish, so the seed alone
#'   does not determine the draws, and the two catenets methods
#'   (\code{"tarnet"}, \code{"dragonnet"}) return different low-order digits
#'   when the thread count changes. Pass
#'   \code{params = list(bcf = list(n_threads = 1L))} for a bit-identical bcf
#'   fit; the thread count the catenets methods see is set outside R, for
#'   instance through \code{OMP_NUM_THREADS}. The caller's R \code{.Random.seed} and any
#'   Python-side RNG state touched by the backend (\code{numpy},
#'   \code{torch}, Python \code{random}) are saved on entry and restored
#'   on exit so seeding does not leak into the user's session. See the
#'   per-method documentation for backend-specific scope.
#' @param preprocess Logical (default \code{TRUE}). When \code{TRUE} the omics
#'   block is library-size normalised (TMM then log-CPM). A
#'   \code{MultiAssayExperiment} is normalised layer by layer, so only its
#'   RNA-seq layers are; a \code{data.frame} that names a mix of RNA-seq and
#'   other layers is refused instead (\code{cateomics_mixed_layer_preprocess}),
#'   because its layers share one matrix. See \code{omics_type}. With
#'   \code{omics = "none"} there is nothing to normalise.
#' @param batch Optional character: column name in colData (SE/MAE) or data
#'   (df) containing batch labels for \code{sva::ComBat} correction on the
#'   log-CPM values. \code{NULL} (default) skips batch
#'   correction.
#' @param feature_method One of \code{"variance"} (default), \code{"none"},
#'   \code{"pca"}, \code{"pathway"}.
#'   \code{"variance"}/\code{"pca"} pick top \code{n_features} by row variance
#'   or principal components; \code{"none"} keeps every assay feature.
#'   \code{"pathway"} runs \code{GSVA::gsva} with \code{gene_sets} (default
#'   Hallmark 50 via msigdbr). The default \code{"variance"} means an SE / MAE
#'   omics assay is included as covariates by default (a \code{data.frame}
#'   without an assay is unaffected).
#' @param n_features Integer, number of features / components / pathways to
#'   keep. \code{NULL} (default) resolves to 100 for a SummarizedExperiment /
#'   data.frame and 50 for a MultiAssayExperiment (per-assay budget).
#' @param gene_sets Named list of gene sets, used only when
#'   \code{feature_method = "pathway"}. \code{NULL} (default) loads Hallmark
#'   50 via msigdbr.
#' @param pathway_method One of \code{"gsva"} (default), \code{"ssgsea"},
#'   \code{"plage"}, \code{"zscore"}; passed to \code{GSVA::gsva}. Scores
#'   for \code{newdata} are computed from the training samples and the new
#'   samples together, keeping the columns of the new samples.
#' @param partition One of \code{"oof"} (default) or
#'   \code{"none"}. When \code{"oof"}, cf_fit runs K-fold OOF + full-data
#'   fit (\code{n_folds + 1} fits per method) and populates
#'   \code{fit$res$tau_hat_oof}, \code{fit$res$fold_id} and
#'   \code{fit$res$nuisance}. When \code{"none"},
#'   only the single full-data fit is run; \code{fit$res$tau_hat_oof} and
#'   \code{fit$res$fold_id} are
#'   \code{NULL} but \code{fit$res$nuisance} is still populated with full-data
#'   in-sample nuisance. Its four fields keep their \code{*_oof} names under
#'   either partition, so that everything downstream reads one shape; under
#'   \code{"none"} the values behind those names are in-sample.
#' @param n_folds Integer K, default \code{5L}. Under
#'   \code{partition = "oof"} it is the stratified K-fold split the
#'   out-of-fold effects and the nuisance are computed on. Under
#'   \code{partition = "none"} there are no outer folds, but it is still the
#'   internal split the Super Learner weights are solved on, so it changes
#'   \code{fit$res$nuisance} and everything downstream that reads it.
#' @param eval_nuisance Learner for the evaluation nuisance (\code{mu0},
#'   \code{mu1}, \code{pi}) needed by DR-based metrics in
#'   \code{\link{cf_evaluate}} and \code{\link{cf_interpret}}.
#'   \code{"super_learner"} (default) fits the members of \code{eval_library}
#'   and combines them with non-negative weights solved on the out-of-fold
#'   predictions. A single learner is also accepted: \code{"rf"},
#'   \code{"lasso"}, \code{"ridge"}, \code{"elastic_net"}, \code{"linear"},
#'   \code{"logistic"} or \code{"xgboost"}. A \code{list(outcome = ,
#'   propensity = )} gives the outcome and propensity slots different
#'   learners, and naming only one slot applies it to both. The weights are
#'   recorded in \code{meta$estimator$eval_sl_weights}.
#' @param eval_library Character vector of Super Learner member names, used
#'   when \code{eval_nuisance = "super_learner"}. \code{NULL} (default) takes
#'   \code{c("mean", "ridge", "elastic_net", "lasso", "rf", "xgboost_d2",
#'   "xgboost_d4")}; \code{"xgboost_d2"} and \code{"xgboost_d4"} are boosted
#'   fits at \code{max_depth} 2 and 4. A \code{list(outcome = ,
#'   propensity = )} gives the two slots different libraries. Naming it
#'   while \code{eval_nuisance} names a single learner has no effect and
#'   warns (\code{cateomics_eval_library_unused}).
#' @param eval_propensity_score Optional known propensity score (scalar in
#'   (0, 1) or length-n vector) used for the evaluation nuisance instead of an
#'   estimated one (RCT context). \code{NULL} (default) estimates it. Giving
#'   this together with \code{eval_nuisance$propensity} is an error. Symmetric
#'   with \code{\link{cf_compare}}.
#' @param keep_kfold_models Logical (default \code{FALSE}). When \code{TRUE}
#'   and \code{partition = "oof"}, the K per-fold fit objects are
#'   retained in \code{fit$fold_cache}.
#' @param bootstrap_ci Logical (default \code{FALSE}); a bootstrap CI for the
#'   ATE and the per-sample effects, resampled by CATEomics the same way for
#'   every method that accepts it: \code{slearner}, \code{tlearner},
#'   \code{xlearner} and \code{drlearner}. Asking for it with any other
#'   method is an error (\code{cateomics_bootstrap_unsupported}, or
#'   \code{cateomics_bootstrap_unsupported_nn} for the neural ones) rather than
#'   a silent no-op. When the method has an interval of its own that the
#'   request would replace, a \code{cateomics_bootstrap_overrides_native}
#'   warning says so.
#' @param device One of \code{"auto"} (default), \code{"cpu"}, \code{"gpu"};
#'   forwarded to the chosen backend, which places its tensors on it.
#'   \code{"auto"} and \code{"gpu"} fall back to the CPU when no GPU is
#'   visible. The device the fit ran on is recorded in
#'   \code{meta$provenance$device_effective}, and
#'   \code{predict(newdata = )} reuses it.
#' @param clip Numeric length-2 vector \code{c(L, U)}; \code{NULL} (the
#'   default) resolves to \code{c(0.025, 0.975)}. Propensity scores below
#'   \code{L} or above \code{U} are replaced by the bound, which changes the
#'   weights and does not drop samples. It reaches the four methods that
#'   estimate a propensity (\code{grf}, \code{bcf}, \code{xlearner},
#'   \code{drlearner}) and is applied the same way in each; naming it for any
#'   other method is reported, since there is no propensity to bound.
#' @param propensity_score Optional scalar or length-\code{n} numeric
#'   vector of known propensity scores (RCT context).
#'   \code{NULL} (default) lets the backend estimate PS internally.
#' @param outcome_type One of \code{"auto"} (default), \code{"continuous"},
#'   or \code{"binary"}. \code{"auto"} detects a 0/1 outcome and
#'   switches every method to a risk-difference / classifier outcome model,
#'   skipping the outcome z-score for NN backends. Force with
#'   \code{"binary"} / \code{"continuous"}.
#'
#'   With a binary outcome the estimates are not constrained to the range the
#'   quantity is interpretable on: \code{tau_hat} can fall outside
#'   \eqn{[-1, 1]} with \code{drlearner}, and \code{mu0_hat} / \code{mu1_hat}
#'   can fall outside \eqn{[0, 1]} with \code{grf} and \code{bcf}.
#' @param propensity_X Character vector of the column names a propensity is
#'   fitted on. \strong{Unset, it is the clinical covariates}: the ones named
#'   in \code{X_baseline}, after factor expansion. The omics block is left to
#'   the outcome models, since a propensity fitted on thousands of
#'   pre-treatment features is unstable and treatment is rarely assigned on
#'   them. A fit with no clinical covariate at all has nothing to fall back to
#'   but every column, which it uses after saying so
#'   (\code{cateomics_propensity_X_fallback}).
#'   It restricts two propensities: the estimator's, for the methods that fit
#'   a separable one (grf / bcf / drlearner / xlearner; naming it for one of
#'   the others emits \code{cateomics_propensity_X_ignored}), and
#'   \strong{the propensity slot of the evaluation nuisance that the proxy
#'   losses and the evaluation metrics read -- for every method, including the
#'   six that have no separable propensity of their own}. The outcome slots
#'   always see every covariate. A known \code{propensity_score} /
#'   \code{eval_propensity_score} fits no propensity on that side, so nothing
#'   is restricted there.
#' @param model_outcome,model_final,model_propensity Base-learner tokens for
#'   the methods that take pluggable learners: \code{"linear"} /
#'   \code{"logistic"} / \code{"lasso"} / \code{"ridge"} /
#'   \code{"elastic_net"} / \code{"rf"} / \code{"xgboost"}. \code{NULL}
#'   keeps each method's own default, which is \code{"ridge"} for
#'   \code{model_outcome}, \code{"logistic"} for \code{model_propensity},
#'   and for \code{model_final} \code{"ridge"} with \code{drlearner} or
#'   whatever \code{model_outcome} resolved to with \code{xlearner}. \code{"lasso"}, \code{"ridge"} and
#'   \code{"elastic_net"} select their penalty by cross-validation over 5
#'   folds of the rows they are given; a slot that would leave fewer than 5
#'   rows (or, for a class label, fewer than 5 of the smallest class) is
#'   refused, and fewer than 10 is reported. One exception: for a continuous
#'   outcome the Python methods give \code{"ridge"} a closed-form
#'   leave-one-out, which forms no folds, so neither applies to it there.
#' @param params Hyperparameters for the fitted method, either flat
#'   (\code{list(num.trees = 4000)}) or keyed by method
#'   (\code{list(grf = list(num.trees = 4000))}).
#'   A method that fits more than one learner takes a block per slot:
#'   \code{outcome} for the S- and T-learner, \code{outcome} and
#'   \code{final} for the X- and DR-learner, as in
#'   \code{list(xlearner = list(outcome = list(alphas = c(0.1, 1))))}. These
#'   slot names are the method's own learners; the \code{outcome} /
#'   \code{propensity} pair of \code{eval_nuisance}, \code{eval_library} and
#'   \code{eval_params} names the evaluation nuisance instead.
#'   The names are the ones the
#'   upstream package takes, so grf's dotted spellings and bcf's
#'   \code{nburn} / \code{nsim} / \code{n_chains}. Unknown names abort, as
#'   does a random seed (\code{seed} / \code{random_seed} /
#'   \code{random_state}): \code{random_state} is the one argument that
#'   drives every draw.
#'   A value a Python-backed method requires to be whole is written as an R
#'   integer (\code{10L}, not \code{10}), because those libraries check the
#'   type; the R-native methods take either.
#' @param propensity_params Hyperparameters for the propensity model
#'   (bcf / xlearner / drlearner).
#' @param eval_params Hyperparameters for the evaluation nuisance models.
#'   With the default Super Learner the block is keyed by library member
#'   (\code{list(xgboost_d4 = list(nrounds = 50))}), optionally inside
#'   \code{list(outcome = , propensity = )}; with a single learner the block
#'   holds its hyperparameters directly. Those two slots are the evaluation
#'   nuisance's outcome and propensity models, not a method's own learners.
#'
#' @return Object of class \code{"cf_model"} with components:
#'   \describe{
#'     \item{method}{Character string identifying the method used.}
#'     \item{res}{Numeric results: \code{tau_hat}, \code{ate},
#'       \code{mu0_hat}, \code{mu1_hat}, \code{summary}, and, where the
#'       method and the partition provide them, \code{tau_hat_oof},
#'       \code{fold_id} (which fold held each row out),
#'       \code{ate_ci}, \code{tau_ci}, \code{ci_method}, \code{ps},
#'       \code{ps_type} and \code{nuisance}.}
#'     \item{model}{The fitted backend container. For the Python methods this
#'       holds the path the model was serialised to and the backend name.}
#'     \item{meta}{Metadata: the \code{provenance} and
#'       \code{estimator} buckets, plus flat entries -- \code{n},
#'       \code{p}, \code{covariate_names}, \code{named_args} (the variables
#'       you named), \code{spec} (the resolved method settings) and the
#'       preprocessing slots the newdata path replays (\code{x_zscore},
#'       \code{y_zscore}, \code{factor_levels}, \code{input_type},
#'       \code{propensity_X}). The seed is
#'       \code{provenance$rng_seed}; the split is
#'       \code{estimator$partition}; the preprocessing \emph{record},
#'       including each omics layer's selection and PCA rotation, is
#'       \code{estimator$preprocess}.}
#'     \item{data}{\code{X_preprocessed}, \code{Y}, \code{T}.}
#'     \item{fold_cache}{Present only with
#'       \code{keep_kfold_models = TRUE}.}
#'   }
#'   See \code{\link{cateomics_result}} for the full structure specification.
#'
#' @seealso
#'   \code{\link{cateomics_result}} for result object structure,
#'   \code{\link{cf_diagnostics}} for post-fit assumption checks
#' @references
#' van der Laan, Polley & Hubbard (2007) Super Learner,
#' \emph{Stat. Appl. Genet. Mol. Biol.} 6(1):Article 25 (the default evaluation
#' nuisance).
#'
#' Robinson & Oshlack (2010) A scaling normalization method for differential
#' expression analysis of RNA-seq data, \emph{Genome Biol.} 11:R25 (TMM).
#'
#' Chen, Chen, Lun, Baldoni & Smyth (2025) edgeR v4, \emph{Nucleic Acids Res.}
#' 53(2):gkaf018 (implementation).
#'
#' Johnson, Li & Rabinovic (2007) Adjusting batch effects in microarray
#' expression data using empirical Bayes methods, \emph{Biostatistics}
#' 8(1):118-127 (ComBat).
#'
#' Hänzelmann, Castelo & Guinney (2013) GSVA: gene set variation analysis for
#' microarray and RNA-Seq data, \emph{BMC Bioinformatics} 14:7.
#'
#' Barbie et al. (2009) Systematic RNA interference reveals that oncogenic
#' KRAS-driven cancers require TBK1, \emph{Nature} 462(7269):108-112 (ssGSEA).
#'
#' @export
#' @examples
#' # Basic example with a GRF causal forest (R-native).
#' set.seed(123)
#' n <- 100
#' demo_data <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X1 + rnorm(n)
#'
#' # GRF causal forest (R-native, ITE-capable).
#' fit <- cf_fit(data = demo_data, Y = "Y", T = "T",
#'               X_baseline = c("X1", "X2"), method = "grf",
#'               partition = "none")
#' predict(fit, type = "ate")        # point ATE
#' head(predict(fit, type = "tau"))  # per-sample CATE
#' summary(fit)                      # ATE + 95% CI
#'
#' \donttest{
#' # SummarizedExperiment input via the named-args API
#' se <- cf_example_data(n = 40, p = 60, n_prog = 8, n_pred = 4)
#' fit_se <- cf_fit(data = se, Y = "outcome", T = "treatment",
#'                  X_baseline = c("age", "sex"), method = "grf",
#'                  n_features = 40)
#' }
cf_fit <- function(data = NULL,
                   method = c("grf",
                              "bcf",
                              "ganite", "drlearner",
                              "cevae",
                              "slearner", "tlearner", "xlearner",
                              "tarnet", "dragonnet"),
                   # ---- The names of the variables ----
                   # `Y` and `T` are required; see `.resolve_cf_api()`.
                   Y          = NULL,    # outcome variable name (character scalar)
                   T          = NULL,    # treatment variable name (character scalar)
                   X_baseline = NULL,    # baseline (non-omics): clinical / demographic
                   omics      = NULL,    # which omics layer(s) / features
                   omics_type = NULL,    # layer type override (rnaseq / other)
                   # PS-specific covariate subset. NULL = use full X for PS;
                   # character vector = use that subset (applied per
                   # PS-using method).
                   propensity_X = NULL,
                   n_bootstrap = 100L,
                   bootstrap_ci = FALSE,
                   random_state = 0L,
                   # ---- preprocessing + feature selection ----
                   preprocess     = TRUE,                # TMM -> log-CPM
                   batch          = NULL,                # ComBat batch column (log values)
                   # Default feature selection = "variance" (top-n by
                   # row variance) so SE/MAE omics assays are actually included
                   # as covariates by default. df without an assay is a no-op.
                   feature_method = c("variance", "none", "pca", "pathway"),
                   # NULL resolves to an input-type-aware default —
                   # 100 for SE / data.frame, 50 for MAE (per-assay budget,
                   # multiple assays). An explicit value overrides both.
                   n_features     = NULL,
                   gene_sets      = NULL,                # for pathway; NULL = Hallmark 50
                   pathway_method = c("gsva", "ssgsea", "plage", "zscore"),
                   partition         = c("oof", "none"),
                   n_folds           = 5L,
                   eval_nuisance     = "super_learner",
                   eval_library      = NULL,      # Super Learner members
                   eval_propensity_score = NULL,  # known eval PS (RCT), symmetric with cf_compare
                   keep_kfold_models = FALSE,
                   # ---- GPU device ----
                   device            = c("auto", "cpu", "gpu"),
                   # ---- PS clip ----
                   # NULL so the dispatcher can tell an explicit value from an
                   # unset one and report a `clip` that cannot reach the chosen
                   # method; the four fitters that bound a propensity carry the
                   # same default, so the resolved value is unchanged.
                   # PS clip applied by PS-using methods (grf/bcf/drlearner/xlearner);
                   # PS-free methods ignore it.
                   clip              = NULL,
                   # ---- known PS specification ----
                   # NULL = estimate (default, observational study) / scalar (RCT 50/50) /
                   # vector length-n (stratified RCT, survey). Mutually exclusive with an
                   # explicit `model_propensity` (cateomics_propensity_args_conflict error).
                   propensity_score  = NULL,
                   # ---- binary outcome support (risk difference / classifier).
                   # "auto" (default) detects a 0/1 outcome and switches every
                   # method to risk-difference / classifier outcome models,
                   # skipping the Y z-score for NN backends. Force with
                   # "binary" / "continuous".
                   outcome_type = c("auto", "continuous", "binary"),
                   # ---- base-learner tokens (which algorithm) ----
                   model_outcome    = NULL,
                   model_final      = NULL,
                   model_propensity = NULL,
                   # ---- hyperparameters (how that algorithm is configured) ----
                   params            = NULL,
                   propensity_params = NULL,
                   eval_params       = NULL) {
  # Every penalised fit this call makes reports its thin split through one
  # line, named for the learners it saw.
  .cf_cv_thin_open()
  on.exit(.cf_cv_thin_close(), add = TRUE)
  # `.resolve_cf_api()` turns the names into the canonical formula the rest of
  # the pipeline parses, plus the `structured` bag the result `meta` records.
  .facade <- .resolve_cf_api(
    Y          = Y,
    T          = T,
    X_baseline = X_baseline,
    omics      = omics,
    omics_type = omics_type,
    data       = data
  )
  formula <- .facade$formula
  .cateomics_named_args <- .facade$structured
  .cateomics_omics_spec <- .facade$omics_spec
  if (is.null(data)) {
    rlang::abort("`data` argument is required.", class = "cateomics_invalid_data")
  }
  # propensity_X-ignored warning. We need `method`
  # in canonical form for the lookup but `match.arg(method)` happens further
  # below; use a lightweight pre-resolve here. The lookup can fail on an
  # unknown method name, which match.arg reports later — swallow it silently.
  .method_pre <- tryCatch(match.arg(method), error = function(e) NA_character_)
  if (!is.na(.method_pre)) {
    .warn_propensity_X_ignored(.method_pre, propensity_X)
  }

  # `feature_method` / `pathway_method` are character options; match.arg gives
  # the canonical first-element default when callers omit the argument.
  feature_method <- match.arg(feature_method)
  pathway_method <- match.arg(pathway_method)

  # Resolve the input-type-aware n_features default. NULL (the
  # default) -> 50 for MAE (per-assay budget, multiple assays) / 100 otherwise
  # (SE / data.frame). An explicit numeric value is honoured as-is.
  n_features <- .cf_resolve_n_features(n_features, data)

  partition <- match.arg(partition)
  device <- match.arg(device)  # GPU device validation
  outcome_type <- match.arg(outcome_type)  # binary outcome routing

  .cf_check_clip(clip)

  # Known PS validation + mutual-exclusion constraint.
  if (!is.null(propensity_score)) {
    if (!is.numeric(propensity_score) ||
        any(propensity_score <= 0 | propensity_score >= 1, na.rm = TRUE)) {
      rlang::abort(
        "`propensity_score` must be numeric with values strictly in (0, 1).",
        class = "cateomics_invalid_propensity_score"
      )
    }
    # For SummarizedExperiment / MultiAssayExperiment, samples live in columns
    # (`ncol`) and features in rows; for data.frame the samples are in rows.
    # Use the sample-count interpretation consistently so `propensity_score`
    # vector length matches the number of subjects, not the number of features.
    n_obs <- if (inherits(data, c("SummarizedExperiment", "MultiAssayExperiment"))) {
      ncol(data)
    } else {
      nrow(data)
    }
    if (length(propensity_score) != 1L &&
        length(propensity_score) != n_obs) {
      rlang::abort(
        sprintf("`propensity_score` length %d does not match number of samples = %d.",
                length(propensity_score), n_obs),
        class = "cateomics_invalid_propensity_score"
      )
    }
    if (!is.null(model_propensity)) {
      rlang::abort(
        c("Both `propensity_score` and `model_propensity` are explicitly specified.",
          i = "For RCT (known PS), use `propensity_score = ...` only.",
          i = "For observational study (estimated PS), use `model_propensity = ...` only."),
        class = "cateomics_propensity_args_conflict"
      )
    }
    # Caveat: drlearner / xlearner internally K-fold
    # cross-fit via econml; FixedPropensityModel.predict_proba(X)
    # cannot recover the row-level subset for vector PS, so vector PS is
    # rejected with an actionable hint pointing at scalar / bcf.
    if (length(propensity_score) > 1L &&
        isTRUE(.method_pre %in% c("drlearner", "xlearner"))) {
      rlang::abort(
        c(sprintf("Vector `propensity_score` is unsupported for method = '%s'.",
                  .method_pre),
          i = "econml cross-fit pipelines cannot map per-sample PS through internal folds.",
          i = "Use a scalar (`propensity_score = 0.5`) or switch to method = 'bcf' for per-sample PS."),
        class = "cateomics_propensity_vector_unsupported"
      )
    }
  }

  # Device resolution: under "auto", probe GPU availability. GPU-unsupported
  # methods = R-native 2 (grf/bcf) + Python CPU-only 4 (slearner / tlearner /
  # xlearner / drlearner: econml = sklearn CPU); device='gpu' on these warns +
  # continues on CPU. GPU-capable = Python NN 4 (cevae/ganite/tarnet/dragonnet),
  # decided backend-side via `torch.cuda.is_available()`.
  gpu_unsupported_methods <- .CF_GPU_UNSUPPORTED
  # Use the pre-resolved scalar method (.method_pre); the formal `method` is
  # still the full default vector until match.arg() below, and passing a vector
  # into the `%in%` device checks would error on `if()` (length > 1).
  device_resolved <- .resolve_device(device, .method_pre, gpu_unsupported_methods)
  if (identical(device, "gpu") && isTRUE(.method_pre %in% gpu_unsupported_methods)) {
    rlang::warn(
      c(
        sprintf("device='gpu' is unsupported for method '%s' (CPU-only backend).", .method_pre),
        i = "Falling back to CPU. Use Python NN methods (cevae/tarnet/dragonnet/ganite) for GPU acceleration."
      ),
      class = "cateomics_gpu_unsupported_warning"
    )
  }

  .cf_check_eval_nuisance(eval_nuisance)
  # The estimator-side and evaluation-side known propensities are separate
  # settings, so one of them alone leaves the other estimated.
  .cf_inform_known_ps_asymmetry(propensity_score, eval_propensity_score)

  # `eval_library` only describes the Super Learner. Passing it beside a single
  # named learner is a setting that cannot reach the fit, which is the class of
  # argument this package refuses to leave silent.
  if (!is.null(eval_library) && !.cf_sl_any_slot_is_sl(eval_nuisance)) {
    rlang::warn(
      c("`eval_library` is not used: `eval_nuisance` names a single learner.",
        i = sprintf("Set `eval_nuisance = \"%s\"` for the library to apply.",
                    .CF_SL_TOKEN)),
      class = "cateomics_eval_library_unused"
    )
  }

  if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L ||
      n_folds != as.integer(n_folds)) {
    rlang::abort(
      sprintf("`n_folds` must be a single integer >= 2 (got: %s).",
              paste(deparse(n_folds), collapse = " ")),
      class = "cateomics_invalid_n_folds"
    )
  }
  n_folds <- as.integer(n_folds)

  if (!is.logical(keep_kfold_models) || length(keep_kfold_models) != 1L) {
    rlang::abort(
      "`keep_kfold_models` must be a single logical value.",
      class = "cateomics_invalid_keep_kfold_models"
    )
  }

  method <- match.arg(method)

  # Resolve the hyperparameter block for this method. `params` may be keyed by
  # method (`list(grf = list(num.trees = 4000))`) or given flat
  # (`list(num.trees = 4000)`); unknown keys and unknown hyperparameters abort
  # here rather than being silently dropped further down.
  .cf_params_check_keys(params, methods = method)
  .cf_params_check_unused(params, methods = method)
  .method_params <- .cf_resolve_params(params, method = method,
                                       methods = method)
  .cf_params_check(method, .method_params)
  propensity_params <- .cf_params_shape(propensity_params, "propensity_params")
  # The propensity block reaches a scikit-learn estimator, whose random_state
  # is in get_params(), so without this the seed is only caught in Python --
  # after the environment boots, and as a ValueError rather than a CATEomics
  # condition.
  .cf_params_reject_rng(propensity_params, "propensity_params")
  eval_params       <- .cf_params_shape(eval_params, "eval_params")
  if (length(propensity_params) &&
      !method %in% c("xlearner", "drlearner", "bcf")) {
    rlang::abort(
      c(sprintf("Method '%s' has no pluggable propensity model, so `propensity_params` would have no effect.",
                method),
        i = if (identical(method, "grf"))
              "grf estimates the propensity with its own forest: use params = list(grf = list(nuisance = list(...)))."
            else "It applies to xlearner / drlearner / bcf."),
      class = "cateomics_propensity_params_unused"
    )
  }
  # bcf fits its propensity in R through the shared base learner, so the names
  # are checked here against that learner. The metalearners and the DR learner
  # fit theirs in Python, where scikit-learn raises on an unknown name itself.
  if (identical(method, "bcf") && length(propensity_params)) {
    .cf_learner_params_check(model_propensity %||% "logistic",
                             propensity_params, "propensity_params",
                             type = "classification")
  }

  # CI policy, decided before any fit work is done. `.cf_check_ci_policy()`
  # raises `cateomics_bootstrap_unsupported_nn` for tarnet/dragonnet/ganite and
  # `cateomics_bootstrap_unsupported` for any other method whose fitter has no
  # `bootstrap_ci` formal (grf/bcf/cevae): those backends never receive the
  # request, so it is refused here rather than dropped on the way down.
  # It also raises `cateomics_bootstrap_overrides_native` for a method that does
  # resample and has a native interval the request replaces.
  .cf_check_ci_policy(method, bootstrap_ci, model_final)

  # ===================================================================
  # cf_fit 7-step flow
  # ===================================================================

  # partition dispatch sits between Step 5 and Step 6. Step 1-5 (log-CPM + batch +
  # feature_method + factor expansion) always run, and both the K-fold loop and
  # the full-data fit operate on the preprocessed X. This guarantees:
  #   - SE/MAE/df + preprocess=TRUE => Step 2 log-CPM is always applied
  #   - K-fold loop and full-data fit run on the same preprocessed X
  #   - meta slots (input_type / preprocess) are
  #     populated in Step 7
  # Step 1-7 run regardless of `partition` (Step 5.5 = the partition='oof'
  # dispatch is the only branch; Step 6-7 are shared by both partitions).

  # Steps 1-5: shared preprocessing (input validation + assay resolution ->
  # log-normalisation -> ComBat -> feature select/aggregation -> factor expand
  # -> assemble X = cbind(X_clinical, t(processed_features))). Extracted to
  # `.cf_preprocess_input()` so cf_fit and cf_compare apply the identical
  # Step 1-5 pipeline to the same SE/MAE/df input (single source of truth — no
  # drift). Binary-outcome resolution + method dispatch stay in cf_fit below.
  # `parsed` / `step3_5` (feature select) / `input_type` /
  # `factor_levels` feed the Step 7 meta-populate path; the normalisation and
  # pathway records are read off `prep` by `.cf_preprocess_record()`.
  prep <- .cf_preprocess_input(
    formula        = formula,
    data           = data,
    omics_spec     = .cateomics_omics_spec,
    preprocess     = preprocess,
    batch          = batch,
    feature_method = feature_method,
    n_features     = n_features,
    gene_sets      = gene_sets,
    pathway_method = pathway_method,
    method         = method
  )
  input_type    <- prep$input_type
  factor_levels <- prep$factor_levels
  X             <- prep$X
  T             <- prep$T
  Y             <- prep$Y
  cov_names     <- prep$cov_names
  # `propensity_X` is checked against the design matrix preprocessing produced,
  # since that is what both propensities are fitted on.
  .cf_check_propensity_X(propensity_X, colnames(X) %||% cov_names)
  # Unnamed, the propensity goes on the clinical covariates alone. Resolved
  # here, after the warning above has already seen what the caller actually
  # wrote: that warning is about an argument not reaching a method, and the
  # default reaches nobody's attention.
  propensity_X <- .cf_default_propensity_X(propensity_X,
                                           prep$clinical_cov_names,
                                           propensity_score,
                                           eval_propensity_score)

  # Resolve binary vs continuous outcome ONCE (before dispatch), inform
  # the user on auto-detection, and thread `binary_outcome` to every backend +
  # the evaluation nuisance so binary outcomes use risk-difference / classifier
  # models and NN backends skip the Y z-score.
  binary_outcome <- .resolve_binary_outcome(Y, outcome_type)
  # Checked here rather than beside the other `params` validation because `T`
  # only becomes the treatment vector after the input steps; up there it is
  # still the column name, and `as.integer()` of that is NA, which made the
  # check pass silently. The same helper runs from `cf_compare`, so a
  # configuration refused by one entry is refused by the other.
  .cf_check_method_splits(
    method, T, Y, binary_outcome,
    n_folds = if (identical(partition, "oof")) n_folds else NULL,
    method_params = .method_params,
    model_outcome = model_outcome, model_propensity = model_propensity,
    model_final = model_final)

  # Wire the known evaluation propensity. `eval_nuisance` stays the user-facing
  # base-learner name; the list spec exists only for the nuisance helpers, and
  # cf_compare builds it from the same helper so the two entries cannot drift.
  eval_nuisance_spec <- .cf_eval_nuisance_spec(eval_nuisance,
                                               eval_propensity_score)

  # ===================================================================
  # Step 5.5: partition dispatch
  # ===================================================================
  # partition='oof' runs `.kfold_oof_predict_matrix()` for per-fold method
  # training + a full-data fit on the preprocessed X. The returned
  # full_inner$fit has the same shape as the partition='none' `fit`, so both
  # merge into the shared Step 6-7 path. oof_ite / nuisance / fold_cache are
  # attached below.
  if (identical(partition, "oof")) {
    # A fold's warning is about this fit, not about this fold, so the caller
    # hears each distinct one once instead of once per fold. A wording that
    # differs between folds is a different warning and still comes through.
    full_inner <- .cf_dedup_conditions(.kfold_oof_predict_matrix(
      X = X, T_vec = T, Y_vec = Y, method = method,
      cov_names         = cov_names,
      n_folds           = n_folds,
      random_state      = random_state,
      keep_kfold_models = keep_kfold_models,
      eval_nuisance     = eval_nuisance_spec,
      eval_library      = eval_library,
      n_bootstrap       = n_bootstrap,
      bootstrap_ci      = bootstrap_ci,
      device            = device_resolved,
      propensity_X      = propensity_X,
      propensity_score  = propensity_score,
      clip              = clip,
      eval_propensity_score = eval_propensity_score,
      binary_outcome    = binary_outcome,
      model_outcome     = model_outcome,
      model_final       = model_final,
      model_propensity  = model_propensity,
      propensity_params = propensity_params,
      eval_params       = eval_params,
      params            = .method_params
    ))
    fit <- full_inner$fit
    oof_ite_value    <- full_inner$oof_ite
    oof_ps_value     <- full_inner$oof_ps
    oof_fold_id      <- full_inner$fold_id
    nuisance_full    <- full_inner$nuisance
    fold_cache_value <- full_inner$fold_cache
  } else {
    # partition = 'none' (Step 5 switch dispatch, in-sample fit)
    oof_ite_value    <- NULL
    oof_ps_value     <- NULL
    oof_fold_id      <- NULL
    fold_cache_value <- NULL
    nuisance_full    <- NULL  # full-data nuisance computed below

  fit <- .cf_dispatch_fit(
    method = method, X = X, T = T, Y = Y, cov_names = cov_names,
    control = list(
      random_state     = as.integer(random_state),
      device           = device_resolved,
      binary_outcome   = binary_outcome,
      bootstrap_ci     = isTRUE(bootstrap_ci),
      n_bootstrap      = as.integer(n_bootstrap),
      clip             = clip,
      propensity_score = propensity_score,
      propensity_X     = propensity_X,
      model_outcome    = model_outcome,
      model_final      = model_final,
      model_propensity = model_propensity,
      propensity_params = propensity_params
    ),
    params = .method_params
  )

    # partition='none' path: fit `eval_nuisance` once on the full data and
    # populate `res$nuisance = list(mu0_oof, mu1_oof, pi_oof, m_oof)` -- the
    # same four elements the OOF path returns, from in-sample predictions.
    nuisance_full <- .compute_full_data_nuisance(
      X, T, Y,
      eval_nuisance         = eval_nuisance_spec,
      random_state          = random_state,
      eval_propensity_score = eval_propensity_score,
      binary_outcome        = binary_outcome,
      eval_params           = eval_params,
      eval_library          = eval_library,
      n_folds               = n_folds,
      propensity_X          = propensity_X
    )
  } # end of partition='none' branch

  # Step 7: Build meta + result. x_zscore / y_zscore are populated per method
  # (cf_fit reads them from fit$res here).
  preprocess_record <- .cf_preprocess_record(prep)

  # Create standardized meta object.
  meta <- create_cf_meta(
    n = length(Y),
    p = ncol(X),
    covariate_names = cov_names,
    method = method,
    bootstrap = bootstrap_ci,  # meta records whether a bootstrap CI was requested
    random_state = random_state,
    # The PCA train-fold rotation, centre, scale and gene order travel in the
    # layer's own record inside this, so `.apply_newdata_preprocessing()` can
    # re-apply the same projection to newdata without recomputing on test
    # samples (data-leakage guard).
    preprocess     = preprocess_record,
    input_type     = input_type,
    factor_levels  = factor_levels,
    # X z-score is applied centrally in `.cf_preprocess_input` (every column,
    # including binary indicators), so the train stats come from
    # `prep$x_zscore` for all methods and
    # cf_predict.R `.apply_newdata_preprocessing` re-applies them to newdata.
    # `y_zscore` stays NN-only: each NN backend standardises Y on the data it
    # receives and inverse-transforms its predictions, returning `res$y_zscore`.
    x_zscore       = prep$x_zscore,
    y_zscore       = if (is.list(fit) && !is.null(fit$res$y_zscore)) fit$res$y_zscore else NULL,
    # The names the caller gave. `.apply_newdata_preprocessing()` reads
    # `named_args$X_baseline` to model.matrix one-hot expand the clinical
    # columns of newdata, so this is the record the newdata path relies on.
    named_args   = .cateomics_named_args,
    propensity_X = propensity_X
  )

  # ===================================================================
  # Step 7: assemble the four-pillar cf_model.
  # ===================================================================
  # Flatten the double-fit nesting: each cf_fit_<method>() returns
  # list(model=, [nuisance siblings], res=list(...)); lift `res` to the
  # top-level `res` pillar and keep the remaining backend container (model +
  # siblings such as grf's forest.Y / forest.W) as the `model` pillar.
  res_pillar   <- fit$res
  model_pillar <- fit
  model_pillar$res <- NULL

  # One canonical place for the interval's provenance, the same here and for a
  # comparison's candidates.
  res_pillar <- .cf_lift_ci_method(res_pillar)
  meta$estimator$ci_method <- res_pillar$ci_method

  # OOF tau + shared nuisance move into the res pillar. `fold_id` goes with
  # them: it says which fold held each row out, so it indexes `tau_hat_oof` and
  # `nuisance` the same way they index the rows. `NULL` under
  # `partition = "none"`, where there is no split.
  res_pillar$tau_hat_oof <- oof_ite_value
  res_pillar$fold_id     <- oof_fold_id

  # The propensity a downstream diagnostic should look at is the held-out one
  # whenever every row has it; the full-model PS (res$ps, set by the per-method
  # fit) is the fallback, and res$ps_type says which of the two this is. A known
  # propensity is known on every row, so the held-out copy is the same number
  # and the label keeps saying where it came from.
  if (!is.null(oof_ps_value) && !anyNA(oof_ps_value) &&
      !identical(res_pillar$ps_type, "known")) {
    res_pillar$ps      <- oof_ps_value
    res_pillar$ps_type <- "oof"
  }
  eval_sl_weights <- nuisance_full$sl_weights
  if (!is.null(nuisance_full)) {
    nuisance_full$sl_weights <- NULL
    res_pillar$nuisance <- nuisance_full
  }

  # The z-score statistics live in meta, so the res copies go.
  res_pillar$x_zscore <- NULL
  res_pillar$y_zscore <- NULL

  # Recompute the canonical IQR summary centrally so every method + partition
  # path produces an identical summary shape.
  res_pillar$summary <- .cf_build_summary(
    ate     = res_pillar$ate,
    ate_ci  = res_pillar$ate_ci,
    tau_hat = res_pillar$tau_hat
  )

  meta$estimator$partition <- partition
  meta$estimator$outcome_type <- if (isTRUE(binary_outcome)) "binary" else "continuous"
  # Which evaluation nuisance produced `res$nuisance`, and for the Super
  # Learner the members and the weights its combination used. The weights are
  # not nuisance vectors, so they live here rather than in `res$nuisance`.
  # A backend that logs its training losses returns them on `res`; the record
  # belongs with the other estimator settings, not among the result vectors.
  if (!is.null(res_pillar$training_trace)) {
    meta$estimator$training_trace <- res_pillar$training_trace
    res_pillar$training_trace <- NULL
  }
  # Only the resolved spec is recorded. The token the caller passed is
  # recoverable from it (one learner per slot, the same name in both when a
  # single token was given), while the reverse is not: a known propensity
  # turns that slot into "known", which the token does not show.
  meta$estimator$eval_nuisance_resolved <-
    .cateomics_resolve_eval_nuisance(eval_nuisance_spec)
  if (.cf_sl_any_slot_is_sl(eval_nuisance_spec)) {
    meta$estimator$eval_library <- .cf_sl_libraries_used(eval_nuisance_spec,
                                                         eval_library)
    meta$estimator$eval_sl_weights <- eval_sl_weights
  }
  # Record the hyperparameters this fit actually ran with, so a result can be
  # reproduced without re-deriving the three-stage resolution.
  # The slots record what actually ran, not the request: the public formals
  # default to NULL and each method's real default lives on its fitter, so
  # storing the argument would leave the default path blank.
  .slots <- .cf_spec_slots(method, model_outcome, model_final,
                           model_propensity, propensity_score, clip)
  meta$spec <- list(
    method            = method,
    clip              = .slots$clip,
    model_outcome     = .slots$model_outcome,
    model_final       = .slots$model_final,
    model_propensity  = .slots$model_propensity,
    params            = .method_params,
    propensity_params = propensity_params,
    eval_params       = eval_params
  )

  result <- structure(
    list(
      method = method,
      res    = res_pillar,   # pillar 1: numeric results
      model  = model_pillar, # pillar 2: fitted backend container
      meta   = meta,         # pillar 3: metadata
      # pillar 4: model-visible training data (cf_interpret surrogate SHAP /
      # cf_predict newdata restore / cf_compare-evaluate-rule nuisance
      # recomputation).
      data   = list(
        X_preprocessed = X,
        Y              = Y,
        T              = T
      )
    ),
    class = "cf_model"
  )
  if (!is.null(fold_cache_value)) {
    result$fold_cache <- fold_cache_value
  }

  # Device record: what the caller asked for (`device`, before resolution) vs
  # the device the backend reports.
  result$meta$provenance <- .cf_record_device(
    result$meta$provenance, device, device_resolved,
    res_pillar$device_effective)
  result$meta$provenance$sample_names_meaningful <-
    isTRUE(prep$sample_names_meaningful)

  # Memory warning.
  obj_bytes <- as.numeric(utils::object.size(result))
  if (obj_bytes > 100 * 1024^2) {
    rlang::warn(
      c(
        sprintf("cf_model object is large (%.1f MB).", obj_bytes / 1024^2),
        i = "The preprocessed X is retained in result$data for downstream cf_evaluate() / cf_interpret() / cf_rule()."
      ),
      class = "cateomics_large_object_size"
    )
  }

  # Diagnostics are on-demand, not auto-attached. cf_fit does not compute
  # E-value / balance / overlap into the result meta; run cf_diagnostics(model)
  # after fitting for overlap, covariate balance (SMD Love), requested-variable
  # overlap, and the E-value.

  # Validate result structure (contract enforcement). Downstream functions
  # (cf_compare/select/evaluate/interpret) rely on fields like tau_hat_oof /
  # ate_ci / tau_ci, so contract violations abort here rather
  # than surfacing only at NPE time in downstream code.
  validate_cateomics_result(result)
  # Checked on the result the caller receives, not inside the dispatch: with
  # `partition = "oof"` the dispatch runs once per fold as well, and those
  # intermediate fits are not what anyone reads.
  .cf_warn_ci_excludes_estimate(result)

  result
}


#' Resolve batch vector from data based on input_type
#'
#' Internal helper for cf_fit Step 3. SE = colData[[batch]] / MAE =
#' colData[complete_primaries, ][[batch]] / df = data[[batch]].
#'
#' @keywords internal
#' @noRd
.resolve_batch_vector <- function(data, batch, input_type, step1) {
  vec <- switch(
    input_type,
    "SE" = {
      cd <- as.data.frame(SummarizedExperiment::colData(data),
                          optional = TRUE)
      cd[[batch]]
    },
    "MAE" = {
      # `optional = TRUE`: a batch column can be named anything.
      cd <- as.data.frame(MultiAssayExperiment::colData(data),
                          optional = TRUE)
      primaries <- step1$mae_complete_primaries
      if (is.null(primaries)) {
        cd[[batch]]
      } else {
        cd[primaries, , drop = FALSE][[batch]]
      }
    },
    "data.frame" = data[[batch]]
  )
  if (is.null(vec)) {
    rlang::abort(
      sprintf("Batch column '%s' not found in colData / data.", batch),
      class = "cateomics_batch_column_missing"
    )
  }
  vec
}

#' Resolve whether the outcome is binary
#'
#' Detects a 0/1 outcome under \code{outcome_type = "auto"} (exactly two
#' distinct non-NA values, both in \code{{0, 1}}), or honours the explicit
#' \code{"binary"} / \code{"continuous"} override. Emits a one-per-session
#' \code{cateomics_binary_outcome_detected} message on auto-detection.
#'
#' @param Y Numeric outcome vector.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"}, \code{"binary"}.
#' @return Logical scalar: \code{TRUE} when the outcome should be treated as binary.
#' @keywords internal
#' @noRd
.resolve_binary_outcome <- function(Y, outcome_type = "auto") {
  y <- Y[!is.na(Y)]
  detected <- length(unique(y)) == 2L && all(y %in% c(0, 1))
  bin <- switch(outcome_type,
    "auto"       = detected,
    "binary"     = TRUE,
    "continuous" = FALSE,
    detected
  )
  if (identical(outcome_type, "auto") && isTRUE(bin)) {
    rlang::inform(
      c("Binary outcome detected (Y in {0, 1}).",
        i = "Using risk-difference / classifier outcome models; NN backends skip the Y z-score.",
        i = "Pass `outcome_type = 'continuous'` to force regression."),
      class = "cateomics_binary_outcome_detected"
    )
  }
  isTRUE(bin)
}

