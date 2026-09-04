# Declare global variables used in ggplot2 NSE to avoid R CMD check NOTEs.
# The cf_plot_* family quotes column names in per-call ggplot2 `aes()`; declaring
# them keeps `R CMD check --as-cran` clean.
utils::globalVariables(c(
  "ITE", "balanced", "ci_lower", "ci_upper", "gate", "group", "ps", "smd", "variable",
  # forest / compare family
  "method", "estimate", "lower", "upper",
  # cf_plot_evaluate (blp calibration)
  "predicted_mean", "observed_mean",
  # cf_plot_evaluate (rate)
  "rate",
  # cf_plot_policy_curve
  "parameter", "policy_value",
  # cf_plot_compare
  "metric", "value"
))

#' @import methods
#' @importFrom stats model.frame as.formula var glm binomial predict rnorm rbinom setNames median reorder prcomp complete.cases plogis sd quantile lm coef confint
#' @importFrom Formula Formula
#' @importFrom rlang abort warn inform .data
#' @importFrom cli cli_alert_success cli_alert_warning cli_alert_info cli_abort cli_h1 cli_h2 cli_text cli_progress_step col_green col_red col_yellow
#' @importFrom glue glue
#' @importFrom utils askYesNo install.packages packageVersion head
#' @importFrom graphics hist legend lines par plot abline
#' @importFrom grDevices rgb
#' @keywords internal
"_PACKAGE"

# NOTE: S3 methods are registered in NAMESPACE via roxygen2 @export tags.
# No .onLoad needed for S3 method registration.
# Python initialization is deliberately NOT performed at package load.
# Python-based methods provision their environment automatically via basilisk
# on first use (in a separate process); no user setup is required.
