#' Render the CALMeR diagnostic report
#'
#' Renders a parameterized Rmarkdown report containing model diagnostics
#' (residual QQ plots, partial effects, calibration), training/scoring
#' summary statistics, and a heatmap of where flags concentrate across
#' n_isoforms x pct_covered.
#'
#' @param fit Object returned by [fit_beinf()].
#' @param scored Data frame returned by [score_reads()].
#' @param out_file Path to the rendered report.
#' @param format Either `"pdf"` (default) or `"html"`.
#'
#' @return Invisibly, the path to the rendered file.
#'
#' @export
calmer_report <- function(fit,
                          scored,
                          out_file,
                          format = c("pdf", "html")) {
  stop("Not yet implemented.")
}
