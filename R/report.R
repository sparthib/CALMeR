#' Render the CALMeR diagnostic report
#'
#' Renders a parameterized Rmarkdown report containing model diagnostics
#' on the training data and (if `scored` is supplied) summaries of where
#' flags landed in the scored data set.
#'
#' Sections:
#' \enumerate{
#'   \item Training summary (predictor distributions, MAPQ histogram).
#'   \item Model fit (formulas, EDFs, pseudo-R^2).
#'   \item Residual diagnostics (global QQ + density vs N(0,1), RQR
#'         density faceted by each predictor's quantile bins).
#'   \item Partial effects on \eqn{E[\mathrm{MAPQ}]}.
#'   \item Coefficient forest plot for all BEINF parameters.
#'   \item Scoring summary — only if `scored` is supplied.
#'   \item Where flags landed — only if `scored` is supplied.
#' }
#'
#' @param fit Object returned by [fit_beinf()].
#' @param scored Optional data frame returned by [score_reads()].
#' @param out_file Path to the rendered report.
#' @param format Either `"pdf"` (default) or `"html"`.
#' @param sample_name Optional string displayed in the report title.
#'
#' @return Invisibly, the absolute path to the rendered file.
#'
#' @examples
#' if (rmarkdown::pandoc_available()) {
#'   ex       <- calmer_example_data()
#'   features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#'   fit <- fit_beinf(
#'     features, n = 100,
#'     formula       = y ~ n_isoforms_sc + tx_length_sc + gc_content_sc + pct_covered_sc,
#'     sigma_formula = ~ 1, nu_formula = ~ 1, tau_formula = ~ 1
#'   )
#'   scored <- score_reads(features, fit, filter_method = "mad", threshold = -3)
#'   out <- file.path(ex$dir, "report.html")
#'   calmer_report(fit, scored, out_file = out, format = "html")
#'   file.exists(out)
#'   unlink(ex$dir, recursive = TRUE)
#' }
#'
#' @export
calmer_report <- function(fit,
                          scored      = NULL,
                          out_file,
                          format      = c("pdf", "html"),
                          sample_name = NULL) {
  format <- match.arg(format)

  if (!inherits(fit, "calmer_fit")) {
    stop("`fit` must be a calmer_fit object returned by fit_beinf().")
  }
  if (!is.null(scored) && !is.data.frame(scored)) {
    stop("`scored` must be a data frame returned by score_reads(), or NULL.")
  }

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render the report.")
  }

  template <- system.file("rmd", "diagnostic_report.Rmd", package = "CALMeR")
  if (!nzchar(template)) {
    stop("Could not find the diagnostic_report.Rmd template inside the ",
         "installed CALMeR package.")
  }

  output_format <- switch(format,
                          pdf  = "pdf_document",
                          html = "html_document")

  out_file <- normalizePath(out_file, mustWork = FALSE)
  rmarkdown::render(
    input         = template,
    output_format = output_format,
    output_file   = out_file,
    params        = list(fit         = fit,
                         scored      = scored,
                         sample_name = sample_name),
    envir         = new.env(parent = globalenv()),
    quiet         = TRUE
  )

  invisible(out_file)
}
