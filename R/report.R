#' Render the CALMeR diagnostic report
#'
#' Renders a parameterized Rmarkdown report containing model diagnostics
#' on the training data and (if `scored` is supplied) summaries of where
#' flags landed in the scored data set.
#'
#' Sections:
#' \enumerate{
#'   \item Training summary (predictor distributions, MAPQ histogram).
#'   \item Model comparison — only if a `calmer_suite` is supplied.
#'         AIC / BIC / global EDF / pseudo-R^2 across all fits in the
#'         suite, side-by-side residual QQ panels.
#'   \item Model fit for the primary model (formulas, EDFs, pseudo-R^2).
#'   \item Residual diagnostics (global QQ + density vs N(0,1), RQR
#'         density faceted by each predictor's quantile bins).
#'   \item Partial effects on \eqn{E[\mathrm{MAPQ}]}.
#'   \item Coefficient forest plot for all BEINF parameters.
#'   \item Scoring summary — only if `scored` is supplied.
#'   \item Where flags landed — only if `scored` is supplied.
#' }
#'
#' @param fit Either a single `calmer_fit` from [fit_beinf()] or a
#'   `calmer_suite` from [fit_beinf_suite()]. When a suite is supplied,
#'   the model-comparison section is rendered and the suite's primary
#'   fit drives the per-model diagnostic sections.
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
#'   clean    <- clean_features(features)
#'   fit <- fit_beinf(
#'     clean, n = 1000,
#'     formula       = y ~ n_isoforms_sc + tx_length_sc + gc_content_sc + pct_covered_sc,
#'     sigma_formula = ~ 1, nu_formula = ~ 1, tau_formula = ~ 1
#'   )
#'   scored <- score_reads(features, fit, filter_method = "mad", threshold = -3)
#'   out <- tempfile(fileext = ".html")
#'   calmer_report(fit, scored, out_file = out, format = "html")
#'   file.exists(out)
#' }
#'
#' @export
calmer_report <- function(fit,
                          scored      = NULL,
                          out_file,
                          format      = c("pdf", "html"),
                          sample_name = NULL) {
  format <- match.arg(format)

  if (!inherits(fit, "calmer_fit") && !inherits(fit, "calmer_suite")) {
    stop("`fit` must be a calmer_fit (or calmer_suite) object.")
  }
  if (!is.null(scored) && !is.data.frame(scored)) {
    stop("`scored` must be a data frame returned by score_reads(), or NULL.")
  }

  suite <- if (inherits(fit, "calmer_suite")) fit else NULL
  primary <- primary_fit(fit)

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
    params        = list(fit         = primary,
                         suite       = suite,
                         scored      = scored,
                         sample_name = sample_name),
    envir         = new.env(parent = globalenv()),
    quiet         = TRUE
  )

  invisible(out_file)
}
