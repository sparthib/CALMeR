#' Render the CALMeR diagnostic report
#'
#' Renders a parameterized Rmarkdown report containing model diagnostics
#' on the training data and (if a scored data set is available) summaries
#' of where flags landed in the scored data set.
#'
#' When `fit` is a `calmer_selection` from [select_beinf()], the model
#' it already picked (minimum BIC among the full model and the
#' backward-selected Strategy A / B models) drives all subsequent
#' sections (residual diagnostics, partial effects, coefficient forest
#' plot) and is used to re-score the data via [score_reads()], so the
#' flag stats in the report match the picked model.
#'
#' When `fit` is a single `calmer_fit`, the fit is used as supplied.
#'
#' Sections:
#' \enumerate{
#'   \item Training summary (predictor distributions, MAPQ histogram).
#'   \item Picked-model fit diagnostics (formulas, EDFs).
#'   \item Residual diagnostics (global QQ + density vs N(0, 1), RQR
#'         density faceted by each predictor's bins).
#'   \item Partial effects on \eqn{E[\mathrm{MAPQ}]}.
#'   \item Coefficient forest plot — 95% Wald intervals for the
#'         non-spline terms in each BEINF distributional parameter
#'         (mu / sigma / nu / tau).
#'   \item Scoring summary (only if `scored` is available).
#'   \item Where flags landed: per-`n_isoforms` flag breakdown and
#'         per-predictor flagged-read distribution.
#' }
#'
#' @param fit Either a single `calmer_fit` from [fit_beinf()] or a
#'   `calmer_selection` from [select_beinf()].
#' @param features Cleaned training data frame from [clean_features()].
#'   When supplied, the report re-scores against the picked fit. If
#'   `NULL`, the report falls back to the pre-computed `scored` argument.
#' @param scored Optional data frame returned by [score_reads()]. Used
#'   only when `features` is `NULL` (legacy path). If both `features`
#'   and `scored` are supplied, `features` wins and the report
#'   re-scores internally.
#' @param out_file Path to the rendered report.
#' @param format Either `"pdf"` (default) or `"html"`.
#' @param sample_name Optional string displayed in the report title.
#' @param score_args Named list of extra arguments to forward to
#'   [score_reads()] when re-scoring internally (e.g.
#'   `list(filter_method = "mad", threshold = -3)`).
#'
#' @return Invisibly, the absolute path to the rendered file, with
#'   attributes `"picked_model"` (the name of the picked model, if a
#'   selection was supplied) and `"scored"` (the data frame the report
#'   rendered against).
#'
#' @examples
#' if (rmarkdown::pandoc_available()) {
#'   ex       <- calmer_example_data()
#'   paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#'   features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                                paralog_table = paralog)
#'   clean    <- clean_features(features)
#'   fit <- fit_beinf(
#'     clean, n = 1000,
#'     formula       = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
#'                         n_isoforms + I(pmax(n_isoforms - 5, 0)) + has_paralog,
#'     sigma_formula = ~ 1, nu_formula = ~ 1, tau_formula = ~ 1
#'   )
#'   out <- tempfile(fileext = ".html")
#'   calmer_report(fit, features = clean, out_file = out, format = "html")
#'   file.exists(out)
#' }
#'
#' @export
calmer_report <- function(fit,
                          features    = NULL,
                          scored      = NULL,
                          out_file,
                          format      = c("pdf", "html"),
                          sample_name = NULL,
                          score_args  = list()) {
  format <- match.arg(format)

  if (!inherits(fit, "calmer_fit") && !inherits(fit, "calmer_selection")) {
    stop("`fit` must be a calmer_fit (or calmer_selection) object.")
  }
  if (!is.null(features) && !is.data.frame(features)) {
    stop("`features` must be a data frame (typically from clean_features()).")
  }
  if (!is.null(scored) && !is.data.frame(scored)) {
    stop("`scored` must be a data frame returned by score_reads(), or NULL.")
  }

  is_selection <- inherits(fit, "calmer_selection")
  rank_table   <- NULL
  picked_name  <- NULL
  picked_fit   <- NULL

  if (is_selection) {
    # select_beinf() already chose the model (min BIC among full / A / B);
    # use its scorable picked fit. The full comparison lives in the
    # selection object; the report renders the picked model's diagnostics.
    picked_fit  <- primary_fit(fit)
    picked_name <- fit$picked
    if (!is.null(features) && is.null(scored)) {
      scored <- do.call(score_reads,
                        c(list(features = features, fit = picked_fit),
                          score_args))
    }
  } else {
    # Single calmer_fit. No ranking; use as supplied.
    picked_fit <- fit
    if (!is.null(features) && is.null(scored)) {
      scored <- do.call(score_reads,
                        c(list(features = features, fit = picked_fit),
                          score_args))
    }
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
    params        = list(fit         = picked_fit,
                         suite       = NULL,
                         scored      = scored,
                         rank_table  = rank_table,
                         picked_name = picked_name,
                         sample_name = sample_name),
    envir         = new.env(parent = globalenv()),
    quiet         = TRUE
  )

  invisible(structure(
    out_file,
    picked_model = picked_name,
    scored       = scored
  ))
}
