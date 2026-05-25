#' Fit a suite of BEINF models for model comparison
#'
#' Convenience wrapper that fits up to 12 BEINF models on the same
#' training data and returns them as a single `calmer_suite` object.
#' The grid crosses three axes:
#' \describe{
#'   \item{Smoother for the scaled predictors}{`pb()` penalized
#'         B-spline vs plain linear.}
#'   \item{Interaction term}{none, `tx_length_sc:n_isoforms` (length by
#'         isoform count), or `tx_length_sc:pct_covered_sc` (length by
#'         coverage). All interactions are linear-by-linear cross terms
#'         on top of the marginal main effects; they test whether the
#'         interaction adds anything beyond the marginals.}
#'   \item{Paralog}{`linear_predictors` (default `has_paralog`)
#'         included vs omitted.}
#' }
#' Across every variant, `piecewise_predictors` (default `n_isoforms`
#' with a knot at 5) enters as the two-piece linear basis; it is not
#' swapped for `pb()` because seven distinct values can't support a
#' useful smooth (see `vignette` / [fit_beinf()] docs).
#'
#' Variant names follow the pattern
#' `<smoother>_<paralog>[_<interaction>]`:
#' \itemize{
#'   \item `pb_paralog`, `pb_paralog_x_niso`, `pb_paralog_x_pcov`
#'   \item `pb_no_paralog`, `pb_no_paralog_x_niso`, `pb_no_paralog_x_pcov`
#'   \item `linear_paralog`, `linear_paralog_x_niso`, `linear_paralog_x_pcov`
#'   \item `linear_no_paralog`, `linear_no_paralog_x_niso`, `linear_no_paralog_x_pcov`
#' }
#' where `_x_niso` adds `tx_length_sc:n_isoforms` and `_x_pcov` adds
#' `tx_length_sc:pct_covered_sc`. Variants whose required predictors
#' are absent from `features` (or from `scaled_predictors` /
#' `piecewise_predictors`) are skipped with a message.
#'
#' Primary fit (the one used by [score_reads()] and [filter_bam()]
#' when the suite is passed in):
#' \itemize{
#'   \item `"user"` if `user_formula` was supplied,
#'   \item otherwise `"pb_paralog"` if it could be fitted,
#'   \item otherwise `"pb_no_paralog"`.
#' }
#'
#' Fitting a dozen BEINF models is several times slower than fitting
#' one. Plan accordingly when calling from an array job.
#'
#' @param features Cleaned training data frame from [clean_features()].
#' @param scaled_predictors,piecewise_predictors,linear_predictors
#'   Forwarded to [fit_beinf()]; see its docstring for the exact
#'   semantics. Defaults match the canonical CALMeR pipeline
#'   (`tx_length`, `gc_content`, `pct_covered` scaled; `n_isoforms`
#'   piecewise at knot 5; `has_paralog` linear).
#' @param user_formula,user_sigma_formula,user_nu_formula,user_tau_formula
#'   Optional gamlss formulas for an additional user-defined fit. If
#'   `user_formula` is `NULL`, no user fit is added. If non-`NULL`, the
#'   sigma/nu/tau formulas default to the same RHS as `user_formula`
#'   (via [fit_beinf()]'s own defaulting).
#' @param inter Integer; number of equal-spaced intervals for each
#'   penalized B-spline term in the `pb_*` variants. Forwarded to the
#'   `inter` argument of [gamlss::pb()]. Default 5.
#' @param min_paralog_in_sample Integer floor. Before fitting the
#'   `*_paralog` variants, the function estimates the expected number
#'   of paralog reads (`has_paralog == TRUE`) in the random training
#'   subset. If that estimate falls below `min_paralog_in_sample`, a
#'   warning is emitted (the variants are still fitted; the warning
#'   just flags that the paralog coefficient may be noisy or
#'   unidentified). Default 100. Set to 0 to disable.
#' @param n,seed,trace Forwarded to [fit_beinf()].
#'
#' @return An object of class `"calmer_suite"`: a named list of
#'   `calmer_fit` objects, with an attribute `"primary"` giving the
#'   name of the fit that downstream functions should use for scoring.
#'
#' @examples
#' \dontrun{
#'   ex       <- calmer_example_data()
#'   paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#'   features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                                paralog_table = paralog)
#'   clean    <- clean_features(features)
#'   suite    <- fit_beinf_suite(clean, n = 5000)
#'   scored   <- score_reads(clean, suite)        # uses primary fit
#'   calmer_report(suite, scored, out_file = "report.html",
#'                 format = "html")
#' }
#'
#' @seealso [fit_beinf()], [calmer_report()], [score_reads()].
#' @export
fit_beinf_suite <- function(features,
                            scaled_predictors     = c("tx_length", "gc_content",
                                                      "pct_covered"),
                            piecewise_predictors  = list(n_isoforms = 5),
                            linear_predictors     = c("has_paralog"),
                            user_formula          = NULL,
                            user_sigma_formula    = NULL,
                            user_nu_formula       = NULL,
                            user_tau_formula      = NULL,
                            inter               = 5,
                            n                     = 50000,
                            seed                  = 42,
                            min_paralog_in_sample = 100,
                            trace                 = FALSE) {
  stopifnot(is.data.frame(features))

  paralog_cols_present <- length(linear_predictors) == 0L ||
                          all(linear_predictors %in% names(features))
  if (!paralog_cols_present) {
    missing <- setdiff(linear_predictors, names(features))
    message(sprintf(
      "fit_beinf_suite: paralog variants skipped (missing predictor%s: %s).",
      if (length(missing) > 1L) "s" else "",
      paste(missing, collapse = ", ")
    ))
  }

  # Surface the paralog-prevalence concern up front. When
  # `has_paralog` is rare in the training population, the random
  # training subset may contain very few TRUE rows, producing a noisy
  # or unidentified linear coefficient on `has_paralog`.
  if (paralog_cols_present && "has_paralog" %in% linear_predictors &&
      "has_paralog" %in% names(features)) {
    primary_mask     <- features$aln_type == "P"
    n_primary        <- sum(primary_mask, na.rm = TRUE)
    n_paralog        <- sum(primary_mask & features$has_paralog,
                            na.rm = TRUE)
    pop_frac         <- if (n_primary > 0) n_paralog / n_primary else 0
    expected_in_samp <- min(n, n_primary) * pop_frac
    if (expected_in_samp < min_paralog_in_sample) {
      warning(sprintf(
        paste0("fit_beinf_suite: expected ~%d paralog reads (has_paralog == TRUE) ",
               "in the %d-read training subset (population fraction = %.4f%%). ",
               "This is below min_paralog_in_sample = %d; the paralog ",
               "coefficient in pb_paralog / linear_paralog may be ",
               "unreliable. Consider raising `n` or treating has_paralog ",
               "as informative only when paralog reads are abundant."),
        round(expected_in_samp), n, 100 * pop_frac,
        min_paralog_in_sample
      ), call. = FALSE)
    }
  }

  build_formulas <- function(scaled, piecewise, linear, smoother, extra = NULL) {
    rhs <- .build_three_class_rhs(
      scaled_predictors    = scaled,
      piecewise_predictors = piecewise,
      linear_predictors    = linear,
      inter              = inter,
      smoother             = smoother
    )
    if (!is.null(extra) && nzchar(extra)) rhs <- paste(rhs, "+", extra)
    list(
      mu    = stats::as.formula(paste("y ~", rhs)),
      sigma = stats::as.formula(paste("~",   rhs)),
      nu    = stats::as.formula(paste("~",   rhs)),
      tau   = stats::as.formula(paste("~",   rhs))
    )
  }

  fit_with <- function(scaled, piecewise, linear, formulas) {
    fit_beinf(features,
              scaled_predictors    = scaled,
              piecewise_predictors = piecewise,
              linear_predictors    = linear,
              formula              = formulas$mu,
              sigma_formula        = formulas$sigma,
              nu_formula           = formulas$nu,
              tau_formula          = formulas$tau,
              n                    = n,
              seed                 = seed,
              trace                = trace)
  }

  # Axes of the variant grid.
  smoothers <- c("pb", "linear")
  interactions <- list(
    list(suffix = "",        extra = NULL),
    list(suffix = "_x_niso", extra = "tx_length_sc:n_isoforms"),
    list(suffix = "_x_pcov", extra = "tx_length_sc:pct_covered_sc")
  )
  paralog_axes <- list(
    list(suffix = "_paralog",    linear = linear_predictors,
         enabled = paralog_cols_present && length(linear_predictors) > 0L),
    list(suffix = "_no_paralog", linear = character(0),
         enabled = TRUE)
  )

  # Defensive: an interaction is only fitted when its constituent
  # predictors are actually in the model (so the formula resolves).
  has_tx_sc      <- "tx_length"   %in% scaled_predictors
  has_pcov_sc    <- "pct_covered" %in% scaled_predictors
  has_n_isoforms <- "n_isoforms"  %in% c(scaled_predictors,
                                         names(piecewise_predictors),
                                         linear_predictors)
  interaction_ok <- function(extra) {
    if (is.null(extra)) return(TRUE)
    if (grepl("tx_length_sc", extra, fixed = TRUE) && !has_tx_sc)
      return(FALSE)
    if (grepl("pct_covered_sc", extra, fixed = TRUE) && !has_pcov_sc)
      return(FALSE)
    if (grepl("n_isoforms", extra, fixed = TRUE) && !has_n_isoforms)
      return(FALSE)
    TRUE
  }

  fits <- list()
  for (sm in smoothers) {
    for (psp in paralog_axes) {
      if (!psp$enabled) next
      for (isp in interactions) {
        if (!interaction_ok(isp$extra)) {
          message(sprintf(
            "  skipping %s%s%s (required predictors not in the model).",
            sm, psp$suffix, isp$suffix
          ))
          next
        }
        name <- paste0(sm, psp$suffix, isp$suffix)
        message(sprintf("  fitting %s ...", name))
        f <- build_formulas(scaled_predictors, piecewise_predictors,
                            psp$linear, smoother = sm, extra = isp$extra)
        fits[[name]] <- fit_with(scaled_predictors, piecewise_predictors,
                                 psp$linear, f)
      }
    }
  }

  if (!is.null(user_formula)) {
    message("  fitting user ...")
    user_linear <- if (paralog_cols_present) linear_predictors else character(0)
    fits$user <- fit_beinf(features,
                           scaled_predictors    = scaled_predictors,
                           piecewise_predictors = piecewise_predictors,
                           linear_predictors    = user_linear,
                           formula              = user_formula,
                           sigma_formula        = user_sigma_formula,
                           nu_formula           = user_nu_formula,
                           tau_formula          = user_tau_formula,
                           n                    = n,
                           seed                 = seed,
                           trace                = trace)
  }

  primary <- if (!is.null(user_formula)) {
    "user"
  } else if ("pb_paralog" %in% names(fits)) {
    "pb_paralog"
  } else {
    "pb_no_paralog"
  }

  structure(fits, class = "calmer_suite", primary = primary)
}

#' @export
print.calmer_suite <- function(x, ...) {
  cat("CALMeR BEINF suite\n")
  cat("  fits:    ", paste(names(x), collapse = ", "), "\n", sep = "")
  cat("  primary: ", attr(x, "primary"), "\n", sep = "")
  invisible(x)
}

#' Return the primary fit from a calmer_suite (or pass through a single fit).
#'
#' Internal-ish helper used by [score_reads()], [calmer_report()] to
#' accept either a single fit or a suite transparently.
#'
#' @param x A `calmer_fit` or `calmer_suite` object.
#' @return A `calmer_fit`.
#' @keywords internal
#' @export
primary_fit <- function(x) {
  if (inherits(x, "calmer_fit"))  return(x)
  if (inherits(x, "calmer_suite")) {
    p <- attr(x, "primary")
    if (is.null(p) || !p %in% names(x)) {
      stop("calmer_suite has no valid primary fit.", call. = FALSE)
    }
    return(x[[p]])
  }
  stop("Expected a calmer_fit or calmer_suite object; got ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}
