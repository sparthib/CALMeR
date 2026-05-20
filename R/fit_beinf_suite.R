#' Fit a suite of BEINF models for model comparison
#'
#' Convenience wrapper that fits up to five BEINF models on the same
#' training data and returns them as a single `calmer_suite` object.
#' Useful for comparing flexible (penalized B-spline) vs. linear
#' parametrisations, with and without paralog context.
#'
#' By default the suite contains five named fits:
#' \describe{
#'   \item{`pb_paralog`}{`pb(., inter = n_knots)` on the four canonical
#'         features **plus a linear term in `max_pct_id`**. The paralog
#'         predictor is intentionally not wrapped in `pb()` because
#'         `max_pct_id` is highly zero-inflated (most transcripts have
#'         no paralog), which destabilises spline fits. Fitted only if
#'         `max_pct_id` is present on `features`.}
#'   \item{`pb_paralog_interact`}{Same as `pb_paralog` plus a
#'         linear-by-linear interaction `tx_length_sc:pct_covered_sc`.
#'         Tests whether the relationship between coverage and read
#'         quality depends on transcript length on top of the marginal
#'         non-linear effects already captured by `pb()`. Same paralog
#'         requirement as `pb_paralog`.}
#'   \item{`pb_no_paralog`}{`pb(., inter = n_knots)` on the four
#'         canonical features.}
#'   \item{`linear_paralog`}{Linear terms on all five predictors
#'         (including `max_pct_id`). Same paralog requirement.}
#'   \item{`linear_no_paralog`}{Linear terms on the four canonical
#'         features.}
#' }
#' If `user_formula` is non-`NULL`, a fifth fit named `"user"` is added
#' using exactly the formulas supplied.
#'
#' Primary fit (the one used by [score_reads()] and [filter_bam()]
#' when the suite is passed in):
#' \itemize{
#'   \item `"user"` if `user_formula` was supplied,
#'   \item otherwise `"pb_paralog"` if it could be fitted,
#'   \item otherwise `"pb_no_paralog"`.
#' }
#'
#' Fitting five to six BEINF models is several times slower than
#' fitting one. Plan accordingly when calling from an array job.
#'
#' @param features Cleaned training data frame from [clean_features()].
#' @param predictors_with_paralog Character vector of predictors for
#'   the `*_paralog` variants. Default
#'   `c("n_isoforms", "tx_length", "gc_content", "pct_covered", "max_pct_id")`.
#' @param predictors_no_paralog Character vector of predictors for
#'   the `*_no_paralog` variants. Default
#'   `c("n_isoforms", "tx_length", "gc_content", "pct_covered")`.
#' @param user_formula,user_sigma_formula,user_nu_formula,user_tau_formula
#'   Optional gamlss formulas for an additional user-defined fit. If
#'   `user_formula` is `NULL`, no user fit is added. If non-`NULL`, the
#'   sigma/nu/tau formulas default to the same RHS as `user_formula`
#'   (via [fit_beinf()]'s own defaulting).
#' @param n_knots Integer passed to `pb(., inter = n_knots)` for the
#'   penalized-spline variants. Default 5.
#' @param min_paralog_in_sample Integer floor. Before fitting the
#'   paralog variants, the function estimates the expected number of
#'   paralog reads (`max_pct_id > 0`) in the random training subset.
#'   If that estimate falls below `min_paralog_in_sample`, a warning is
#'   emitted (the variants are still fitted; the warning just flags
#'   that the paralog coefficient may be noisy or unidentified).
#'   Default 100. Set to 0 to disable.
#' @param n,seed,trace Forwarded to [fit_beinf()].
#'
#' @return An object of class `"calmer_suite"`: a named list of
#'   `calmer_fit` objects, with an attribute `"primary"` giving the
#'   name of the fit that downstream functions should use for scoring.
#'
#' @examples
#' \dontrun{
#'   ex       <- calmer_example_data()
#'   features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                                paralog_table = my_paralog_tsv)
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
                            predictors_with_paralog = c("n_isoforms", "tx_length",
                                                         "gc_content", "pct_covered",
                                                         "max_pct_id"),
                            predictors_no_paralog   = c("n_isoforms", "tx_length",
                                                         "gc_content", "pct_covered"),
                            user_formula            = NULL,
                            user_sigma_formula      = NULL,
                            user_nu_formula         = NULL,
                            user_tau_formula        = NULL,
                            n_knots                 = 5,
                            n                       = 50000,
                            seed                    = 42,
                            min_paralog_in_sample   = 100,
                            trace                   = FALSE) {
  stopifnot(is.data.frame(features))

  paralog_ok <- all(predictors_with_paralog %in% names(features))
  if (!paralog_ok) {
    missing <- setdiff(predictors_with_paralog, names(features))
    message(sprintf(
      "fit_beinf_suite: paralog variants skipped (missing predictor%s: %s).",
      if (length(missing) > 1) "s" else "",
      paste(missing, collapse = ", ")
    ))
  }

  # When fitting paralog variants, fit_beinf() will sample n primaries
  # uniformly. If paralog reads are rare in the population, the sample
  # may contain very few of them, producing a noisy or unidentified
  # paralog coefficient. Surface this up front rather than letting the
  # user discover it from a near-zero coefficient in the report.
  if (paralog_ok && "max_pct_id" %in% names(features)) {
    primary_mask     <- features$aln_type == "P"
    n_primary        <- sum(primary_mask, na.rm = TRUE)
    n_paralog        <- sum(primary_mask & features$max_pct_id > 0,
                            na.rm = TRUE)
    pop_frac         <- if (n_primary > 0) n_paralog / n_primary else 0
    expected_in_samp <- min(n, n_primary) * pop_frac
    if (expected_in_samp < min_paralog_in_sample) {
      warning(sprintf(
        paste0("fit_beinf_suite: expected ~%d paralog reads (max_pct_id > 0) in ",
               "the %d-read training subset (population fraction = %.4f%%). ",
               "This is below min_paralog_in_sample = %d; the paralog ",
               "coefficient in pb_paralog / linear_paralog may be ",
               "unreliable. Consider raising `n` or treating max_pct_id ",
               "as informative only when paralog reads are abundant."),
        round(expected_in_samp), n, 100 * pop_frac,
        min_paralog_in_sample
      ), call. = FALSE)
    }
  }

  build_linear_formula <- function(preds, response = "y") {
    rhs <- paste(paste0(preds, "_sc"), collapse = " + ")
    list(
      mu    = stats::as.formula(paste(response, "~", rhs)),
      sigma = stats::as.formula(paste("~", rhs)),
      nu    = stats::as.formula(paste("~", rhs)),
      tau   = stats::as.formula(paste("~", rhs))
    )
  }

  fits <- list()

  if (paralog_ok) {
    # Hybrid: pb() on the canonical four predictors + linear max_pct_id.
    # `max_pct_id` is intentionally not in pb() because the column is
    # heavily zero-inflated (most transcripts have no paralog), and a
    # penalized spline on a near-degenerate predictor fits poorly.
    message("  fitting pb_paralog (hybrid: pb on canonical 4 + linear max_pct_id) ...")
    np_sc       <- paste0(predictors_no_paralog, "_sc")
    paralog_sc  <- setdiff(paste0(predictors_with_paralog, "_sc"), np_sc)
    pb_terms    <- paste0("pb(", np_sc, ", inter = ", n_knots, ")",
                          collapse = " + ")
    hybrid_rhs  <- if (length(paralog_sc) > 0L) {
      paste(pb_terms, "+", paste(paralog_sc, collapse = " + "))
    } else {
      pb_terms
    }
    hf <- list(
      mu    = stats::as.formula(paste("y ~", hybrid_rhs)),
      sigma = stats::as.formula(paste("~",   hybrid_rhs)),
      nu    = stats::as.formula(paste("~",   hybrid_rhs)),
      tau   = stats::as.formula(paste("~",   hybrid_rhs))
    )
    fits$pb_paralog <- fit_beinf(features,
                                 predictors    = predictors_with_paralog,
                                 formula       = hf$mu,
                                 sigma_formula = hf$sigma,
                                 nu_formula    = hf$nu,
                                 tau_formula   = hf$tau,
                                 n             = n,
                                 seed          = seed,
                                 trace         = trace)

    # Same hybrid backbone as pb_paralog, plus a linear-by-linear
    # interaction between tx_length and pct_covered. Tests whether the
    # coverage-vs-quality relationship shifts with transcript length on
    # top of what the marginal pb() splines already explain.
    if (all(c("tx_length_sc", "pct_covered_sc") %in%
            c(np_sc, paralog_sc))) {
      message("  fitting pb_paralog_interact (pb_paralog + tx_length:pct_covered) ...")
      interact_rhs <- paste(hybrid_rhs, "+ tx_length_sc:pct_covered_sc")
      hif <- list(
        mu    = stats::as.formula(paste("y ~", interact_rhs)),
        sigma = stats::as.formula(paste("~",   interact_rhs)),
        nu    = stats::as.formula(paste("~",   interact_rhs)),
        tau   = stats::as.formula(paste("~",   interact_rhs))
      )
      fits$pb_paralog_interact <- fit_beinf(features,
                                            predictors    = predictors_with_paralog,
                                            formula       = hif$mu,
                                            sigma_formula = hif$sigma,
                                            nu_formula    = hif$nu,
                                            tau_formula   = hif$tau,
                                            n             = n,
                                            seed          = seed,
                                            trace         = trace)
    } else {
      message("  skipping pb_paralog_interact (tx_length/pct_covered not in predictor set).")
    }
  }

  message("  fitting pb_no_paralog ...")
  fits$pb_no_paralog <- fit_beinf(features,
                                  predictors = predictors_no_paralog,
                                  n_knots    = n_knots,
                                  n          = n,
                                  seed       = seed,
                                  trace      = trace)

  if (paralog_ok) {
    message("  fitting linear_paralog ...")
    lf <- build_linear_formula(predictors_with_paralog)
    fits$linear_paralog <- fit_beinf(features,
                                     predictors    = predictors_with_paralog,
                                     formula       = lf$mu,
                                     sigma_formula = lf$sigma,
                                     nu_formula    = lf$nu,
                                     tau_formula   = lf$tau,
                                     n             = n,
                                     seed          = seed,
                                     trace         = trace)
  }

  message("  fitting linear_no_paralog ...")
  lf <- build_linear_formula(predictors_no_paralog)
  fits$linear_no_paralog <- fit_beinf(features,
                                      predictors    = predictors_no_paralog,
                                      formula       = lf$mu,
                                      sigma_formula = lf$sigma,
                                      nu_formula    = lf$nu,
                                      tau_formula   = lf$tau,
                                      n             = n,
                                      seed          = seed,
                                      trace         = trace)

  if (!is.null(user_formula)) {
    message("  fitting user ...")
    # fit_beinf() requires `predictors` to map terms onto column names
    # for scaling. Best-effort guess: union of the canonical four with
    # paralog if present in features.
    user_preds <- if (paralog_ok) predictors_with_paralog else predictors_no_paralog
    fits$user <- fit_beinf(features,
                           predictors    = user_preds,
                           formula       = user_formula,
                           sigma_formula = user_sigma_formula,
                           nu_formula    = user_nu_formula,
                           tau_formula   = user_tau_formula,
                           n             = n,
                           seed          = seed,
                           trace         = trace)
  }

  primary <- if (!is.null(user_formula)) {
    "user"
  } else if (paralog_ok) {
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
