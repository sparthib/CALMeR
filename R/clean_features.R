#' Clean a feature table before fitting
#'
#' Applies the CALMeR canonical preprocessing pipeline to a feature
#' table from [extract_features()], in this order:
#' \enumerate{
#'   \item Drop secondary/supplementary if `primary_only = TRUE`.
#'   \item Drop rows with `NA` in `MAPQ` or any column listed in
#'         `required`.
#'   \item `n_isoforms`: lump at `max_n_isoforms` (default 7) by
#'         clamping (or drop above the cap if `mode = "drop"`).
#'   \item `tx_length`:  natural-log in place if `log_tx_length = TRUE`.
#'   \item `pct_covered`: winsorize to the supplied or data-driven
#'         quantile bounds if `winsorize_pct_covered = TRUE`.
#' }
#'
#' This function performs only the variable-specific transforms ("step
#' 1" of the spec). Z-scoring of the smoothed predictors is applied
#' downstream in [fit_beinf()] (and mirrored in [score_reads()]) using
#' parameters computed on the training data, so it is not done here.
#'
#' Defaults are opinionated to match the BEINF model used downstream:
#' `tx_length` is log-transformed (its dynamic range otherwise
#' dominates the z-score), `pct_covered` is winsorized (the 0-100 bound
#' concentrates mass near 100 with rare extreme low values), and
#' `n_isoforms` is lumped at 7 because the upper tail is sparse and the
#' variable enters the model as a piecewise linear term, not as a
#' smooth.
#'
#' When winsorizing data-drivenly, the bounds used are attached as
#' attribute `"pct_covered_winsor_bounds"` on the returned data frame.
#' Pass those same bounds back via `pct_covered_winsor_bounds` when
#' cleaning a test/scoring set to apply the training cutoffs rather
#' than recomputing on the new data.
#'
#' @param features Data frame from [extract_features()].
#' @param primary_only Logical; keep only primary alignments. Default `TRUE`.
#' @param max_n_isoforms Integer cap for `n_isoforms` (lump at this
#'   value). Default `7L`. Set `NULL` to disable.
#' @param log_tx_length Logical; if `TRUE` (default), overwrite
#'   `tx_length` in place with `log(tx_length)`. Requires all values
#'   to be strictly positive.
#' @param winsorize_pct_covered Logical; if `TRUE` (default), clip
#'   `pct_covered` to `pct_covered_winsor_bounds` (or to the
#'   data-driven quantiles at `pct_covered_winsor_probs` when bounds
#'   are `NULL`).
#' @param pct_covered_winsor_probs Length-2 numeric in (0, 1); the
#'   quantiles to use when winsorize bounds are not supplied.
#'   Default `c(0.01, 0.99)`.
#' @param pct_covered_winsor_bounds Optional length-2 numeric; if
#'   supplied, used as the lower/upper clip bounds for `pct_covered`
#'   without recomputing quantiles. Use this to apply training bounds
#'   to a scoring set.
#' @param mode Either `"clamp"` (default) or `"drop"`. Controls how
#'   `max_n_isoforms` is applied: `"clamp"` lumps values above the cap
#'   down to the cap; `"drop"` removes those rows.
#' @param required Character vector of additional columns that must
#'   be non-NA. Defaults to the four canonical predictors.
#' @param verbose Logical; print row counts (and clamp counts, if any)
#'   after each step.
#'
#' @return A data frame with predictor columns transformed in place
#'   per the steps above. When `winsorize_pct_covered = TRUE`, the
#'   returned data frame has attribute `"pct_covered_winsor_bounds"`
#'   set to the length-2 numeric used.
#'
#' @examples
#' ex       <- calmer_example_data()
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#' clean    <- clean_features(features)
#' nrow(features); nrow(clean)
#' attr(clean, "pct_covered_winsor_bounds")
#'
#' @importFrom stats quantile complete.cases
#' @export
clean_features <- function(features,
                           primary_only              = TRUE,
                           max_n_isoforms            = 7L,
                           log_tx_length             = TRUE,
                           winsorize_pct_covered     = TRUE,
                           pct_covered_winsor_probs  = c(0.01, 0.99),
                           pct_covered_winsor_bounds = NULL,
                           mode                      = c("clamp", "drop"),
                           required                  = c("n_isoforms", "tx_length",
                                                         "gc_content", "pct_covered"),
                           verbose                   = FALSE) {
  stopifnot(is.data.frame(features))
  mode <- match.arg(mode)

  log_step <- function(label, df) {
    if (verbose) message(sprintf("  %-32s n = %d", label, nrow(df)))
  }
  log_clamp <- function(label, n) {
    if (verbose) message(sprintf("  %-32s clamped = %d", label, n))
  }
  log_note <- function(label, info) {
    if (verbose) message(sprintf("  %-32s %s", label, info))
  }

  log_step("input", features)

  if (primary_only) {
    if (!"aln_type" %in% names(features)) {
      stop("primary_only = TRUE but `aln_type` column is missing.")
    }
    features <- features[features$aln_type == "P", , drop = FALSE]
    log_step("primary only", features)
  }

  cols_for_na <- intersect(c("MAPQ", required), names(features))
  features <- features[stats::complete.cases(features[, cols_for_na, drop = FALSE]),
                       , drop = FALSE]
  log_step("non-NA predictors", features)

  if (!is.null(max_n_isoforms) && "n_isoforms" %in% names(features)) {
    if (mode == "clamp") {
      n_hit <- sum(features$n_isoforms > max_n_isoforms, na.rm = TRUE)
      features$n_isoforms <- pmin(features$n_isoforms, max_n_isoforms)
      log_clamp(sprintf("n_isoforms <= %d", max_n_isoforms), n_hit)
    } else {
      features <- features[features$n_isoforms <= max_n_isoforms, , drop = FALSE]
      log_step(sprintf("n_isoforms <= %d", max_n_isoforms), features)
    }
  }

  if (isTRUE(log_tx_length) && "tx_length" %in% names(features)) {
    if (any(features$tx_length <= 0, na.rm = TRUE)) {
      stop("log_tx_length = TRUE but tx_length has values <= 0. ",
           "Drop or impute non-positive lengths upstream, or pass ",
           "log_tx_length = FALSE.", call. = FALSE)
    }
    features$tx_length <- log(features$tx_length)
    log_note("tx_length", "natural-log transformed in place")
  }

  bounds_used <- NULL
  if (isTRUE(winsorize_pct_covered) && "pct_covered" %in% names(features)) {
    if (is.null(pct_covered_winsor_bounds)) {
      stopifnot(length(pct_covered_winsor_probs) == 2L,
                all(pct_covered_winsor_probs > 0),
                all(pct_covered_winsor_probs < 1),
                pct_covered_winsor_probs[1] < pct_covered_winsor_probs[2])
      bounds_used <- unname(stats::quantile(
        features$pct_covered,
        probs = pct_covered_winsor_probs,
        na.rm = TRUE,
        names = FALSE
      ))
    } else {
      stopifnot(length(pct_covered_winsor_bounds) == 2L,
                pct_covered_winsor_bounds[1] <= pct_covered_winsor_bounds[2])
      bounds_used <- pct_covered_winsor_bounds
    }
    n_hit <- sum(features$pct_covered < bounds_used[1] |
                 features$pct_covered > bounds_used[2], na.rm = TRUE)
    features$pct_covered <- pmax(pmin(features$pct_covered, bounds_used[2]),
                                 bounds_used[1])
    log_clamp(sprintf("pct_covered winsorize [%.3g, %.3g]",
                      bounds_used[1], bounds_used[2]), n_hit)
  }

  if (!is.null(bounds_used)) {
    attr(features, "pct_covered_winsor_bounds") <- bounds_used
  }
  features
}
