#' Clean a feature table before fitting
#'
#' Prepares a feature table for [fit_beinf()] and [score_reads()]. By
#' default the range-style arguments **clamp** values to their bounds
#' (so no rows are dropped due to a single out-of-range predictor); set
#' `mode = "drop"` to fall back to the legacy row-dropping behaviour.
#'
#' Filters and transforms (in order):
#' \enumerate{
#'   \item Drop secondary/supplementary if `primary_only = TRUE`.
#'   \item Drop rows with `NA` in `MAPQ` or any column listed in
#'         `required` (NA cannot be clamped).
#'   \item `n_isoforms`: clamp to `[1, max_n_isoforms]` (or drop above
#'         the cap if `mode = "drop"`).
#'   \item `gc_content`: clamp to `gc_range` (or drop outside).
#'   \item `tx_length`:  clamp to `tx_length_range` (or drop outside).
#' }
#'
#' Using `mode = "clamp"` (the default) is recommended when the same
#' table will be used for both [fit_beinf()] (training) and
#' [score_reads()] (scoring). Out-of-bound reads get scored against the
#' fitted model at the boundary value rather than via spline
#' extrapolation past the training domain.
#'
#' @param features Data frame from [extract_features()].
#' @param primary_only Logical; keep only primary alignments. Default `TRUE`.
#' @param max_n_isoforms Integer upper bound for `n_isoforms`. Default
#'   `10`. Set `NULL` to disable.
#' @param gc_range Numeric length-2; inclusive `c(low, high)` percent G+C.
#'   Default `c(30, 70)`. Set `NULL` to disable.
#' @param tx_length_range Numeric length-2; inclusive `c(min, max)`
#'   transcript length in nucleotides. Default `c(100, 10000)`. Set
#'   `NULL` to disable.
#' @param mode Either `"clamp"` (default) or `"drop"`. Controls how
#'   `max_n_isoforms`, `gc_range`, and `tx_length_range` are applied.
#' @param required Character vector of additional columns that must be
#'   non-NA. Defaults to the four canonical predictors.
#' @param verbose Logical; print row counts (and clamp counts, if any)
#'   after each step.
#'
#' @return A data frame with the same rows as `features` minus any
#'   dropped by `primary_only` / NA / `mode = "drop"` filtering, and
#'   with predictor columns clamped if `mode = "clamp"`.
#'
#' @examples
#' ex       <- calmer_example_data()
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#' clean    <- clean_features(features)         # clamp (default)
#' nrow(features); nrow(clean)
#' clean_drop <- clean_features(features, mode = "drop")  # legacy
#'
#' @export
clean_features <- function(features,
                           primary_only    = TRUE,
                           max_n_isoforms  = 10,
                           gc_range        = c(30, 70),
                           tx_length_range = c(100, 10000),
                           mode            = c("clamp", "drop"),
                           required        = c("n_isoforms", "tx_length",
                                               "gc_content", "pct_covered"),
                           verbose         = FALSE) {
  stopifnot(is.data.frame(features))
  mode <- match.arg(mode)

  log_step <- function(label, df) {
    if (verbose) message(sprintf("  %-32s n = %d", label, nrow(df)))
  }
  log_clamp <- function(label, n) {
    if (verbose) message(sprintf("  %-32s clamped = %d", label, n))
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

  if (!is.null(gc_range) && "gc_content" %in% names(features)) {
    stopifnot(length(gc_range) == 2L, gc_range[1] <= gc_range[2])
    if (mode == "clamp") {
      n_hit <- sum(features$gc_content < gc_range[1] |
                   features$gc_content > gc_range[2], na.rm = TRUE)
      features$gc_content <- pmax(pmin(features$gc_content, gc_range[2]),
                                  gc_range[1])
      log_clamp(sprintf("gc_content in [%g, %g]", gc_range[1], gc_range[2]),
                n_hit)
    } else {
      features <- features[features$gc_content >= gc_range[1] &
                           features$gc_content <= gc_range[2], , drop = FALSE]
      log_step(sprintf("gc_content in [%g, %g]", gc_range[1], gc_range[2]),
               features)
    }
  }

  if (!is.null(tx_length_range) && "tx_length" %in% names(features)) {
    stopifnot(length(tx_length_range) == 2L,
              tx_length_range[1] <= tx_length_range[2])
    if (mode == "clamp") {
      n_hit <- sum(features$tx_length < tx_length_range[1] |
                   features$tx_length > tx_length_range[2], na.rm = TRUE)
      features$tx_length <- pmax(pmin(features$tx_length, tx_length_range[2]),
                                 tx_length_range[1])
      log_clamp(sprintf("tx_length in [%g, %g]",
                        tx_length_range[1], tx_length_range[2]),
                n_hit)
    } else {
      features <- features[features$tx_length >= tx_length_range[1] &
                           features$tx_length <= tx_length_range[2],
                           , drop = FALSE]
      log_step(sprintf("tx_length in [%g, %g]",
                       tx_length_range[1], tx_length_range[2]),
               features)
    }
  }

  features
}
