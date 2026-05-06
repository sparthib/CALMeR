#' Clean a feature table before fitting
#'
#' Applies a set of standard filters commonly used before fitting a CALMeR
#' BEINF model. Each filter can be relaxed or disabled by passing `NULL`.
#' Defaults reflect choices made on long-read RNA-seq data; if your study
#' has different characteristics, override them.
#'
#' Filters (in order):
#' \enumerate{
#'   \item Drop secondary/supplementary if `primary_only = TRUE`.
#'   \item Drop rows with `NA` in `MAPQ` or any column listed in
#'         `required`.
#'   \item Cap `n_isoforms` at `max_n_isoforms` (drops rows above the cap).
#'   \item Restrict `gc_content` to the half-open range `gc_range`.
#'   \item Restrict `tx_length` to the half-open range `tx_length_range`.
#' }
#'
#' @param features Data frame from [extract_features()].
#' @param primary_only Logical; keep only primary alignments. Default `TRUE`.
#' @param max_n_isoforms Integer; drop rows with `n_isoforms` above this.
#'   Default `10`. Set `NULL` to disable.
#' @param gc_range Numeric length-2; inclusive `c(low, high)` percent G+C.
#'   Default `c(30, 70)`. Set `NULL` to disable.
#' @param tx_length_range Numeric length-2; inclusive `c(min, max)`
#'   transcript length in nucleotides. Default `c(100, 10000)`. Set
#'   `NULL` to disable.
#' @param required Character vector of additional columns that must be
#'   non-NA. Defaults to the four canonical predictors.
#' @param verbose Logical; print row counts after each filter.
#'
#' @return A filtered data frame.
#'
#' @examples
#' ex       <- calmer_example_data()
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#' clean    <- clean_features(features)
#' nrow(features); nrow(clean)
#' unlink(ex$dir, recursive = TRUE)
#'
#' @export
clean_features <- function(features,
                           primary_only    = TRUE,
                           max_n_isoforms  = 10,
                           gc_range        = c(30, 70),
                           tx_length_range = c(100, 10000),
                           required        = c("n_isoforms", "tx_length",
                                               "gc_content", "pct_covered"),
                           verbose         = FALSE) {
  stopifnot(is.data.frame(features))

  log_step <- function(label, df) {
    if (verbose) message(sprintf("  %-28s n = %d", label, nrow(df)))
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
    features <- features[features$n_isoforms <= max_n_isoforms, , drop = FALSE]
    log_step(sprintf("n_isoforms <= %d", max_n_isoforms), features)
  }

  if (!is.null(gc_range) && "gc_content" %in% names(features)) {
    stopifnot(length(gc_range) == 2L, gc_range[1] <= gc_range[2])
    features <- features[features$gc_content >= gc_range[1] &
                         features$gc_content <= gc_range[2], , drop = FALSE]
    log_step(sprintf("gc_content in [%g, %g]", gc_range[1], gc_range[2]),
             features)
  }

  if (!is.null(tx_length_range) && "tx_length" %in% names(features)) {
    stopifnot(length(tx_length_range) == 2L,
              tx_length_range[1] <= tx_length_range[2])
    features <- features[features$tx_length >= tx_length_range[1] &
                         features$tx_length <= tx_length_range[2],
                         , drop = FALSE]
    log_step(sprintf("tx_length in [%g, %g]",
                     tx_length_range[1], tx_length_range[2]),
             features)
  }

  features
}
