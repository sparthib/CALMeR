#' Write a filtered BAM and a TSV of flagged reads
#'
#' Removes all alignments belonging to flagged reads from `bam` and writes
#' the result to `out_bam`. Also writes a TSV listing the flagged reads
#' along with their feature values and randomized quantile residual.
#'
#' All alignments (primary + secondary) of a flagged read are removed
#' from the output BAM. The flag itself is determined per primary
#' alignment in [score_reads()]; if a read's primary alignment was
#' flagged, the whole read is dropped.
#'
#' @param bam Path to the input BAM file.
#' @param scored Data frame returned by [score_reads()].
#' @param out_bam Path for the filtered BAM. An index is created
#'   alongside it.
#' @param out_tsv Path for the flagged-reads TSV. Gzipped automatically
#'   if the path ends in `.gz`.
#' @param read_id_col Name of the read-ID column in `scored`. Default
#'   `"read_id"`.
#'
#' @return Invisibly, a list with `n_alignments_in`, `n_flagged_reads`,
#'   `n_alignments_kept`, `n_alignments_removed`, `out_bam`, `out_tsv`.
#'
#' @examples
#' \donttest{
#' filter_bam(bam, scored,
#'            out_bam = "filtered.bam",
#'            out_tsv = "flagged.tsv.gz")
#' }
#'
#' @importFrom Rsamtools filterBam ScanBamParam indexBam
#' @importFrom S4Vectors FilterRules
#' @importFrom readr write_tsv
#' @export
filter_bam <- function(bam,
                       scored,
                       out_bam,
                       out_tsv,
                       read_id_col = "read_id") {
  stopifnot(
    file.exists(bam),
    is.data.frame(scored),
    "flag_low"      %in% names(scored),
    read_id_col     %in% names(scored)
  )

  flagged <- scored[scored$flag_low %in% TRUE, , drop = FALSE]
  flagged_ids <- unique(flagged[[read_id_col]])
  n_flagged_reads <- length(flagged_ids)

  if (!file.exists(paste0(bam, ".bai"))) {
    Rsamtools::indexBam(bam)
  }

  filter_func <- function(df) !(df$qname %in% flagged_ids)
  rules <- S4Vectors::FilterRules(list(not_flagged = filter_func))

  Rsamtools::filterBam(
    file        = bam,
    destination = out_bam,
    param       = Rsamtools::ScanBamParam(what = "qname"),
    filter      = rules,
    indexDestination = TRUE
  )

  drop_cols <- c(grep("_sc$", names(flagged), value = TRUE), "y")
  out_df    <- flagged[, setdiff(names(flagged), drop_cols), drop = FALSE]
  readr::write_tsv(out_df, out_tsv)

  n_in   <- nrow(scored)
  n_kept <- n_in - nrow(flagged)
  message(sprintf(
    "filter_bam: %d alignments in, %d flagged reads (%d alignments removed), %d alignments kept.",
    n_in, n_flagged_reads, nrow(flagged), n_kept
  ))

  invisible(list(
    n_alignments_in      = n_in,
    n_flagged_reads      = n_flagged_reads,
    n_alignments_kept    = n_kept,
    n_alignments_removed = nrow(flagged),
    out_bam              = out_bam,
    out_tsv              = out_tsv
  ))
}
