#' Write a filtered BAM and a TSV of flagged reads
#'
#' Removes reads flagged by [score_reads()] from the input BAM and writes a
#' TSV recording the flagged read IDs and their feature values + RQR.
#'
#' @param bam Path to input BAM file.
#' @param scored Data frame returned by [score_reads()].
#' @param out_bam Path for the output filtered BAM.
#' @param out_tsv Path for the output flagged-reads TSV (gzipped if path
#'   ends in `.gz`).
#' @param read_id_col Name of the read-ID column in `scored` (default
#'   `"read_id"`).
#'
#' @return Invisibly, a list with `n_input`, `n_flagged`, `n_kept`,
#'   `out_bam`, `out_tsv`.
#'
#' @export
filter_bam <- function(bam,
                       scored,
                       out_bam,
                       out_tsv,
                       read_id_col = "read_id") {
  stop("Not yet implemented.")
}
