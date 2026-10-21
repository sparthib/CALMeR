#' Quantify transcript counts from CALMeR-filtered reads
#'
#' Produces a per-transcript count table from a scored read table by
#' counting the reads CALMeR retains. By default this counts each
#' transcript's **primary alignments that were not flagged** by
#' [score_reads()] --- i.e. the "CALMeR-unflagged" quantification
#' strategy: keep every primary read CALMeR did not flag (this preserves
#' context-consistent low-MAPQ reads and drops only the context-anomalous
#' ones), then count reads per transcript.
#'
#' This is a light-weight, primary-alignment counting quantifier. It does
#' not perform expectation-maximization; each retained primary read
#' contributes a count of 1 to its assigned `transcript_id`. For
#' EM-based abundance estimation, pass the filtered BAM from
#' [filter_bam()] to a dedicated quantifier (e.g. \code{oarfish},
#' \code{salmon}).
#'
#' @param scored Data frame returned by [score_reads()]. Must contain
#'   `transcript_id` and `flag_low`; `aln_type` is used when present to
#'   restrict to primary alignments.
#' @param transcript_col Name of the transcript-ID column. Default
#'   `"transcript_id"`.
#' @param primary_only Logical; if `TRUE` (default) count only primary
#'   alignments (`aln_type == "P"`). Ignored with a warning if no
#'   `aln_type` column is present.
#' @param include_unmatched Logical; if `TRUE`, transcripts with a count
#'   of zero are still returned when `all_transcripts` is supplied.
#'   Default `FALSE`.
#' @param all_transcripts Optional character vector of transcript IDs to
#'   report counts for (e.g. the full annotation). Transcripts not seen
#'   in `scored` are reported with count 0 when `include_unmatched = TRUE`.
#'
#' @return A data frame with columns `transcript_id` and `count`, one row
#'   per transcript, sorted by descending count. The result carries an
#'   attribute `"quantify_summary"` with `n_reads_total`,
#'   `n_reads_flagged`, `n_reads_counted`, and `n_transcripts`.
#'
#' @examples
#' ex       <- calmer_example_data()
#' paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                              paralog_table = paralog)
#' clean    <- clean_features(features)
#' fit <- fit_beinf(
#'   clean, n = 1000,
#'   formula = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
#'                 n_isoforms + I(pmax(n_isoforms - 5, 0)) + has_paralog
#' )
#' scored <- score_reads(clean, fit, filter_method = "mad", threshold = -3)
#' counts <- quantify_counts(scored)
#' head(counts)
#'
#' @export
quantify_counts <- function(scored,
                            transcript_col    = "transcript_id",
                            primary_only      = TRUE,
                            include_unmatched = FALSE,
                            all_transcripts   = NULL) {
  stopifnot(
    is.data.frame(scored),
    transcript_col %in% names(scored),
    "flag_low"     %in% names(scored)
  )

  keep <- scored

  # restrict to primary alignments when requested and possible
  if (primary_only) {
    if ("aln_type" %in% names(keep)) {
      keep <- keep[keep$aln_type == "P", , drop = FALSE]
    } else {
      warning("quantify_counts: no 'aln_type' column; counting all ",
              "alignments (primary_only ignored).", call. = FALSE)
    }
  }

  n_reads_total <- nrow(keep)

  # drop CALMeR-flagged reads (flag_low == TRUE); NA flags are retained
  flagged        <- keep$flag_low %in% TRUE
  n_reads_flagged <- sum(flagged)
  keep            <- keep[!flagged, , drop = FALSE]
  n_reads_counted <- nrow(keep)

  tx <- keep[[transcript_col]]
  tx <- tx[!is.na(tx)]
  tab <- table(tx)
  counts <- data.frame(
    transcript_id = names(tab),
    count         = as.integer(tab),
    stringsAsFactors = FALSE
  )

  # optionally pad with zero-count transcripts from a supplied universe
  if (!is.null(all_transcripts) && include_unmatched) {
    missing <- setdiff(unique(all_transcripts), counts$transcript_id)
    if (length(missing)) {
      counts <- rbind(
        counts,
        data.frame(transcript_id = missing, count = 0L,
                   stringsAsFactors = FALSE)
      )
    }
  }

  counts <- counts[order(-counts$count, counts$transcript_id), , drop = FALSE]
  rownames(counts) <- NULL

  names(counts)[1] <- transcript_col

  attr(counts, "quantify_summary") <- list(
    n_reads_total   = n_reads_total,
    n_reads_flagged = n_reads_flagged,
    n_reads_counted = n_reads_counted,
    n_transcripts   = nrow(counts)
  )
  counts
}
