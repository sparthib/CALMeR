#' CALMeR: Context-Aware Long-read MAPQ Rescuer
#'
#' Filters BAM files from long-read transcriptome alignments by flagging reads
#' whose MAPQ is anomalously low given transcript-level structural features
#' (number of isoforms, transcript length, GC content, percent covered).
#'
#' @section Workflow:
#' \enumerate{
#'   \item [extract_features()] — read BAM + GTF + transcriptome FASTA, build
#'         per-read feature table.
#'   \item [clean_features()] — apply standard pre-fit filters
#'         (primary alignments, n_isoforms cap, GC range, length range).
#'   \item [fit_beinf()] — fit BEINF regression on a sampled subset.
#'   \item [score_reads()] — compute randomized quantile residuals on the full
#'         feature table and flag reads.
#'   \item [filter_bam()] — write filtered BAM and a TSV of flagged reads.
#'   \item [calmer_report()] — render a diagnostic report.
#' }
#'
#' @keywords internal
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point geom_line
#'   geom_density geom_abline labs theme_minimal facet_wrap
#' @importFrom patchwork wrap_plots
"_PACKAGE"
