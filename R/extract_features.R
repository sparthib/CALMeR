#' Extract per-read features from a BAM file
#'
#' Reads a transcriptome-aligned BAM file and joins per-read alignment
#' features (MAPQ, alignment length, alignment type) with per-transcript
#' annotation features derived from a GTF and transcriptome FASTA
#' (transcript length, isoform count, exon count, GC content).
#'
#' @param bam Path to BAM file.
#' @param gtf Path to GTF annotation file.
#' @param fasta Path to transcriptome FASTA.
#' @param gene_types Character vector of gene_type values to retain
#'   (default: `c("lncRNA", "protein_coding")`).
#' @param spike_in_seqnames Character vector of seqnames treated as
#'   spike-in chromosomes (always retained regardless of gene_type).
#' @param include_secondary Logical; include secondary alignments (default
#'   `FALSE`, primary only).
#'
#' @return A data.frame with one row per alignment containing
#'   `transcript_id`, `MAPQ`, `alignment_length`, `aln_type`,
#'   `tx_length`, `n_isoforms`, `n_exons`, `gc_content`, `pct_covered`.
#'
#' @export
extract_features <- function(bam,
                             gtf,
                             fasta,
                             gene_types        = c("lncRNA", "protein_coding"),
                             spike_in_seqnames = c("SIRVomeERCCome", "chrIS"),
                             include_secondary = FALSE) {
  stop("Not yet implemented.")
}
