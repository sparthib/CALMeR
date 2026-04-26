#' Extract per-read features from a BAM file
#'
#' Reads a transcriptome-aligned BAM file and joins per-read alignment
#' features with per-transcript annotation features derived from a GTF and
#' transcriptome FASTA.
#'
#' Per-read features extracted from the BAM:
#'   `read_id`, `transcript_id`, `MAPQ`, `alignment_length`, `aln_type`
#'   (`"P"` for primary, `"S"` for secondary).
#'
#' Per-transcript features derived from the GTF and FASTA:
#'   `tx_length` (sum of exon widths), `n_exons`, `n_isoforms`
#'   (transcripts per gene), `gc_content` (percent G+C in the transcript
#'   sequence).
#'
#' The two are joined on `transcript_id`. If a read does not match by
#'   full transcript ID (including version suffix), a fallback match is
#'   attempted with the version stripped.
#'
#' @param bam Path to BAM file. An index is created if missing.
#' @param gtf Path to GTF annotation file.
#' @param fasta Path to transcriptome FASTA used as the alignment reference.
#' @param gene_types Optional character vector of `gene_type` values to
#'   retain (e.g. `c("lncRNA", "protein_coding")`). If `NULL` (default),
#'   all gene types are kept.
#' @param extra_seqnames Optional character vector of seqnames that should
#'   be retained regardless of `gene_types` (useful for spike-in
#'   chromosomes such as `"SIRVomeERCCome"` or `"chrIS"`).
#' @param include_secondary Logical; if `FALSE` (default), secondary
#'   alignments are excluded.
#'
#' @return A data.frame with one row per alignment containing
#'   `read_id`, `transcript_id`, `MAPQ`, `alignment_length`, `aln_type`,
#'   `tx_length`, `n_exons`, `n_isoforms`, `gc_content`, `pct_covered`.
#'
#' @examples
#' \dontrun{
#' features <- extract_features(
#'   bam   = "sample.sorted.bam",
#'   gtf   = "gencode.v44.gtf",
#'   fasta = "transcripts.fa"
#' )
#' }
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag indexBam
#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace
#' @importFrom rtracklayer import
#' @importFrom Biostrings readDNAStringSet letterFrequency
#' @importFrom dplyr group_by summarise filter mutate left_join select
#'   distinct n_distinct n
#' @importFrom rlang .data
#' @export
extract_features <- function(bam,
                             gtf,
                             fasta,
                             gene_types        = NULL,
                             extra_seqnames    = NULL,
                             include_secondary = FALSE) {
  stopifnot(
    file.exists(bam),
    file.exists(gtf),
    file.exists(fasta),
    is.logical(include_secondary), length(include_secondary) == 1L
  )

  tx_lookup <- .build_tx_lookup(gtf, fasta, gene_types, extra_seqnames)
  reads_df  <- .read_bam_features(bam, include_secondary)
  reads_df  <- .join_tx_features(reads_df, tx_lookup)

  reads_df$pct_covered <- ifelse(
    !is.na(reads_df$tx_length) & reads_df$tx_length > 0,
    pmin(reads_df$alignment_length / reads_df$tx_length * 100, 100),
    NA_real_
  )

  reads_df[, c("read_id", "transcript_id", "MAPQ", "alignment_length",
               "aln_type", "tx_length", "n_exons", "n_isoforms",
               "gc_content", "pct_covered")]
}

# Internal: build per-transcript lookup from GTF + FASTA.
.build_tx_lookup <- function(gtf, fasta, gene_types, extra_seqnames) {
  gtf_obj <- rtracklayer::import(gtf)
  exons   <- gtf_obj[gtf_obj$type == "exon"]

  tx_info <- as.data.frame(exons) |>
    dplyr::group_by(.data$transcript_id, .data$gene_id,
                    .data$gene_type, .data$seqnames) |>
    dplyr::summarise(
      tx_length = sum(.data$width),
      n_exons   = dplyr::n(),
      .groups   = "drop"
    )

  if (!is.null(gene_types) || !is.null(extra_seqnames)) {
    keep <- rep(FALSE, nrow(tx_info))
    if (!is.null(gene_types)) {
      keep <- keep | (tx_info$gene_type %in% gene_types)
    }
    if (!is.null(extra_seqnames)) {
      keep <- keep | (as.character(tx_info$seqnames) %in% extra_seqnames)
    }
    tx_info <- tx_info[keep, , drop = FALSE]
  }

  isoform_counts <- tx_info |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(
      n_isoforms = dplyr::n_distinct(.data$transcript_id),
      .groups    = "drop"
    )

  tx_info <- dplyr::left_join(tx_info, isoform_counts, by = "gene_id")

  txome <- Biostrings::readDNAStringSet(fasta)
  gc_freq <- Biostrings::letterFrequency(txome, letters = c("G", "C"),
                                         as.prob = TRUE)
  transcript_gc <- data.frame(
    transcript_id = names(txome),
    gc_content    = rowSums(gc_freq) * 100,
    stringsAsFactors = FALSE
  )

  dplyr::left_join(tx_info, transcript_gc, by = "transcript_id")
}

# Internal: scan BAM and build per-read frame.
.read_bam_features <- function(bam, include_secondary) {
  if (!file.exists(paste0(bam, ".bai"))) {
    Rsamtools::indexBam(bam)
  }

  param <- Rsamtools::ScanBamParam(
    what = c("qname", "rname", "mapq", "cigar", "flag"),
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery          = FALSE,
      isSupplementaryAlignment = FALSE,
      isSecondaryAlignment     = if (include_secondary) NA else FALSE
    )
  )
  bam_data <- Rsamtools::scanBam(bam, param = param)[[1]]

  aln_length <- GenomicAlignments::cigarWidthAlongQuerySpace(
    bam_data$cigar, after.soft.clipping = TRUE
  )

  data.frame(
    read_id          = bam_data$qname,
    transcript_id    = as.character(bam_data$rname),
    MAPQ             = bam_data$mapq,
    alignment_length = aln_length,
    aln_type         = ifelse(bitwAnd(bam_data$flag, 0x100) > 0, "S", "P"),
    stringsAsFactors = FALSE
  )
}

# Internal: join reads with tx features, with no-version fallback.
.join_tx_features <- function(reads_df, tx_lookup) {
  out <- dplyr::left_join(
    reads_df,
    dplyr::select(tx_lookup, "transcript_id", "tx_length", "n_exons",
                  "n_isoforms", "gc_content"),
    by = "transcript_id"
  )

  unmatched <- is.na(out$tx_length)
  if (any(unmatched)) {
    tx_lookup_nv <- tx_lookup |>
      dplyr::mutate(
        transcript_id_nv = sub("\\.\\d+$", "", .data$transcript_id)
      ) |>
      dplyr::distinct(.data$transcript_id_nv, .keep_all = TRUE) |>
      dplyr::select("transcript_id_nv", "tx_length", "n_exons",
                    "n_isoforms", "gc_content")

    out$transcript_id_nv <- sub("\\.\\d+$", "", out$transcript_id)
    fallback <- dplyr::left_join(
      data.frame(transcript_id_nv = out$transcript_id_nv[unmatched],
                 stringsAsFactors = FALSE),
      tx_lookup_nv,
      by = "transcript_id_nv"
    )

    out$tx_length[unmatched]  <- fallback$tx_length
    out$n_exons[unmatched]    <- fallback$n_exons
    out$n_isoforms[unmatched] <- fallback$n_isoforms
    out$gc_content[unmatched] <- fallback$gc_content
    out$transcript_id_nv      <- NULL
  }

  out
}
