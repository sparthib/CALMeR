#' Build a small example dataset for CALMeR examples
#'
#' Creates a temporary directory containing a tiny BAM, GTF, FASTA, and PAF
#' suitable for exercising the CALMeR pipeline in examples and tests
#' without external data. Five transcripts on three genes, ~200 synthetic
#' reads with MAPQ values spanning the full 0-60 range, varying transcript
#' lengths, varying GC content, and varying alignment coverage.
#'
#' Cleanup is the caller's responsibility: `unlink(result$dir, recursive
#' = TRUE)`.
#'
#' @return A named list with character paths:
#' \describe{
#'   \item{`dir`}{The tempdir holding all files (delete to clean up).}
#'   \item{`bam`}{Tiny BAM (sorted, indexed).}
#'   \item{`gtf`}{Matching GTF.}
#'   \item{`fasta`}{Matching transcriptome FASTA.}
#'   \item{`paf`}{Tiny PAF with one cross-gene paralog hit.}
#' }
#'
#' @examples
#' ex <- calmer_example_data()
#' file.exists(ex$bam)
#' unlink(ex$dir, recursive = TRUE)
#'
#' @importFrom Rsamtools asBam
#' @keywords internal
#' @export
calmer_example_data <- function() {
  td <- tempfile("calmer_ex_")
  dir.create(td)

  tx_table <- data.frame(
    tx       = c("tx1", "tx2", "tx3", "tx4", "tx5"),
    gene     = c("gene1", "gene1", "gene1", "gene2", "gene3"),
    length   = c(200L, 400L, 800L, 300L, 600L),
    gc_pct   = c(30, 40, 50, 60, 70),
    stringsAsFactors = FALSE
  )

  # GTF: one exon per transcript (covers full length).
  gtf_path <- file.path(td, "ex.gtf")
  gtf_lines <- vapply(seq_len(nrow(tx_table)), function(i) {
    paste("chr1", "ex", "exon", 1, tx_table$length[i], ".", "+", ".",
          sprintf(
            'gene_id "%s"; transcript_id "%s"; gene_type "protein_coding";',
            tx_table$gene[i], tx_table$tx[i]
          ),
          sep = "\t")
  }, character(1))
  writeLines(gtf_lines, gtf_path)

  # FASTA: synthetic sequences with the target GC fraction.
  fa_path <- file.path(td, "ex.fa")
  fa_lines <- character(0)
  for (i in seq_len(nrow(tx_table))) {
    n_gc <- round(tx_table$length[i] * tx_table$gc_pct[i] / 100)
    n_at <- tx_table$length[i] - n_gc
    bases <- c(rep("G", n_gc %/% 2), rep("C", n_gc - n_gc %/% 2),
               rep("A", n_at %/% 2), rep("T", n_at - n_at %/% 2))
    fa_lines <- c(fa_lines, paste0(">", tx_table$tx[i]),
                  paste(bases, collapse = ""))
  }
  writeLines(fa_lines, fa_path)

  # SAM: 5 transcripts x 8 MAPQ levels x 5 alignment-length levels = 200 reads.
  sam_path <- file.path(td, "ex.sam")
  header <- c(
    "@HD\tVN:1.6\tSO:unsorted",
    sprintf("@SQ\tSN:%s\tLN:%d", tx_table$tx, tx_table$length)
  )
  mapq_levels   <- c(0L, 1L, 10L, 20L, 30L, 40L, 50L, 60L)
  aln_pct_levels <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  read_id <- 0L
  read_lines <- character(0)
  for (i in seq_len(nrow(tx_table))) {
    L <- tx_table$length[i]
    for (mapq in mapq_levels) {
      for (pct in aln_pct_levels) {
        read_id <- read_id + 1L
        aln_len <- max(20L, as.integer(L * pct))
        seq_str  <- paste(rep("A", aln_len), collapse = "")
        qual_str <- paste(rep("I", aln_len), collapse = "")
        read_lines <- c(read_lines, paste(
          paste0("read", read_id),  # qname
          "0",                       # flag (primary, mapped, +)
          tx_table$tx[i],            # rname
          "1",                       # pos
          mapq,                      # mapq
          paste0(aln_len, "M"),      # cigar
          "*", "0", "0",
          seq_str, qual_str,
          sep = "\t"
        ))
      }
    }
  }
  writeLines(c(header, read_lines), sam_path)

  bam_prefix <- sub("\\.sam$", "", sam_path)
  bam_path <- Rsamtools::asBam(sam_path, bam_prefix,
                               indexDestination = TRUE, overwrite = TRUE)

  # Tiny PAF: tx1 has one cross-gene paralog hit to tx4 (different gene),
  # 90% identity, full-length alignment.
  paf_path <- file.path(td, "ex.paf")
  writeLines(paste(
    "tx1", 400, 0, 400, "+", "tx4", 400, 0, 400,
    360, 400, 60, sep = "\t"
  ), paf_path)

  list(
    dir   = td,
    bam   = as.character(bam_path),
    gtf   = gtf_path,
    fasta = fa_path,
    paf   = paf_path
  )
}
