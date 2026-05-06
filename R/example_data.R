#' Paths to CALMeR's bundled example dataset
#'
#' CALMeR ships with a small example bundle in `inst/extdata/`: a BAM
#' subsampled from a real PCR cDNA long-read sample (~2000 reads), a
#' matching GTF and transcriptome FASTA pruned to only the transcripts
#' those reads align to, and a `minimap2` self-alignment PAF used by
#' [make_paralog_table()]. Returns the resolved paths.
#'
#' These files live inside the installed package, so they should be
#' treated as read-only. Use [tempfile()] / [tempdir()] for any output
#' files an example writes (filtered BAM, report, etc.).
#'
#' @return A named list of character paths:
#' \describe{
#'   \item{`bam`}{Sorted, indexed transcriptome BAM (~2000 alignments).}
#'   \item{`gtf`}{Matching GTF (transcripts present in `bam` only).}
#'   \item{`fasta`}{Matching transcriptome FASTA.}
#'   \item{`paf`}{`minimap2` self-alignment of `fasta`, ready for
#'     [make_paralog_table()].}
#' }
#'
#' @examples
#' ex <- calmer_example_data()
#' file.exists(ex$bam) && file.exists(ex$gtf) &&
#'   file.exists(ex$fasta) && file.exists(ex$paf)
#'
#' @export
calmer_example_data <- function() {
  list(
    bam   = system.file("extdata", "ex.bam",     package = "CALMeR",
                        mustWork = TRUE),
    gtf   = system.file("extdata", "ex.gtf.gz",  package = "CALMeR",
                        mustWork = TRUE),
    fasta = system.file("extdata", "ex.fa.gz",   package = "CALMeR",
                        mustWork = TRUE),
    paf   = system.file("extdata", "ex.paf.gz", package = "CALMeR",
                        mustWork = TRUE)
  )
}
