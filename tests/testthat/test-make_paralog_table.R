test_that("make_paralog_table parses a precomputed PAF and returns one binary row per transcript", {
  td <- tempfile("calmer_paralog_test_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  gtf_path <- file.path(td, "test.gtf")
  paf_path <- file.path(td, "test.paf")
  fa_path  <- file.path(td, "test.fa")

  # Minimal GTF: three transcripts on two genes (gene1 has tx1 and tx3,
  # gene2 has tx2). All length 100 so cross-gene alignments span 100%
  # of the query in either direction and clear the query_cov filter.
  # rtracklayer requires the standard 9-column layout.
  gtf_lines <- c(
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx1";', sep = "\t"),
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene2"; transcript_id "tx2";', sep = "\t"),
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx3";', sep = "\t")
  )
  writeLines(gtf_lines, gtf_path)

  # FASTA only checked for existence by make_paralog_table (no parsing here).
  writeLines(c(">tx1", "ACGT", ">tx2", "ACGT", ">tx3", "ACGT"), fa_path)

  # Realistic minimap2 self-alignment output (no -X): both directions of
  # the tx1/tx2 pair plus the diagonal self-hits. tx1 <-> tx2 are
  # cross-gene paralogs; the qname != tname / cross-gene filters drop
  # the self-hits so they cannot masquerade as paralog calls. tx3 has
  # no cross-gene PAF rows and should backfill to FALSE.
  paf_lines <- c(
    paste("tx1", 100, 0, 100, "+", "tx2", 100, 0, 100,  90, 100, 60,
          sep = "\t"),
    paste("tx2", 100, 0, 100, "+", "tx1", 100, 0, 100,  90, 100, 60,
          sep = "\t"),
    paste("tx1", 100, 0, 100, "+", "tx1", 100, 0, 100, 100, 100, 60,
          sep = "\t"),
    paste("tx2", 100, 0, 100, "+", "tx2", 100, 0, 100, 100, 100, 60,
          sep = "\t")
  )
  writeLines(paf_lines, paf_path)

  result <- make_paralog_table(
    fasta         = fa_path,
    gtf           = gtf_path,
    paf_path      = paf_path,
    aln_min_len   = 50,
    query_cov_min = 70
  )

  # One row per annotated transcript with the three canonical columns.
  expect_equal(nrow(result), 3L)
  expect_setequal(result$transcript_id, c("tx1", "tx2", "tx3"))
  expect_named(result, c("transcript_id", "gene_id", "has_paralog"))
  expect_type(result$has_paralog, "logical")

  # tx1 <-> tx2 are cross-gene; both directions are emitted so both
  # flag TRUE. tx3 has no qualifying hit and backfills to FALSE.
  expect_true(result$has_paralog[result$transcript_id == "tx1"])
  expect_true(result$has_paralog[result$transcript_id == "tx2"])
  expect_false(result$has_paralog[result$transcript_id == "tx3"])
})

test_that("make_paralog_table excludes within-gene hits", {
  td <- tempfile("calmer_paralog_within_gene_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  gtf_path <- file.path(td, "test.gtf")
  paf_path <- file.path(td, "test.paf")
  fa_path  <- file.path(td, "test.fa")

  # tx1 and tx3 share gene1 -> a tx1 -> tx3 hit must not count as a paralog.
  gtf_lines <- c(
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx1";', sep = "\t"),
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx3";', sep = "\t")
  )
  writeLines(gtf_lines, gtf_path)
  writeLines(c(">tx1", "ACGT", ">tx3", "ACGT"), fa_path)

  paf_lines <- paste(
    "tx1", 100, 0, 100, "+", "tx3", 100, 0, 100, 99, 100, 60,
    sep = "\t"
  )
  writeLines(paf_lines, paf_path)

  result <- make_paralog_table(
    fasta    = fa_path,
    gtf      = gtf_path,
    paf_path = paf_path
  )

  # No cross-gene hits for either transcript.
  expect_true(all(!result$has_paralog))
})
