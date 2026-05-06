test_that("make_paralog_table parses a precomputed PAF and returns one row per transcript", {
  td <- tempfile("calmer_paralog_test_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  gtf_path <- file.path(td, "test.gtf")
  paf_path <- file.path(td, "test.paf")
  fa_path  <- file.path(td, "test.fa")

  # Minimal GTF: three transcripts on two genes (gene1 has tx1 and tx3,
  # gene2 has tx2). rtracklayer requires the standard 9-column layout.
  gtf_lines <- c(
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx1";', sep = "\t"),
    paste("chr1", "test", "transcript", 1, 200, ".", "+", ".",
          'gene_id "gene2"; transcript_id "tx2";', sep = "\t"),
    paste("chr1", "test", "transcript", 1, 100, ".", "+", ".",
          'gene_id "gene1"; transcript_id "tx3";', sep = "\t")
  )
  writeLines(gtf_lines, gtf_path)

  # FASTA only checked for existence by make_paralog_table (no parsing here).
  writeLines(c(">tx1", "ACGT", ">tx2", "ACGT", ">tx3", "ACGT"), fa_path)

  # Two PAF rows: tx1 -> tx2 and tx2 -> tx1, both cross-gene, 90% identity,
  # full query coverage. Tx3 has no PAF rows -> should backfill to 0/FALSE.
  paf_lines <- c(
    paste("tx1", 100, 0, 100, "+", "tx2", 200, 0, 100, 90, 100, 60,
          sep = "\t"),
    paste("tx2", 200, 0, 100, "+", "tx1", 100, 0, 100, 90, 100, 60,
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

  # One row per annotated transcript.
  expect_equal(nrow(result), 3L)
  expect_setequal(result$transcript_id, c("tx1", "tx2", "tx3"))

  # tx1 has a 90% cross-gene hit (tx2): paralog_50/80/90 TRUE, 95 FALSE.
  tx1 <- result[result$transcript_id == "tx1", ]
  expect_equal(tx1$max_pct_id, 90)
  expect_equal(tx1$n_cross_hits, 1L)
  expect_true(tx1$has_paralog_50)
  expect_true(tx1$has_paralog_80)
  expect_true(tx1$has_paralog_90)
  expect_false(tx1$has_paralog_95)

  # tx3 has no PAF rows -> backfilled to zero / FALSE.
  tx3 <- result[result$transcript_id == "tx3", ]
  expect_equal(tx3$max_pct_id, 0)
  expect_equal(tx3$n_cross_hits, 0L)
  expect_false(tx3$has_paralog_50)
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
  expect_true(all(result$max_pct_id == 0))
  expect_true(all(result$n_cross_hits == 0L))
  expect_true(all(!result$has_paralog_50))
})
