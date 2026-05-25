#' Compute a per-transcript paralog status table via minimap2 self-alignment
#'
#' For each annotated transcript, identifies whether it has at least
#' one cross-gene transcript with significant sequence similarity, by
#' running a minimap2 self-alignment of the transcriptome FASTA.
#' Returns one row per transcript suitable for joining into
#' [extract_features()] as a binary context predictor.
#'
#' Cross-gene definition: the target transcript must belong to a
#' different `gene_id` than the query. A transcript is flagged with
#' `has_paralog = TRUE` if at least one qualifying cross-gene
#' alignment survives the `aln_min_len` / `query_cov_min` filters and
#' the `asm20` preset's effective identity floor (~80%); otherwise
#' `FALSE`. No continuous identity score is carried through: with the
#' canonical pipeline, the alignment-stage filters and `asm20` preset
#' together encode the paralog cutoff, so the relevant signal is just
#' "qualifying hit, yes or no."
#'
#' If `paf_path` already exists and is non-empty, alignment is skipped and
#' the existing PAF is parsed. If it does not exist (or is empty), minimap2
#' is invoked via [base::system2()] with the supplied threads, preset, and
#' `-N` setting.
#'
#' @param fasta Path to transcriptome FASTA (used as both query and target
#'   in the self-alignment).
#' @param gtf Path to GTF annotation file. Used to map `transcript_id` to
#'   `gene_id` for cross-gene filtering.
#' @param paf_path Path to a PAF file. If `NULL` (default), a temporary
#'   file is created and removed at exit. Provide an explicit path to
#'   cache the alignment for reuse. A precomputed PAF must contain both
#'   mapping directions of each pair (i.e. not be generated with
#'   minimap2's `-X` flag); otherwise per-transcript scores are asymmetric.
#' @param aln_min_len Minimum alignment length in bp (default 50).
#' @param query_cov_min Minimum query coverage as a percent (default 70).
#' @param minimap2_path Path to the minimap2 executable. Default
#'   `"minimap2"` (resolved via `PATH`). Tilde expansion is applied.
#' @param threads Threads passed to minimap2 (default 8).
#' @param preset minimap2 alignment preset (default `"asm20"`, which
#'   targets pairs with <20% divergence).
#' @param max_secondary Value passed to minimap2's `-N` (default 100).
#'
#' @return A data frame with one row per annotated transcript. Columns:
#' \describe{
#'   \item{transcript_id, gene_id}{Identifiers from the GTF.}
#'   \item{has_paralog}{Logical; `TRUE` if the transcript has at
#'     least one qualifying cross-gene alignment, `FALSE` otherwise.}
#' }
#'
#' @section minimap2 dependency:
#' This function calls the [minimap2](https://github.com/lh3/minimap2)
#' executable when computing a fresh PAF; minimap2 is declared in
#' `SystemRequirements`. If you supply a precomputed PAF (or a precomputed
#' paralog TSV consumed directly by [extract_features()]), the external
#' tool is not required at runtime.
#'
#' @examples
#' ex      <- calmer_example_data()
#' paralog <- make_paralog_table(
#'   fasta    = ex$fasta,
#'   gtf      = ex$gtf,
#'   paf_path = ex$paf
#' )
#' head(paralog)
#'
#' @importFrom rtracklayer import
#' @importFrom data.table fread
#' @importFrom dplyr filter select distinct left_join mutate
#' @importFrom rlang .data
#' @export
make_paralog_table <- function(fasta,
                               gtf,
                               paf_path      = NULL,
                               aln_min_len   = 50L,
                               query_cov_min = 70,
                               minimap2_path = "minimap2",
                               threads       = 8,
                               preset        = "asm20",
                               max_secondary = 100L) {
  stopifnot(
    file.exists(fasta),
    file.exists(gtf),
    is.numeric(aln_min_len),   length(aln_min_len)   == 1L, aln_min_len   > 0,
    is.numeric(query_cov_min), length(query_cov_min) == 1L,
    query_cov_min >= 0, query_cov_min <= 100,
    is.numeric(threads), length(threads) == 1L, threads > 0
  )

  minimap2_path <- path.expand(minimap2_path)

  if (is.null(paf_path)) {
    paf_path <- tempfile(fileext = ".paf")
    on.exit(unlink(paf_path), add = TRUE)
  }

  paf_ok <- file.exists(paf_path) && file.size(paf_path) > 0
  if (!paf_ok) {
    if (!.minimap2_available(minimap2_path)) {
      stop("PAF not found (or empty) at: ", paf_path,
           "\nand minimap2 is not available at: ", minimap2_path,
           "\nInstall minimap2 (https://github.com/lh3/minimap2) ",
           "or supply a precomputed PAF via paf_path.",
           call. = FALSE)
    }
    .run_minimap2_self_aln(
      fasta         = fasta,
      paf_out       = paf_path,
      threads       = threads,
      preset        = preset,
      max_secondary = max_secondary,
      minimap2_path = minimap2_path
    )
  }

  # transcript_id -> gene_id from GTF.
  gtf_obj <- rtracklayer::import(gtf)
  tx2gene <- as.data.frame(gtf_obj) |>
    dplyr::filter(.data$type == "transcript",
                  !is.na(.data$transcript_id)) |>
    dplyr::select("transcript_id", "gene_id") |>
    dplyr::distinct()

  # Parse PAF (12 fixed cols; PAF can have a variable number of optional
  # tag cols per row, so use fill = Inf so fread doesn't lock onto the
  # column count from the first few rows).
  paf <- data.table::fread(
    paf_path, select = seq_len(12L), header = FALSE, sep = "\t",
    col.names = c("qname", "qlen", "qstart", "qend", "strand",
                  "tname", "tlen", "tstart", "tend",
                  "nmatch", "alen", "mapq"),
    fill = Inf
  )
  paf$pct_id    <- paf$nmatch / paf$alen * 100
  paf$query_cov <- (paf$qend - paf$qstart) / paf$qlen * 100

  # Cross-gene + length/coverage filter.
  paf_cg <- paf[paf$qname != paf$tname &
                  paf$alen >= aln_min_len &
                  paf$query_cov >= query_cov_min, ]
  paf_cg$q_gene <- tx2gene$gene_id[match(paf_cg$qname, tx2gene$transcript_id)]
  paf_cg$t_gene <- tx2gene$gene_id[match(paf_cg$tname, tx2gene$transcript_id)]
  paf_cg <- paf_cg[is.na(paf_cg$q_gene) | is.na(paf_cg$t_gene) |
                     paf_cg$q_gene != paf_cg$t_gene, ]

  # Set of query transcripts with at least one qualifying cross-gene
  # alignment. Every other annotated transcript backfills to FALSE.
  hit_ids <- unique(paf_cg$qname)

  out <- data.frame(
    transcript_id    = unique(tx2gene$transcript_id),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(tx2gene, by = "transcript_id") |>
    dplyr::mutate(has_paralog = .data$transcript_id %in% hit_ids) |>
    dplyr::select("transcript_id", "gene_id", "has_paralog")

  out
}

# Internal: probe minimap2 availability without erroring.
.minimap2_available <- function(path) {
  out <- tryCatch(
    suppressWarnings(
      system2(path, "--version", stdout = TRUE, stderr = TRUE)
    ),
    error = function(e) NULL
  )
  length(out) > 0
}

# Internal: run minimap2 self-alignment, write PAF.
.run_minimap2_self_aln <- function(fasta, paf_out, threads, preset,
                                   max_secondary, minimap2_path) {
  dir.create(dirname(paf_out), recursive = TRUE, showWarnings = FALSE)
  log_out <- sub("\\.paf$", ".log", paf_out)
  # NB: -D (skip self mappings) rather than -X. -X is equivalent to
  # -D + --dual=no, and --dual=no drops one direction of every pair, so
  # per-query aggregation in make_paralog_table() would only see one
  # member of each paralog pair and backfill the other to
  # has_paralog = FALSE. -D keeps both directions while still
  # suppressing the trivial diagonal at the source. The qname != tname
  # filter at the parse step is kept as a safety net for precomputed
  # PAFs from other sources.
  status <- system2(
    command = minimap2_path,
    args    = c("-D",                          # skip self-hits (keep both dirs)
                "-x", preset,                  # preset (default: asm20)
                "--secondary=yes",
                "-N", as.character(max_secondary),
                "-c",                          # CIGAR for accurate identity
                "-t", as.character(threads),
                shQuote(fasta),                # query
                shQuote(fasta)),               # target -> self-align
    stdout  = paf_out,
    stderr  = log_out
  )
  if (!identical(as.integer(status), 0L)) {
    stop("minimap2 failed (exit ", status, "). See log: ", log_out,
         call. = FALSE)
  }
  invisible(paf_out)
}
