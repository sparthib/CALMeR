# CALMeR

**C**ontext-**A**ware **L**ong-read **MAPQ** **R**escuer.

CALMeR filters BAM files from long-read transcriptome alignments by flagging
reads whose MAPQ is anomalously low given transcript-level structural features
(number of isoforms, transcript length, GC content, percent covered, and
cross-gene paralog presence). It fits a Beta-Inflated (BEINF) regression
model on a sampled subset of reads, computes randomized quantile residuals
on the full BAM, and removes reads with anomalously low residuals.

## Status

Early development. Function stubs are in place; core implementations are
being ported from analysis scripts.

## Workflow

Single-model workflow (one BEINF fit, one diagnostic report):

```r
paralog  <- make_paralog_table(fasta, gtf, paf_path = "self.paf")
features <- extract_features(bam, gtf, fasta, paralog_table = paralog)
clean    <- clean_features(features)
fit      <- fit_beinf(clean, n = 50000, seed = 42, inter = 5)
scored   <- score_reads(clean, fit, filter_method = "mad", threshold = -3)
filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")
calmer_report(fit, features = clean,
              out_file = "report.html", format = "html")
```

Model-comparison workflow (fits a suite of up to 12 default models, ranks
them, picks the best, and drives the report and scoring against it):

```r
paralog  <- make_paralog_table(fasta, gtf, paf_path = "self.paf")
features <- extract_features(bam, gtf, fasta, paralog_table = paralog)
clean    <- clean_features(features)

# Fits up to 12 variants on the grid:
#   {pb, linear} smoother for the scaled predictors
#   x { none, + tx_length_sc:n_isoforms, + tx_length_sc:pct_covered_sc }
#   x { with has_paralog, without }
# Variant names follow `<smoother>_<paralog>[_<interaction>]`, e.g.
# `pb_paralog_x_pcov` or `linear_no_paralog`. See ?fit_beinf_suite.
suite <- fit_beinf_suite(clean, n = 50000, seed = 42, inter = 5)

# calmer_report() ranks the suite by AIC + BIC + Nagelkerke pseudo-R^2,
# picks the model with the smallest average rank, re-scores against the
# picked model, and renders the diagnostic report. The picked model
# name and the scored data are returned as attributes on the path.
out    <- calmer_report(suite, features = clean,
                        out_file = "report.html", format = "html")
picked <- attr(out, "picked_model")
scored <- attr(out, "scored")

filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")
```

## Installation

```r
# devtools::install_github("sparthib/CALMeR")
```

(Not yet installable — package is in scaffolding.)
