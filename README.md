# CALMeR

**C**ontext-**A**ware **L**ong-read **MAPQ** **R**escuer.

CALMeR filters BAM files from long-read transcriptome alignments by flagging
reads whose MAPQ is anomalously low given transcript-level structural features
(number of isoforms, transcript length, GC content, percent covered). It fits
a Beta-Inflated regression model on a sampled subset of reads, computes
randomized quantile residuals on the full BAM, and removes reads with low
residuals.

## Status

Early development. Function stubs are in place; core implementations are
being ported from analysis scripts.

## Workflow

Single-model workflow (one BEINF fit, one diagnostic report):

```r
features <- extract_features(bam, gtf, fasta)
clean    <- clean_features(features)
fit      <- fit_beinf(clean, n = 50000, seed = 42, n_knots = 5)
scored   <- score_reads(clean, fit, filter_method = "mad", threshold = -3)
filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")
calmer_report(fit, scored, out_file = "report.html", format = "html")
```

Model-comparison workflow (fits a suite of default models, picks one
to drive scoring):

```r
features <- extract_features(bam, gtf, fasta, paralog_table = my_paralog_tsv)
clean    <- clean_features(features)

# Fits four defaults:
#   pb_paralog        : pb() on canonical 4 + linear max_pct_id (hybrid)
#   pb_no_paralog     : pb() on canonical 4
#   linear_paralog    : linear on all 5
#   linear_no_paralog : linear on canonical 4
# A fifth "user" fit is added if `user_formula` is supplied.
suite <- fit_beinf_suite(clean, n = 50000, seed = 42, n_knots = 5)

# Scoring + filtering use the suite's PRIMARY fit:
#   - "user"          if user_formula was supplied,
#   - "pb_paralog"    if max_pct_id is on `features`,
#   - "pb_no_paralog" otherwise.
scored <- score_reads(clean, suite, filter_method = "mad", threshold = -3)
filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")

# Report adds a model-comparison section (AIC / BIC / EDF / pseudo-R^2
# table + side-by-side QQ panels). Per-model diagnostic sections run
# on the primary fit.
calmer_report(suite, scored, out_file = "report.html", format = "html")
```

## Installation

```r
# devtools::install_github("sparthib/CALMeR")
```

(Not yet installable — package is in scaffolding.)
