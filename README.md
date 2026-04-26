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

```r
features <- extract_features(bam, gtf, fasta)
fit      <- fit_beinf(features, n = 50000, seed = 42)
scored   <- score_reads(features, fit, filter_method = "mad", threshold = -3)
filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")
calmer_report(fit, scored, out_file = "report.pdf")
```

## Installation

```r
# devtools::install_github("sparthib/CALMeR")
```

(Not yet installable — package is in scaffolding.)
