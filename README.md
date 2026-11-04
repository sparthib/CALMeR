# CALMeR

**C**ontext-**A**ware **L**ong-read **MAPQ** **R**escuer.

CALMeR filters BAM files from long-read transcriptome alignments by flagging
reads whose MAPQ is anomalously low given transcript-level structural features
(number of isoforms, transcript length, GC content, percent covered, and
cross-gene paralog presence). It fits a Beta-Inflated (BEINF) regression
model on a sampled subset of reads, computes randomized quantile residuals
on the full BAM, and removes reads with anomalously low residuals.

## Status

Core pipeline implemented: feature extraction, feature cleaning, BEINF
model fitting, backward GAMLSS model selection, read scoring/flagging,
BAM filtering, per-transcript quantification, and a diagnostic report.

## Workflow

Single-model workflow (one BEINF fit, one diagnostic report):

```r
paralog  <- make_paralog_table(fasta, gtf, paf_path = "self.paf")
features <- extract_features(bam, gtf, fasta, paralog_table = paralog)
clean    <- clean_features(features)
fit      <- fit_beinf(clean, n = 50000, seed = 42, inter = 5)
scored   <- score_reads(clean, fit, filter_method = "mad",
                        threshold = -2, dev_threshold = 2)
filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")

# per-transcript counts from CALMeR-retained (unflagged) reads
counts   <- quantify_counts(scored)

calmer_report(fit, features = clean,
              out_file = "report.html", format = "html")
```

Model-selection workflow (backward GAMLSS selection, then report and
scoring against the picked model):

```r
paralog  <- make_paralog_table(fasta, gtf, paf_path = "self.paf")
features <- extract_features(bam, gtf, fasta, paralog_table = paralog)
clean    <- clean_features(features)

# Backward GAMLSS selection from the FULL BEINF model (pb smooths for the
# scaled predictors, the two-piece isoform basis, the paralog indicator,
# and the transcript-feature interactions). At each step the term whose
# removal most improves the GAIC is dropped, subject to marginality;
# k = log(n) => BIC by default (pass k = 2 for AIC). Both GAMLSS
# strategies are run and the picked (minimum-BIC) model is kept:
#   Strategy A (stepGAICAll.A): per-parameter sweeps (mu, sigma, nu, tau)
#   Strategy B (stepGAICAll.B): one interleaved search over all parameters
# Optionally add held-out scoring with evaluate = c("gaic", "holdout")
# (validation deviance via gamlss::gamlssVGD) or K-fold CV (cv_folds).
sel <- select_beinf(clean, strategies = c("A", "B"), evaluate = "gaic",
                    n = 50000, seed = 42, inter = 5, ncpus = 4)

# calmer_report() accepts the selection object directly: it re-scores
# against the picked model and renders the diagnostic report. The picked
# model name and the scored data are returned as attributes on the path.
out    <- calmer_report(sel, features = clean,
                        out_file = "report.html", format = "html")
picked <- attr(out, "picked_model")
scored <- attr(out, "scored")

filter_bam(bam, scored, out_bam = "filtered.bam", out_tsv = "flagged.tsv.gz")
counts <- quantify_counts(scored)
```

## Installation

```r
# devtools::install_github("sparthib/CALMeR")
```

(Not yet installable — package is in scaffolding.)
