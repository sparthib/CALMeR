#' Score and flag reads using a fitted BEINF model
#'
#' For every read in `features`, predicts the four BEINF parameters
#' (mu, sigma, nu, tau) under the fitted model, then computes a
#' randomized quantile residual (RQR) and flags reads with anomalously
#' low residuals.
#'
#' Preprocessing must be consistent between training and scoring data.
#' The expected pattern is: clean the training data with
#' [clean_features()], fit with [fit_beinf()] or [fit_beinf_suite()],
#' then clean the scoring data with the **same** transforms — including
#' the training winsorize bounds stored on the fit
#' (`fit$winsor_bounds`):
#' \preformatted{
#'   clean_train <- clean_features(features)
#'   fit         <- fit_beinf(clean_train)
#'   clean_test  <- clean_features(
#'     features,
#'     pct_covered_winsor_bounds = fit$winsor_bounds
#'   )
#'   scored      <- score_reads(clean_test, fit)
#' }
#' Inside `score_reads`, only the z-score step is re-applied (using
#' `fit$scale_params` on `scaled_predictors`); the log / winsorize /
#' lump-at-7 transforms are expected to already be in place on the
#' input.
#'
#' For continuous response values (MAPQ in 1..59), the RQR is the
#' standard normal quantile of the model CDF at that value. For point
#' masses (MAPQ = 0 or MAPQ = 60), the residual is randomized over the
#' corresponding CDF interval (Dunn & Smyth 1996), so the overall
#' residual distribution is approximately N(0, 1) under a correctly
#' specified model.
#'
#' Flag rule is asymmetric and applied per region:
#' \itemize{
#'   \item **MAPQ in (0, 60)** (continuous part) — flag if `rqr` falls
#'         in the lower tail (`rqr < flag_threshold`, where the
#'         threshold is set by `filter_method` + `threshold` on the
#'         continuous-only RQRs).
#'   \item **MAPQ = 0** (lower point mass) — flag if the model is
#'         surprised by the zero, i.e. if the binomial deviance for
#'         "is the read at the zero mass?" is large:
#'         `dev_zero = sqrt(-2 log p_hat_0) > dev_threshold`.
#'         Equivalently, flag if
#'         `p_hat_0 < exp(-dev_threshold^2 / 2)`.
#'         This catches well-covered, single-isoform, no-paralog reads
#'         that landed at MAPQ = 0 despite the model predicting a tiny
#'         probability of zero-mass outcome.
#'   \item **MAPQ = 60** (upper point mass) — never flagged. MAPQ = 60
#'         is the maximum possible value, so it cannot be "unusually
#'         low" given any context.
#' }
#'
#' Reads with `NA` in any required predictor are kept in the output but
#' assigned `rqr = NA` and `flag_low = FALSE` (i.e. retained by default).
#' Reads with at least one predictor outside the training-data range
#' are marked `out_of_range = TRUE` and not scored, since the model's
#' cubic-spline extrapolation is unreliable beyond the training domain.
#'
#' @param features Data frame from [clean_features()] (or
#'   [extract_features()] if no cleaning is needed). Must contain all
#'   predictors used by the fit, transformed the same way as the
#'   training data.
#' @param fit Either a single `calmer_fit` returned by [fit_beinf()],
#'   or a `calmer_suite` returned by [fit_beinf_suite()]. If a suite is
#'   supplied, the suite's primary fit is used.
#' @param filter_method One of `"mad"` (default) or `"zscore"`. Controls
#'   the RQR-tail rule for continuous-MAPQ reads.
#'   For `"mad"`: flag if `rqr < median(rqr) + threshold * mad(rqr)`,
#'   computed on the continuous-MAPQ subset only.
#'   For `"zscore"`: flag if `(rqr - mean(rqr)) / sd(rqr) < threshold`.
#' @param threshold Numeric scalar for the RQR rule; default `-3`.
#' @param dev_threshold Numeric scalar for the MAPQ = 0 deviance rule;
#'   default `3` (mirrors `abs(threshold)`). Flag MAPQ = 0 reads when
#'   `sqrt(-2 log p_hat_0) > dev_threshold`, i.e. when the model gave
#'   them less than `exp(-dev_threshold^2 / 2)` probability of being
#'   at the zero mass. With `dev_threshold = 3` that probability cutoff
#'   is ~ 0.011.
#' @param seed Random seed for the uniform draws on point masses.
#' @param primary_only Logical; if `TRUE` (default) score only primary
#'   alignments. Secondary alignments inherit `flag_low = FALSE`.
#'
#' @return A data frame with the input columns plus `mu_p`, `sigma_p`,
#'   `nu_p`, `tau_p`, `p0_hat`, `p1_hat`, `rqr`, `dev_zero`,
#'   `flag_low`, `flag_source` (one of `"rqr_low"`, `"p0_dev_high"`,
#'   `"none"`), `out_of_range`, and an attribute `"flag_threshold"`
#'   recording the RQR cutoff used.
#'
#' @examples
#' ex       <- calmer_example_data()
#' paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                              paralog_table = paralog)
#' clean    <- clean_features(features)
#' fit <- fit_beinf(
#'   clean, n = 1000,
#'   formula       = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
#'                       n_isoforms + I(pmax(n_isoforms - 5, 0)) + has_paralog,
#'   sigma_formula = ~ 1, nu_formula = ~ 1, tau_formula = ~ 1
#' )
#' scored <- score_reads(clean, fit,
#'                       filter_method = "mad", threshold = -3)
#' table(scored$flag_low)
#'
#' @importFrom stats predict median mad sd complete.cases runif qnorm
#' @export
score_reads <- function(features,
                        fit,
                        filter_method = c("mad", "zscore"),
                        threshold     = -3,
                        dev_threshold = 3,
                        seed          = 42,
                        primary_only  = TRUE) {
  filter_method <- match.arg(filter_method)
  stopifnot(is.numeric(dev_threshold), length(dev_threshold) == 1L,
            dev_threshold > 0)

  fit <- primary_fit(fit)
  if (!inherits(fit, "calmer_fit")) {
    stop("`fit` must be a calmer_fit (or calmer_suite) object.")
  }
  stopifnot(
    is.data.frame(features),
    "MAPQ" %in% names(features),
    all(fit$predictors %in% names(features)),
    is.numeric(threshold), length(threshold) == 1L
  )

  scaled_predictors    <- fit$scaled_predictors
  piecewise_predictors <- fit$piecewise_predictors
  linear_predictors    <- fit$linear_predictors
  all_predictors       <- fit$predictors
  scale_params         <- fit$scale_params

  # Re-apply z-scoring to the scaled predictors using training params.
  # Piecewise and linear predictors enter raw.
  for (p in scaled_predictors) {
    features[[paste0(p, "_sc")]] <- apply_scale(features[[p]],
                                                scale_params[[p]])
  }
  features$y <- features$MAPQ / 60

  features$mu_p        <- NA_real_
  features$sigma_p     <- NA_real_
  features$nu_p        <- NA_real_
  features$tau_p       <- NA_real_
  features$p0_hat      <- NA_real_
  features$p1_hat      <- NA_real_
  features$rqr         <- NA_real_
  features$dev_zero    <- NA_real_
  features$flag_low    <- FALSE
  features$flag_source <- "none"

  scorable <- stats::complete.cases(
    features[, c("MAPQ", all_predictors), drop = FALSE]
  )
  if (primary_only && "aln_type" %in% names(features)) {
    scorable <- scorable & (features$aln_type == "P")
  }

  # Out-of-range check against the training-data range, for every
  # predictor (scaled, piecewise, and linear). Logical predictors
  # naturally pass the [0, 1] check; piecewise/lumped predictors are
  # only out-of-range if cleaning was inconsistent with training.
  features$out_of_range <- FALSE
  for (p in all_predictors) {
    rng     <- range(fit$training_data[[p]], na.rm = TRUE)
    too_low <- !is.na(features[[p]]) & features[[p]] < rng[1]
    too_hi  <- !is.na(features[[p]]) & features[[p]] > rng[2]
    out     <- too_low | too_hi
    if (any(out)) {
      message(sprintf(
        "score_reads: %d reads have %s outside training range [%g, %g]; marking unscored.",
        sum(out), p, rng[1], rng[2]
      ))
    }
    features$out_of_range <- features$out_of_range | out
  }
  scorable <- scorable & !features$out_of_range

  if (!any(scorable)) {
    warning("No scorable reads (all rows have NA or out-of-range predictors).")
    attr(features, "flag_threshold") <- NA_real_
    return(features)
  }

  # newdata for predict.gamlss: needs every column referenced in the
  # model formula. pb(<x>_sc) wants <x>_sc; bare piecewise/linear
  # terms (and the I(pmax(...)) basis) want the raw predictor names.
  scaled_sc    <- paste0(scaled_predictors, "_sc")
  newdata_cols <- unique(c(all_predictors, scaled_sc, "y"))
  newdata      <- features[scorable, newdata_cols, drop = FALSE]

  mu_p    <- stats::predict(fit$fit, newdata = newdata, what = "mu",
                            type = "response", data = fit$training_data)
  sigma_p <- stats::predict(fit$fit, newdata = newdata, what = "sigma",
                            type = "response", data = fit$training_data)
  nu_p    <- stats::predict(fit$fit, newdata = newdata, what = "nu",
                            type = "response", data = fit$training_data)
  tau_p   <- stats::predict(fit$fit, newdata = newdata, what = "tau",
                            type = "response", data = fit$training_data)

  rqr <- compute_rqr_beinf(
    y     = newdata$y,
    mu    = mu_p,
    sigma = sigma_p,
    nu    = nu_p,
    tau   = tau_p,
    seed  = seed
  )

  # Predicted point-mass probabilities (BEINF parametrisation: nu and
  # tau are the unnormalised mass weights for y = 0 and y = 1).
  p0_hat <- nu_p  / (1 + nu_p + tau_p)
  p1_hat <- tau_p / (1 + nu_p + tau_p)

  # Binomial deviance residual for "is the read at the zero mass?",
  # restricted to MAPQ = 0 observations (where its sign is unambiguous).
  y_scored <- newdata$y
  at_zero  <- y_scored == 0
  dev_zero <- ifelse(at_zero & p0_hat > 0 & p0_hat < 1,
                     sqrt(-2 * log(p0_hat)),
                     NA_real_)

  features$mu_p[scorable]     <- mu_p
  features$sigma_p[scorable]  <- sigma_p
  features$nu_p[scorable]     <- nu_p
  features$tau_p[scorable]    <- tau_p
  features$p0_hat[scorable]   <- p0_hat
  features$p1_hat[scorable]   <- p1_hat
  features$rqr[scorable]      <- rqr
  features$dev_zero[scorable] <- dev_zero

  # RQR-tail threshold computed on the continuous-MAPQ subset only, so
  # the median/MAD aren't distorted by the point-mass residuals.
  in_continuous <- y_scored > 0 & y_scored < 1
  if (any(in_continuous)) {
    rqr_cont <- rqr[in_continuous]
    flag_thresh <- if (filter_method == "mad") {
      stats::median(rqr_cont) + threshold * stats::mad(rqr_cont)
    } else {
      threshold * stats::sd(rqr_cont) + mean(rqr_cont)
    }
  } else {
    flag_thresh <- NA_real_
  }

  flag_rqr <- !is.na(flag_thresh) & in_continuous & rqr < flag_thresh
  flag_p0  <- at_zero & !is.na(dev_zero) & dev_zero > dev_threshold
  # MAPQ = 60 (y == 1) is never flagged: the read achieved the maximum
  # possible MAPQ and so cannot be "unusually low" given any context.

  flag_low <- flag_rqr | flag_p0
  flag_src <- ifelse(flag_rqr, "rqr_low",
                     ifelse(flag_p0, "p0_dev_high", "none"))

  features$flag_low[scorable]    <- flag_low
  features$flag_source[scorable] <- flag_src

  attr(features, "flag_threshold")     <- flag_thresh
  attr(features, "filter_method")      <- filter_method
  attr(features, "threshold")          <- threshold
  attr(features, "dev_threshold")      <- dev_threshold
  attr(features, "p0_threshold")       <- exp(-dev_threshold^2 / 2)
  features
}
