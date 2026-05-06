#' Score and flag reads using a fitted BEINF model
#'
#' For every read in `features`, predicts the four BEINF parameters
#' (mu, sigma, nu, tau) under the fitted model, then computes a
#' randomized quantile residual (RQR) and flags reads with anomalously
#' low residuals.
#'
#' For continuous response values, the RQR is the standard normal
#' quantile of the model CDF at that value. For point masses (MAPQ = 0
#' or MAPQ = 60), the residual is randomized over the corresponding CDF
#' interval, following Dunn & Smyth (1996), so that the overall residual
#' distribution is approximately N(0, 1) under a correctly specified
#' model.
#'
#' Reads with `NA` in any required predictor are kept in the output but
#' assigned `rqr = NA` and `flag_low = FALSE` (i.e. retained by default).
#'
#' @param features Data frame from [extract_features()] (typically the
#'   full set, not just the training subset).
#' @param fit Object returned by [fit_beinf()].
#' @param filter_method One of `"mad"` (default) or `"zscore"`.
#'   For `"mad"`: flag if `rqr < median(rqr) + threshold * mad(rqr)`.
#'   For `"zscore"`: flag if `(rqr - mean(rqr)) / sd(rqr) < threshold`.
#' @param threshold Numeric scalar; default `-3`.
#' @param seed Random seed for the uniform draws on point masses.
#' @param primary_only Logical; if `TRUE` (default) score only primary
#'   alignments. Secondary alignments inherit `flag_low = FALSE`.
#'
#' Reads with at least one predictor outside the training-data range
#'   are marked as `out_of_range = TRUE` and not scored, since the
#'   model's cubic-spline extrapolation is unreliable beyond the
#'   data it saw during fitting.
#'
#' @return A data frame with the input columns plus `mu_p`, `sigma_p`,
#'   `nu_p`, `tau_p`, `rqr`, `flag_low`, `out_of_range`, and an
#'   attribute `"flag_threshold"` recording the absolute residual
#'   threshold used.
#'
#' @examples
#' ex       <- calmer_example_data()
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#' clean    <- clean_features(features)
#' fit <- fit_beinf(
#'   clean, n = 1000,
#'   formula       = y ~ n_isoforms_sc + tx_length_sc + gc_content_sc + pct_covered_sc,
#'   sigma_formula = ~ 1, nu_formula = ~ 1, tau_formula = ~ 1
#' )
#' scored <- score_reads(features, fit,
#'                       filter_method = "mad", threshold = -3)
#' table(scored$flag_low)
#'
#' @importFrom stats predict median mad sd complete.cases runif qnorm
#' @export
score_reads <- function(features,
                        fit,
                        filter_method = c("mad", "zscore"),
                        threshold     = -3,
                        seed          = 42,
                        primary_only  = TRUE) {
  filter_method <- match.arg(filter_method)

  if (!inherits(fit, "calmer_fit")) {
    stop("`fit` must be a calmer_fit object returned by fit_beinf().")
  }
  stopifnot(
    is.data.frame(features),
    "MAPQ" %in% names(features),
    all(fit$predictors %in% names(features)),
    is.numeric(threshold), length(threshold) == 1L
  )

  predictors    <- fit$predictors
  predictors_sc <- paste0(predictors, "_sc")
  scale_params  <- fit$scale_params

  for (p in predictors) {
    features[[paste0(p, "_sc")]] <- apply_scale(features[[p]],
                                                scale_params[[p]])
  }
  features$y <- features$MAPQ / 60

  features$mu_p     <- NA_real_
  features$sigma_p  <- NA_real_
  features$nu_p     <- NA_real_
  features$tau_p    <- NA_real_
  features$rqr      <- NA_real_
  features$flag_low <- FALSE

  scorable <- stats::complete.cases(
    features[, c("MAPQ", predictors, predictors_sc), drop = FALSE]
  )
  if (primary_only && "aln_type" %in% names(features)) {
    scorable <- scorable & (features$aln_type == "P")
  }

  features$out_of_range <- FALSE
  for (p in predictors) {
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

  newdata <- features[scorable,
                      unique(c(predictors, predictors_sc, "y")),
                      drop = FALSE]

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

  features$mu_p[scorable]    <- mu_p
  features$sigma_p[scorable] <- sigma_p
  features$nu_p[scorable]    <- nu_p
  features$tau_p[scorable]   <- tau_p
  features$rqr[scorable]     <- rqr

  flag_thresh <- if (filter_method == "mad") {
    stats::median(rqr) + threshold * stats::mad(rqr)
  } else {
    threshold * stats::sd(rqr) + mean(rqr)
  }

  features$flag_low[scorable] <- rqr < flag_thresh

  attr(features, "flag_threshold") <- flag_thresh
  attr(features, "filter_method")  <- filter_method
  attr(features, "threshold")      <- threshold
  features
}
