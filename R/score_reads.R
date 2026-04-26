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
#' @return A data frame with the input columns plus `mu_p`, `sigma_p`,
#'   `nu_p`, `tau_p`, `rqr`, `flag_low`, and an integer attribute
#'   `"flag_threshold"` recording the absolute residual threshold used.
#'
#' @examples
#' \dontrun{
#' features <- extract_features(bam, gtf, fasta)
#' clean    <- clean_features(features)
#' fit      <- fit_beinf(clean)
#' scored   <- score_reads(features, fit,
#'                         filter_method = "mad", threshold = -3)
#' }
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
  if (!any(scorable)) {
    warning("No scorable reads (all rows have NA predictors).")
    attr(features, "flag_threshold") <- NA_real_
    return(features)
  }

  newdata <- features[scorable,
                      unique(c(predictors, predictors_sc, "y")),
                      drop = FALSE]

  mu_p    <- stats::predict(fit$fit, newdata = newdata, what = "mu",
                            type = "response")
  sigma_p <- stats::predict(fit$fit, newdata = newdata, what = "sigma",
                            type = "response")
  nu_p    <- stats::predict(fit$fit, newdata = newdata, what = "nu",
                            type = "response")
  tau_p   <- stats::predict(fit$fit, newdata = newdata, what = "tau",
                            type = "response")

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
