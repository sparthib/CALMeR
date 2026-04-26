#' Score and flag reads using a fitted BEINF model
#'
#' Computes randomized quantile residuals (RQR) for every read in `features`
#' using a model fitted by [fit_beinf()], then flags reads with anomalously
#' low residuals using either a MAD-based or z-score-based threshold.
#'
#' @param features Data frame from [extract_features()] (the full read set).
#' @param fit Object returned by [fit_beinf()].
#' @param filter_method Either `"mad"` (default) or `"zscore"`.
#' @param threshold Numeric scalar; threshold applied to the residual scale.
#'   For `"mad"`: flag if `rqr < median(rqr) + threshold * mad(rqr)`.
#'   For `"zscore"`: flag if `(rqr - mean(rqr)) / sd(rqr) < threshold`.
#'   Default `-3`.
#' @param seed Random seed for the RQR uniform draws on point masses.
#'
#' @return A data.frame with the original `features` columns plus
#'   `mu_p`, `sigma_p`, `nu_p`, `tau_p`, `rqr`, and `flag_low`.
#'
#' @export
score_reads <- function(features,
                        fit,
                        filter_method = c("mad", "zscore"),
                        threshold     = -3,
                        seed          = 42) {
  stop("Not yet implemented.")
}
