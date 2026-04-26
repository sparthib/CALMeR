#' Fit a BEINF regression on a sampled subset of reads
#'
#' Samples `n` primary-alignment reads (with optional pre-filtering on
#' transcript features) and fits a Beta-Inflated (BEINF) regression of
#' `MAPQ / 60` against the supplied predictors. Default uses penalized
#' B-splines (`pb()`) for all numeric predictors.
#'
#' @param features Data frame from [extract_features()].
#' @param n Number of reads to sample for training (default 50000).
#' @param seed Random seed for reproducible sampling (default 42).
#' @param predictors Character vector of predictor column names. Default
#'   `c("n_isoforms", "tx_length", "gc_content", "pct_covered")`.
#' @param formula Optional gamlss formula for `mu`. If `NULL`, a default
#'   `pb()`-based formula is built from `predictors`.
#' @param sigma_formula,nu_formula,tau_formula Optional formulas for the
#'   other BEINF parameters. If `NULL`, the same default formula is used.
#' @param prefilter Named list of pre-filtering predicates applied before
#'   sampling. See vignette for examples.
#'
#' @return A list with components:
#' \describe{
#'   \item{fit}{Fitted gamlss BEINF model.}
#'   \item{scale_params}{Scaling parameters used to standardize predictors.}
#'   \item{training_data}{The cleaned training data frame used for the fit.}
#'   \item{call}{The matched call.}
#' }
#'
#' @export
fit_beinf <- function(features,
                      n              = 50000,
                      seed           = 42,
                      predictors     = c("n_isoforms", "tx_length",
                                         "gc_content", "pct_covered"),
                      formula        = NULL,
                      sigma_formula  = NULL,
                      nu_formula     = NULL,
                      tau_formula    = NULL,
                      prefilter      = list()) {
  stop("Not yet implemented.")
}
