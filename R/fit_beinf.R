#' Fit a BEINF regression on a sampled subset of reads
#'
#' Samples `n` primary-alignment reads and fits a Beta-Inflated (BEINF)
#' regression of `MAPQ / 60` against the supplied predictors. By default,
#' all numeric predictors enter through penalized B-splines (`pb()`) on
#' pre-scaled versions; users can override any of the four BEINF parameter
#' formulas (mu, sigma, nu, tau) to fit a different functional form.
#'
#' Filtering on annotation features (e.g. dropping rare `n_isoforms`,
#' extreme GC content, very short or very long transcripts) is the
#' caller's responsibility — pass an already-filtered `features` data
#' frame.
#'
#' @param features Data frame from [extract_features()].
#' @param n Number of primary-alignment reads to sample for training
#'   (default 50000). If `features` has fewer primary alignments than
#'   `n`, all are used.
#' @param seed Random seed for the subset (default 42).
#' @param predictors Character vector of column names in `features` to
#'   use as predictors. Default is the four canonical features:
#'   `c("n_isoforms", "tx_length", "gc_content", "pct_covered")`.
#' @param formula Optional gamlss formula for `mu`. If `NULL`, a default
#'   `pb()`-based formula is built from `predictors` after pre-scaling.
#'   In a custom formula, refer to scaled predictors as
#'   `<predictor>_sc` (e.g. `n_isoforms_sc`).
#' @param sigma_formula,nu_formula,tau_formula Optional formulas for the
#'   other BEINF distributional parameters. If `NULL`, the same default
#'   formula is used for each.
#' @param trace Logical; if `TRUE`, print gamlss convergence messages.
#'   Default `FALSE`.
#'
#' @return A list of class `"calmer_fit"` with:
#' \describe{
#'   \item{fit}{Fitted gamlss BEINF model.}
#'   \item{scale_params}{Named list of `list(center, scale)` per
#'     predictor, used to standardize new data identically.}
#'   \item{predictors}{Character vector of predictor names.}
#'   \item{training_data}{Cleaned + scaled training data.}
#'   \item{formula}{The four formulas actually used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \donttest{
#' features <- extract_features(bam, gtf, fasta)
#' fit      <- fit_beinf(features, n = 50000, seed = 42)
#' }
#'
#' @importFrom stats as.formula sd
#' @importFrom gamlss pb
#' @export
fit_beinf <- function(features,
                      n             = 50000,
                      seed          = 42,
                      predictors    = c("n_isoforms", "tx_length",
                                        "gc_content", "pct_covered"),
                      formula       = NULL,
                      sigma_formula = NULL,
                      nu_formula    = NULL,
                      tau_formula   = NULL,
                      trace         = FALSE) {
  call <- match.call()

  stopifnot(
    is.data.frame(features),
    "aln_type" %in% names(features),
    "MAPQ"     %in% names(features),
    is.numeric(n), length(n) == 1L, n > 0,
    is.numeric(seed), length(seed) == 1L,
    all(predictors %in% names(features))
  )

  primary <- features[features$aln_type == "P", , drop = FALSE]
  if (nrow(primary) == 0L) {
    stop("No primary alignments (aln_type == 'P') in features.")
  }

  if (nrow(primary) > n) {
    # Bioconductor flags set.seed() in package code; with_seed() restores
    # the user's RNG state when the block exits.
    primary <- withr::with_seed(
      seed,
      primary[sample(nrow(primary), size = n), , drop = FALSE]
    )
  }

  primary <- primary[stats::complete.cases(primary[, c("MAPQ", predictors)]),
                     , drop = FALSE]

  primary$y <- primary$MAPQ / 60

  scale_params <- stats::setNames(
    lapply(predictors, function(p) compute_scale_params(primary[[p]])),
    predictors
  )

  for (p in predictors) {
    primary[[paste0(p, "_sc")]] <- apply_scale(primary[[p]], scale_params[[p]])
  }

  predictors_sc <- paste0(predictors, "_sc")
  default_rhs <- paste(paste0("pb(", predictors_sc, ", inter = 5)"),
                       collapse = " + ")
  if (is.null(formula))       formula       <- stats::as.formula(paste("y ~", default_rhs))
  if (is.null(sigma_formula)) sigma_formula <- stats::as.formula(paste("~",   default_rhs))
  if (is.null(nu_formula))    nu_formula    <- stats::as.formula(paste("~",   default_rhs))
  if (is.null(tau_formula))   tau_formula   <- stats::as.formula(paste("~",   default_rhs))

  fit <- gamlss::gamlss(
    formula       = formula,
    sigma.formula = sigma_formula,
    nu.formula    = nu_formula,
    tau.formula   = tau_formula,
    family        = gamlss.dist::BEINF(),
    data          = primary,
    control       = gamlss::gamlss.control(trace = trace)
  )

  structure(
    list(
      fit           = fit,
      scale_params  = scale_params,
      predictors    = predictors,
      training_data = primary,
      formula       = list(mu    = formula,
                           sigma = sigma_formula,
                           nu    = nu_formula,
                           tau   = tau_formula),
      call          = call
    ),
    class = "calmer_fit"
  )
}

#' @export
print.calmer_fit <- function(x, ...) {
  cat("CALMeR BEINF fit\n")
  cat("  predictors:      ", paste(x$predictors, collapse = ", "), "\n", sep = "")
  cat("  training reads:  ", nrow(x$training_data), "\n", sep = "")
  cat("  mu formula:      ", deparse1(x$formula$mu), "\n",    sep = "")
  cat("  sigma formula:   ", deparse1(x$formula$sigma), "\n", sep = "")
  cat("  nu formula:      ", deparse1(x$formula$nu), "\n",    sep = "")
  cat("  tau formula:     ", deparse1(x$formula$tau), "\n",   sep = "")
  invisible(x)
}
