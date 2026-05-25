#' Fit a BEINF regression on a sampled subset of reads
#'
#' Samples `n` primary-alignment reads and fits a Beta-Inflated (BEINF)
#' regression of `MAPQ / 60` against the supplied predictors. Predictors
#' are split into three classes that enter the formula differently:
#' \describe{
#'   \item{`scaled_predictors`}{Z-scored on the training subset and
#'         entered through penalized B-splines (`pb(<x>_sc, inter = ...)`)
#'         in the default formula, where the number of intervals is
#'         set by the `inter` argument. Default: `tx_length`,
#'         `gc_content`, `pct_covered`.}
#'   \item{`piecewise_predictors`}{A named list mapping predictor name
#'         to a knot value. The predictor enters as two linear basis
#'         terms with a kink at the knot:
#'         `<p> + I(pmax(<p> - knot, 0))`. Not z-scored. Default:
#'         `list(n_isoforms = 5)`.}
#'   \item{`linear_predictors`}{Plain linear terms (1 df). Not z-scored;
#'         expected to be 0/1 indicators or other already-bounded
#'         covariates. Default: `has_paralog`.}
#' }
#'
#' Filtering on annotation features (e.g. lumping rare `n_isoforms`,
#' log-transforming `tx_length`, winsorizing `pct_covered`) is the
#' caller's responsibility — pass an already-cleaned `features` data
#' frame, typically from [clean_features()]. If the data frame carries
#' a `pct_covered_winsor_bounds` attribute, it is preserved on the fit
#' so [score_reads()] can apply the same bounds to a scoring set.
#'
#' @param features Data frame from [extract_features()] /
#'   [clean_features()].
#' @param n Number of primary-alignment reads to sample for training
#'   (default 50000). If `features` has fewer primary alignments than
#'   `n`, all are used.
#' @param seed Random seed for the subset (default 42).
#' @param scaled_predictors Character vector of predictor columns to
#'   z-score and (in the default formula) wrap in `pb()`. Default
#'   `c("tx_length", "gc_content", "pct_covered")`.
#' @param piecewise_predictors Named list mapping predictor name to
#'   knot value. Each entry enters as `<p> + I(pmax(<p> - knot, 0))`.
#'   Default `list(n_isoforms = 5)`. Pass `list()` to disable.
#' @param linear_predictors Character vector of predictor columns to
#'   enter as plain linear terms (no z-scoring). Default
#'   `c("has_paralog")`. Pass `character(0)` to disable.
#' @param formula Optional gamlss formula for `mu`. If `NULL`, a default
#'   formula is built combining `pb(<x>_sc)` for the scaled set, the
#'   two-piece basis for each piecewise predictor, and raw terms for
#'   the linear set. In a custom formula, refer to scaled predictors
#'   as `<predictor>_sc` (e.g. `tx_length_sc`).
#' @param sigma_formula,nu_formula,tau_formula Optional formulas for the
#'   other BEINF distributional parameters. If `NULL`, each defaults to
#'   the same RHS as `formula` (the default-built one if `formula` is
#'   `NULL`, otherwise the user-supplied one).
#' @param inter Integer; number of equal-spaced intervals for each
#'   penalized B-spline term in the default formula. Forwarded to the
#'   `inter` argument of [gamlss::pb()]. Ignored when a custom
#'   `formula` is supplied. Default `5`.
#' @param trace Logical; if `TRUE`, print gamlss convergence messages.
#'   Default `FALSE`.
#'
#' @return A list of class `"calmer_fit"` with:
#' \describe{
#'   \item{fit}{Fitted gamlss BEINF model.}
#'   \item{scale_params}{Named list of `list(center, scale)` per
#'     scaled predictor, used to standardize new data identically.}
#'   \item{scaled_predictors, piecewise_predictors, linear_predictors}{
#'     The three predictor specs, mirrored from the arguments.}
#'   \item{predictors}{Character vector of all predictor names (union
#'     of the three sets) for backwards-compatible access.}
#'   \item{winsor_bounds}{Length-2 numeric copy of
#'     `attr(features, "pct_covered_winsor_bounds")` if present on the
#'     input, otherwise `NULL`. Used by [score_reads()] to apply
#'     identical winsorize cutoffs to a scoring set.}
#'   \item{training_data}{Cleaned + scaled training subset.}
#'   \item{formula}{The four formulas actually used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' ex       <- calmer_example_data()
#' paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                              paralog_table = paralog)
#' clean    <- clean_features(features)
#' # On the small bundled dataset, override the default pb() formula with
#' # linear terms so the fit converges quickly.
#' fit <- fit_beinf(
#'   clean, n = 1000,
#'   formula       = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
#'                       n_isoforms + I(pmax(n_isoforms - 5, 0)) +
#'                       has_paralog,
#'   sigma_formula = ~ 1,
#'   nu_formula    = ~ 1,
#'   tau_formula   = ~ 1
#' )
#' fit
#'
#' @importFrom stats as.formula sd
#' @importFrom gamlss pb
#' @export
fit_beinf <- function(features,
                      n                    = 50000,
                      seed                 = 42,
                      scaled_predictors    = c("tx_length", "gc_content",
                                               "pct_covered"),
                      piecewise_predictors = list(n_isoforms = 5),
                      linear_predictors    = c("has_paralog"),
                      formula              = NULL,
                      sigma_formula        = NULL,
                      nu_formula           = NULL,
                      tau_formula          = NULL,
                      inter              = 5,
                      trace                = FALSE) {
  call <- match.call()

  all_predictors <- c(scaled_predictors,
                      names(piecewise_predictors),
                      linear_predictors)

  stopifnot(
    is.data.frame(features),
    "aln_type" %in% names(features),
    "MAPQ"     %in% names(features),
    is.numeric(n), length(n) == 1L, n > 0,
    is.numeric(seed), length(seed) == 1L,
    is.numeric(inter), length(inter) == 1L, inter >= 1,
    is.character(scaled_predictors),
    is.list(piecewise_predictors),
    length(piecewise_predictors) == 0L ||
      !is.null(names(piecewise_predictors)),
    is.character(linear_predictors),
    all(all_predictors %in% names(features))
  )
  if (length(piecewise_predictors) > 0L) {
    stopifnot(all(vapply(piecewise_predictors, is.numeric, logical(1))),
              all(vapply(piecewise_predictors, length,    integer(1)) == 1L))
  }

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

  primary <- primary[stats::complete.cases(primary[, c("MAPQ", all_predictors)]),
                     , drop = FALSE]

  primary$y <- primary$MAPQ / 60

  # Z-score only the scaled predictors. Piecewise and linear predictors
  # enter raw (the piecewise basis is computed inside the formula via
  # I(pmax(...))).
  scale_params <- stats::setNames(
    lapply(scaled_predictors, function(p) compute_scale_params(primary[[p]])),
    scaled_predictors
  )
  for (p in scaled_predictors) {
    primary[[paste0(p, "_sc")]] <- apply_scale(primary[[p]],
                                               scale_params[[p]])
  }

  # Build the default RHS from the three predictor classes if the user
  # did not supply one.
  default_rhs <- .build_three_class_rhs(
    scaled_predictors    = scaled_predictors,
    piecewise_predictors = piecewise_predictors,
    linear_predictors    = linear_predictors,
    inter              = inter
  )

  if (is.null(formula)) {
    formula <- stats::as.formula(paste("y ~", default_rhs))
    rhs_for_others <- default_rhs
  } else {
    # Defaults for sigma/nu/tau mirror the user's mu formula RHS.
    rhs_for_others <- paste(deparse(formula[[3L]]), collapse = " ")
  }
  if (is.null(sigma_formula)) sigma_formula <- stats::as.formula(paste("~", rhs_for_others))
  if (is.null(nu_formula))    nu_formula    <- stats::as.formula(paste("~", rhs_for_others))
  if (is.null(tau_formula))   tau_formula   <- stats::as.formula(paste("~", rhs_for_others))

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
      fit                  = fit,
      scale_params         = scale_params,
      scaled_predictors    = scaled_predictors,
      piecewise_predictors = piecewise_predictors,
      linear_predictors    = linear_predictors,
      predictors           = all_predictors,
      winsor_bounds        = attr(features, "pct_covered_winsor_bounds"),
      training_data        = primary,
      formula              = list(mu    = formula,
                                  sigma = sigma_formula,
                                  nu    = nu_formula,
                                  tau   = tau_formula),
      call                 = call
    ),
    class = "calmer_fit"
  )
}

# Internal: build a `pb()` + piecewise + linear RHS string from the
# three predictor specs. Shared by fit_beinf() and fit_beinf_suite().
.build_three_class_rhs <- function(scaled_predictors,
                                   piecewise_predictors,
                                   linear_predictors,
                                   inter,
                                   smoother = c("pb", "linear")) {
  smoother <- match.arg(smoother)
  rhs_terms <- character(0)
  if (length(scaled_predictors) > 0L) {
    scaled_sc <- paste0(scaled_predictors, "_sc")
    rhs_terms <- c(rhs_terms,
                   if (smoother == "pb") {
                     paste0("pb(", scaled_sc, ", inter = ", inter, ")")
                   } else {
                     scaled_sc
                   })
  }
  for (p in names(piecewise_predictors)) {
    knot <- piecewise_predictors[[p]]
    rhs_terms <- c(rhs_terms,
                   p,
                   sprintf("I(pmax(%s - %g, 0))", p, knot))
  }
  if (length(linear_predictors) > 0L) {
    rhs_terms <- c(rhs_terms, linear_predictors)
  }
  if (length(rhs_terms) == 0L) {
    stop("No predictors supplied: scaled/piecewise/linear are all empty.",
         call. = FALSE)
  }
  paste(rhs_terms, collapse = " + ")
}

#' @export
print.calmer_fit <- function(x, ...) {
  cat("CALMeR BEINF fit\n")
  cat("  scaled:    ", paste(x$scaled_predictors, collapse = ", "),
      "\n", sep = "")
  if (length(x$piecewise_predictors) > 0L) {
    knots <- vapply(names(x$piecewise_predictors),
                    function(p) sprintf("%s@%g", p, x$piecewise_predictors[[p]]),
                    character(1))
    cat("  piecewise: ", paste(knots, collapse = ", "), "\n", sep = "")
  }
  if (length(x$linear_predictors) > 0L) {
    cat("  linear:    ", paste(x$linear_predictors, collapse = ", "),
        "\n", sep = "")
  }
  cat("  training reads:  ", nrow(x$training_data), "\n", sep = "")
  cat("  mu formula:      ", deparse1(x$formula$mu), "\n",    sep = "")
  cat("  sigma formula:   ", deparse1(x$formula$sigma), "\n", sep = "")
  cat("  nu formula:      ", deparse1(x$formula$nu), "\n",    sep = "")
  cat("  tau formula:     ", deparse1(x$formula$tau), "\n",   sep = "")
  invisible(x)
}
