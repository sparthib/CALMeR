#' Backward GAMLSS model selection for the BEINF MAPQ model
#'
#' Runs principled backward term selection for the CALMeR Beta-Inflated
#' (BEINF) model of `MAPQ / 60`, starting from the \emph{full} model
#' (penalized B-spline smooths for the scaled predictors, the two-piece
#' isoform basis, the paralog indicator, and the transcript-feature
#' interactions) and dropping terms that do not improve a generalized
#' AIC (GAIC). Two GAMLSS strategies are available:
#' \describe{
#'   \item{Strategy A (`stepGAICAll.A`)}{Selects terms for each
#'         distribution parameter in turn (\eqn{\mu}, then \eqn{\sigma},
#'         \eqn{\nu}, \eqn{\tau}), holding the others fixed at each
#'         sweep.}
#'   \item{Strategy B (`stepGAICAll.B`)}{A single interleaved search
#'         over all parameters at once.}
#' }
#' In both, at each step the term whose removal most improves GAIC is
#' dropped, subject to marginality (an interaction is retained whenever
#' a parent term is in the model). The penalty is `k = log(n)` by
#' default, i.e. the Bayesian Information Criterion (BIC); pass `k = 2`
#' for AIC.
#'
#' The training frame is built identically to [fit_beinf()] (sample `n`
#' primary-alignment reads at `seed`, drop incomplete cases, z-score the
#' scaled predictors). The two-piece isoform basis and all interaction
#' terms are materialized as explicit columns (`<p>_hi`, etc.) so the
#' stepwise scope is well-formed.
#'
#' @section Evaluation routes:
#' Beyond the training-penalty GAIC used by the stepwise search itself,
#' the candidate models (full, Strategy-A-selected, Strategy-B-selected)
#' can be scored on held-out data:
#' \describe{
#'   \item{`"gaic"`}{Always computed: AIC, BIC, residual global
#'         deviance, effective df, Nagelkerke (Cragg-Uhler) pseudo-R2 on
#'         the training subset.}
#'   \item{`"holdout"`}{A single train/validation split
#'         (`holdout_frac` held out); validation global deviance via
#'         [gamlss::gamlssVGD()]. This is the large-\eqn{n} route
#'         recommended by Stasinopoulos et al. (2017).}
#'   \item{`"cv"`}{`cv_folds`-fold cross-validated global deviance via
#'         [gamlss::gamlssCV()]. A supporting robustness check; much
#'         more expensive (K refits per candidate).}
#' }
#'
#' @param features Cleaned training data frame from [clean_features()].
#' @param scaled_predictors,piecewise_predictors,linear_predictors
#'   Predictor specs, same semantics as [fit_beinf()]. Defaults match
#'   the canonical CALMeR pipeline.
#' @param interactions Character vector of interaction terms to include
#'   in the full start model, expressed against the materialized columns
#'   (`<scaled>_sc`, `<piecewise>`, `<piecewise>_hi`, `<linear>`). If
#'   `NULL` (default), the three manuscript interactions are used:
#'   length x coverage, length x isoform, and paralog x isoform (the
#'   latter two each expanded over the two-piece isoform basis).
#' @param strategies Character vector; which backward strategies to run.
#'   Any of `"A"`, `"B"`. Default both.
#' @param evaluate Character vector of evaluation routes. `"gaic"` is
#'   always included. Add `"holdout"` and/or `"cv"` to score candidates
#'   on held-out data. Default `"gaic"`.
#' @param select_sigma Logical; if `FALSE` (default) `sigma` is held at
#'   `~ 1` (matching the manuscript) and excluded from selection. If
#'   `TRUE`, `sigma` starts from the full RHS and is selected like the
#'   other parameters.
#' @param holdout_frac Fraction of the training subset held out for the
#'   `"holdout"` route. Default 0.2.
#' @param cv_folds Number of folds for the `"cv"` route. Default 5.
#' @param n,seed,inter Forwarded to the training-frame construction;
#'   same defaults as [fit_beinf()] (50000, 42, 5).
#' @param k GAIC penalty per parameter. If `NULL` (default), `log(n)`
#'   where `n` is the number of training reads actually used (BIC).
#' @param ncpus Integer; cores for the parallel stepwise search
#'   (`parallel = "multicore"` when `> 1`). Default 1.
#' @param trace Logical; print gamlss / stepwise progress. Default
#'   `FALSE`.
#'
#' @return An object of class `"calmer_selection"`: a list with
#' \describe{
#'   \item{full}{The fitted full-model gamlss object.}
#'   \item{selected}{Named list of selected gamlss fits, one per
#'     strategy run (`A`, `B`).}
#'   \item{formulas}{Named list of the selected four-parameter formulas
#'     per strategy, plus `full`.}
#'   \item{comparison}{Data frame, one row per candidate model, with
#'     AIC / BIC / deviance / df / Nagelkerke (and validation deviance /
#'     CV deviance when those routes are requested), plus a `picked`
#'     flag (min BIC among candidates).}
#'   \item{scale_params}{Per-scaled-predictor center/scale, as in
#'     [fit_beinf()].}
#'   \item{training_data}{The cleaned + scaled + basis-materialized
#'     training subset.}
#'   \item{k, call}{The GAIC penalty used and the matched call.}
#' }
#'
#' @references Stasinopoulos, M. D., Rigby, R. A., Heller, G. Z.,
#'   Voudouris, V., & De Bastiani, F. (2017). \emph{Flexible Regression
#'   and Smoothing: Using GAMLSS in R}. Chapman & Hall/CRC.
#'
#' @seealso [fit_beinf()], [score_reads()].
#'
#' @examples
#' \dontrun{
#'   ex       <- calmer_example_data()
#'   paralog  <- make_paralog_table(ex$fasta, ex$gtf, paf_path = ex$paf)
#'   features <- extract_features(ex$bam, ex$gtf, ex$fasta,
#'                                paralog_table = paralog)
#'   clean    <- clean_features(features)
#'   sel <- select_beinf(clean, n = 50000,
#'                        strategies = c("A", "B"),
#'                        evaluate   = c("gaic", "holdout"))
#'   sel
#'   sel$comparison
#' }
#'
#' @importFrom stats as.formula sd complete.cases AIC BIC
#' @export
select_beinf <- function(features,
                         scaled_predictors    = c("tx_length", "gc_content",
                                                  "pct_covered"),
                         piecewise_predictors = list(n_isoforms = 5),
                         linear_predictors    = c("has_paralog"),
                         interactions         = NULL,
                         strategies           = c("A", "B"),
                         evaluate             = c("gaic"),
                         select_sigma         = FALSE,
                         holdout_frac         = 0.2,
                         cv_folds             = 5,
                         n                    = 50000,
                         seed                 = 42,
                         inter                = 5,
                         k                    = NULL,
                         ncpus                = 1,
                         trace                = FALSE) {
  call <- match.call()

  strategies <- match.arg(strategies, c("A", "B"), several.ok = TRUE)
  evaluate   <- unique(c("gaic", intersect(evaluate, c("holdout", "cv"))))

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
    length(piecewise_predictors) == 0L || !is.null(names(piecewise_predictors)),
    is.character(linear_predictors),
    all(all_predictors %in% names(features)),
    is.numeric(holdout_frac), holdout_frac > 0, holdout_frac < 1,
    is.numeric(cv_folds), cv_folds >= 2
  )

  # ----- training frame: replicate fit_beinf() prep exactly -----
  primary <- features[features$aln_type == "P", , drop = FALSE]
  if (nrow(primary) == 0L) stop("No primary alignments (aln_type == 'P').")
  if (nrow(primary) > n) {
    primary <- withr::with_seed(
      seed, primary[sample(nrow(primary), size = n), , drop = FALSE])
  }
  primary <- primary[stats::complete.cases(primary[, c("MAPQ", all_predictors)]),
                     , drop = FALSE]
  primary$y <- primary$MAPQ / 60

  winsor_bounds <- attr(features, "pct_covered_winsor_bounds")

  scale_params <- stats::setNames(
    lapply(scaled_predictors, function(p) compute_scale_params(primary[[p]])),
    scaled_predictors)
  for (p in scaled_predictors) {
    primary[[paste0(p, "_sc")]] <- apply_scale(primary[[p]], scale_params[[p]])
  }
  # Materialize the two-piece isoform basis so interactions referencing
  # it (and the stepwise scope) are well-formed.
  for (p in names(piecewise_predictors)) {
    knot <- piecewise_predictors[[p]]
    primary[[paste0(p, "_hi")]] <- pmax(primary[[p]] - knot, 0)
  }

  n_train <- nrow(primary)
  if (is.null(k)) k <- log(n_train)

  # ----- build the full-model RHS from materialized columns -----
  if (is.null(interactions)) {
    interactions <- .default_beinf_interactions(
      scaled_predictors, piecewise_predictors, linear_predictors)
  }
  main_terms <- c(
    paste0("pb(", paste0(scaled_predictors, "_sc"), ", inter = ", inter, ")"),
    unlist(lapply(names(piecewise_predictors),
                  function(p) c(p, paste0(p, "_hi"))), use.names = FALSE),
    linear_predictors
  )
  full_rhs <- paste(c(main_terms, interactions), collapse = " + ")
  full_formula  <- stats::as.formula(paste("y ~", full_rhs))
  sigma_formula <- if (select_sigma) full_formula else stats::as.formula("~ 1")
  upper_mu      <- stats::as.formula(paste("~", full_rhs))
  lower         <- stats::as.formula("~ 1")

  # ----- fit the full start model -----
  if (trace) message("select_beinf: fitting full start model (n = ", n_train, ") ...")
  full_fit <- gamlss::gamlss(
    formula       = full_formula,
    sigma.formula = sigma_formula,
    nu.formula    = stats::as.formula(paste("~", full_rhs)),
    tau.formula   = stats::as.formula(paste("~", full_rhs)),
    family        = gamlss.dist::BEINF(),
    data          = primary,
    control       = gamlss::gamlss.control(trace = trace)
  )

  par_arg <- if (ncpus > 1L) "multicore" else "no"

  # ----- run the requested backward strategies -----
  # Each strategy is wrapped so a crash in one (e.g. a parameter equation
  # that gamlss deems too complex for the data at small n) records NULL
  # and does not abort the whole selection.
  selected <- list()
  if ("A" %in% strategies) {
    if (trace) message("select_beinf: Strategy A (stepGAICAll.A, backward) ...")
    selected$A <- tryCatch(
      gamlss::stepGAICAll.A(
        full_fit,
        scope      = list(lower = lower, upper = upper_mu),
        nu.scope   = list(lower = lower, upper = upper_mu),
        tau.scope  = list(lower = lower, upper = upper_mu),
        sigma.try  = select_sigma,
        mu.try = TRUE, nu.try = TRUE, tau.try = TRUE,
        direction  = "backward",
        k          = k,
        parallel   = par_arg, ncpus = ncpus,
        trace      = trace),
      error = function(e) {
        warning("select_beinf: Strategy A failed (", conditionMessage(e),
                "); dropping it.", call. = FALSE); NULL })
  }
  if ("B" %in% strategies) {
    if (trace) message("select_beinf: Strategy B (stepGAICAll.B, backward) ...")
    selected$B <- tryCatch(
      gamlss::stepGAICAll.B(
        full_fit,
        scope     = list(lower = lower, upper = upper_mu),
        direction = "backward",
        k         = k,
        parallel  = par_arg, ncpus = ncpus,
        trace     = trace),
      error = function(e) {
        warning("select_beinf: Strategy B failed (", conditionMessage(e),
                "); dropping it.", call. = FALSE); NULL })
  }
  # Drop any strategy that failed (NULL) before assembling candidates.
  selected <- selected[!vapply(selected, is.null, logical(1))]

  # ----- assemble candidate set + comparison table -----
  candidates <- c(list(full = full_fit), selected)
  rand <- NULL
  if ("holdout" %in% evaluate) {
    rand <- withr::with_seed(
      seed + 1L,
      sample(c(1L, 2L), n_train, replace = TRUE,
             prob = c(1 - holdout_frac, holdout_frac)))
  }
  comparison <- do.call(rbind, lapply(names(candidates), function(nm) {
    .select_metrics_row(candidates[[nm]], nm, primary, rand,
                        cv_folds, seed, evaluate, trace)
  }))
  # pick = min BIC among candidates that have one
  comparison$picked <- FALSE
  okb <- !is.na(comparison$BIC)
  if (any(okb)) comparison$picked[which(okb)[which.min(comparison$BIC[okb])]] <- TRUE

  formulas <- lapply(candidates, .beinf_formula_list)

  # Wrap the picked model as a scorable calmer_fit (same structure as
  # fit_beinf()'s return) so score_reads() / calmer_report() can consume
  # a calmer_selection directly.
  picked_name <- comparison$model[comparison$picked][1]
  picked_fit  <- .wrap_calmer_fit(
    candidates[[picked_name]],
    scaled_predictors    = scaled_predictors,
    piecewise_predictors = piecewise_predictors,
    linear_predictors    = linear_predictors,
    scale_params         = scale_params,
    winsor_bounds        = winsor_bounds,
    training_data        = primary,
    call                 = call
  )

  structure(
    list(
      fit           = picked_fit,
      picked        = picked_name,
      full          = full_fit,
      selected      = selected,
      formulas      = formulas,
      comparison    = comparison,
      scale_params  = scale_params,
      training_data = primary,
      k             = k,
      call          = call
    ),
    class = "calmer_selection"
  )
}

# Internal: wrap a fitted gamlss BEINF model in the calmer_fit structure
# that score_reads() / calmer_report() expect. Mirrors fit_beinf()'s
# return list exactly.
.wrap_calmer_fit <- function(gfit, scaled_predictors, piecewise_predictors,
                             linear_predictors, scale_params, winsor_bounds,
                             training_data, call) {
  fl <- .beinf_formula_list(gfit)
  structure(
    list(
      fit                  = gfit,
      scale_params         = scale_params,
      scaled_predictors    = scaled_predictors,
      piecewise_predictors = piecewise_predictors,
      linear_predictors    = linear_predictors,
      predictors           = c(scaled_predictors,
                               names(piecewise_predictors),
                               linear_predictors),
      winsor_bounds        = winsor_bounds,
      training_data        = training_data,
      formula              = fl,
      call                 = call
    ),
    class = "calmer_fit"
  )
}

#' Resolve a scorable fit from a calmer_fit or calmer_selection
#'
#' Internal-ish helper used by [score_reads()] and [calmer_report()] to
#' accept either a single [fit_beinf()] fit or the picked model of a
#' [select_beinf()] selection transparently.
#'
#' @param x A `calmer_fit` or `calmer_selection` object.
#' @return A `calmer_fit`.
#' @keywords internal
#' @export
primary_fit <- function(x) {
  if (inherits(x, "calmer_fit")) return(x)
  if (inherits(x, "calmer_selection")) {
    if (is.null(x$fit) || !inherits(x$fit, "calmer_fit")) {
      stop("calmer_selection has no scorable picked fit.", call. = FALSE)
    }
    return(x$fit)
  }
  stop("Expected a calmer_fit or calmer_selection object; got ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}

# Internal: the manuscript default interaction set, against materialized
# columns. length x coverage; length x isoform (2-piece); paralog x
# isoform (2-piece). Only emitted when the constituent predictors exist.
.default_beinf_interactions <- function(scaled_predictors,
                                        piecewise_predictors,
                                        linear_predictors) {
  ints <- character(0)
  has_tx   <- "tx_length"   %in% scaled_predictors
  has_pcov <- "pct_covered" %in% scaled_predictors
  has_niso <- "n_isoforms"  %in% names(piecewise_predictors)
  has_para <- "has_paralog" %in% linear_predictors
  if (has_tx && has_pcov) ints <- c(ints, "tx_length_sc:pct_covered_sc")
  if (has_tx && has_niso)
    ints <- c(ints, "tx_length_sc:n_isoforms", "tx_length_sc:n_isoforms_hi")
  if (has_para && has_niso)
    ints <- c(ints, "has_paralog:n_isoforms", "has_paralog:n_isoforms_hi")
  ints
}

# Internal: pull the 4 parameter formulas out of a fitted gamlss object.
# gamlss stores these as $<par>.formula slots; the stepGAICAll return
# object carries the selected-model formulas in the same slots.
.beinf_formula_list <- function(m) {
  if (is.null(m)) return(NULL)
  list(mu = m$mu.formula, sigma = m$sigma.formula,
       nu = m$nu.formula, tau = m$tau.formula)
}

# Internal: one comparison-table row for a candidate fit, with optional
# held-out (VGD) and CV deviance.
.select_metrics_row <- function(m, nm, data, rand, cv_folds, seed,
                                evaluate, trace) {
  if (is.null(m)) {
    return(data.frame(model = nm, AIC = NA_real_, BIC = NA_real_,
                      deviance = NA_real_, df = NA_real_,
                      Nagelkerke = NA_real_, VGD = NA_real_, CV = NA_real_,
                      stringsAsFactors = FALSE))
  }
  nag <- tryCatch(suppressWarnings(gamlss::Rsq(m, type = "Cragg Uhler")),
                  error = function(e) NA_real_)
  fl  <- .beinf_formula_list(m)

  vgd_val <- NA_real_
  if ("holdout" %in% evaluate && !is.null(rand)) {
    vgd_val <- tryCatch({
      vg <- gamlss::gamlssVGD(
        formula = fl$mu, sigma.formula = fl$sigma,
        nu.formula = fl$nu, tau.formula = fl$tau,
        data = data, family = gamlss.dist::BEINF(), rand = rand,
        control = gamlss::gamlss.control(trace = trace))
      as.numeric(vg$oVGD)
    }, error = function(e) { if (trace) message("  VGD failed (", nm, "): ", conditionMessage(e)); NA_real_ })
  }

  cv_val <- NA_real_
  if ("cv" %in% evaluate) {
    cv_val <- tryCatch({
      cvo <- gamlss::gamlssCV(
        formula = fl$mu, sigma.formula = fl$sigma,
        nu.formula = fl$nu, tau.formula = fl$tau,
        data = data, family = gamlss.dist::BEINF(),
        K.fold = cv_folds, set.seed = seed,
        control = gamlss::gamlss.control(trace = trace))
      as.numeric(gamlss::CV(cvo))
    }, error = function(e) { if (trace) message("  CV failed (", nm, "): ", conditionMessage(e)); NA_real_ })
  }

  data.frame(
    model      = nm,
    AIC        = tryCatch(stats::AIC(m),   error = function(e) NA_real_),
    BIC        = tryCatch(stats::BIC(m),   error = function(e) NA_real_),
    deviance   = tryCatch(m$G.deviance,    error = function(e) NA_real_),
    df         = tryCatch(m$df.fit,        error = function(e) NA_real_),
    Nagelkerke = as.numeric(nag),
    VGD        = vgd_val,
    CV         = cv_val,
    stringsAsFactors = FALSE
  )
}

#' @export
print.calmer_selection <- function(x, ...) {
  cat("CALMeR BEINF backward selection\n")
  cat("  strategies: ", paste(names(x$selected), collapse = ", "), "\n", sep = "")
  cat("  GAIC penalty k = ", format(x$k, digits = 5),
      if (isTRUE(all.equal(x$k, 2))) "  (AIC)" else "  (~ BIC)", "\n", sep = "")
  cat("  training reads: ", nrow(x$training_data), "\n", sep = "")
  picked <- x$comparison$model[x$comparison$picked]
  cat("  picked (min BIC): ", if (length(picked)) picked else "NA", "\n", sep = "")
  cat("\n")
  cmp <- x$comparison
  num <- vapply(cmp, is.numeric, logical(1))
  cmp[num] <- lapply(cmp[num], function(v) round(v, 3))
  print(cmp, row.names = FALSE)
  invisible(x)
}
