# Internal helpers — not exported.

# Compute scaling parameters (center + scale) for a numeric vector.
compute_scale_params <- function(x) {
  list(center = mean(x, na.rm = TRUE),
       scale  = stats::sd(x, na.rm = TRUE))
}

# Apply pre-computed scaling parameters.
apply_scale <- function(x, params) {
  (x - params$center) / params$scale
}

# Compute randomized quantile residuals for a BEINF fit on new data.
# y: response vector (already in [0, 1])
# mu, sigma, nu, tau: predicted parameter vectors (response scale)
# seed: optional integer for reproducibility
compute_rqr_beinf <- function(y, mu, sigma, nu, tau, seed = NULL) {
  p1        <- tau / (1 + nu + tau)
  F_y       <- gamlss.dist::pBEINF(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
  F_y_minus <- ifelse(y == 0, 0,
                      ifelse(y == 1, 1 - p1, F_y))
  # Bioconductor flags set.seed() in package code; with_seed() restores
  # the user's RNG state when the block exits.
  u <- if (is.null(seed)) {
    stats::runif(length(y), min = F_y_minus, max = F_y)
  } else {
    withr::with_seed(seed,
      stats::runif(length(y), min = F_y_minus, max = F_y)
    )
  }
  stats::qnorm(u)
}

# Build a default pb()-based gamlss formula from a vector of scaled
# predictor names. `linear` names enter as plain linear terms.
build_default_formula <- function(predictors_sc, linear = character(0),
                                   response = "y") {
  rhs_terms <- c(
    paste0("pb(", predictors_sc, ")"),
    linear
  )
  stats::as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
}

# Rank a calmer_suite by AIC, BIC, and Nagelkerke pseudo-R^2 (Cragg-Uhler)
# and pick the model with the best (= smallest) average rank. Returns
# a data.frame with one row per fit and columns:
#   Model, AIC, BIC, df, Nagelkerke,
#   rank_AIC, rank_BIC, rank_Nagelkerke, avg_rank, picked
# Fits that error out on any metric get NA for that column and NA
# avg_rank; they are excluded from the pick.
.rank_suite <- function(suite) {
  if (!inherits(suite, "calmer_suite") || length(suite) == 0L) {
    stop(".rank_suite needs a non-empty calmer_suite.", call. = FALSE)
  }

  metrics_one <- function(m, nm) {
    if (is.null(m) || is.null(m$fit)) {
      return(data.frame(
        Model = nm, AIC = NA_real_, BIC = NA_real_,
        df = NA_real_, Nagelkerke = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    nag <- tryCatch(
      suppressWarnings(gamlss::Rsq(m$fit, type = "Cragg Uhler")),
      error = function(e) NA_real_
    )
    data.frame(
      Model      = nm,
      AIC        = tryCatch(stats::AIC(m$fit), error = function(e) NA_real_),
      BIC        = tryCatch(stats::BIC(m$fit), error = function(e) NA_real_),
      df         = tryCatch(m$fit$df.fit,      error = function(e) NA_real_),
      Nagelkerke = as.numeric(nag),
      stringsAsFactors = FALSE
    )
  }

  tbl <- do.call(rbind, lapply(names(suite),
                               function(nm) metrics_one(suite[[nm]], nm)))
  if (is.null(tbl) || nrow(tbl) == 0L) {
    stop(".rank_suite: no usable metrics from any suite fit.", call. = FALSE)
  }

  # Lower AIC/BIC = better; higher Nagelkerke = better. Rank 1 = best.
  # Use ties.method = "average" so tied fits share the rank position.
  safe_rank <- function(x, decreasing = FALSE) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    r <- rep(NA_real_, length(x))
    ok <- !is.na(x)
    r[ok] <- if (decreasing) {
      rank(-x[ok], ties.method = "average")
    } else {
      rank( x[ok], ties.method = "average")
    }
    r
  }
  tbl$rank_AIC        <- safe_rank(tbl$AIC)
  tbl$rank_BIC        <- safe_rank(tbl$BIC)
  tbl$rank_Nagelkerke <- safe_rank(tbl$Nagelkerke, decreasing = TRUE)

  tbl$avg_rank <- rowMeans(
    tbl[, c("rank_AIC", "rank_BIC", "rank_Nagelkerke")],
    na.rm = FALSE
  )

  ok <- !is.na(tbl$avg_rank)
  tbl$picked <- FALSE
  if (any(ok)) {
    pick_idx <- which(ok)[which.min(tbl$avg_rank[ok])]
    tbl$picked[pick_idx] <- TRUE
  }
  tbl
}

# Extract a long-format coefficient table from a gamlss BEINF fit,
# one row per (parameter, term). Returns:
#   parameter ("mu" / "sigma" / "nu" / "tau"), term, estimate, se,
#   lower, upper (95% Wald CI), is_pb (TRUE if the term is a pb()
#   smoother). The forest-plot chunk filters out is_pb rows.
.extract_coef_table <- function(model, ci_z = 1.96) {
  s <- tryCatch(
    suppressWarnings(suppressMessages(
      summary(model, save = TRUE, type = "qr")
    )),
    error = function(e) NULL
  )
  if (is.null(s)) return(data.frame())

  params <- c("mu", "sigma", "nu", "tau")
  rows <- lapply(params, function(p) {
    tbl <- s[[paste0(p, ".coef.table")]]
    if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
    df <- data.frame(
      parameter = p,
      term      = rownames(tbl),
      estimate  = tbl[, "Estimate"],
      se        = tbl[, "Std. Error"],
      stringsAsFactors = FALSE
    )
    df$lower <- df$estimate - ci_z * df$se
    df$upper <- df$estimate + ci_z * df$se
    df$is_pb <- grepl("^pb\\(", df$term)
    df
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}
