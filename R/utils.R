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
