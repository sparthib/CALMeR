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
  if (!is.null(seed)) set.seed(seed)
  p1        <- tau / (1 + nu + tau)
  F_y       <- gamlss.dist::pBEINF(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
  F_y_minus <- ifelse(y == 0, 0,
                      ifelse(y == 1, 1 - p1, F_y))
  u <- stats::runif(length(y), min = F_y_minus, max = F_y)
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
