# Unit tests for the formula-building helper. These are pure-string
# manipulations and run in milliseconds.

test_that(".build_three_class_rhs builds a pb + piecewise + linear RHS by default", {
  rhs <- CALMeR:::.build_three_class_rhs(
    scaled_predictors    = c("tx_length", "gc_content", "pct_covered"),
    piecewise_predictors = list(n_isoforms = 5),
    linear_predictors    = c("has_paralog"),
    inter              = 5,
    smoother             = "pb"
  )
  expect_equal(
    rhs,
    paste(
      "pb(tx_length_sc, inter = 5)",
      "pb(gc_content_sc, inter = 5)",
      "pb(pct_covered_sc, inter = 5)",
      "n_isoforms",
      "I(pmax(n_isoforms - 5, 0))",
      "has_paralog",
      sep = " + "
    )
  )
})

test_that(".build_three_class_rhs respects the linear smoother", {
  rhs <- CALMeR:::.build_three_class_rhs(
    scaled_predictors    = c("tx_length", "gc_content", "pct_covered"),
    piecewise_predictors = list(n_isoforms = 5),
    linear_predictors    = c("has_paralog"),
    inter              = 5,
    smoother             = "linear"
  )
  # Scaled predictors enter as plain `<x>_sc` terms (no pb()).
  expect_equal(
    rhs,
    paste(
      "tx_length_sc", "gc_content_sc", "pct_covered_sc",
      "n_isoforms",   "I(pmax(n_isoforms - 5, 0))",
      "has_paralog",
      sep = " + "
    )
  )
})

test_that(".build_three_class_rhs omits the linear group when empty", {
  rhs <- CALMeR:::.build_three_class_rhs(
    scaled_predictors    = c("tx_length"),
    piecewise_predictors = list(n_isoforms = 5),
    linear_predictors    = character(0),
    inter              = 5,
    smoother             = "pb"
  )
  expect_false(grepl("has_paralog", rhs, fixed = TRUE))
  expect_match(rhs, "pb\\(tx_length_sc, inter = 5\\)")
  expect_match(rhs, "I\\(pmax\\(n_isoforms - 5, 0\\)\\)")
})

test_that(".build_three_class_rhs handles a non-default knot value", {
  rhs <- CALMeR:::.build_three_class_rhs(
    scaled_predictors    = character(0),
    piecewise_predictors = list(n_isoforms = 3),
    linear_predictors    = character(0),
    inter              = 5,
    smoother             = "pb"
  )
  expect_equal(rhs, "n_isoforms + I(pmax(n_isoforms - 3, 0))")
})

test_that(".build_three_class_rhs errors when all three classes are empty", {
  expect_error(
    CALMeR:::.build_three_class_rhs(
      scaled_predictors    = character(0),
      piecewise_predictors = list(),
      linear_predictors    = character(0),
      inter              = 5,
      smoother             = "pb"
    ),
    "No predictors"
  )
})

# Smoke test: synthetic data + intercept-only sigma/nu/tau formulas so
# the BEINF fit converges quickly. Goal is to verify the new return
# shape and that score_reads can consume it end-to-end, not to make
# inference claims.
make_synthetic <- function(n = 400, seed = 42) {
  withr::with_seed(seed, {
    # MAPQ: a mix of zeros (point mass at 0), sixties (point mass at 1
    # after y = MAPQ/60), and a uniform continuous tail so BEINF has
    # support in all three regions.
    n_zero <- 40
    n_sixty <- floor(0.6 * n)
    n_cont <- n - n_zero - n_sixty
    mapq <- c(rep(0L, n_zero),
              rep(60L, n_sixty),
              sample(1:59, n_cont, replace = TRUE))
    df <- data.frame(
      aln_type    = rep("P", n),
      MAPQ        = sample(mapq),
      n_isoforms  = sample(1:7, n, replace = TRUE),
      # tx_length pre-logged, gc/pct in their natural ranges.
      tx_length   = stats::runif(n, log(100), log(10000)),
      gc_content  = stats::runif(n, 30, 70),
      pct_covered = stats::runif(n, 50, 100),
      has_paralog = stats::runif(n) < 0.3,
      stringsAsFactors = FALSE
    )
    attr(df, "pct_covered_winsor_bounds") <- c(50, 100)
    df
  })
}

test_that("fit_beinf returns the three-class structure and preserves winsor_bounds", {
  df  <- make_synthetic(400)
  fit <- fit_beinf(
    df, n = 300, seed = 1,
    formula       = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
                        n_isoforms + I(pmax(n_isoforms - 5, 0)) + has_paralog,
    sigma_formula = ~ 1,
    nu_formula    = ~ 1,
    tau_formula   = ~ 1
  )

  expect_s3_class(fit, "calmer_fit")
  expect_equal(fit$scaled_predictors,
               c("tx_length", "gc_content", "pct_covered"))
  expect_equal(names(fit$piecewise_predictors), "n_isoforms")
  expect_equal(fit$piecewise_predictors$n_isoforms, 5)
  expect_equal(fit$linear_predictors, "has_paralog")
  expect_setequal(
    fit$predictors,
    c("tx_length", "gc_content", "pct_covered", "n_isoforms", "has_paralog")
  )
  # Only the scaled set is z-scored, so scale_params has exactly those names.
  expect_named(fit$scale_params,
               c("tx_length", "gc_content", "pct_covered"))
  # Winsor bounds carried through from the input attribute.
  expect_equal(fit$winsor_bounds, c(50, 100))
})

test_that("score_reads round-trips on training data with the new fit shape", {
  df  <- make_synthetic(400)
  fit <- fit_beinf(
    df, n = 300, seed = 1,
    formula       = y ~ tx_length_sc + gc_content_sc + pct_covered_sc +
                        n_isoforms + I(pmax(n_isoforms - 5, 0)) + has_paralog,
    sigma_formula = ~ 1,
    nu_formula    = ~ 1,
    tau_formula   = ~ 1
  )

  scored <- score_reads(df, fit, filter_method = "mad", threshold = -3)

  # New output columns are present and have values for scorable rows.
  for (col in c("mu_p", "sigma_p", "nu_p", "tau_p",
                "rqr", "flag_low", "flag_source", "out_of_range")) {
    expect_true(col %in% names(scored), info = col)
  }
  expect_true(any(!is.na(scored$rqr)))
  expect_type(scored$flag_low, "logical")
})
