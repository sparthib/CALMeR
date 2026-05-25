# Minimal features table used across tests. clean_features() filters on
# aln_type and complete cases by default, so each row has both.
make_features <- function() {
  data.frame(
    aln_type    = c("P", "P", "P", "P", "P", "P", "P", "P", "S"),
    MAPQ        = c(60,  60,  60,  60,  60,  60,  60,  60,  0),
    n_isoforms  = c(1,   2,   5,   7,   8,   12,  3,   4,   2),
    tx_length   = c(100, 500, 1e3, 2e3, 5e3, 1e4, 8e3, 1.5e3, 700),
    gc_content  = c(20,  35,  50,  55,  60,  65,  70,  40,   45),
    pct_covered = c(0.5, 5,   50,  80,  95,  98,  99,  100,  60),
    stringsAsFactors = FALSE
  )
}

test_that("clean_features defaults: lumps n_isoforms at 7, logs tx_length, winsorizes pct_covered, no gc/length clamps", {
  feat  <- make_features()
  clean <- clean_features(feat)

  # primary-only filter dropped the one secondary row.
  expect_equal(nrow(clean), sum(feat$aln_type == "P"))

  # n_isoforms lumped at 7 (default cap), no other clamps.
  expect_lte(max(clean$n_isoforms), 7)
  expect_true(any(clean$n_isoforms == 7))   # the 8 and 12 rows pulled down to 7

  # tx_length log-transformed in place. Smallest input was 100 -> log(100).
  expect_equal(min(clean$tx_length), log(100))
  expect_lt(max(clean$tx_length), log(1e4) + 1e-9)

  # gc_content untouched: no default clamp.
  primary_gc <- feat$gc_content[feat$aln_type == "P"]
  expect_equal(sort(clean$gc_content), sort(primary_gc))

  # pct_covered winsorized to data-driven 1st/99th and exposed via attr.
  bounds <- attr(clean, "pct_covered_winsor_bounds")
  expect_length(bounds, 2L)
  expect_true(bounds[1] <= bounds[2])
  expect_gte(min(clean$pct_covered), bounds[1])
  expect_lte(max(clean$pct_covered), bounds[2])
})

test_that("clean_features uses caller-supplied winsorize bounds without recomputing", {
  feat <- make_features()

  fit_bounds <- c(5, 95)
  clean_test <- clean_features(feat, pct_covered_winsor_bounds = fit_bounds)

  # Bounds are used as-given and round-tripped on the result.
  expect_equal(attr(clean_test, "pct_covered_winsor_bounds"), fit_bounds)
  expect_gte(min(clean_test$pct_covered), 5)
  expect_lte(max(clean_test$pct_covered), 95)
})

test_that("clean_features can disable the log transform", {
  feat  <- make_features()
  clean <- clean_features(feat, log_tx_length = FALSE)

  primary_tx <- feat$tx_length[feat$aln_type == "P"]
  # n_isoforms cap pulls no rows out, so the set is the same.
  expect_equal(sort(clean$tx_length), sort(primary_tx))
})

test_that("clean_features errors clearly when log_tx_length = TRUE meets a non-positive value", {
  feat <- make_features()
  feat$tx_length[1] <- 0
  expect_error(clean_features(feat, log_tx_length = TRUE), "tx_length")
})

