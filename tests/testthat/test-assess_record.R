
test_that("assess_record detects gaps and classifies segments correctly", {

  # Synthetic: 2 years of hourly data with a 100-day gap in the middle
  P_annual <- 365.25 * 86400
  t1 <- seq(0, by = 3600, length.out = 200 * 24)  # 200 days
  t2 <- seq(300 * 86400, by = 3600, length.out = 500 * 24)  # 500 days starting day 300
  t_all <- c(t1, t2)
  vals <- sin(2 * pi * t_all / P_annual) + 10
  z <- zoo::zoo(vals, order.by = t_all)

 a <- assess_record(z, period = P_annual)

  # Should be a record_assessment
  expect_s3_class(a, "record_assessment")

  # Should detect 1 gap
  expect_equal(nrow(a$gaps), 1)

  # Gap should be ~100 days = 100*86400 seconds
  expect_gt(a$gaps$duration_s[1], 90 * 86400)
  expect_lt(a$gaps$duration_s[1], 110 * 86400)

  # Gap is ~27% of annual period — should exceed default 17% tolerance
  expect_false(a$gaps$tolerable[1])

  # Should have 2 segments
  expect_equal(nrow(a$segments), 2)

  # First segment ~200 days (0.55 cycles) — insufficient
  expect_equal(a$segments$verdict[1], "insufficient")

  # Second segment ~500 days (1.37 cycles) — marginal (>= 1.0 but < 2.0)
  expect_equal(a$segments$verdict[2], "marginal")

  # Overall verdict should be marginal
  expect_equal(a$overall_verdict, "marginal")
})


test_that("assess_record passes clean multi-year record", {

  P_annual <- 365.25 * 86400
  t_clean <- seq(0, by = 3600, length.out = 3 * 365 * 24)  # 3 years hourly
  vals <- sin(2 * pi * t_clean / P_annual) + 10
  z <- zoo::zoo(vals, order.by = t_clean)

  a <- assess_record(z, period = P_annual, min_cycles = 2.0)

  expect_equal(nrow(a$gaps), 0)
  expect_equal(nrow(a$segments), 1)
  expect_equal(a$segments$verdict[1], "sufficient")
  expect_equal(a$overall_verdict, "sufficient")
})


test_that("assess_record works for diel period", {

  # 30 days of hourly data, no gaps
  t_diel <- seq(0, by = 3600, length.out = 30 * 24)
  vals <- sin(2 * pi * t_diel / 86400) + 15
  z <- zoo::zoo(vals, order.by = t_diel)

  a <- assess_record(z, period = 86400, min_cycles = 5.0)

  expect_equal(a$overall_verdict, "sufficient")
  expect_gt(a$segments$cycles[1], 29)  # ~30 diel cycles
})


test_that("assess_record rejects record shorter than min_cycles", {

  P_annual <- 365.25 * 86400
  # Only 100 days — well under 2 annual cycles
  t_short <- seq(0, by = 3600, length.out = 100 * 24)
  vals <- rnorm(length(t_short))
  z <- zoo::zoo(vals, order.by = t_short)

  a <- assess_record(z, period = P_annual, min_cycles = 2.0)

  expect_equal(a$overall_verdict, "insufficient")
})


test_that("assess_record respects custom max_gap_frac", {

  P_annual <- 365.25 * 86400

  # 2 years with a 80-day gap
  t1 <- seq(0, by = 3600, length.out = 300 * 24)
  t2 <- seq(380 * 86400, by = 3600, length.out = 400 * 24)
  z <- zoo::zoo(c(rnorm(length(t1)), rnorm(length(t2))), order.by = c(t1, t2))

  # Default 0.17 (62 days) — gap of 80 days should NOT be tolerable
  a1 <- assess_record(z, period = P_annual, max_gap_frac = 0.17)
  expect_false(a1$gaps$tolerable[1])

  # Relaxed 0.25 (91 days) — gap of 80 days SHOULD be tolerable
  a2 <- assess_record(z, period = P_annual, max_gap_frac = 0.25)
  expect_true(a2$gaps$tolerable[1])
})



# === Diel gap tolerance tests =================================================

test_that("assess_record flags intolerable diel gap (10-hour gap)", {

 # 14 days of hourly data with a 10-hour gap at day 7.
  # 10 hours = 42% of diel period, well above the 17% default tolerance.
  # The gap should be intolerable, splitting the record into two segments,
  # but each segment (~7 days) should be independently sufficient.

  P_diel <- 86400
  t1 <- seq(0, by = 3600, length.out = 7 * 24)               # days 0-7
  t2 <- seq(7 * 86400 + 10 * 3600, by = 3600, length.out = 7 * 24)  # resumes 10h later
  vals <- sin(2 * pi * c(t1, t2) / P_diel) + 15
  z <- zoo::zoo(vals, order.by = c(t1, t2))

  a <- assess_record(z, period = P_diel, min_cycles = 5.0)

  # Should detect 1 gap
  expect_equal(nrow(a$gaps), 1)

  # Gap should be ~10-11 hours (depends on hourly stepping)
  expect_gt(a$gaps$duration_s[1], 9 * 3600)
  expect_lt(a$gaps$duration_s[1], 12 * 3600)

  # Gap should exceed 17% diel tolerance
  expect_gt(a$gaps$frac_of_period[1], 0.17)
  expect_false(a$gaps$tolerable[1])

  # Should split into 2 segments, each ~7 days = ~7 cycles
  expect_equal(nrow(a$segments), 2)
  expect_gt(a$segments$cycles[1], 6)
  expect_gt(a$segments$cycles[2], 6)

  # Both segments sufficient (>5 diel cycles each)
  expect_equal(a$segments$verdict[1], "sufficient")
  expect_equal(a$segments$verdict[2], "sufficient")
  expect_equal(a$overall_verdict, "sufficient")
})


test_that("assess_record tolerates short diel gap (3-hour gap)", {

  # 14 days of hourly data with a 3-hour gap at day 7.
  # 3 hours = 12.5% of diel period, below the 17% default tolerance.
  # The gap should be tolerable — no record splitting.

  P_diel <- 86400
  t1 <- seq(0, by = 3600, length.out = 7 * 24)                # days 0-7
  t2 <- seq(7 * 86400 + 3 * 3600, by = 3600, length.out = 7 * 24)  # resumes 3h later
  vals <- sin(2 * pi * c(t1, t2) / P_diel) + 15
  z <- zoo::zoo(vals, order.by = c(t1, t2))

  a <- assess_record(z, period = P_diel, min_cycles = 5.0)

  # Should still detect the gap (> 3x median dt)
  expect_equal(nrow(a$gaps), 1)

  # Gap should be tolerable (< 17% of diel period)
  expect_lt(a$gaps$frac_of_period[1], 0.17)
  expect_true(a$gaps$tolerable[1])

  # Record splits at the gap structurally, but gap is tolerable
  # Both segments are sufficient regardless
  expect_equal(a$overall_verdict, "sufficient")
})


# === min_segment_duration tests ===============================================

test_that("min_segment_duration downgrades short segments that pass cycle count", {

  # 4 days of hourly data — passes min_cycles = 2.0 for diel (4 cycles > 2)
  # but fails a 7-day min_segment_duration requirement.

  P_diel <- 86400
  t <- seq(0, by = 3600, length.out = 4 * 24)
  vals <- sin(2 * pi * t / P_diel) + 15
  z <- zoo::zoo(vals, order.by = t)

  # Without min_segment_duration: sufficient (4 cycles > 2)
  a1 <- assess_record(z, period = P_diel, min_cycles = 2.0)
  expect_equal(a1$overall_verdict, "sufficient")

  # With min_segment_duration = 7 days: marginal (4 days < 7, but >= 1 cycle)
  a2 <- assess_record(z, period = P_diel, min_cycles = 2.0,
                       min_segment_duration = 7 * 86400)
  expect_equal(a2$segments$verdict[1], "marginal")
  expect_equal(a2$overall_verdict, "marginal")
})


test_that("min_segment_duration does not affect already-long segments", {

  # 30 days of clean hourly data — easily passes both min_cycles and
  # min_segment_duration.

  P_diel <- 86400
  t <- seq(0, by = 3600, length.out = 30 * 24)
  vals <- sin(2 * pi * t / P_diel) + 15
  z <- zoo::zoo(vals, order.by = t)

  a <- assess_record(z, period = P_diel, min_cycles = 5.0,
                      min_segment_duration = 7 * 86400)

  expect_equal(a$overall_verdict, "sufficient")
  expect_gt(a$segments$cycles[1], 29)
  expect_equal(a$params$min_segment_duration, 7 * 86400)
})


test_that("min_segment_duration is stored in params and NULL by default", {

  P_diel <- 86400
  t <- seq(0, by = 3600, length.out = 10 * 24)
  z <- zoo::zoo(sin(2 * pi * t / P_diel), order.by = t)

  # Default: NULL
  a1 <- assess_record(z, period = P_diel)
  expect_null(a1$params$min_segment_duration)

  # Explicit value: stored
  a2 <- assess_record(z, period = P_diel, min_segment_duration = 5 * 86400)
  expect_equal(a2$params$min_segment_duration, 5 * 86400)
})
