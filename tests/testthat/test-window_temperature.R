test_that("window_temperature drops windows spanning temporal gaps", {

  # 21 days of hourly data with a 3-day gap at days 10-13.
  # 7-day windows that overlap the gap should be dropped because
  # their temporal coverage falls below 90%.
  # Windows fully within days 0-10 or days 13-21 should survive.

  t1 <- seq(0, by = 3600, length.out = 10 * 24)
  t2 <- seq(13 * 86400, by = 3600, length.out = 8 * 24)
  vals <- sin(2 * pi * c(t1, t2) / 86400) + 15
  z <- zoo::zoo(vals, order.by = c(t1, t2))

  wins <- window_temperature(z, width = 7 * 86400, step = 86400, min_coverage = 0.9)

  # Windows starting at days 0-3 and day 13 should pass (fully within clean data)
  # Windows starting at days 4-12 should be dropped (span the gap)
  expect_equal(nrow(wins), 5)

  # All retained windows should start before or after the gap
  starts_days <- wins$window_start / 86400
  expect_true(all(starts_days <= 3 | starts_days >= 13))
})


test_that("window_temperature retains all windows for gap-free record", {

  # 14 days clean hourly data, 7-day windows with 1-day step
  t <- seq(0, by = 3600, length.out = 14 * 24)
  z <- zoo::zoo(sin(2 * pi * t / 86400) + 15, order.by = t)

  wins <- window_temperature(z, width = 7 * 86400, step = 86400, min_coverage = 0.9)

  # Should produce windows starting at days 0 through 7
  expect_equal(nrow(wins), 7)
})


test_that("window_temperature drops windows with NA values", {

  # 14 days clean hourly, but inject NAs in days 3-4 (48 of 168 obs in any
  # 7-day window = 29% missing → below 90%)
  t <- seq(0, by = 3600, length.out = 14 * 24)
  vals <- sin(2 * pi * t / 86400) + 15
  # Set days 3-4 to NA (indices 73:120)
  vals[73:120] <- NA
  z <- zoo::zoo(vals, order.by = t)

  wins <- window_temperature(z, width = 7 * 86400, step = 86400, min_coverage = 0.9)

  # Windows spanning days 3-4 should be dropped
  # Windows starting at day 0: covers 0-7, has 48 NAs → 120/168 = 71% → fail
  # Windows starting at day 5+: no NAs → pass
  expect_lt(nrow(wins), 7)
  # Windows after the NA region should survive
  expect_true(any(wins$window_start >= 5 * 86400))
})


test_that("window_temperature detects gap-as-missing-index vs NA-at-index", {

  # Compare: same effective gap, represented two different ways:
  # A) Missing index entries (the bug case)
  # B) NA values at expected index entries
  # Both should be detected with the temporal coverage check.

  # Scenario: 10 days hourly, 2-day gap at days 4-6

  # (A) Missing index entries
  t1 <- seq(0, by = 3600, length.out = 4 * 24)
  t2 <- seq(6 * 86400, by = 3600, length.out = 4 * 24)
  z_gap <- zoo::zoo(sin(2 * pi * c(t1, t2) / 86400), order.by = c(t1, t2))

  # (B) NAs at expected positions
  t_full <- seq(0, by = 3600, length.out = 10 * 24)
  vals_na <- sin(2 * pi * t_full / 86400)
  vals_na[(4 * 24 + 1):(6 * 24)] <- NA
  z_na <- zoo::zoo(vals_na, order.by = t_full)

  wins_gap <- window_temperature(z_gap, width = 7 * 86400, step = 86400, min_coverage = 0.9)
  wins_na  <- window_temperature(z_na,  width = 7 * 86400, step = 86400, min_coverage = 0.9)

  # Both representations should drop gap-spanning windows
  expect_equal(nrow(wins_gap), nrow(wins_na))
})

