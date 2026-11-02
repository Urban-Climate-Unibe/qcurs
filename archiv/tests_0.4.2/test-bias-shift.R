# Tests of T_QC_9_bias_shift (0.4.0-0.4.2), archived with the level. They use
# shift_fixture() from test-clock-offset.R.

test_that("T9 finds a +2.5 K shift from day 40 on and leaves fronts and weather days alone", {
  f <- shift_fixture()
  X <- zoo::coredata(f$x)
  d <- as.integer(f$day - f$day[1]) + 1
  # the planted shift: S2 reads 2.5 K high from day 40 on
  X[d >= 40, 2] <- X[d >= 40, 2] + 2.5
  # a one-day front at S3 (day 20): +2.5 K for a single day
  X[d == 20, 3] <- X[d == 20, 3] + 2.5
  # a week on which the network diverges (days 60-66): S1 +3 K, S5 -3 K; with the
  # shifted S2 that is three of six stations off their baselines -> weather, not judged
  X[d %in% 60:66, 1] <- X[d %in% 60:66, 1] + 3
  X[d %in% 60:66, 5] <- X[d %in% 60:66, 5] - 3
  x <- xts::xts(X, order.by = f$ti)
  r <- T_QC_9_bias_shift(x, verbose = FALSE)
  i <- r$qc_info$t9_bias_shift
  FLG <- zoo::coredata(r$qc_data_flagged)
  # every value of S2 from day 40 on is objected to, nothing before; the
  # weather days in between are not judged and stay unflagged
  expect_true(all(FLG[d >= 40 & !(d %in% 60:66), 2] == 9))
  expect_true(all(FLG[d < 40 | d %in% 60:66, 2] == 0))
  # the other stations stay clean: the front is too short, the weather days are unjudged
  expect_true(all(FLG[, -2] == 0))
  expect_equal(i$days_unjudged, format(as.Date("2025-04-01") + 59:65))
  # the episode table names the shift, running across the unjudged days
  expect_equal(nrow(i$episodes), 1)
  expect_equal(i$episodes$station, "S2")
  expect_equal(i$episodes$from, "2025-05-10")
  expect_equal(i$episodes$to, "2025-06-29")
  expect_equal(i$episodes$days, 51 - 7)
  expect_equal(i$episodes$shift_K, 2.5, tolerance = 0.1)
  # the baselines carry the site offsets (against the median of the others), not the shift
  expect_equal(unname(i$baseline["S2"]), 0 - 0.8, tolerance = 0.15)
  expect_equal(unname(i$baseline["S4"]), 0.4 - 0.8, tolerance = 0.15)
  # a station that read high from day one has it in the baseline and is not objected to
  X2 <- zoo::coredata(f$x); X2[, 4] <- X2[, 4] + 3
  r2 <- T_QC_9_bias_shift(xts::xts(X2, order.by = f$ti), verbose = FALSE)
  expect_equal(r2$qc_info$t9_bias_shift$n_flagged, 0)
  expect_equal(unname(r2$qc_info$t9_bias_shift$baseline["S4"]), 3.4 - 0.8, tolerance = 0.2)
})

test_that("T9 refuses a record shorter than its baseline and skips a tiny network", {
  f <- shift_fixture(days = 20)
  r <- T_QC_9_bias_shift(f$x, verbose = FALSE)
  i <- r$qc_info$t9_bias_shift
  expect_equal(i$n_flagged, 0)
  expect_length(i$never_evaluated, 6)
  expect_true(all(is.na(i$baseline)))
  expect_true(all(i$n_judged_by_station == 0))
  r2 <- T_QC_9_bias_shift(f$x[, 1:2], verbose = FALSE)
  expect_true(r2$qc_info$t9_bias_shift$skipped)
  expect_error(T_QC_9_bias_shift(f$x, bias_tol = 0), "bias_tol")
  expect_error(T_QC_9_bias_shift(f$x, ref_days = 2), "ref_days")
})

test_that("T9 and T10 run clean on the Biel network and leave their records", {
  ch <- biel_chain()
  i9  <- ch$res$qc_info$t9_bias_shift
  i10 <- ch$res$qc_info$t10_clock_offset
  expect_equal(i9$n_flagged, 0)
  expect_equal(i10$n_flagged, 0)
  expect_true(all(!is.na(i9$baseline)))
  expect_true(all(i9$n_judged_by_station > 0))
  expect_true(all(i10$n_judged_by_station > 0))
  expect_equal(nrow(i9$episodes), 0)
  expect_equal(nrow(i10$episodes), 0)
})
