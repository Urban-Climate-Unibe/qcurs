# Synthetic six-station network with a diurnal cycle and a synoptic swing:
# one station shifts by +2.5 K from day 40 on (T9), one runs an hour early
# from day 8 on (T10). Healthy stations, a one-day front, a week on which
# the network diverges and a site delay of 20 minutes must stay unobjected.
# The site offsets are not in station order, so that the station with the
# site delay (S3) is never the median of the others.

shift_fixture <- function(days = 90, seed = 21) {
  set.seed(seed)
  ti <- seq(as.POSIXct("2025-04-01", tz = "UTC"), by = "10 min", length.out = days * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  t_num <- as.numeric(ti - ti[1], units = "days")
  # diurnal cycle plus a slow synoptic swing, so daily means move with the weather
  base <- 12 + 5 * sin(2 * pi * (mins - 540) / 1440) + 4 * sin(2 * pi * t_num / 11) + 2 * sin(2 * pi * t_num / 29)
  off <- c(-0.4, 0, 1.6, 0.4, 0.8, 1.2)
  TT <- sapply(1:6, function(k) base + rnorm(length(ti), 0, 0.3) + off[k])
  colnames(TT) <- paste0("S", 1:6)
  list(x = xts::xts(TT, order.by = ti), ti = ti, day = as.Date(ti), base = base)
}

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
  expect_length(i$refused, 6)
  expect_true(all(is.na(i$baseline)))
  expect_true(all(i$n_judged_by_station == 0))
  r2 <- T_QC_9_bias_shift(f$x[, 1:2], verbose = FALSE)
  expect_true(r2$qc_info$t9_bias_shift$skipped)
  expect_error(T_QC_9_bias_shift(f$x, bias_tol = 0), "bias_tol")
  expect_error(T_QC_9_bias_shift(f$x, ref_days = 2), "ref_days")
})

test_that("T10 finds a clock an hour behind and leaves a 20-minute site delay alone", {
  f <- shift_fixture(days = 35)
  X <- zoo::coredata(f$x)
  d <- as.integer(f$day - f$day[1]) + 1
  n <- nrow(X)
  # S2's stamps run an hour early from day 8 on: at stamp t it shows the
  # value of t + 6 steps (its clock is an hour behind)
  late <- which(d >= 8)
  src <- pmin(late + 6, n)
  X[late, 2] <- f$base[src] + rnorm(length(late), 0, 0.3)
  # S3 answers 20 minutes late throughout: a site, not a clock
  X[, 3] <- c(rep(NA, 2), f$base[1:(n - 2)]) + rnorm(n, 0, 0.3) + 1.6
  x <- xts::xts(X, order.by = f$ti)
  r <- T_QC_10_clock_offset(x, verbose = FALSE)
  i <- r$qc_info$t10_clock_offset
  FLG <- zoo::coredata(r$qc_data_flagged)
  # the lag table: +6 at S2 from window 2 on, about -2 at S3, near 0 elsewhere.
  # In a six-station network the hour at S2 leaks into everybody's yardstick
  # (the median of five), so the healthy lags scatter by a step or so
  expect_true(all(round(i$lags[2:5, "S2"]) == 6))
  expect_true(all(abs(i$lags[, c("S1", "S4", "S5", "S6")]) <= 2, na.rm = TRUE))
  expect_true(all(round(i$lags[, "S3"]) %in% -4:-1, na.rm = TRUE))
  # every value of S2 in the offset windows is objected to, nothing before day 8
  expect_true(all(FLG[d >= 8, 2] == 10))
  expect_true(all(FLG[d < 8, 2] == 0))
  expect_true(all(FLG[, -2] == 0, na.rm = TRUE))
  # the episode carries the lag in minutes
  expect_equal(nrow(i$episodes), 1)
  expect_equal(i$episodes$station, "S2")
  expect_equal(i$episodes$lag_min, 60, tolerance = 0.1)
  expect_equal(i$episodes$from, "2025-04-08")
  expect_gt(i$episodes$r, 0.9)
  # one offset window alone is not enough
  X3 <- zoo::coredata(f$x)
  w2 <- which(d %in% 8:14)
  X3[w2, 4] <- f$base[pmin(w2 + 6, n)] + rnorm(length(w2), 0, 0.3) + 0.4
  r3 <- T_QC_10_clock_offset(xts::xts(X3, order.by = f$ti), verbose = FALSE)
  expect_gte(r3$qc_info$t10_clock_offset$lags[2, "S4"], 4)
  expect_equal(r3$qc_info$t10_clock_offset$n_flagged, 0)
})

test_that("T10 skips a tiny network and refuses bad parameters", {
  f <- shift_fixture(days = 10)
  r <- T_QC_10_clock_offset(f$x[, 1:2], verbose = FALSE)
  expect_true(r$qc_info$t10_clock_offset$skipped)
  expect_error(T_QC_10_clock_offset(f$x, min_lag = 20), "min_lag")
  expect_error(T_QC_10_clock_offset(f$x, min_lag = 0), "min_lag")
  expect_error(T_QC_10_clock_offset(f$x, window_days = 0.5), "window_days")
  expect_error(T_QC_10_clock_offset(f$x, min_cor = 1), "min_cor")
  # a clean network: every window judged, no lag beyond one step, nothing flagged
  r2 <- T_QC_10_clock_offset(f$x, verbose = FALSE)
  expect_equal(r2$qc_info$t10_clock_offset$n_flagged, 0)
  expect_true(all(abs(r2$qc_info$t10_clock_offset$lags[1, ]) <= 1.5))
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
