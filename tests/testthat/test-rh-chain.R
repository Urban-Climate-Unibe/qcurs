# Synthetic four-station network with a healthy diurnal cycle, one planted
# error per RH level, run through the whole humidity chain.

rh_fixture <- function() {
  set.seed(11)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 20 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  tcyc <- 18 + 5 * sin(2 * pi * (mins - 540) / 1440)             # diurnal temperature
  TT <- sapply(1:4, function(k) tcyc + rnorm(length(ti), 0, 0.3) + (k - 2) * 0.4)
  # healthy negative coupling, with nights around 89: the network must NOT
  # brush the 95 percent saturation criterion outside the fog day, because one
  # healthy saturated step per segment acquits a drifter by construction
  RH <- 90 - 2.5 * (TT - 13) + matrix(rnorm(length(ti) * 4, 0, 1.5), ncol = 4)
  colnames(TT) <- colnames(RH) <- paste0("S", 1:4)
  day <- as.Date(ti)

  # E2: an impossible value and a tolerance-zone value at S1
  RH[10, 1] <- 130;  RH[20, 1] <- 103
  # E3: isolated +30 spike at S2
  RH[500, 2] <- RH[500, 2] + 30
  # E4: fog day 8 (whole network saturated) with S3 drifting at ~80
  fog <- which(day == as.Date("2025-06-08"))
  # on a real fog day the TEMPERATURE is flat too - otherwise RH5 would rightly
  # call the pinned humidity "decoupled"; the flat T triggers its trange skip
  TT[fog, ] <- 12 + matrix(rnorm(length(fog) * 4, 0, 0.3), ncol = 4)
  # the healthy stations sit NEAR-CONSTANT at 97: sd below sd_tol, so without
  # the saturation exception RH4 would flag the whole fog day
  RH[fog, ] <- 97 + matrix(rnorm(length(fog) * 4, 0, 0.05), ncol = 4)
  RH[fog, 3] <- 80 + rnorm(length(fog), 0, 0.5)
  # E5: stuck block at S4 (60.0 for 40 steps, well below saturation)
  RH[2000:2039, 4] <- 60.0
  # E6: decoupled days 14-16 at S2 (constant-ish, uncorrelated with its T);
  # three days, because a single day is a front, not a sensor (RH5 needs two)
  dec <- which(day >= as.Date("2025-06-14") & day <= as.Date("2025-06-16"))
  RH[dec, 2] <- 55 + rnorm(length(dec), 0, 1.5)
  # E7: dewpoint offset at S1, day 18, midday: a DRY bias of -30 percentage
  # points. A wet bias saturates against the 100 percent ceiling and barely
  # moves the dewpoint; a dry bias shifts it by several Kelvin.
  off <- which(day == as.Date("2025-06-18"))[61:96]
  RH[off, 1] <- pmax(RH[off, 1] - 30, 20)

  md <- data.frame(ID = paste0("S", 1:4),
                   LON = 7.44 + (0:3) * 0.004, LAT = 46.95 + (0:3) * 0.002,
                   Landuse = "Sealed Areas", stringsAsFactors = FALSE)
  list(rh = xts::xts(RH, order.by = ti), tt = xts::xts(TT, order.by = ti),
       md = md, fog = fog, dec = dec, off = off)
}

test_that("the humidity chain finds every planted error class", {
  f <- rh_fixture()
  r <- RH_QC_2_range(f$rh, verbose = FALSE)
  r <- RH_QC_3_spike(r, verbose = FALSE)
  r <- RH_QC_4_stuck_values(r, verbose = FALSE)
  r <- RH_QC_5_decoupling(r, temperature = f$tt, verbose = FALSE)
  r <- RH_QC_6_saturation_drift(r, verbose = FALSE)
  r <- RH_QC_7_dewpoint_consistency(r, temperature = f$tt, metadata = f$md, verbose = FALSE)
  FLG <- zoo::coredata(r$qc_data_flagged)

  expect_equal(unname(FLG[10, 1]), 2)                        # impossible value flagged
  expect_equal(unname(zoo::coredata(r$qc_data)[20, 1]), 103) # tolerance value kept AS MEASURED
  expect_equal(unname(FLG[20, 1]), 0)                        # ... and unflagged: QC never alters data
  expect_equal(unname(FLG[500, 2]), 3)                       # spike
  expect_gte(mean(FLG[f$fog, 3] == 6), 0.95)                 # saturation drifter
  expect_true(all(FLG[2000:2039, 4] == 4))                   # stuck block
  expect_gte(mean(FLG[f$dec, 2] == 5), 0.95)                 # decoupled day
  expect_gte(mean(FLG[f$off, 1] == 7), 0.8)                  # dewpoint offset
  # the healthy remainder stays essentially untouched: total flags stay close
  # to the planted budget (fog drift + stuck + decoupled day + offset + slack)
  expect_lt(sum(FLG %in% 2:7, na.rm = TRUE),
            length(f$fog) + 40 + length(f$dec) + length(f$off) + 20)
  # every level left its record
  expect_true(all(c("dataset_humidity", "rh2_range", "rh3_spike", "rh4_stuck_values",
                    "rh5_decoupling", "rh6_saturation_drift", "rh7_dewpoint_consistency")
                  %in% names(r$qc_info)))
})

test_that("the saturation exception keeps fog days unflagged in RH4", {
  f <- rh_fixture()
  r5 <- RH_QC_4_stuck_values(f$rh, verbose = FALSE)
  # the fog day is near-constant at 97 but sits ABOVE sat_max: exempt
  expect_true(all(zoo::coredata(r5$qc_data_flagged)[f$fog, c(1, 2, 4)] == 0))
})

test_that("RH5 and RH7 accept the temperature chain result and refuse misalignment", {
  f <- rh_fixture()
  t_res <- T_QC_1_gross_error(f$tt, verbose = FALSE)
  r <- RH_QC_5_decoupling(f$rh, temperature = t_res, verbose = FALSE)
  expect_true("rh5_decoupling" %in% names(r$qc_info))
  # RH5 matches by stamp and name like RH7: a shorter temperature series is
  # judged where it overlaps, and a foreign one skips the level
  r6 <- RH_QC_5_decoupling(f$rh, temperature = f$tt[1:(5 * 144), ], verbose = FALSE)
  expect_equal(r6$qc_info$rh5_decoupling$n_time_with_temperature, 5 * 144)
  expect_true(all(r6$qc_info$rh5_decoupling$n_judged_by_station <= 5 * 144))
  foreign <- xts::xts(zoo::coredata(f$tt), order.by = zoo::index(f$tt) + 3600 * 24 * 400)
  expect_true(RH_QC_5_decoupling(f$rh, temperature = foreign, verbose = FALSE)$qc_info$rh5_decoupling$skipped)
  expect_warning(
    r7 <- RH_QC_7_dewpoint_consistency(f$rh, temperature = f$tt, metadata = f$md[1:2, ], verbose = FALSE),
    "No metadata for")
  expect_true(r7$qc_info$rh7_dewpoint_consistency$skipped)   # <3 stations: honest skip
  # RH7 matches the temperature by stamp and name: a shorter, reordered
  # temperature series still yields dewpoints where it overlaps
  tt2 <- f$tt[1:(10 * 144), c("S4", "S3", "S2", "S1")]
  r7b <- RH_QC_7_dewpoint_consistency(f$rh, temperature = tt2, metadata = f$md, verbose = FALSE)
  expect_false(isTRUE(r7b$qc_info$rh7_dewpoint_consistency$skipped))
  expect_equal(r7b$qc_info$rh7_dewpoint_consistency$n_time_with_temperature, 10 * 144)
  expect_true(all(zoo::coredata(r7b$qc_data_flagged)[-(1:(10 * 144)), ] != 7, na.rm = TRUE))
  # no temperature at all: skip, data untouched
  r7c <- RH_QC_7_dewpoint_consistency(f$rh, metadata = f$md, verbose = FALSE)
  expect_true(r7c$qc_info$rh7_dewpoint_consistency$skipped)
  expect_identical(zoo::coredata(r7c$qc_data), zoo::coredata(f$rh))
})

test_that("RH5 tells a front from a sensor: network-wide days and single days are not objected to", {
  f <- rh_fixture()
  RH <- zoo::coredata(f$rh); day <- as.Date(zoo::index(f$rh))
  # a front on day 10: every station looks decoupled at once
  front <- which(day == as.Date("2025-06-10"))
  RH[front, ] <- 60 + matrix(rnorm(length(front) * 4, 0, 1.5), ncol = 4)
  # a single decoupled day at S3, no neighbour affected
  lone <- which(day == as.Date("2025-06-12"))
  RH[lone, 3] <- 55 + rnorm(length(lone), 0, 1.5)
  rh <- xts::xts(RH, order.by = zoo::index(f$rh))
  r <- RH_QC_5_decoupling(rh, temperature = f$tt, verbose = FALSE)
  FLG <- zoo::coredata(r$qc_data_flagged)
  expect_true(all(FLG[front, ] == 0))                        # weather day: nobody judged
  expect_equal(r$qc_info$rh5_decoupling$days_unjudged, "2025-06-10")
  expect_true(all(FLG[lone, 3] == 0))                        # one day alone: not enough
  expect_gte(mean(FLG[f$dec, 2] == 5), 0.95)                 # three days alone: objected to
  # the coverage count excludes the weather day
  expect_equal(unname(r$qc_info$rh5_decoupling$n_judged_by_station["S3"]) %% 144, 0)
  # both rules off: the plain daily test objects to the front and the single day
  r0 <- RH_QC_5_decoupling(rh, temperature = f$tt, min_consecutive_days = 1,
                           max_network_frac = 1, verbose = FALSE)
  F0 <- zoo::coredata(r0$qc_data_flagged)
  expect_true(all(F0[front, ] == 5))
  expect_true(all(F0[lone, 3] == 5))
  expect_length(r0$qc_info$rh5_decoupling$days_unjudged, 0)
})

test_that("RH1 inherits only where time stamp AND logger match, and skips otherwise", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 6)
  rh <- xts::xts(cbind(Log_1 = 55:60, Log_2 = 70:75, Log_3 = 40:45), order.by = ti)
  tf <- xts::xts(cbind(Log_2 = c(0, 0, 4, 0, 0, 0), Log_1 = c(0, 1, 0, 0, 0, 0)), order.by = ti)
  # fewer loggers, other column order: matched by name, Log_3 untouched
  r <- RH_QC_1_inherit_temperature(rh, tf, verbose = FALSE)
  f <- zoo::coredata(r$qc_data_flagged)
  expect_equal(unname(f[2, "Log_1"]), 1); expect_equal(unname(f[3, "Log_2"]), 1)
  expect_true(all(f[, "Log_3"] == 0))
  expect_equal(r$qc_info$rh1_inherit_temperature$loggers_unmatched, "Log_3")
  # partial time overlap: matched by exact stamp, not by row position
  tf2 <- xts::xts(cbind(Log_1 = c(1, 0, 0, 0)), order.by = ti[3:6])
  r2 <- RH_QC_1_inherit_temperature(rh, tf2, verbose = FALSE)
  expect_equal(unname(zoo::coredata(r2$qc_data_flagged)[3, "Log_1"]), 1)
  expect_equal(r2$qc_info$rh1_inherit_temperature$n_time_matched, 4)
  # unusable input: skip with reason, data untouched, chain continues
  for (bad in list(NULL, zoo::coredata(tf),                                   # nothing / no names+index
                   xts::xts(cbind(Log_1 = 18:23), order.by = ti),            # values, not flags
                   xts::xts(cbind(Log_9 = rep(0, 6)), order.by = ti),        # no common logger
                   xts::xts(zoo::coredata(tf), order.by = ti + 86400))) {    # no common stamp
    rb <- RH_QC_1_inherit_temperature(rh, bad, verbose = FALSE)
    expect_false(rb$qc_info$rh1_inherit_temperature$inherited)
    expect_identical(zoo::coredata(rb$qc_data), zoo::coredata(rh))
  }
})

test_that("RH6 judges per segment and acquits a station that reaches saturation once", {
  set.seed(1)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 40 * 144)
  RH <- matrix(70 + rnorm(length(ti) * 5, 0, 3), ncol = 5, dimnames = list(NULL, paste0("S", 1:5)))
  day <- as.Date(ti)
  fog1 <- day == as.Date("2025-06-08"); fog2 <- day == as.Date("2025-07-05")   # segment 0 and 1
  RH[fog1 | fog2, ] <- 97
  RH[fog1 | fog2, 3] <- 80                              # S3 drifts on both fog days
  RH[fog2, 4] <- 80; RH[which(fog2)[1:5], 4] <- 96      # S4 drifts but touches saturation 5 times
  r <- RH_QC_6_saturation_drift(xts::xts(RH, order.by = ti), verbose = FALSE)
  n <- r$qc_info$rh6_saturation_drift$n_flagged_by_station
  expect_equal(unname(n[["S3"]]), 2 * 144)              # both segments, evidence steps only
  expect_equal(unname(n[["S4"]]), 0)                    # acquitted for the whole segment
  expect_equal(sum(n[c("S1", "S2", "S5")]), 0)
})

test_that("RH4 exempts a plateau at the sensor's own ceiling, but not a stuck block below it", {
  set.seed(3)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 20 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  RH <- 70 - 20 * sin(2 * pi * (mins - 540) / 1440) + matrix(rnorm(length(ti) * 4, 0, 1), ncol = 4)
  colnames(RH) <- paste0("S", 1:4); day <- as.Date(ti)
  # S1 reads dry at the top: a low-cost element whose ceiling is 82, not 100
  RH[, 1] <- RH[, 1] * 0.9
  # a fog day: the healthy stations sit flat at 97, S1 flat at its ceiling
  fog <- which(day == as.Date("2025-06-08"))
  RH[fog, ] <- 97 + matrix(rnorm(length(fog) * 4, 0, 0.05), ncol = 4)
  RH[fog, 1] <- 82 + rnorm(length(fog), 0, 0.05)
  # and a genuinely stuck block at S1, well below its ceiling
  RH[1500:1539, 1] <- 60
  x <- xts::xts(RH, order.by = ti)
  r <- RH_QC_4_stuck_values(x, verbose = FALSE)
  FLG <- zoo::coredata(r$qc_data_flagged); i <- r$qc_info$rh4_stuck_values
  expect_true(all(FLG[fog, ] == 0))                       # the plateau at the own ceiling is saturation
  expect_true(all(FLG[1500:1539, 1] == 4))                # the stuck block is not
  expect_lt(i$sat_ceiling[["S1"]], 82); expect_gt(i$sat_ceiling[["S2"]], 87)   # each at its own top
  # without the station rule the ceiling plateau is a stuck window again
  r0 <- RH_QC_4_stuck_values(x, sat_quantile = NULL, verbose = FALSE)
  expect_gte(mean(zoo::coredata(r0$qc_data_flagged)[fog, 1] == 4), 0.9)
  expect_true(all(r0$qc_info$rh4_stuck_values$sat_ceiling == 95))
  expect_error(RH_QC_4_stuck_values(x, sat_quantile = 1), "sat_quantile")
})

test_that("RH6 skips itself audibly when the network never saturates", {
  set.seed(4)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 40 * 144)
  RH <- matrix(70 + rnorm(length(ti) * 5, 0, 3), ncol = 5, dimnames = list(NULL, paste0("S", 1:5)))
  expect_message(r <- RH_QC_6_saturation_drift(xts::xts(RH, order.by = ti)), "never reached network_min")
  i <- r$qc_info$rh6_saturation_drift
  expect_true(i$skipped); expect_equal(i$n_saturated_steps, 0L); expect_lt(i$network_max, 95)
  expect_identical(zoo::coredata(r$qc_data), RH)
  # the level leaves no row in the summary, like every level that could not judge
  expect_false("rh6_saturation_drift" %in% qc_summary(r)$level)
})
