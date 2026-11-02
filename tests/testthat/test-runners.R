runner_fixture <- function() {
  set.seed(21)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 10 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  tcyc <- 18 + 5 * sin(2 * pi * (mins - 540) / 1440)
  TT <- sapply(1:4, function(k) tcyc + rnorm(length(ti), 0, 0.3) + (k - 2) * 0.4)
  RH <- 90 - 2.5 * (TT - 13) + matrix(rnorm(length(ti) * 4, 0, 1.5), ncol = 4)
  colnames(TT) <- colnames(RH) <- paste0("S", 1:4)
  TT[100, 1] <- -45                                   # a gross temperature error
  # midday, where the baseline is low enough that +30 stays below 100 -
  # planted at night it would exceed the range and be caught by level 2 instead
  RH[360, 2] <- RH[360, 2] + 30
  RH[500:502, 3] <- NA                                # a short plain gap
  md <- data.frame(ID = paste0("S", 1:4),
                   LON = 7.44 + (0:3) * 0.004, LAT = 46.95 + (0:3) * 0.002,
                   Landuse = "Sealed Areas", stringsAsFactors = FALSE)
  list(tt = xts::xts(TT, order.by = ti), rh = xts::xts(RH, order.by = ti), md = md)
}

test_that("run_qc_temperature equals the manual chain and records skips", {
  f <- runner_fixture()
  r_auto <- run_qc_temperature(f$tt, metadata = f$md, verbose = FALSE)
  r_man  <- T_QC_1_gross_error(f$tt, verbose = FALSE)
  r_man  <- T_QC_2_out_of_range(r_man, verbose = FALSE)
  r_man  <- T_QC_3_time_consistency(r_man, verbose = FALSE)
  r_man  <- T_QC_4_stuck_values(r_man, verbose = FALSE)
  r_man  <- T_QC_5_climatic_outliers(r_man, verbose = FALSE)
  r_man  <- T_QC_6_spatial_consistency(r_man, metadata = f$md, verbose = FALSE)
  r_man  <- T_QC_7_spatiotemporal_consistency(r_man, metadata = f$md, verbose = FALSE)
  r_man  <- T_QC_8_diurnal_range(r_man, verbose = FALSE)
  r_man  <- T_QC_9_clock_offset(r_man, verbose = FALSE)
  expect_identical(zoo::coredata(r_auto$qc_data_flagged), zoo::coredata(r_man$qc_data_flagged))
  expect_identical(r_auto$qc_info$run_qc_temperature$levels_skipped, "t10")
  expect_identical(r_auto$qc_info$run_qc_temperature$levels, 1:10)
  # without metadata the spatial levels are skipped and say so in the record
  r_no <- run_qc_temperature(f$tt, verbose = FALSE)
  expect_true(all(c("t6", "t7") %in% r_no$qc_info$run_qc_temperature$levels_skipped))
  expect_false("t6_spatial_consistency" %in% names(r_no$qc_info))
})

test_that("params overrides reach their level", {
  f <- runner_fixture()
  b <- list(winter = list(min_val = -30, max_val = 20), spring = list(min_val = -25, max_val = 35),
            summer = list(min_val = -10, max_val = 42), autumn = list(min_val = -15, max_val = 38))
  r <- run_qc_temperature(f$tt, params = list(t2 = list(season_thresholds = b)), verbose = FALSE)
  expect_equal(r$qc_info$t2_out_of_range$season_thresholds$summer$min_val, -10)
  expect_error(run_qc_temperature(f$tt, params = list(nope = list())), "params keys")
})

test_that("levels selects or excludes levels, by number or by name", {
  f <- runner_fixture()
  # exclude one: its record is absent, it is listed as skipped, the rest runs
  r <- run_qc_temperature(f$tt, metadata = f$md, levels = -9, verbose = FALSE)
  expect_false("t9_clock_offset" %in% names(r$qc_info))
  expect_true(all(c("t8_diurnal_range", "t6_spatial_consistency") %in% names(r$qc_info)))
  expect_identical(r$qc_info$run_qc_temperature$levels, c(1:8, 10L))
  expect_true(all(c("t9", "t10") %in% r$qc_info$run_qc_temperature$levels_skipped))
  # select a few: only those run, in the chain order
  r2 <- run_qc_temperature(f$tt, metadata = f$md, levels = c(6, 1:3), verbose = FALSE)
  expect_equal(grep("^t[0-9]+_", names(r2$qc_info), value = TRUE),
               c("t1_gross_error", "t2_out_of_range", "t3_time_consistency", "t6_spatial_consistency"))
  expect_identical(r2$qc_info$run_qc_temperature$levels, c(1L, 2L, 3L, 6L))
  # names work too, in either case and with the sign
  r3 <- run_qc_temperature(f$tt, levels = c("T1", "t4"), verbose = FALSE)
  expect_equal(grep("^t[0-9]+_", names(r3$qc_info), value = TRUE), c("t1_gross_error", "t4_stuck_values"))
  r4 <- run_qc_temperature(f$tt, levels = "-T9", verbose = FALSE)
  expect_false("t9_clock_offset" %in% names(r4$qc_info))
  # level 10 still needs interpolate, level 6 still needs metadata
  r5 <- run_qc_temperature(f$tt, levels = c(1, 6, 10), verbose = FALSE)
  expect_equal(grep("^t[0-9]+_", names(r5$qc_info), value = TRUE), "t1_gross_error")
  expect_true(all(c("t6", "t10") %in% r5$qc_info$run_qc_temperature$levels_skipped))
  # bad input is refused
  expect_error(run_qc_temperature(f$tt, levels = c(1, -2)), "mix")
  expect_error(run_qc_temperature(f$tt, levels = 11), "between 1 and 10")
  expect_error(run_qc_temperature(f$tt, levels = "T9x"), "names like")
  expect_error(run_qc_temperature(f$tt, levels = 0), "between 1 and 10")
  # the humidity runner does the same
  t_res <- run_qc_temperature(f$tt, verbose = FALSE)
  rh <- run_qc_humidity(f$rh, temperature = t_res, levels = -3, verbose = FALSE)
  expect_false("rh3_spike" %in% names(rh$qc_info))
  expect_true("rh3" %in% rh$qc_info$run_qc_humidity$levels_skipped)
  expect_identical(rh$qc_info$run_qc_humidity$levels, c(1L, 2L, 4:8))
})

test_that("run_qc_humidity wires the temperature result into levels 1, 5 and 7", {
  f <- runner_fixture()
  t_res <- run_qc_temperature(f$tt, metadata = f$md, verbose = FALSE)
  r <- run_qc_humidity(f$rh, temperature = t_res, metadata = f$md,
                       interpolate = TRUE, verbose = FALSE)
  # the inherited gross error, the spike, and the filled gap all carry codes
  FLG <- zoo::coredata(r$qc_data_flagged)
  expect_equal(unname(FLG[100, 1]), 51)               # inherited (1), then refilled by 8
  expect_equal(unname(FLG[360, 2]), 53)               # spike (3), then refilled by 8
  expect_true(all(FLG[500:502, 3] == 50))             # plain gap filled
  expect_gt(r$qc_info$rh8_interpolate$n_refilled, 0)
  # the combined record documents BOTH chains
  expect_true(all(c("dataset_temperature", "dataset_humidity",
                    "t1_gross_error", "rh1_inherit_temperature") %in% names(r$qc_info)))
  expect_identical(r$qc_info$run_qc_humidity$levels_skipped, character(0))
})

test_that("run_qc_humidity skips audibly with a bare series or nothing", {
  f <- runner_fixture()
  # bare series: inheritance impossible, decoupling and dewpoint still run
  r1 <- run_qc_humidity(f$rh, temperature = f$tt, metadata = f$md, verbose = FALSE)
  expect_identical(r1$qc_info$run_qc_humidity$levels_skipped, "rh8")
  expect_false(r1$qc_info$rh1_inherit_temperature$inherited)      # ran, inherited nothing
  expect_true("rh5_decoupling" %in% names(r1$qc_info))
  # nothing: levels 1, 6, 7 skip themselves, the rest runs
  r0 <- run_qc_humidity(f$rh, verbose = FALSE)
  expect_identical(r0$qc_info$run_qc_humidity$levels_skipped, "rh8")
  expect_true(r0$qc_info$rh5_decoupling$skipped)                    # ran, skipped itself
  expect_true(r0$qc_info$rh7_dewpoint_consistency$skipped)          # ran, skipped itself
  expect_false(r0$qc_info$rh1_inherit_temperature$inherited)
  expect_true("rh6_saturation_drift" %in% names(r0$qc_info))
})
