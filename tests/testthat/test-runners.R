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
  r_man  <- T_QC_4_temporal_persistence(r_man, verbose = FALSE)
  r_man  <- T_QC_5_climatic_outliers(r_man, verbose = FALSE)
  r_man  <- T_QC_6_spatial_consistency(r_man, metadata = f$md, verbose = FALSE)
  r_man  <- T_QC_7_spatiotemporal_consistency(r_man, metadata = f$md, verbose = FALSE)
  r_man  <- T_QC_8_diurnal_range(r_man, verbose = FALSE)
  expect_identical(zoo::coredata(r_auto$qc_data_flagged), zoo::coredata(r_man$qc_data_flagged))
  expect_identical(r_auto$qc_info$run_qc_temperature$levels_skipped, "t9")
  # without metadata the spatial levels are skipped and say so in the record
  r_no <- run_qc_temperature(f$tt, verbose = FALSE)
  expect_true(all(c("t6", "t7") %in% r_no$qc_info$run_qc_temperature$levels_skipped))
  expect_false("t6_spatial_consistency" %in% names(r_no$qc_info))
})

test_that("params overrides reach their level", {
  f <- runner_fixture()
  r <- run_qc_temperature(f$tt, params = list(t2 = list(summer_min = -10)), verbose = FALSE)
  expect_equal(r$qc_info$t2_out_of_range$season_thresholds$summer$min_val, -10)
  expect_error(run_qc_temperature(f$tt, params = list(nope = list())), "params keys")
})

test_that("run_qc_humidity wires the temperature result into levels 1, 6 and 7", {
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
  expect_identical(r1$qc_info$run_qc_humidity$levels_skipped, c("rh1", "rh8"))
  expect_true("rh6_decoupling" %in% names(r1$qc_info))
  # nothing: levels 1, 6, 7 skipped, the rest runs
  r0 <- run_qc_humidity(f$rh, verbose = FALSE)
  expect_true(all(c("rh1", "rh6", "rh7", "rh8") %in% r0$qc_info$run_qc_humidity$levels_skipped))
  expect_true("rh4_saturation_drift" %in% names(r0$qc_info))
})
