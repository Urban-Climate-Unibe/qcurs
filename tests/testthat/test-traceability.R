# Every level leaves the same three counters: what it flagged, per station,
# and how many cells it actually judged per station. A cell with flag 0 is
# only "checked" by the levels whose coverage reached it.

test_that("every level records its coverage, and flags never exceed it", {
  ch <- biel_chain()
  X0 <- zoo::coredata(ch$x)
  n_valid <- colSums(!is.na(X0))
  keys <- c("t1_gross_error", "t2_out_of_range", "t3_time_consistency", "t4_stuck_values",
            "t5_climatic_outliers", "t6_spatial_consistency", "t7_spatiotemporal_consistency",
            "t8_diurnal_range")
  for (k in keys) {
    rec <- ch$res$qc_info[[k]]
    expect_true(all(c("n_flagged_by_station", "n_judged_by_station") %in% names(rec)), info = k)
    expect_true(all(rec$n_flagged_by_station <= rec$n_judged_by_station), info = k)
    expect_true(all(rec$n_judged_by_station <= n_valid[names(rec$n_judged_by_station)]), info = k)
  }
  # the first two levels see every value; the refusing levels see none
  expect_equal(unname(ch$res$qc_info$t1_gross_error$n_judged_by_station), unname(n_valid))
  expect_true(all(ch$res$qc_info$t5_climatic_outliers$n_judged_by_station == 0))
  expect_true(all(ch$res$qc_info$t7_spatiotemporal_consistency$n_judged_by_station == 0))
  # level 6 reached exactly the stations it does not list as never evaluated
  i6 <- ch$res$qc_info$t6_spatial_consistency
  expect_identical(names(i6$n_judged_by_station)[i6$n_judged_by_station == 0], i6$never_evaluated)
})

test_that("a second run of the same level warns and keeps the first flags", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 20)
  v  <- rep(20, 20); v[10] <- 35
  x  <- xts::xts(matrix(v, ncol = 1, dimnames = list(NULL, "S1")), order.by = ti)
  r1 <- T_QC_3_time_consistency(x, verbose = FALSE)
  expect_warning(r2 <- T_QC_3_time_consistency(r1, verbose = FALSE), "already ran")
  expect_equal(unname(zoo::coredata(r2$qc_data_flagged)[10, 1]), 3)
})

test_that("level 6 skips a too-small network instead of stopping the chain", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 144)
  x  <- xts::xts(matrix(20 + rnorm(288, 0, 0.2), ncol = 2, dimnames = list(NULL, c("A", "B"))), order.by = ti)
  md <- data.frame(ID = c("A", "B"), LON = c(7.44, 7.45), LAT = c(46.95, 46.95), Landuse = "Sealed Areas")
  r  <- T_QC_6_spatial_consistency(x, metadata = md, verbose = FALSE)
  expect_true(r$qc_info$t6_spatial_consistency$skipped)
  expect_identical(zoo::coredata(r$qc_data), zoo::coredata(x))
  # the runner therefore gets through with metadata for a two-station network
  expect_true("t8_diurnal_range" %in% names(run_qc_temperature(x, metadata = md, verbose = FALSE)$qc_info))
})

test_that("RH1 records which temperature codes it inherited", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 6)
  rh <- xts::xts(cbind(Log_1 = 55:60), order.by = ti)
  tf <- xts::xts(cbind(Log_1 = c(0, 1, 4, 4, 50, 0)), order.by = ti)
  r  <- RH_QC_1_inherit_temperature(rh, tf, verbose = FALSE)
  tab <- r$qc_info$rh1_inherit_temperature$n_flagged_by_t_code
  expect_equal(unname(tab[c("1", "4")]), c(1L, 2L))
  expect_true(is.integer(tab))
  expect_false("50" %in% names(tab))                                # gap fill not inherited
  expect_equal(unname(r$qc_info$rh1_inherit_temperature$n_judged_by_station), 6L)
})

test_that("levels 6 and 7 survive metadata that match no station at all", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 144)
  x  <- xts::xts(matrix(20 + rnorm(432, 0, 0.2), ncol = 3, dimnames = list(NULL, c("A", "B", "C"))), order.by = ti)
  md <- data.frame(ID = c("Q1", "Q2"), LON = c(7.44, 7.45), LAT = c(46.95, 46.95), Landuse = "Sealed Areas")
  expect_warning(r6 <- T_QC_6_spatial_consistency(x, metadata = md, verbose = FALSE), "No metadata for: A, B, C")
  expect_true(r6$qc_info$t6_spatial_consistency$skipped)
  expect_warning(r7 <- T_QC_7_spatiotemporal_consistency(x, metadata = md, verbose = FALSE), "No metadata for: A, B, C")
  expect_equal(r7$qc_info$t7_spatiotemporal_consistency$skipped, paste0(c("A", "B", "C"), "(no metadata)"))
  expect_true(all(r7$qc_info$t7_spatiotemporal_consistency$n_judged_by_station == 0))
  # the humidity twin of level 6 skips the same way
  expect_warning(r_rh7 <- RH_QC_7_dewpoint_consistency(x, temperature = x, metadata = md, verbose = FALSE), "No metadata")
  expect_true(r_rh7$qc_info$rh7_dewpoint_consistency$skipped)
})
