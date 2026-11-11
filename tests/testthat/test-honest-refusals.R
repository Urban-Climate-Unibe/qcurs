test_that("level 5 refuses a single-season campaign instead of inventing a climatology", {
  i5 <- biel_chain()$res$qc_info$t5_climatic_outliers
  expect_equal(i5$n_flagged, 0)
  expect_gte(length(i5$refused), 21)          # 7 stations x at least 3 months
})

test_that("level 7 skips every Biel station and says so", {
  i7 <- biel_chain()$res$qc_info$t7_spatiotemporal_consistency
  expect_equal(i7$n_flagged, 0)
  expect_length(i7$never_evaluated, 7)
  # by name, with the neighbours each one found kept apart
  expect_setequal(i7$never_evaluated, colnames(biel_chain()$x))
  expect_true(all(i7$n_neighbours < 5))
})

test_that("level 5 lists a station-month that reported but never had a full day", {
  # 40 days in one month and the next, 30 values a day: no day is dense enough
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 40 * 144)
  v <- 18 + rnorm(length(ti), 0, 1)
  v[as.numeric(format(ti, "%H")) >= 5] <- NA
  x <- xts::xts(matrix(v, ncol = 1, dimnames = list(NULL, "S1")), order.by = ti)
  i5 <- T_QC_5_climatic_outliers(x, verbose = FALSE)$qc_info$t5_climatic_outliers
  expect_equal(i5$refused, c("S1/06(days=0,y=0)", "S1/07(days=0,y=0)"))
})

test_that("level 6 names the stations it could never evaluate", {
  i6 <- biel_chain()$res$qc_info$t6_spatial_consistency
  expect_true(all(c("Log_B1", "Log_B7") %in% i6$never_evaluated))
})

test_that("level 10 reports what it added and never bridges a long block", {
  ch  <- biel_chain()
  i10 <- ch$res$qc_info$t10_interpolate
  FLG <- zoo::coredata(ch$res$qc_data_flagged)
  expect_gt(i10$n_gap_filled, 0)
  expect_equal(sum(FLG == 50, na.rm = TRUE), i10$n_gap_filled)
  expect_false(any(FLG[planted_cells("E4", ch$res$qc_data_flagged)] == 54))
})
