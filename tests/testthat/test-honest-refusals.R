test_that("level 5 refuses a single-season campaign instead of inventing a climatology", {
  i5 <- biel_chain()$res$qc_info$t5_climatic_outliers
  expect_equal(i5$n_flagged, 0)
  expect_gte(length(i5$refused), 21)          # 7 stations x at least 3 months
})

test_that("level 7 skips every Biel station and says so", {
  i7 <- biel_chain()$res$qc_info$t7_spatiotemporal
  expect_equal(i7$n_flagged, 0)
  expect_length(i7$skipped, 7)
})

test_that("level 6 names the stations it could never evaluate", {
  i6 <- biel_chain()$res$qc_info$t6_spatial_consistency
  expect_true(all(c("Log_B1", "Log_B7") %in% i6$never_evaluated))
})

test_that("level 9 reports what it added and never bridges a long block", {
  ch  <- biel_chain()
  i9  <- ch$res$qc_info$t9_interpolate
  FLG <- zoo::coredata(ch$res$qc_data_flagged)
  expect_gt(i9$n_gap_filled, 0)
  expect_equal(sum(FLG == 50, na.rm = TRUE), i9$n_gap_filled)
  expect_false(any(FLG[planted_cells("E4", ch$res$qc_data_flagged)] == 54))
})
