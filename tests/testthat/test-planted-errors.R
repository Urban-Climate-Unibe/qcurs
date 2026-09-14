test_that("every planted error reaches its intended detector", {
  ch  <- biel_chain()
  FLG <- zoo::coredata(ch$res$qc_data_flagged)
  flg <- ch$res$qc_data_flagged

  # Codes 51/52/53 = the level's own code plus the interpolation base (50):
  # T_QC_9 refills the single-value gaps these levels left behind.
  expect_true(all(FLG[planted_cells("E1", flg)] == 51))   # fault code -45
  expect_true(all(FLG[planted_cells("E2", flg)] == 52))   # 55 degC in summer
  expect_true(all(FLG[planted_cells("E3", flg)] == 53))   # isolated +15 K spike
  expect_true(all(FLG[planted_cells("E4", flg)] ==  4))   # 40-point constant block

  # E5 is spatial: only judged where both compatible neighbours report a value.
  e5 <- planted_cells("E5", flg)
  X  <- zoo::coredata(ch$x)
  judged <- !is.na(X[e5[, 1], "Log_B2"]) & !is.na(X[e5[, 1], "Log_B6"])
  expect_gte(sum(judged), 15)
  expect_true(all(FLG[e5[judged, , drop = FALSE]] %in% c(6, 56)))

  # E6 is the indoor episode at a station level 6 cannot see.
  expect_gte(mean(FLG[planted_cells("E6", flg)] == 8), 0.95)
})

test_that("the real Biel minimum survives the whole chain", {
  ch <- biel_chain()
  expect_equal(unname(zoo::coredata(ch$res$qc_data_flagged)[planted_cells("N1", ch$res$qc_data_flagged)]), 0)
})

test_that("the chain leaves the unplanted record untouched", {
  ch  <- biel_chain()
  FLG <- zoo::coredata(ch$res$qc_data_flagged)
  lin <- unlist(lapply(c("E1", "E2", "E3", "E4", "E5", "E6"), function(id) {
    cc <- planted_cells(id, ch$res$qc_data_flagged); cc[, 1] + (cc[, 2] - 1) * nrow(FLG) }))
  rest <- FLG[-lin]
  expect_lt(sum(rest %in% 1:8, na.rm = TRUE), 0.0015 * sum(!is.na(rest)))
})
