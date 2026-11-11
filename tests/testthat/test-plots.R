# The four plot functions run on the Biel chain into a null device and hand
# back the numbers they drew.

test_that("the plot functions draw the Biel result and return their tables", {
  ch <- biel_chain()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  tl <- qc_plot_timeline(ch$res)
  expect_named(tl, c("day", "level", "n"))
  expect_equal(sum(tl$n), sum(zoo::coredata(ch$res$qc_data_flagged) %in% c(1:9, 51:59)))
  tl1 <- qc_plot_timeline(ch$res, station = "Log_B1")
  expect_equal(sum(tl1$n[tl1$level == "t8_diurnal_range"]), 720)
  lv <- qc_plot_levels(ch$res)
  expect_named(lv, c("level", "flagged", "judged", "share"))
  expect_equal(lv$flagged[lv$level == "t4_stuck_values"], 40)
  st <- qc_plot_levels(ch$res, by = "station")
  expect_equal(dim(st), c(6, 7))                                   # six levels objected (1,2,3,4,6,8), seven stations
  expect_gt(st["t8_diurnal_range", "Log_B1"], 1)
  cal <- qc_plot_calendar(ch$res)
  expect_equal(dim(cal), c(94, 7))
  expect_equal(unname(cal["2025-08-16", "Log_B1"]), 100)           # an indoor day: every value objected to
  expect_true(is.na(cal["2025-05-31", "Log_B7"]) || cal["2025-05-31", "Log_B7"] >= 0)
  cal6 <- qc_plot_calendar(ch$res, level = "t6_spatial_consistency")
  expect_equal(sum(cal6[, "Log_B1"], na.rm = TRUE), 0)             # never reached by level 6
  expect_error(qc_plot_calendar(ch$res, level = "nope"), "level must be")
  sp <- qc_plot_station(ch$res, "Log_B4", raw = ch$x, from = "2025-08-05/2025-08-07")
  expect_named(sp, c("time", "value", "level"))
  expect_true(all(sp$level == "t6_spatial_consistency"))
  expect_false(anyNA(sp$value))
  sp2 <- qc_plot_station(ch$res, "Log_B4", from = "2025-08-05/2025-08-07")   # no raw: ticks, values NA
  expect_true(all(is.na(sp2$value)))
  expect_error(qc_plot_station(ch$res, "nope"), "not in the result")
  expect_error(qc_plot_timeline(ch$x), "chain result")
})

test_that("the plot functions survive a clean result, an emptied station and a lone `to`", {
  set.seed(9)
  n <- 5 * 144
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = n)
  h <- as.numeric(ti - ti[1], units = "hours")
  X <- sapply(1:4, function(j) 18 + 6 * sin(2 * pi * (h - 9) / 24) + 0.3 * j + rnorm(n, 0, 0.15))
  colnames(X) <- paste0("S", 1:4)
  x <- xts::xts(X, order.by = ti)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  clean <- run_qc_temperature(x, verbose = FALSE)
  expect_equal(sum(zoo::coredata(clean$qc_data_flagged) > 0), 0)
  expect_equal(dim(qc_plot_levels(clean, by = "station")), c(0, 4))
  # a result only the interpolation went through has no level to show
  expect_equal(nrow(qc_plot_levels(T_QC_10_interpolate(x, verbose = FALSE))), 0)
  # a station that only ever sent its fault code: every value objected to
  xb <- x
  xb[, "S2"] <- 999
  res <- run_qc_temperature(xb, verbose = FALSE)
  sp <- qc_plot_station(res, "S2")
  expect_equal(nrow(sp), n)
  # `to` without `from` cuts the window
  sp_to <- qc_plot_station(res, "S2", to = "2025-06-02")
  expect_equal(nrow(sp_to), 2 * 144)
})
