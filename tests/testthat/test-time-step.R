# The values a day needs and the lags of level 9 are given in time, not in
# steps: the same call must work on 10-, 15- and 30-minute data.

step_net <- function(days, st, step_min, start = "2024-01-01", seed = 1) {
  set.seed(seed)
  n <- days * 1440 / step_min
  ti <- seq(as.POSIXct(start, tz = "UTC"), by = step_min * 60, length.out = n)
  h <- as.numeric(ti - ti[1], units = "hours")
  X <- sapply(1:st, function(j) 12 + 8 * sin(2 * pi * (h / 24 - 110) / 365) + 6 * sin(2 * pi * (h - 9) / 24) +
                4 * sin(2 * pi * h / 216) + 0.3 * j + rnorm(n, 0, 0.15))
  colnames(X) <- paste0("S", 1:st)
  xts::xts(X, order.by = ti)
}

test_that("levels 5 and 8 and RH5 judge 30-minute data with their defaults", {
  x <- step_net(730, 4, 30)
  d <- format(zoo::index(x), "%Y-%m-%d")
  # S4 indoors for 20 days
  ind <- d >= "2025-03-01" & d <= "2025-03-20"
  x[ind, "S4"] <- 21 + rnorm(sum(ind), 0, 0.1)
  i5 <- T_QC_5_climatic_outliers(x, verbose = FALSE)$qc_info$t5_climatic_outliers
  expect_equal(i5$min_obs_per_day, 24)
  expect_true(all(i5$n_judged_by_station > 0))
  r8 <- T_QC_8_diurnal_range(x, verbose = FALSE)
  expect_equal(r8$qc_info$t8_diurnal_range$min_obs_per_day, 24)
  expect_true(all(zoo::coredata(r8$qc_data_flagged)[ind, "S4"] == 8))
  q5 <- RH_QC_5_decoupling(130 - 4 * x, temperature = x, verbose = FALSE)$qc_info$rh5_decoupling
  expect_equal(q5$decouple_min_n, 20)
  expect_true(all(q5$n_judged_by_station > 0))
})

test_that("level 9 finds an hour on 15-minute data and keeps its numbers on 10-minute data", {
  x <- step_net(28, 6, 15)
  n <- nrow(x)
  # S3 one hour early: four steps
  v <- as.numeric(zoo::coredata(x[, "S3"]))
  x[, "S3"] <- c(v[5:n], rep(NA, 4))
  i <- T_QC_9_clock_offset(x, verbose = FALSE)$qc_info$t9_clock_offset
  expect_equal(c(i$min_lag, i$max_lag), c(4, 8))
  expect_gt(i$n_flagged_by_station[["S3"]], 0)
  expect_equal(sum(i$n_flagged_by_station[-3]), 0)
  expect_equal(i$episodes$lag_min, 60)
  # the 10-minute defaults are what they were
  x10 <- step_net(14, 4, 10)
  i10 <- T_QC_9_clock_offset(x10, verbose = FALSE)$qc_info$t9_clock_offset
  expect_equal(c(i10$min_lag, i10$max_lag), c(6, 12))
  expect_equal(T_QC_8_diurnal_range(x10, verbose = FALSE)$qc_info$t8_diurnal_range$min_obs_per_day, 72)
  expect_equal(T_QC_5_climatic_outliers(x10, verbose = FALSE)$qc_info$t5_climatic_outliers$min_obs_per_day, 72)
  expect_equal(RH_QC_5_decoupling(130 - 4 * x10, temperature = x10, verbose = FALSE)$qc_info$rh5_decoupling$decouple_min_n, 60)
})
