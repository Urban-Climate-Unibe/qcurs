# Levels 5 and 7 legitimately refuse on the Biel campaign, so their positive
# detection paths are proven on minimal synthetic records instead.

test_that("level 5 catches a climatic outlier on a two-year record", {
  set.seed(7)
  ti <- seq(as.POSIXct("2023-01-01", tz = "UTC"), by = "1 hour", length.out = 2 * 365 * 24)
  doy <- as.numeric(format(ti, "%j")); hod <- as.numeric(format(ti, "%H"))
  v <- 10 + 10 * sin(2 * pi * (doy - 100) / 365) + 4 * sin(2 * pi * (hod - 9) / 24) +
       stats::rnorm(length(ti), 0, 1.5)
  v[5000] <- v[5000] + 40
  x <- xts::xts(matrix(v, ncol = 1, dimnames = list(NULL, "S1")), order.by = ti)
  f <- zoo::coredata(T_QC_5_climatic_outliers(x, verbose = FALSE)$qc_data_flagged)
  expect_equal(unname(f[5000, 1]), 5)
  expect_lte(sum(f == 5, na.rm = TRUE), 6)
})

test_that("level 7 catches a spatiotemporal spike on a dense network", {
  set.seed(8)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 30 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  common <- 18 + 6 * sin(2 * pi * (mins - 540) / 1440)
  X <- sapply(1:6, function(k) common + stats::rnorm(length(ti), 0, 0.4) + (k - 3) * 0.3)
  colnames(X) <- paste0("S", 1:6)
  X[2000, 3] <- X[2000, 3] + 10
  md <- data.frame(ID = paste0("S", 1:6), LON = 7.44 + (0:5) * 0.003, LAT = 46.95 + (0:5) * 0.001)
  f <- zoo::coredata(T_QC_7_spatiotemporal_consistency(
         xts::xts(X, order.by = ti), metadata = md, verbose = FALSE)$qc_data_flagged)
  expect_equal(unname(f[2000, 3]), 7)
  expect_equal(sum(f == 7, na.rm = TRUE), 1)
})
