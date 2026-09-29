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

test_that("qc_find_spikes flags an isolated spike, keeps a front, refuses thin context", {
  base <- c(18.0, 18.2, 18.1, 18.3, 18.4, 18.2, 18.5, 18.3, 18.1, 18.0)
  spike <- base; spike[5] <- 30
  sp <- qc_find_spikes(spike, dt = 3, threshold = 6)
  expect_equal(which(sp$hit), 5L)
  expect_equal(which(sp$judged), 2:9)                    # ends never judged
  front <- c(18.0, 18.2, 18.1, 18.3, 30.0, 30.2, 30.1, 30.3, 30.0, 30.2)
  expect_false(any(qc_find_spikes(front, dt = 3, threshold = 6)$hit))
  thin <- c(18.0, 18.2, NA, NA, 30, NA, NA, 18.5, 18.3, 18.1)
  th <- qc_find_spikes(thin, dt = 3, threshold = 6)
  expect_false(any(th$hit)); expect_false(th$judged[5])   # refused, not acquitted
  # the two levels reach the same verdict through the same code
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 10)
  x <- xts::xts(matrix(spike, ncol = 1, dimnames = list(NULL, "S1")), order.by = ti)
  t3  <- T_QC_3_time_consistency(x, verbose = FALSE)
  rh3 <- RH_QC_3_spike(x, threshold = 6, verbose = FALSE)
  expect_equal(unname(zoo::coredata(t3$qc_data_flagged)[5, 1]), 3)
  expect_equal(unname(zoo::coredata(rh3$qc_data_flagged)[5, 1]), 3)
})

test_that("qc_window_points converts durations and refuses what cannot work", {
  expect_equal(qc_window_points("6 hours", 600), 37L)
  expect_equal(qc_window_points("6 h", 600), 37L)
  expect_equal(qc_window_points(as.difftime(90, units = "mins"), 600), 10L)
  expect_equal(qc_window_points("2 days", 600), 289L)
  expect_warning(w <- qc_window_points("45 mins", 600), "using 4 steps")
  expect_equal(w, 5L)
  expect_error(qc_window_points("3 weeks!", 600), "should be one of")
  expect_error(qc_window_points("6h", 600), "must look like")
  expect_error(qc_window_points("10 mins", 600), "need at least 2")
})

test_that("qc_find_stuck flags a constant block, respects the guard and the exemption", {
  set.seed(4)
  v <- 20 + sin(seq_len(200) / 15) + rnorm(200, 0, 0.05)
  v[80:125] <- 24                                        # 46 stuck points, window 37
  st <- qc_find_stuck(v, width = 37, sd_tol = 0, min_valid = 19)
  expect_equal(which(st$hit), 80:125)
  expect_true(all(st$judged))                            # a full series: every point judged
  # a block shorter than the window is invisible
  v2 <- v; v2[80:125] <- 20 + sin(80:125 / 15); v2[100:130] <- 24
  expect_false(any(qc_find_stuck(v2, 37, 0, 19)$hit))
  # two identical survivors in an otherwise empty window do not count
  v3 <- v; v3[60:100] <- NA; v3[101:102] <- 21.5
  s3 <- qc_find_stuck(v3, 37, 0, 19)
  expect_false(any(s3$hit))
  expect_false(any(s3$judged[60:100]))                   # gaps are never "judged"
  # the saturation exemption: the same constant block at 97 is not judged at all
  v4 <- v; v4[80:125] <- 97
  expect_true(any(qc_find_stuck(v4, 37, 0, 19)$hit))
  s4 <- qc_find_stuck(v4, 37, 0, 19, exempt_above = 95)
  expect_false(any(s4$hit)); expect_false(s4$judged[100])
  # a series with fewer values than one window IS judged when a window still
  # holds min_valid of them (the guard counts min_valid, not width)
  v5 <- c(rep(20, 25), rep(NA, 40))
  expect_equal(which(qc_find_stuck(v5, 37, 0, 19)$hit), 1:25)
})

test_that("qc_fill_gaps fills short gaps, codes refills by level, leaves long gaps and ends", {
  X <- matrix(c(1, NA, 3, NA, NA, NA, NA, NA, NA, 10, 11, NA), ncol = 1, dimnames = list(NULL, "S1"))
  F <- matrix(c(0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA), ncol = 1, dimnames = list(NULL, "S1"))
  r <- qc_fill_gaps(X, F, maxgap = 5, refill_flagged = TRUE)
  expect_equal(unname(r$X[2, 1]), 2); expect_equal(unname(r$previous_flag[2, 1]), 53)   # spike refilled
  expect_true(all(is.na(r$X[4:9, 1])))                                  # 6-gap > maxgap
  expect_true(is.na(r$X[12, 1]))                                        # end never extrapolated
  expect_equal(c(r$n_gap, r$n_ref), c(0, 1))
  r2 <- qc_fill_gaps(X, F, maxgap = 5, refill_flagged = FALSE)
  expect_true(is.na(r2$X[2, 1])); expect_equal(unname(r2$previous_flag[2, 1]), 3)  # protected
})
