test_that("qc_as_xts converts the Biel table without losing an observation", {
  wide <- biel_wide()
  x <- qc_as_xts(wide, time_col = "time", verbose = FALSE)
  expect_s3_class(x, "xts")
  expect_equal(dim(x), c(13308L, 7L))
  # round trip against the tidy source: every (time, station) value survives
  tidy <- biel_tidy()
  back <- data.frame(time    = rep(format(zoo::index(x), "%Y-%m-%d %H:%M:%S"), times = ncol(x)),
                     station = rep(colnames(x), each = nrow(x)),
                     temperature = as.vector(zoo::coredata(x)), stringsAsFactors = FALSE)
  ki <- paste(tidy$time, tidy$station); ko <- paste(back$time, back$station)
  expect_setequal(ki, ko)
  expect_equal(tidy$temperature[order(ki)], back$temperature[order(ko)])
})

test_that("qc_as_xts takes the time column by exact name or by number", {
  tt <- c("2025-06-01 00:10:00", "2025-06-01 00:20:00")
  d  <- data.frame(time = tt, Log_1 = c(1, 2), Log_2 = c(3, 4))
  expect_equal(dim(qc_as_xts(d, time_col = "time", verbose = FALSE)), c(2L, 2L))
  expect_equal(dim(qc_as_xts(d, time_col = 1, verbose = FALSE)), c(2L, 2L))
  # a POSIXct column is taken as it is
  p <- d; p$time <- as.POSIXct(tt, tz = "UTC")
  expect_equal(dim(qc_as_xts(p, time_col = "time", verbose = FALSE)), c(2L, 2L))
  # an entirely empty logger column (read.csv makes it logical) is accepted
  e <- d; e$Log_9 <- NA
  expect_equal(ncol(qc_as_xts(e, time_col = "time", verbose = FALSE)), 3L)
})

test_that("qc_as_xts refuses everything it would otherwise have to guess", {
  tt <- c("2025-06-01 00:10:00", "2025-06-01 00:20:00")
  d  <- data.frame(time = tt, Log_1 = c(1, 2), Log_2 = c(3, 4))
  expect_error(qc_as_xts(d), "time_col is required")
  expect_error(qc_as_xts(d, time_col = "Time"), "not found")        # exact, case-sensitive
  expect_error(qc_as_xts(d, time_col = 9), "column number")
  # a POSIXct time column with a gap is refused as well, not only text stamps
  p <- d; p$time <- as.POSIXct(tt, tz = "UTC"); p$time[2] <- NA
  expect_error(qc_as_xts(p, time_col = "time"), "missing")
  emp <- d; names(emp)[3] <- ""
  expect_error(qc_as_xts(emp, time_col = "time"), "non-empty")
  expect_error(qc_as_xts(d, time_col = c(1, 2)), "ONE column")
  expect_error(qc_as_xts("file.csv", time_col = "time"), "data frame")
  dup <- d; names(dup) <- c("time", "Log_1", "Log_1")
  expect_error(qc_as_xts(dup, time_col = "time"), "unique")
  txt <- data.frame(time = tt, Log_1 = c("1,5", "2,5"))
  expect_error(qc_as_xts(txt, time_col = "time"), "not numeric")
})

test_that("qc_as_xts parses with exactly the given format and catches truncation", {
  ch <- data.frame(time = c("01.06.2025 00:10", "01.06.2025 00:20"), Log_1 = c(1, 2))
  expect_error(qc_as_xts(ch, time_col = "time"), "not matching")    # default is ISO
  x <- qc_as_xts(ch, time_col = "time", time_format = "%d.%m.%Y %H:%M", verbose = FALSE)
  expect_equal(format(zoo::index(x)[2]), "2025-06-01 00:20:00")
  # a date-only format on date-time strings collapses the day onto midnight
  d <- data.frame(time = c("2025-06-01 00:10:00", "2025-06-01 00:20:00"), Log_1 = c(1, 2))
  expect_error(qc_as_xts(d, time_col = "time", time_format = "%Y-%m-%d"), "duplicated")
})

test_that("the conversion report reaches qc_info and names what was odd", {
  tt <- c("2025-06-01 00:20:00", "2025-06-01 00:10:00")          # not in time order
  d  <- data.frame(time = tt, Log_1 = c(1, 2), Log_9 = c(NA, NA))  # one empty logger
  expect_message(x <- qc_as_xts(d, time_col = "time"), "empty loggers")
  rec <- xts::xtsAttributes(x)$qc_conversion
  expect_equal(rec$empty_loggers, "Log_9")
  expect_true(rec$resorted)
  # the first level moves the report into qc_info and strips it from the data
  res <- T_QC_1_gross_error(x, verbose = FALSE)
  expect_equal(res$qc_info$conversion_temperature$time_col, "time")
  expect_null(xts::xtsAttributes(res$qc_data)$qc_conversion)
  # a humidity run files its own report under its own name
  rh <- RH_QC_2_range(qc_as_xts(d, time_col = "time", verbose = FALSE), verbose = FALSE)
  expect_true("conversion_humidity" %in% names(rh$qc_info))
})

test_that("qc_prepare_input validates the pair and records the dataset", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 20)
  x  <- xts::xts(matrix(1:40, 20, 2, dimnames = list(NULL, c("A", "B"))), order.by = ti)
  p  <- qc_prepare_input(x)
  expect_named(p, c("qc_data", "qc_data_flagged", "qc_info"))
  expect_equal(p$qc_info$dataset_temperature$time_step_sec, 600)
  expect_equal(p$qc_info$dataset_temperature$stations, c("A", "B"))
  expect_equal(p$qc_info$dataset_temperature$variable, "temperature")
  expect_error(qc_prepare_input(data.frame(a = 1)), "data frame")
  bad <- p; bad$qc_data_flagged <- bad$qc_data_flagged[, 1]
  expect_error(qc_prepare_input(bad), "identical dimensions")
  dupn <- x; colnames(dupn) <- c("A", "A")
  expect_error(qc_prepare_input(dupn), "duplicated station names")
})

test_that("T_QC_4 refuses a window too short to ever reach min_non_na", {
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 50)
  x  <- xts::xts(matrix(20, 50, 1, dimnames = list(NULL, "S")), order.by = ti)
  expect_error(T_QC_4_stuck_values(x, window_size = "30 mins"), "min_non_na")
  expect_silent(T_QC_4_stuck_values(x, window_size = "40 mins", verbose = FALSE))
})
