test_that("qc_as_xts pivots the tidy layout without losing an observation", {
  tidy <- biel_tidy()
  x <- qc_as_xts(tidy, verbose = FALSE)
  expect_s3_class(x, "xts")
  expect_equal(dim(x), c(13308L, 7L))
  back <- data.frame(
    time    = rep(format(zoo::index(x), "%Y-%m-%d %H:%M:%S"), times = ncol(x)),
    station = rep(colnames(x), each = nrow(x)),
    temperature = as.vector(zoo::coredata(x)), stringsAsFactors = FALSE)
  ki <- paste(tidy$time, tidy$station); ko <- paste(back$time, back$station)
  expect_setequal(ki, ko)
  expect_equal(tidy$temperature[order(ki)], back$temperature[order(ko)])
})

test_that("qc_as_xts refuses the ambiguous and the destructive cases", {
  tt <- c("2025-06-01 00:10", "2025-06-01 00:20")
  d  <- data.frame(time = tt, Log_1 = c(1, 2), Log_2 = c(3, 4))
  expect_s3_class(qc_as_xts(d, verbose = FALSE), "xts")

  dup <- d; names(dup) <- c("time", "Log_1", "Log_1")
  expect_error(qc_as_xts(dup, verbose = FALSE), "Duplicated column names")

  gap <- data.frame(time = c(tt, NA), Log_1 = c(1, 2, NA))
  expect_error(qc_as_xts(gap, verbose = FALSE), "missing or empty")

  expect_error(qc_as_xts(d, time_col = 9, verbose = FALSE), "out of range")
  expect_error(qc_as_xts(d, time_col = "nope", verbose = FALSE), "not found")

  lg <- data.frame(time = rep(tt, each = 2), station = rep(c("A", "B"), 2), temperature = 1:4)
  expect_equal(dim(qc_as_xts(lg, verbose = FALSE)), c(2L, 2L))
  expect_error(qc_as_xts(rbind(lg, lg[1, ]), verbose = FALSE), "repeated")

  # a date-only format must never collapse real time stamps onto midnight
  x <- qc_as_xts(d, verbose = FALSE)
  expect_equal(length(unique(zoo::index(x))), 2L)
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
