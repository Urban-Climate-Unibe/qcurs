# The elevation rule of the neighbour geometry: a station on the slope is not
# a neighbour of the valley floor.

test_that("ELEV in the metadata keeps a station on the slope out of the consensus", {
  set.seed(4)
  ti <- seq(as.POSIXct("2025-07-01", tz = "UTC"), by = "10 min", length.out = 10 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  base <- 20 + 6 * sin(2 * pi * (mins - 540) / 1440)
  # six stations on the valley floor within a kilometre, one 380 m higher 2 km away
  md <- data.frame(ID = paste0("S", 1:7),
                   LON = 7.62 + c(0, 0.004, 0.008, 0, 0.004, 0.008, 0.026),
                   LAT = 46.75 + c(0, 0, 0, 0.003, 0.003, 0.003, 0.004),
                   Landuse = "Vegetated Areas", ELEV = c(rep(565, 6), 944), stringsAsFactors = FALSE)
  X <- sapply(1:7, function(k) base + rnorm(length(ti), 0, 0.2)); colnames(X) <- md$ID
  # the high site: 3.5 K colder on afternoons, as lapse rate and shade would make it
  X[, 7] <- X[, 7] - 3.5 * pmax(0, sin(2 * pi * (mins - 480) / 1440))
  x <- xts::xts(X, order.by = ti)
  with_elev <- T_QC_6_spatial_consistency(x, metadata = md, verbose = FALSE)
  no_elev   <- T_QC_6_spatial_consistency(x, metadata = md[, setdiff(names(md), "ELEV")], verbose = FALSE)
  i1 <- with_elev$qc_info$t6_spatial_consistency; i0 <- no_elev$qc_info$t6_spatial_consistency
  expect_true(i1$elevation_used); expect_false(i0$elevation_used)
  expect_gt(i0$n_flagged_by_station[["S7"]], 100)          # judged against the valley: objected to
  expect_equal(i1$n_flagged_by_station[["S7"]], 0)          # with elevation: not a neighbour of anyone
  expect_equal(i1$never_evaluated, "S7")
  expect_equal(sum(i1$n_flagged_by_station[1:6]), 0)        # the floor stays clean either way
  expect_equal(sum(i0$n_flagged_by_station[1:6]), 0)
  # a wide tolerance brings the station back in
  wide <- T_QC_6_spatial_consistency(x, metadata = md, max_elev_diff = 500, verbose = FALSE)
  expect_gt(wide$qc_info$t6_spatial_consistency$n_judged_by_station[["S7"]], 0)
  expect_error(T_QC_6_spatial_consistency(x, metadata = md, max_elev_diff = 0), "max_elev_diff")
})
