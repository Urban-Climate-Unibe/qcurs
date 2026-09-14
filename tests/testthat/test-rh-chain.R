# Synthetic four-station network with a healthy diurnal cycle, one planted
# error per RH level, run through the whole humidity chain.

rh_fixture <- function() {
  set.seed(11)
  ti <- seq(as.POSIXct("2025-06-01", tz = "UTC"), by = "10 min", length.out = 20 * 144)
  mins <- as.numeric(format(ti, "%H")) * 60 + as.numeric(format(ti, "%M"))
  tcyc <- 18 + 5 * sin(2 * pi * (mins - 540) / 1440)             # diurnal temperature
  TT <- sapply(1:4, function(k) tcyc + rnorm(length(ti), 0, 0.3) + (k - 2) * 0.4)
  # healthy negative coupling, with nights around 89: the network must NOT
  # brush the 95 percent saturation criterion outside the fog day, because one
  # healthy saturated step per segment acquits a drifter by construction
  RH <- 90 - 2.5 * (TT - 13) + matrix(rnorm(length(ti) * 4, 0, 1.5), ncol = 4)
  colnames(TT) <- colnames(RH) <- paste0("S", 1:4)
  day <- as.Date(ti)

  # E2: impossible + clampable values at S1
  RH[10, 1] <- 130;  RH[20, 1] <- 103
  # E3: isolated +30 spike at S2
  RH[500, 2] <- RH[500, 2] + 30
  # E4: fog day 8 (whole network saturated) with S3 drifting at ~80
  fog <- which(day == as.Date("2025-06-08"))
  # on a real fog day the TEMPERATURE is flat too - otherwise RH6 would rightly
  # call the pinned humidity "decoupled"; the flat T triggers its trange skip
  TT[fog, ] <- 12 + matrix(rnorm(length(fog) * 4, 0, 0.3), ncol = 4)
  # the healthy stations sit NEAR-CONSTANT at 97: sd < stuck_sd, so without
  # the saturation exception RH5 would flag the whole fog day
  RH[fog, ] <- 97 + matrix(rnorm(length(fog) * 4, 0, 0.05), ncol = 4)
  RH[fog, 3] <- 80 + rnorm(length(fog), 0, 0.5)
  # E5: stuck block at S4 (60.0 for 40 steps, well below saturation)
  RH[2000:2039, 4] <- 60.0
  # E6: decoupled day 15 at S2 (constant-ish, uncorrelated with its T)
  dec <- which(day == as.Date("2025-06-15"))
  RH[dec, 2] <- 55 + rnorm(length(dec), 0, 1.5)
  # E7: dewpoint offset at S1, day 18, midday: a DRY bias of -30 percentage
  # points. A wet bias saturates against the 100 percent ceiling and barely
  # moves the dewpoint; a dry bias shifts it by several Kelvin.
  off <- which(day == as.Date("2025-06-18"))[61:96]
  RH[off, 1] <- pmax(RH[off, 1] - 30, 20)

  md <- data.frame(ID = paste0("S", 1:4),
                   LON = 7.44 + (0:3) * 0.004, LAT = 46.95 + (0:3) * 0.002,
                   Landuse = "Sealed Areas", stringsAsFactors = FALSE)
  list(rh = xts::xts(RH, order.by = ti), tt = xts::xts(TT, order.by = ti),
       md = md, fog = fog, dec = dec, off = off)
}

test_that("the humidity chain finds every planted error class", {
  f <- rh_fixture()
  r <- RH_QC_2_range(f$rh, verbose = FALSE)
  r <- RH_QC_3_spike(r, verbose = FALSE)
  r <- RH_QC_4_saturation_drift(r, verbose = FALSE)
  r <- RH_QC_5_persistence(r, verbose = FALSE)
  r <- RH_QC_6_decoupling(r, temperature = f$tt, verbose = FALSE)
  r <- RH_QC_7_dewpoint_consistency(r, temperature = f$tt, metadata = f$md, verbose = FALSE)
  FLG <- zoo::coredata(r$qc_data_flagged)

  expect_equal(unname(FLG[10, 1]), 2)                        # impossible value flagged
  expect_equal(unname(zoo::coredata(r$qc_data)[20, 1]), 100) # tolerance value clamped, kept
  expect_gte(r$qc_info$rh2_range$n_clamped, 1)               # the plant, plus any natural 100.x
  expect_equal(unname(FLG[500, 2]), 3)                       # spike
  expect_gte(mean(FLG[f$fog, 3] == 4), 0.95)                 # saturation drifter
  expect_true(all(FLG[2000:2039, 4] == 5))                   # stuck block
  expect_gte(mean(FLG[f$dec, 2] == 6), 0.95)                 # decoupled day
  expect_gte(mean(FLG[f$off, 1] == 7), 0.8)                  # dewpoint offset
  # the healthy remainder stays essentially untouched: total flags stay close
  # to the planted budget (fog drift + stuck + decoupled day + offset + slack)
  expect_lt(sum(FLG %in% 2:7, na.rm = TRUE),
            length(f$fog) + 40 + length(f$dec) + length(f$off) + 20)
  # every level left its record
  expect_true(all(c("dataset_humidity", "rh2_range", "rh3_spike", "rh4_saturation_drift",
                    "rh5_persistence", "rh6_decoupling", "rh7_dewpoint_consistency")
                  %in% names(r$qc_info)))
})

test_that("the saturation exception keeps fog days unflagged in RH5", {
  f <- rh_fixture()
  r5 <- RH_QC_5_persistence(f$rh, verbose = FALSE)
  # the fog day is near-constant at 97 but sits ABOVE stuck_sat_max: exempt
  expect_true(all(zoo::coredata(r5$qc_data_flagged)[f$fog, c(1, 2, 4)] == 0))
})

test_that("RH6 and RH7 accept the temperature chain result and refuse misalignment", {
  f <- rh_fixture()
  t_res <- T_QC_1_gross_error(f$tt, verbose = FALSE)
  r <- RH_QC_6_decoupling(f$rh, temperature = t_res, verbose = FALSE)
  expect_true("rh6_decoupling" %in% names(r$qc_info))
  expect_error(RH_QC_6_decoupling(f$rh, temperature = f$tt[1:100, ]), "Row mismatch")
  r7 <- RH_QC_7_dewpoint_consistency(f$rh, temperature = f$tt, metadata = f$md[1:2, ], verbose = FALSE)
  expect_true(r7$qc_info$rh7_dewpoint_consistency$skipped)   # <3 stations: honest skip
})
