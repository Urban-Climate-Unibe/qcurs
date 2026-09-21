# ------------------------------------------------------------------------------
# make_testdata.R - builds the QC test dataset from the clean Biel summer 2025
# campaign by planting DOCUMENTED errors, one per detectable failure class.
#
# Output is TIDY (long): one row per observation, one column per variable
#   time | station | temperature
# Missing observations are kept as explicit NA rows, so the time grid of the
# campaign survives the round trip and a gap stays visibly a gap.
# qc_as_xts() pivots this to the wide layout the QC chain needs.
#
# The manifest planted_errors.csv is tidy too - one row per planted cell:
#   id | time | station | what | expect_code
#
# Planted classes and their intended detector:
#   E1  logger fault code      -45 degC          -> level 1 (gross error)
#   E2  impossible summer heat  55 degC          -> level 2 (out of range)
#   E3  isolated spike          value + 15 K     -> level 3 (spike); the gap it
#                                                  leaves is refilled by level 9
#                                                  with code 50 + 3 = 53
#   E4  stuck sensor            40 x 24.00 degC  -> level 4 (stuck values)
#   E5  night offset            +12 K, 8 hours   -> level 6 (spatial) - planted
#                                                  at Log_B4, the ONLY Biel
#                                                  station with >= 2 compatible
#                                                  neighbours after the fixed
#                                                  k-cap; at night so the sum
#                                                  stays below the seasonal max
#   E6  indoor episode          5 days 21+-0.2   -> level 8 (diurnal collapse) -
#                                                  planted at Log_B1 (Water,
#                                                  ZERO compatible neighbours),
#                                                  the station class level 6
#                                                  cannot see; noise sd 0.2 so
#                                                  level 4 does not claim it
#
# Levels 5 and 7 are demonstrated on the campaign as HONEST REFUSALS (no
# multi-year climatology / no station with 5 neighbours within 2.5 km); their
# positive detection paths are proven on synthetic data in test_all_levels.R.
# The real coldest Biel reading (Log_B7, 7.18 degC) is kept as a NEGATIVE
# control: it must survive the whole chain unflagged.
# ------------------------------------------------------------------------------

Sys.setlocale("LC_ALL", "C.UTF-8")           # the container needs this for UTF-8 reads
suppressMessages({library(xts); library(zoo)})
options(xts_check_TZ = FALSE)

# run this from the package root (the folder containing R/ and testdata/)
base <- Sys.getenv("QC_PKG_ROOT", ".")
if (!dir.exists(file.path(base, "testdata"))) stop("Run from the package root: testdata/ not found.")

# ---- read the clean campaign (delivered wide) -------------------------------
raw <- read.csv(file.path(base, "inst/extdata/Biel_summer_2025_raw.csv"),
                sep = ";", dec = ".", fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
names(raw)[1] <- "time"
stopifnot(nrow(raw) == 13308, ncol(raw) == 8)          # verify the read BEFORE using it
tv <- as.POSIXct(raw$time, format = "%d.%m.%Y %H:%M", tz = "UTC")
stopifnot(!anyNA(tv))
stations <- setdiff(names(raw), "time")

# helper: row index of a given time stamp
at <- function(stamp) which(tv == as.POSIXct(stamp, tz = "UTC"))
# helper: first row >= stamp where station s has a value and a full +/-3 context
with_context <- function(stamp, s) {
  i0 <- at(stamp); v <- raw[[s]]
  for (i in i0:(nrow(raw) - 3)) if (all(!is.na(v[(i - 3):(i + 3)]))) return(i)
  stop("no context found")
}
# helper: record a plant in the manifest, as station/time pairs (never row numbers,
# because row numbers do not survive a reshape - the keys do)
plants <- list()
note <- function(id, s, rows, what, code)
  plants[[id]] <<- data.frame(id = id, time = format(tv[rows], "%Y-%m-%d %H:%M:%S"),
                              station = s, what = what, expect_code = code,
                              stringsAsFactors = FALSE)

# ---- E1: logger fault code at Log_B2 ----------------------------------------
i <- at("2025-06-10 12:00"); raw$Log_B2[i] <- -45
note("E1", "Log_B2", i, "fault code -45", 1)

# ---- E2: impossible summer heat at Log_B5 -----------------------------------
i <- at("2025-07-05 15:00"); raw$Log_B5[i] <- 55
note("E2", "Log_B5", i, "55 degC in summer", 2)

# ---- E3: isolated spike at Log_B3 (needs temporal context) ------------------
i <- with_context("2025-06-20 03:00", "Log_B3"); raw$Log_B3[i] <- raw$Log_B3[i] + 15
note("E3", "Log_B3", i, "+15 K isolated spike", 3)

# ---- E4: stuck sensor at Log_B6 (40 identical readings = 6.5 h) -------------
i <- at("2025-07-15 00:00"); rows4 <- i:(i + 39); raw$Log_B6[rows4] <- 24.00
note("E4", "Log_B6", rows4, "constant 24.00 x 40", 4)

# ---- E5: +12 K night offset at Log_B4 (the only spatially covered station) --
i <- at("2025-08-05 22:00"); rows5 <- i:(i + 47)
ok5 <- rows5[!is.na(raw$Log_B4[rows5])]                 # only shift existing readings
raw$Log_B4[ok5] <- raw$Log_B4[ok5] + 12
note("E5", "Log_B4", ok5, "+12 K offset (night, 8 h)", 6)

# ---- E6: indoor episode at Log_B1 (5 whole days, no spatial support) --------
set.seed(42)                                            # reproducible noise
i <- at("2025-08-15 00:00"); rows6 <- i:(i + 5 * 144 - 1)
raw$Log_B1[rows6] <- round(21 + rnorm(length(rows6), 0, 0.2), 2)
note("E6", "Log_B1", rows6, "indoor 21 +- 0.2, 5 days", 8)

# ---- negative control: the real Biel minimum must survive -------------------
j <- which.min(raw$Log_B7); stopifnot(abs(raw$Log_B7[j] - 7.18) < 0.01)
note("N1", "Log_B7", j, "real minimum 7.18 - must NOT be flagged", 0)

# ---- reshape wide -> tidy and write -----------------------------------------
# one row per (time, station); NA observations are KEPT so the grid stays complete
to_tidy <- function(df, stamps) {
  out <- data.frame(
    time        = rep(format(stamps, "%Y-%m-%d %H:%M:%S"), times = length(stations)),
    station     = rep(stations, each = nrow(df)),
    temperature = unlist(df[stations], use.names = FALSE),
    stringsAsFactors = FALSE)
  out[order(out$time, out$station), ]                   # sorted by observation
}
tidy <- to_tidy(raw, tv)
stopifnot(nrow(tidy) == 13308 * 7)                      # every cell became a row
write.table(tidy, file.path(base, "inst/extdata/Biel_test_tidy.csv"),
            sep = ";", dec = ".", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8", na = "")

manifest <- do.call(rbind, plants)
write.table(manifest, file.path(base, "inst/extdata/planted_errors.csv"),
            sep = ";", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

message(sprintf("Tidy test dataset: %d observations (%d time steps x %d stations), %.1f%% valid.",
                nrow(tidy), nrow(raw), length(stations), 100 * mean(!is.na(tidy$temperature))))
message(sprintf("Manifest: %d planted cells across 6 error classes + 1 negative control.",
                sum(manifest$expect_code > 0)))
