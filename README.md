# qcurs

Quality control for urban low-cost temperature and humidity networks
(URBNET/URS campaigns).

## Install

```r
# install.packages("remotes")
remotes::install_local("qcurs_0.3.0.tar.gz")
```

## Use

```r
library(qcurs)

# 1. read the file YOURSELF and check the result - qc_as_xts() deliberately
#    does not read files, because CSV failures are silent
raw <- read.csv("Biel_summer_2025_raw.csv", sep = ";", dec = ".")
stopifnot(nrow(raw) == 13308)

# 2. convert: a WIDE table (one time column, one numeric column per logger),
#    the time column named exactly, the format given explicitly
x <- qc_as_xts(raw, time_col = "time", time_format = "%d.%m.%Y %H:%M")

# 3a. the runner fixes the order (every level stays individually callable)
res    <- run_qc_temperature(x, metadata = meta)
rh_res <- run_qc_humidity(rh, temperature = res, metadata = meta)

# 3b. or run the chain by hand; every level takes and returns the same list
res <- T_QC_1_gross_error(x)
res <- T_QC_2_out_of_range(res)
res <- T_QC_3_time_consistency(res)
res <- T_QC_4_stuck_values(res)
res <- T_QC_5_climatic_outliers(res)
res <- T_QC_6_spatial_consistency(res, metadata = meta)
res <- T_QC_7_spatiotemporal_consistency(res, metadata = meta)
res <- T_QC_8_diurnal_range(res)
res <- T_QC_9_interpolate(res)          # OPTIONAL - this one ADDS values

res$qc_data           # cleaned series: objected cells are NA
res$qc_data_flagged   # flag matrix: which level objected to which cell
res$qc_info           # what every level did, judged, refused, skipped
```

Every level can also be called on a bare `xts` on its own; it then starts a
fresh chain (`qc_data_flagged` all 0 where a value exists).

## Flag codes

| code | temperature | humidity |
|------|-------------|----------|
| `NA` | no observation | no observation |
| `0` | not objected by any level that reached the cell | same |
| `1` | gross error | inherited from the temperature run (`n_flagged_by_t_code` says from which level) |
| `2` | seasonal range | physical range |
| `3` | isolated spike | isolated spike |
| `4` | stuck values | saturation drift |
| `5` | climatic outlier | stuck values |
| `6` | spatial consensus | decoupled from own temperature |
| `7` | spatiotemporal | dewpoint consensus |
| `8` | diurnal range collapse | - |
| `50` | plain gap filled by the interpolation | same |
| `50 + N` | removed by level N, refilled by the interpolation | same |

A cell is flagged once: the first level that objects writes its code, later
levels never overwrite it. The order of the levels is the order of the
records in `qc_info`.

## Reading a result

`0` does not mean "checked by every level". Levels 5, 6, 7 and 8 refuse where
the data cannot support the test (single-season campaign, no compatible
neighbours, too few days). Every level therefore records, per station, how
many cells it actually judged:

```r
res$qc_info$t6_spatial_consistency$n_judged_by_station   # 0 = never reached
res$qc_info$t6_spatial_consistency$never_evaluated        # the same stations by name
res$qc_info$t5_climatic_outliers$refused                  # station-months without a verdict
res$qc_info$t7_spatiotemporal_consistency$skipped         # stations without a full neighbourhood
res$qc_info$t8_diurnal_range$skipped_stations             # stations without a baseline
```

`n_flagged_by_station <= n_judged_by_station <= valid values` holds for every
level. A level that runs twice on the same chain warns: its earlier flags
stay, its report is replaced.

## Test data

`inst/extdata/` ships the Biel summer 2025 campaign with six documented error
classes planted into it (tidy layout, manifest in `planted_errors.csv`; the
test helper pivots it to the wide layout `qc_as_xts()` expects).
`data-raw/make_testdata.R` regenerates it. The test suite asserts that every
class reaches its intended detector, that levels 5 and 7 refuse honestly on a
single-season campaign, that the real record stays untouched, and that every
level's coverage record is consistent with its flags.

## Divergences from the published pipeline

This is a rebuild of the chain in Amini et al. (2026), *Scientific Data* 13:658,
with the review findings corrected. Each divergence is documented in the
docstring of the level it affects: neighbour selection order (level 6), the
removed sigma gate (level 6), the restored NA tolerance (level 4), the
seasonal bound provenance and sign error (level 2), position-based expansion
(level 4), and the honest refusals (levels 5 and 7). Level 8 is new.
