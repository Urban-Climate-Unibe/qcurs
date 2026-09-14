# qcurbnet

Quality control for urban low-cost temperature networks (URBNET/URS campaigns).

## Install

```r
# install.packages("remotes")
remotes::install_local("qcurbnet_0.1.0.tar.gz")
```

## Use

```r
library(qcurbnet)

# 1. read the file YOURSELF and check the result - qc_as_xts() deliberately
#    does not read files, because CSV failures are silent
raw <- read.csv("Biel_test_tidy.csv", sep = ";", dec = ".")
stopifnot(nrow(raw) == 93156)

# 2. convert (wide or tidy input, both work)
x <- qc_as_xts(raw)

# 3a. or let the runner fix the order (every level stays individually callable)
res    <- run_qc_temperature(x, metadata = meta)
rh_res <- run_qc_humidity(rh, temperature = res, metadata = meta)

# 3b. run the chain manually; every level takes and returns the same list
res <- T_QC_1_gross_error(x)
res <- T_QC_2_out_of_range(res)
res <- T_QC_3_time_consistency(res)
res <- T_QC_4_temporal_persistence(res)
res <- T_QC_5_climatic_outliers(res)
res <- T_QC_6_spatial_consistency(res, metadata = meta)
res <- T_QC_7_spatiotemporal_consistency(res, metadata = meta)
res <- T_QC_8_diurnal_range(res)
res <- T_QC_9_interpolate(res)          # OPTIONAL - this one ADDS values

res$qc_data           # cleaned series
res$qc_data_flagged   # flag matrix
res$qc_info           # what every level did, refused and could not reach
```

## Flag codes

| code | meaning |
|------|---------|
| `0` | checked, unobjected |
| `NA` | no observation |
| `1`-`8` | objected by that level |
| `50` | interpolated plain gap (level 9) |
| `50 + N` | value removed by level N and refilled by level 9 |

## Test data

`inst/extdata/` ships the Biel summer 2025 campaign with six documented error
classes planted into it (tidy layout, manifest in `planted_errors.csv`).
`data-raw/make_testdata.R` regenerates it. The test suite asserts that every
class reaches its intended detector, that levels 5 and 7 refuse honestly on a
single-season campaign, and that the real record stays untouched.

## Divergences from the published pipeline

This is a rebuild of the chain in Amini et al. (2026), *Scientific Data* 13:658,
with the review findings corrected. Each divergence is documented in the
docstring of the level it affects: neighbour selection order (level 6), the
removed sigma gate (level 6), the restored NA tolerance (level 4), the
seasonal bound provenance and sign error (level 2), position-based expansion
(level 4), and the honest refusals (levels 5 and 7). Level 8 is new.
