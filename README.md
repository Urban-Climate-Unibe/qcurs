# qcurs

## Author

Patrick Kallabis

Research Associate

Institute of Geography & Oeschger Center for Climate Change Research, University of Bern

## What is this about

Quality control for the URS/URBNET low-cost temperature and humidity networks. R package, built on `xts`. Temperature runs through nine levels, humidity through eight; every level reads the same list, writes the same list, and leaves a record of what it did.

The chain follows Amini et al. (2026), *Scientific Data* 13:658, with the review findings fixed (see "Divergences" at the end).

## Install

From GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Urban-Climate-Unibe/qcurs")
library(qcurs)
```

While the repository is private, `remotes` needs a GitHub token: create one under GitHub \> Settings \> Developer settings \> Personal access tokens (scope `repo`) and put it in `~/.Renviron` as `GITHUB_PAT=...`, then restart R.

From a downloaded tarball or zip:

``` r
remotes::install_local("qcurs_0.3.0.tar.gz")
```

Needs `xts`, `zoo`, `geosphere`; `remotes` installs them.

## Input

**Data**: one table per variable, wide layout. One time column, one numeric column per logger, column name = logger ID. Missing values as `NA`.

```         
time;Log_B1;Log_B2;Log_B3
31.05.2025 22:00;17.19;18.08;17.13
31.05.2025 22:10;17.05;17.88;NA
```

**Metadata** (levels 6 and 7 only): a data frame with `ID`, `LAT`, `LON`, `Landuse`. Landuse classes as in the URS metadata: `Sealed Areas`, `Vegetated Areas`, `Forest`, `Water`. Only stations present in both the data and the metadata are used; the rest are reported, not dropped in silence.

Read the files yourself and check them - `qc_as_xts()` does not read files on purpose, a silent CSV misread is the most expensive mistake in this pipeline:

``` r
raw <- read.csv("Biel_summer_2025_raw.csv", sep = ";", dec = ".")
stopifnot(nrow(raw) == 13308)
x <- qc_as_xts(raw, time_col = "time", time_format = "%d.%m.%Y %H:%M")
```

`time_col` is required, by name or column number. `time_format` as in `strptime()`; the default is `"%Y-%m-%d %H:%M:%S"`, time zone `UTC`.

## Run

The runner keeps the order and wires the pieces together:

``` r
t_res  <- run_qc_temperature(x, metadata = meta)
rh_res <- run_qc_humidity(rh, temperature = t_res, metadata = meta)
```

Pass the whole `t_res` to the humidity runner, not just the series: level 1 inherits the temperature flags from it, levels 5 and 7 use its cleaned data.

Per-level parameters go in as a list keyed by level:

``` r
t_res <- run_qc_temperature(x, metadata = meta,
                            params = list(t2 = list(summer_min = -10),
                                          t4 = list(window_size = "1 hour")))
```

Or call the levels yourself, in any order, on a bare `xts` or on the list the previous level returned:

``` r
res <- T_QC_1_gross_error(x)
res <- T_QC_2_out_of_range(res)
res <- T_QC_3_time_consistency(res)
res <- T_QC_4_stuck_values(res)
res <- T_QC_5_climatic_outliers(res)
res <- T_QC_6_spatial_consistency(res, metadata = meta)
res <- T_QC_7_spatiotemporal_consistency(res, metadata = meta)
res <- T_QC_8_diurnal_range(res)
res <- T_QC_9_interpolate(res)   # only if you want filled gaps - see below
```

Every level prints one line with its tally. `verbose = FALSE` switches that off.

There is an R Markdown template with the whole workflow, one chunk per level, parameters written out, summary tables, plot and export:

``` r
rmarkdown::draft("qc_biel.Rmd", template = "qc-workflow", package = "qcurs")
```

## Result

Three elements, same for both chains:

``` r
res$qc_data          # the series, objected cells set to NA
res$qc_data_flagged  # same shape: which level objected to which cell
res$qc_info          # one record per level, in the order they ran
```

### Flag codes

| code | temperature | humidity |
|------------------------|------------------------|------------------------|
| `NA` | no observation | no observation |
| `0` | not objected by any level that reached the cell | same |
| `1` | gross error, fault code | inherited from the temperature run |
| `2` | outside the seasonal range | outside 0-105 % |
| `3` | isolated spike | isolated spike |
| `4` | stuck sensor | stuck sensor |
| `5` | climatic outlier (monthly IQR) | not following its own temperature |
| `6` | off the neighbour consensus | stays below saturation while the network is saturated |
| `7` | extreme in space and time at once | dewpoint off the neighbour consensus |
| `8` | diurnal range collapsed (indoors) | \- |
| `50` | plain gap, filled by the interpolation | same |
| `50+N` | removed by level N, refilled | same |

The first level that objects writes its code; nothing overwrites it later. The code is the level number and cannot be changed.

### Reading `qc_info`

`0` is not "checked by everything". Levels 5-8 refuse where the data cannot carry the test: a single-season campaign has no monthly climatology, a station without a compatible neighbour in 3 km has no consensus. So every level records, per station, what it flagged and what it actually judged:

``` r
i6 <- res$qc_info$t6_spatial_consistency
i6$n_flagged_by_station   # what it objected to
i6$n_judged_by_station    # what it looked at - 0 means the test never reached the station
i6$never_evaluated        # those stations by name
```

The other refusals by name: `t5_climatic_outliers$refused` (station-months), `t7_spatiotemporal_consistency$skipped` (stations), `t8_diurnal_range$skipped_stations`, `rh1_inherit_temperature$loggers_unmatched`, `rh5_decoupling$skipped_stations`, `rh7_dewpoint_consistency$never_evaluated`. A level that could not run at all (no metadata, no temperature, too few stations) leaves `skipped = TRUE` and a `reason`.

`rh1_inherit_temperature$n_flagged_by_t_code` says which temperature levels the humidity inherited from - the humidity flag itself is always 1.

Every record also holds the parameters the level ran with. Together with `qc_data_flagged` that is the full audit trail: which cell, which level, which settings, and where a level was blind.

### Interpolation

Levels T9 and RH8 ADD values. They are off in the runners (`interpolate = TRUE` turns them on) because the method paper says QC flags and never alters the data. If you publish an interpolated series, say so and publish the codes.

## Levels and defaults

| level | test | default |
|------------------------|------------------------|------------------------|
| T1 | outside physical bounds | \> 60 or \<= -40 °C |
| T2 | outside the seasonal bounds | Bern-derived, summer floor -5.85 °C |
| T3 | off the window median AND the past AND the future | 6 K, ±3 steps |
| T4 | no variation in a window with enough valid values | 6 h, sd 0, ≥ 50 % valid, ≥ 5 values |
| T5 | beyond Q1/Q3 ± 4·IQR of the station-month | ≥ 300 values, more than one year |
| T6 | off the weighted neighbour consensus | max(6σ, 3 K), 3 km, 5 compatible neighbours, ≥ 2 valid |
| T7 | 99 % extreme against all 5 neighbours AND both time neighbours | 2500 m |
| T8 | diurnal range below 0.4 of the station median, ≥ 2 days in a row | ≥ 72 values/day, ≥ 14 days |
| T9 | linear fill | gaps ≤ 5 steps |
| RH1 | inherit temperature flags (not the plain gap fill 50) |  |
| RH2 | outside the physical range | 0-105 % |
| RH3 | as T3 | 20 %, ±3 steps |
| RH4 | as T4, except when the window median is ≥ 95 % | 6 h, sd 0.1 |
| RH5 | daily T-RH correlation ≥ -0.3 on days with ≥ 3 K amplitude | ≥ 60 pairs/day |
| RH6 | below 90 % while the network median is ≥ 95 % | 30-day segments, ≥ 12 hits |
| RH7 | as T6, on the dewpoint | same as T6 |
| RH8 | as T9 |  |

Windows are given as text (`"6 hours"`, `"90 mins"`, `"2 days"`) or as a `difftime`; the level converts them with the time step of the data.

## Files

Three kinds of file in `R/`:

- `T_QC_<n>_*.R`, `RH_QC_<n>_*.R` - one level each: parameter checks, the call to the shared test, the flag code, the report.
- `qc_find_*.R` - the tests themselves, written once and used by both chains: spikes (T3/RH3), stuck values (T4/RH4), neighbour consensus (T6/RH7).
- `qc_*.R` - everything around the tests: `qc_as_xts` and `qc_prepare_input` (input), `qc_window_points` (windows), `qc_distance_matrix` and `qc_neighbour_weights` (geometry), `qc_match_grid` (temperature onto the humidity grid), `qc_fill_gaps` (interpolation).
- `run_qc.R` - the two runners.

Exported: the 17 levels, the 2 runners, `qc_as_xts`, `qc_prepare_input`. Everything else is internal.

## Tests

`inst/extdata/` holds the Biel summer 2025 campaign (7 loggers, 10 min, June-August) with six planted error classes; `planted_errors.csv` says where. `tests/testthat/` checks that each class reaches its level, that levels 5 and 7 refuse on a single season, that the real minima survive, and that every level's coverage counts are consistent with its flags. Levels 5 and 7 prove their positive path on synthetic data.

## Divergences from the published pipeline

Documented in the docstring of the level concerned: the k-nearest cap is applied AFTER the landuse masking and there is no sigma gate (T6, RH7); the NA tolerance of the stuck test is real (T4); the seasonal bounds are configurable and their provenance is written down (T2); levels 5 and 7 refuse instead of computing nonsense on short campaigns; level 8 is new.
