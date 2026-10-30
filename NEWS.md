# qcurs 0.4.1 (2026-10-30)

## Method

* Humidity levels checked against two summers of the Thun network with
  MeteoSwiss Thun as the reference. The sensors read dry at high humidity:
  where MeteoSwiss stands at 99 percent they stand at 73 to 90, each at its
  own ceiling, the devices installed in July 2026 five to eight points above
  the old ones. Three consequences for the chain:
* `RH_QC_4_stuck_values()`: the saturation exception now also applies at the
  station's own ceiling (`sat_quantile`, the 90th percentile of its valid
  values; the threshold used is recorded per station as `sat_ceiling`).
  Before, 95 percent of the 353 stuck flags of 2025 sat on dew nights at the
  forest stations - the sensor flat at 85 while MeteoSwiss read 100. A sensor
  stuck at one value for the whole record has that value as its ceiling and
  is left to level 5.
* `RH_QC_6_saturation_drift()` skips itself and says so when the network
  median never reaches `network_min`, instead of reporting a clean run with
  nothing judged; `n_saturated_steps` and `network_max` are in the record.
* `RH_QC_7_dewpoint_consistency()`: floor 4 K instead of 3. Healthy low-cost
  sensors differ by 2 to 3 K of dewpoint at saturation, and the 3 K floor
  objected to the wettest sensor on humid nights and the driest two on hot
  afternoons - the ones MeteoSwiss agreed with (97 of 116 flags in 2025).
  The one real fault, a sensor at 70 to 92 percent on dry August days, is
  still found. Temperature level 6 keeps its 3 K.

# qcurs 0.4.0 (2026-10-30)

## Method

* Two new temperature levels for the faults the single-value tests cannot
  see. `T_QC_9_bias_shift()` (code 9): the daily-mean anomaly of a station
  against the median of the other stations, smoothed over five days, is
  compared with its baseline from the first 30 judged days; a departure of
  more than 1.5 K for five judged days in a row is objected to, and a day on
  which more than a quarter of the stations depart is weather and not
  judged. It finds the change, not an offset from day one; the `episodes`
  table (station, from, to, days, shift) says where. On two summers of Thun
  no healthy station departed by more than 1.15 K; a planted +2 K was found
  with gaps at the threshold, a planted +3 K to 99 percent.
  `T_QC_10_clock_offset()` (code 10): per station and week the lag at which
  the station correlates best with the median of the other stations; a lag
  that rounds to an hour or more, with the same sign over two weeks, is
  objected to, and the `episodes` table gives the lag in minutes. Forest and
  lake sites run 20 to 55 minutes behind the network and the open reference
  site up to 30 minutes ahead, in both summers, which is why the threshold
  is an hour and not the 30 minutes first planned; a planted hour was found
  in every full week. Both levels need no metadata.

## Structure

* BREAKING: the interpolation is now the last level, `T_QC_11_interpolate()`
  with the record `t11_interpolate` and the runner key `t11`; the file
  `T_QC_9_interpolate.R` was renamed. Flag codes and behaviour are unchanged.
  Level 1 of the humidity chain accepts the codes 0-11 and 50-61.
* `qc_row_median()` (internal): row medians by one sort instead of one
  `median()` per row, so the median of the other stations costs a second,
  not a minute.

# qcurs 0.3.7 (2026-10-30)

## Method

* Levels T6 and RH7 take the elevation into account when the metadata carry
  an `ELEV` column (metres): two stations are neighbours only within
  `max_elev_diff` (100 m) of each other. On Thun 2025 the published chain
  and qcurs alike objected to Goldiwil (944 m against a network at 565 m,
  2.6 K colder on afternoons by lapse rate) - 151 values that were the site,
  not the sensor. Without the column the rule is off and the record says so
  (`elevation_used`). The Thun acceptance run otherwise confirmed the chain:
  the one real fault of 2025 (Thunerallmend stuck, 11 243 values) is found
  cell for cell as in URS V1, and the four dead or stuck sensors of 2026 are
  caught completely with the device changes matching to the hour.

# qcurs 0.3.6 (2026-10-30)

## Structure

* Four plot functions, base graphics, one colour per level, each returning
  the numbers it drew: `qc_plot_timeline()` (flags per day, one line per
  level), `qc_plot_levels()` (bars per level with the share of the judged
  values, or stacked bars per station sorted by share), `qc_plot_calendar()`
  (station-by-day share of objected values: weather days are vertical
  stripes, a broken sensor a horizontal bar, gaps grey) and
  `qc_plot_station()` (one series with its flags, raw behind it when
  given). The template uses them.

# qcurs 0.3.5 (2026-10-29)

## Method

* `T_QC_8_diurnal_range()` is now specifically the test for a logger that was
  taken down and kept logging indoors. A collapsed diurnal range alone also
  describes a station alone under high fog (Grenchen, 20/21 January 2026:
  2.2 K against 7.7 K in the network, 288 values flagged); what separates
  the two is the weather: the fog station still takes every front with the
  network (daily means r = 0.95), the logger in a heated room sits at 21. A
  collapsed run is therefore objected to only when the swing of its daily
  means stays below `range_frac` of the network's swing over the same days.
  Runs during a calm spell (network swing below `min_network_swing`, 2 K)
  are not decidable and are listed in `runs_undecided`.
  `min_consecutive_days` defaults to 3 (was 2). The planted Biel indoor
  episode is found as before.

# qcurs 0.3.4 (2026-10-29)

## Method

* `T_QC_5_climatic_outliers()` judges the DAILY EXTREMES, as the method paper
  intends, instead of the 10-minute values: per station and pooled calendar
  month, a day whose maximum exceeds Q3 + m * IQR of the daily maxima or whose
  minimum falls below Q1 - m * IQR of the daily minima is objected to, on the
  cells beyond the bound. Outer tails only. `ext_lim_factor` defaults to 2.5:
  the published 4 leaves the bounds 25 K beyond the quartiles on 10-minute
  values and on daily extremes alike and never bites; 2.5 objected to four
  values on six years of three MeteoSwiss stations (two cold December
  nights at Thun). `min_n` (values) is replaced by
  `min_days` (60) and `min_obs_per_day` (72); the applied bounds are recorded
  in `qc_info$t5_climatic_outliers$bounds`. Single-season campaigns are
  refused as before.
* `landuse_mode = "none"` for `T_QC_6_spatial_consistency()` and
  `RH_QC_7_dewpoint_consistency()`: distance only, every pair compatible. For
  small networks where the landuse mask leaves hardly any station with
  neighbours - Biel goes from 1 to 6 evaluable stations. The price: on the
  clean Biel summer "none" objects to 0.75 percent of the values, mostly late
  afternoon at the vegetated and the water site, which genuinely differ from
  the sealed neighbours by more than 3 K then. Default stays "graded".

# qcurs 0.3.3 (2026-10-29)

## Structure

* `qc_summary()` is part of the package: one row per level and station with
  flagged, judged and judged as a share of the observed values; it used to
  live in the R Markdown template only. The template calls it now.
* The vignette, a copy of the template, moved to `archiv/`; `R CMD check` is
  clean of WARNINGs.
* DESCRIPTION names both chains.

# qcurs 0.3.2 (2026-10-29)

## Method

* `T_QC_2_out_of_range()`: the autumn minimum is -19.90 (was -8.42). Same
  sign error as the summer bound: the published value is the seasonal extreme
  plus 6 instead of minus 6 (the Bern-Zollikofen autumn minimum 1995-2023 is
  -13.9). Healthy MeteoSwiss stations in Bern, Thun and Interlaken fell below
  -8.42 in November 2024 and 2025.
* `T_QC_3_time_consistency()`: `threshold` default 4 K (was 6). On six years
  of 10-minute data from three stations the largest genuine isolated
  deviation is below 4 K; 6 K caught nothing that 4 K does not.
* `qc_prepare_input()` records the value `resolution` of the dataset (the
  smallest step between distinct values); `T_QC_4_stuck_values()` reports it
  and warns when the window is shorter than 6 hours at 0.05 K or coarser,
  because at that resolution fog produces genuine plateaus of that length
  (2020-2025 MeteoSwiss series: 1 hour 2000-4200 false flags, 6 hours none).

# qcurs 0.3.1 (2026-10-29)

## Method

* `RH_QC_5_decoupling()` separates weather from a sensor. The daily T-RH
  correlation also breaks down on days with an air mass change (warm front
  with rain, cold advection); on healthy MeteoSwiss stations the plain daily
  test objected to 5 percent of all days, 13 percent in winter, and on a
  36-station network 90 percent of its objections fell on 13 days with 3 to
  30 stations at once. Two new rules, both after the pattern of `T_QC_8`:
  `max_network_frac` (default 0.25) - a day on which more than that share of
  the judged stations (and at least two) look decoupled is a weather day and
  is judged at no station, listed in `days_unjudged`; `min_consecutive_days`
  (default 2) - a station is objected to only after that many judged days in
  a row. `min_consecutive_days = 1, max_network_frac = 1` gives the 0.3.0
  behaviour. Temperature flags are unchanged.

## Structure

* The 0.2.0 humidity levels and their Rd files were still in `R/` and `man/`
  and compiled into the namespace (`RH_QC_5_stuck_values` twice). Moved to
  `archiv/`, excluded from the build; `load_packages.R` archived too.

# qcurs 0.3.0 (2026-09-29)

Review of the whole temperature and humidity chain. Flags and cleaned data on
the shipped Biel test set are bit-identical to 0.2.0; what changed is what the
result RECORDS and how the edges behave.

## Traceability

* Every level records `n_judged_by_station`: how many cells of each station
  the test actually reached. `n_flagged <= n_judged <= valid values` holds
  everywhere and is tested. A flag of 0 means "not objected by any level that
  reached the cell" - the README says so now.
* `RH_QC_1_inherit_temperature()` records `n_flagged_by_t_code`: which
  temperature levels were inherited how often (the humidity flag is always 1).
* A level that runs a second time on the same chain warns; its earlier flags
  stay, its report is replaced (`qc_prepare_input(level = ...)`).
* `T_QC_8_diurnal_range()` lists `days_unjudged` (no network yardstick).
* `T_QC_7_spatiotemporal_consistency()` warns about stations without a
  metadata row and lists them in `skipped` as `<id>(no metadata)`; they were
  silently dropped before.

## Method

* `T_QC_8_diurnal_range()` judges the diurnal range against the NETWORK's
  median range of the same day instead of the station's own median: a fog
  week collapses every station together and stays invisible, an indoor sensor
  collapses alone and shows. `min_reference_days` is gone, `min_stations`
  (default 3) is new; with fewer stations the level skips itself.
* `qc_find_spatial_outliers()` (levels T6 and RH7) takes its objections one
  at a time per time step: remove the worst offender, rebuild consensus and
  spread, test again. A bad neighbour can no longer drag a good station out
  with it or hide behind the spread it inflates.
* `T_QC_7_spatiotemporal_consistency()` judges every station against the
  matrix as it entered the level; the result no longer depends on the column
  order.
* `T_QC_2_out_of_range()`: the `summer_min` shortcut is gone, the bounds list
  is the parameter. Default summer minimum stays -5.85 (the published 7.85 is
  a sign error).

## Fixes

* `qc_find_stuck()`: the population guard compared the whole series against
  the window width instead of `min_valid`; a short series of identical values
  was never judged.
* `T_QC_6_spatial_consistency()` skips itself with a reason when fewer than
  `min_neighbours + 1` stations match (it stopped the chain before), and both
  spatial temperature levels survive metadata that match no station at all
  (`subscript out of bounds` before).
* `RH_QC_5_decoupling()` matches the temperature by exact time stamp and
  logger name through `qc_match_grid()`, like levels 1 and 7; a partial
  overlap is judged where it exists, no overlap skips the level. The
  row-count check is gone.
* Calendar days in `T_QC_8` and `RH_QC_5` come from `format(index, "%Y-%m-%d")`,
  which honours the index time zone on every R version (`as.Date()` only from
  R 4.3 on).

## Template

* `inst/rmarkdown/templates/qc-workflow`: an R Markdown template of the
  whole workflow (`rmarkdown::draft(..., template = "qc-workflow",
  package = "qcurs")`) - read, convert, one chunk per level with the
  parameters written out, summary tables, plot, export; humidity part
  switched on by setting `rh_file`.

## Structure

* `qc_distance_matrix()` (new, internal): one distance matrix with the LAT/LON
  checks, used by `qc_neighbour_weights()` and `T_QC_7`; `geosphere::distm()`
  replaces the double loop.
* `qc_find_spikes()`, `qc_find_stuck()`, `qc_find_spatial_outliers()` return
  `list(hit, judged)`.
* `T_QC_8` and `RH_QC_5` split the rows by day once instead of comparing the
  whole day vector per day (T8 on 52560 x 40: 2.1 s -> 0.5 s).
* `run_qc_humidity()`: levels 5 and 7 always run and skip themselves;
  `levels_skipped` lists only `rh8`.

## Breaking

* Humidity levels 4-6 are renumbered: 4 = stuck values (was 5), 5 = decoupling
  (was 6), 6 = saturation drift (was 4). Function names, `qc_info` keys and
  flag codes follow the new numbers; a humidity flag matrix from 0.2.0 reads
  differently in these three codes.

* The `qc_info` key of level 7 is `t7_spatiotemporal_consistency`
  (was `t7_spatiotemporal`).
* `RH_QC_5_decoupling()` no longer errors on a row mismatch.
* The three internal helpers return lists.
