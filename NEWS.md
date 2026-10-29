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
