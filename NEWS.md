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
* `T_QC_8_diurnal_range()` lists `skipped_stations` (no baseline).
* `T_QC_7_spatiotemporal_consistency()` warns about stations without a
  metadata row and lists them in `skipped` as `<id>(no metadata)`; they were
  silently dropped before.

## Fixes

* `qc_find_stuck()`: the population guard compared the whole series against
  the window width instead of `min_valid`; a short series of identical values
  was never judged.
* `T_QC_6_spatial_consistency()` skips itself with a reason when fewer than
  `min_neighbours + 1` stations match (it stopped the chain before), and both
  spatial temperature levels survive metadata that match no station at all
  (`subscript out of bounds` before).
* `RH_QC_6_decoupling()` matches the temperature by exact time stamp and
  logger name through `qc_match_grid()`, like levels 1 and 7; a partial
  overlap is judged where it exists, no overlap skips the level. The
  row-count check is gone.
* Calendar days in `T_QC_8` and `RH_QC_6` come from `format(index, "%Y-%m-%d")`,
  which honours the index time zone on every R version (`as.Date()` only from
  R 4.3 on).

## Structure

* `qc_distance_matrix()` (new, internal): one distance matrix with the LAT/LON
  checks, used by `qc_neighbour_weights()` and `T_QC_7`; `geosphere::distm()`
  replaces the double loop.
* `qc_find_spikes()`, `qc_find_stuck()`, `qc_find_spatial_outliers()` return
  `list(hit, judged)`.
* `T_QC_8` and `RH_QC_6` split the rows by day once instead of comparing the
  whole day vector per day (T8 on 52560 x 40: 2.1 s -> 0.5 s).
* `run_qc_humidity()`: levels 6 and 7 always run and skip themselves;
  `levels_skipped` lists only `rh8`.

## Breaking

* The `qc_info` key of level 7 is `t7_spatiotemporal_consistency`
  (was `t7_spatiotemporal`).
* `RH_QC_6_decoupling()` no longer errors on a row mismatch.
* The three internal helpers return lists.
