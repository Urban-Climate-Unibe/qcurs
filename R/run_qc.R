#' Run the whole temperature QC chain in the fixed order
#'
#' Convenience wrapper around levels T_QC_1 to T_QC_9 (and optionally 10). It
#' exists to prevent the two classic application mistakes: running the levels
#' in the wrong order, and silently forgetting the spatial levels when the
#' metadata are at hand. Every level remains individually callable; this
#' function adds nothing the manual chain does not have.
#'
#' Level 8 runs before levels 6 and 7: it judges whole days of a station, and
#' the neighbour comparisons would thin a logger indoors out before it gets
#' there. Level 9 runs last: a daily course distorted by sun on the sensor
#' shows the same lag as a wrong clock, and level 6 has to take those hours
#' out first. Numbers and flag codes stay as they are.
#'
#' Without `metadata` the spatial levels 6 and 7 skip themselves, with the
#' reason in their record - useful for a single isolated station, where they
#' could never run anyway. The
#' interpolation is off by default, because it alters the data (see
#' `T_QC_10_interpolate`). `levels` limits the run to a subset of the levels,
#' with the sign convention of R indexing.
#'
#' @param input xts of temperature (or an already started chain list).
#' @param metadata Data frame with ID, LAT, LON, Landuse for levels 6 and 7,
#'   or NULL to skip both.
#' @param interpolate Also run the optional level 10 gap interpolation.
#' @param levels The levels allowed to run: numbers or names like `"T3"`,
#'   positive to select, negative to exclude - `-9` runs everything but
#'   level 9, `c(1:3, 6)` only those four. Levels 6 and 7 still need
#'   `metadata`, level 10 still needs `interpolate = TRUE`.
#' @param params Named list of per-level parameter overrides, keyed t1..t10,
#'   e.g. `list(t3 = list(threshold = 5), t6 = list(k_sigma = 5))`. Every
#'   entry is passed to that level as arguments, exactly as given.
#' @param verbose Passed to every level.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info; the
#'   levels run and skipped are listed in `qc_info$run_qc_temperature$levels`
#'   and `$levels_skipped`.
#'
#' @examples
#' \dontrun{
#' res <- run_qc_temperature(x, metadata = meta)
#' res <- run_qc_temperature(x, metadata = meta, levels = -9)
#' res <- run_qc_temperature(x, metadata = meta, interpolate = TRUE,
#'                           params = list(t3 = list(threshold = 5)))
#' }
#'
#' @import xts
#' @import zoo
#' @export
run_qc_temperature <- function(input,
                               metadata = NULL,
                               interpolate = FALSE,
                               levels = 1:10,
                               params = list(),
                               verbose = TRUE) {
  # the overrides must be a named list keyed by level
  if (!is.list(params)) stop("params must be a named list, e.g. list(t3 = list(threshold = 5)).")
  ok_keys <- paste0("t", 1:10)
  if (length(params) && (is.null(names(params)) || !all(names(params) %in% ok_keys)))
    stop(sprintf("params keys must be among %s.", paste(ok_keys, collapse = ", ")))
  # the levels this run may touch; the rest is left out and recorded as skipped
  run <- qc_pick_levels(levels, 10, "t")
  # the chain list right away, so the record below can be written even when
  # no level runs
  input <- qc_prepare_input(input, what = "temperature")
  
  # one level call: base arguments, overridden by the caller's params entry
  step <- function(fun, key, base = list()) {
    # a level outside `levels` is left out
    if (!(as.integer(sub("^t", "", key)) %in% run)) return(input)
    args <- qc_level_args(c(list(verbose = verbose), base), params[[key]], key)
    do.call(fun, c(list(input = input), args))
  }
  
  # the fixed order of the chain
  input <- step(T_QC_1_gross_error,           "t1")
  input <- step(T_QC_2_out_of_range,          "t2")
  input <- step(T_QC_3_time_consistency,      "t3")
  input <- step(T_QC_4_stuck_values,          "t4")
  input <- step(T_QC_5_climatic_outliers,     "t5")
  # the indoor test comes before the neighbour comparisons: those would thin
  # such a logger out and leave level 8 too little to judge
  input <- step(T_QC_8_diurnal_range,         "t8")
  # the spatial levels always run: without metadata each skips itself and says why
  input <- step(T_QC_6_spatial_consistency,        "t6", list(metadata = metadata))
  input <- step(T_QC_7_spatiotemporal_consistency, "t7", list(metadata = metadata))
  # the clock test comes last: sun on the sensor gives the same lag as a wrong
  # clock, and level 6 has to take those hours out first
  input <- step(T_QC_9_clock_offset,          "t9")
  # the interpolation alters the data: only on explicit request
  if (isTRUE(interpolate)) input <- step(T_QC_10_interpolate, "t10")
  
  # record what ran and what was left out, so the result says so even without the console
  skipped <- sort(unique(c(setdiff(seq_len(10), run), if (!isTRUE(interpolate)) 10)))
  input$qc_info$run_qc_temperature <- list(levels         = setdiff(run, skipped),
                                           levels_skipped = if (length(skipped)) paste0("t", skipped) else character(0))
  input
}

#' Run the whole humidity QC chain in the fixed order
#'
#' Convenience wrapper around levels RH_QC_1 to RH_QC_7 (and optionally 8).
#' Every level remains individually callable; this function only fixes the
#' order and wires the temperature result into the three levels that use it
#' (1, 5 and 7).
#'
#' `temperature` is ideally the RESULT LIST of the temperature chain: level 1
#' then inherits its flags, and levels 5 and 7 use its CLEANED series. A bare
#' temperature xts also works for levels 5 and 7, but carries no flags - level
#' 1 then runs, inherits nothing, and says so. Levels 5 and 7 always run and
#' skip themselves, with the reason in their record, when temperature or
#' metadata are missing.
#'
#' @param input xts of relative humidity (or an already started chain list).
#' @param temperature The temperature chain result list (preferred), or a bare
#'   temperature xts, or NULL.
#' @param metadata Data frame with ID, LAT, LON, Landuse for level 7, or NULL
#'   to skip it.
#' @param interpolate Also run the optional level 8 gap interpolation.
#' @param levels The levels allowed to run, as in `run_qc_temperature()`:
#'   `-5` runs everything but level 5, `c(1:4, 7)` only those five. Level 8
#'   still needs `interpolate = TRUE`.
#' @param params Named list of per-level parameter overrides, keyed rh1..rh8,
#'   e.g. `list(rh3 = list(threshold = 15))`.
#' @param verbose Passed to every level.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info; the
#'   levels run and skipped are listed in `qc_info$run_qc_humidity$levels`
#'   and `$levels_skipped`.
#'
#' @examples
#' \dontrun{
#' t_res  <- run_qc_temperature(t_xts, metadata = meta)
#' rh_res <- run_qc_humidity(rh_xts, temperature = t_res, metadata = meta)
#' }
#'
#' @import xts
#' @import zoo
#' @export
run_qc_humidity <- function(input,
                            temperature = NULL,
                            metadata = NULL,
                            interpolate = FALSE,
                            levels = 1:8,
                            params = list(),
                            verbose = TRUE) {
  # the overrides must be a named list keyed by level
  if (!is.list(params)) stop("params must be a named list, e.g. list(rh3 = list(threshold = 15)).")
  ok_keys <- paste0("rh", 1:8)
  if (length(params) && (is.null(names(params)) || !all(names(params) %in% ok_keys)))
    stop(sprintf("params keys must be among %s.", paste(ok_keys, collapse = ", ")))
  # the levels this run may touch; the rest is left out and recorded as skipped
  run <- qc_pick_levels(levels, 8, "rh")
  # the chain list right away, so the record below can be written even when
  # no level runs
  input <- qc_prepare_input(input, what = "humidity")
  
  # does the temperature argument carry flags (i.e. is it a chain result)?
  has_flags <- is.list(temperature) && !is.data.frame(temperature) &&
    !inherits(temperature, "xts") && "qc_data_flagged" %in% names(temperature)

  # one level call: base arguments, overridden by the caller's params entry
  step <- function(fun, key, base = list()) {
    # a level outside `levels` is left out
    if (!(as.integer(sub("^rh", "", key)) %in% run)) return(input)
    args <- qc_level_args(c(list(verbose = verbose), base), params[[key]], key)
    do.call(fun, c(list(input = input), args))
  }
  
  # level 1 always runs: it inherits from a chain result (flags + records); given
  # a bare temperature series it gets NULL and reports itself that nothing can
  # be inherited. The level, not the runner, says why.
  input <- step(RH_QC_1_inherit_temperature, "rh1",
                list(temperature_flags = if (has_flags) temperature else NULL))
  input <- step(RH_QC_2_range,            "rh2")
  input <- step(RH_QC_3_spike,            "rh3")
  input <- step(RH_QC_4_stuck_values,     "rh4")
  # levels 5 and 7 always run: given what is there, each skips itself and says why
  input <- step(RH_QC_5_decoupling, "rh5", list(temperature = temperature))
  input <- step(RH_QC_6_saturation_drift, "rh6")
  input <- step(RH_QC_7_dewpoint_consistency, "rh7",
                list(temperature = temperature, metadata = metadata))
  # the interpolation alters the data: only on explicit request
  if (isTRUE(interpolate)) input <- step(RH_QC_8_interpolate, "rh8")

  # record what ran and what was left out, so the result says so even without the console
  skipped <- sort(unique(c(setdiff(seq_len(8), run), if (!isTRUE(interpolate)) 8)))
  input$qc_info$run_qc_humidity <- list(levels         = setdiff(run, skipped),
                                        levels_skipped = if (length(skipped)) paste0("rh", skipped) else character(0))
  input
}

# the arguments of one level call: the runner's base arguments with the
# caller's params entry on top. An entry replaces the base value as given,
# NULL and whole tables included
qc_level_args <- function(base, over, key) {
  if (is.null(over)) return(base)
  if (!is.list(over) || is.data.frame(over) || (length(over) && is.null(names(over))))
    stop(sprintf("params$%s must be a named list of arguments.", key))
  for (nm in names(over)) base[nm] <- list(over[[nm]])
  base
}

# the levels a runner may run, from `levels`: numbers or names like "T3"/"t3",
# positive to select, negative to exclude, as in R indexing
qc_pick_levels <- function(levels, n, prefix) {
  if (is.character(levels)) {
    v <- suppressWarnings(as.numeric(sub(paste0("^(-?)", prefix, "([0-9]+)$"), "\\1\\2", levels, ignore.case = TRUE)))
    if (anyNA(v))
      stop(sprintf("levels must be numbers or names like \"%s3\" / \"-%s3\".", toupper(prefix), toupper(prefix)))
    levels <- v
  }
  if (!is.numeric(levels) || !length(levels) || anyNA(levels) || any(levels != round(levels)) ||
      any(levels == 0) || any(abs(levels) > n))
    stop(sprintf("levels must be whole numbers between 1 and %d, or their negatives.", n))
  if (any(levels > 0) && any(levels < 0)) stop("levels cannot mix positive and negative entries.")
  if (all(levels < 0)) setdiff(seq_len(n), -levels) else sort(unique(as.integer(levels)))
}
