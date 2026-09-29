#' Run the whole temperature QC chain in the fixed order
#'
#' Convenience wrapper around levels T_QC_1 to T_QC_8 (and optionally 9). It
#' exists to prevent the two classic application mistakes: running the levels
#' in the wrong order, and silently forgetting the spatial levels when the
#' metadata are at hand. Every level remains individually callable; this
#' function adds nothing the manual chain does not have.
#'
#' Without `metadata` the spatial levels 6 and 7 are SKIPPED audibly - useful
#' for a single isolated station, where they could never run anyway. The
#' interpolation is off by default, because it alters the data (see
#' `T_QC_9_interpolate`).
#'
#' @param input xts of temperature (or an already started chain list).
#' @param metadata Data frame with ID, LAT, LON, Landuse for levels 6 and 7,
#'   or NULL to skip both.
#' @param interpolate Also run the optional level 9 gap interpolation.
#' @param params Named list of per-level parameter overrides, keyed t1..t9,
#'   e.g. `list(t2 = list(summer_min = -10), t6 = list(k_sigma = 5))`. Every
#'   entry is passed to that level as arguments.
#' @param verbose Passed to every level.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info; skipped
#'   levels are listed in `qc_info$run_qc_temperature$levels_skipped`.
#'
#' @examples
#' \dontrun{
#' res <- run_qc_temperature(x, metadata = meta)
#' res <- run_qc_temperature(x, metadata = meta, interpolate = TRUE,
#'                           params = list(t2 = list(summer_min = -10)))
#' }
#'
#' @import xts
#' @import zoo
#' @export
run_qc_temperature <- function(input,
                               metadata = NULL,
                               interpolate = FALSE,
                               params = list(),
                               verbose = TRUE) {
  # the overrides must be a named list keyed by level
  if (!is.list(params)) stop("params must be a named list, e.g. list(t2 = list(summer_min = -10)).")
  ok_keys <- paste0("t", 1:9)
  if (length(params) && (is.null(names(params)) || !all(names(params) %in% ok_keys)))
    stop(sprintf("params keys must be among %s.", paste(ok_keys, collapse = ", ")))
  
  # one level call: base arguments, overridden by the caller's params entry
  step <- function(fun, key, base = list()) {
    args <- utils::modifyList(c(list(verbose = verbose), base),
                              if (is.null(params[[key]])) list() else params[[key]])
    do.call(fun, c(list(input = input), args))
  }
  
  # the fixed order of the chain
  input <- step(T_QC_1_gross_error,           "t1")
  input <- step(T_QC_2_out_of_range,          "t2")
  input <- step(T_QC_3_time_consistency,      "t3")
  input <- step(T_QC_4_stuck_values,          "t4")
  input <- step(T_QC_5_climatic_outliers,     "t5")
  # the spatial levels need the metadata; without them, skip audibly
  if (!is.null(metadata)) {
    input <- step(T_QC_6_spatial_consistency,        "t6", list(metadata = metadata))
    input <- step(T_QC_7_spatiotemporal_consistency, "t7", list(metadata = metadata))
  } else if (isTRUE(verbose)) {
    message("run_qc_temperature: no metadata supplied - levels 6 (spatial) and 7 (spatiotemporal) skipped.")
  }
  input <- step(T_QC_8_diurnal_range,         "t8")
  # the interpolation alters the data: only on explicit request
  if (isTRUE(interpolate)) input <- step(T_QC_9_interpolate, "t9")
  
  # record what was skipped, so the result says so even without the console
  skipped <- c(if (is.null(metadata)) c("t6", "t7"), if (!isTRUE(interpolate)) "t9")
  input$qc_info$run_qc_temperature <- list(levels_skipped = if (length(skipped)) skipped else character(0))
  input
}

#' Run the whole humidity QC chain in the fixed order
#'
#' Convenience wrapper around levels RH_QC_1 to RH_QC_7 (and optionally 8).
#' Every level remains individually callable; this function only fixes the
#' order and wires the temperature result into the three levels that use it.
#'
#' `temperature` is ideally the RESULT LIST of the temperature chain: level 1
#' then inherits its flags, and levels 6 and 7 use its CLEANED series. A bare
#' temperature xts also works for levels 6 and 7, but carries no flags - level
#' 1 then runs, inherits nothing, and says so. Levels 6 and 7 always run and
#' skip themselves, with the reason in their record, when temperature or
#' metadata are missing.
#'
#' @param input xts of relative humidity (or an already started chain list).
#' @param temperature The temperature chain result list (preferred), or a bare
#'   temperature xts, or NULL.
#' @param metadata Data frame with ID, LAT, LON, Landuse for level 7, or NULL
#'   to skip it.
#' @param interpolate Also run the optional level 8 gap interpolation.
#' @param params Named list of per-level parameter overrides, keyed rh1..rh8,
#'   e.g. `list(rh3 = list(threshold = 15))`.
#' @param verbose Passed to every level.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info; skipped
#'   levels are listed in `qc_info$run_qc_humidity$levels_skipped`.
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
                            params = list(),
                            verbose = TRUE) {
  # the overrides must be a named list keyed by level
  if (!is.list(params)) stop("params must be a named list, e.g. list(rh3 = list(threshold = 15)).")
  ok_keys <- paste0("rh", 1:8)
  if (length(params) && (is.null(names(params)) || !all(names(params) %in% ok_keys)))
    stop(sprintf("params keys must be among %s.", paste(ok_keys, collapse = ", ")))
  
  # does the temperature argument carry flags (i.e. is it a chain result)?
  has_flags <- is.list(temperature) && !is.data.frame(temperature) &&
    !inherits(temperature, "xts") && "qc_data_flagged" %in% names(temperature)

  # one level call: base arguments, overridden by the caller's params entry
  step <- function(fun, key, base = list()) {
    args <- utils::modifyList(c(list(verbose = verbose), base),
                              if (is.null(params[[key]])) list() else params[[key]])
    do.call(fun, c(list(input = input), args))
  }
  
  # level 1 always runs: it inherits from a chain result (flags + records); given
  # a bare temperature series it gets NULL and reports itself that nothing can
  # be inherited. The level, not the runner, says why.
  input <- step(RH_QC_1_inherit_temperature, "rh1",
                list(temperature_flags = if (has_flags) temperature else NULL))
  input <- step(RH_QC_2_range,            "rh2")
  input <- step(RH_QC_3_spike,            "rh3")
  input <- step(RH_QC_4_saturation_drift, "rh4")
  input <- step(RH_QC_5_stuck_values,     "rh5")
  # levels 6 and 7 always run: given what is there, each skips itself and says why
  input <- step(RH_QC_6_decoupling, "rh6", list(temperature = temperature))
  input <- step(RH_QC_7_dewpoint_consistency, "rh7",
                list(temperature = temperature, metadata = metadata))
  # the interpolation alters the data: only on explicit request
  if (isTRUE(interpolate)) input <- step(RH_QC_8_interpolate, "rh8")

  # record what was skipped, so the result says so even without the console
  skipped <- c(if (!isTRUE(interpolate)) "rh8")
  input$qc_info$run_qc_humidity <- list(levels_skipped = if (length(skipped)) skipped else character(0))
  input
}
