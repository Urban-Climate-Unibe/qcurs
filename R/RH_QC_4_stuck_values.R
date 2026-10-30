#' Humidity QC Level 4: stuck values, with a saturation exception
#'
#' Same search as temperature level 4, via the shared `qc_find_stuck()`: a
#' window with no variation marks a stuck sensor - EXCEPT when the window
#' median sits at saturation, because in fog and continuous rain the humidity
#' genuinely stops moving; without this exception the test flags every winter
#' systematically. A window is judged only if enough of it holds valid values.
#'
#' Saturation is not 95 percent for every sensor. The low-cost elements of
#' the Thun network stand at 73 to 90 percent - each at its own ceiling -
#' while MeteoSwiss Thun reads 99, and in a dew night a forest station sits
#' flat at 85 for ten hours: the sensor at its ceiling, not a stuck sensor.
#' So a window is also exempt when its median reaches the station's own
#' `sat_quantile` quantile of all its valid values (the 90th: the plateau
#' nights ARE the top of the record, and a stricter quantile cuts through
#' them); the ceiling that resulted is recorded per station
#' (`sat_ceiling`). A sensor stuck at ONE
#' value for the whole record has that value as its ceiling and escapes
#' here; level 5 catches it, because such a sensor no longer follows its own
#' temperature. `sat_quantile = NULL` leaves only the fixed `sat_max`.
#'
#' The window is also the detection delay (see `T_QC_4_stuck_values()`).
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param window_size Window length, either a `difftime` or a string
#'   "<n> <unit>" with a unit `as.difftime()` understands (secs, mins, hours,
#'   days, weeks).
#' @param sd_tol Variation at or below this counts as "no variation". Unlike
#'   the temperature twin the default is not 0: capacitive elements quantise
#'   more coarsely.
#' @param sat_max Windows whose median is at or above this are exempt (the
#'   saturation exception).
#' @param sat_quantile Windows whose median reaches this quantile of the
#'   station's own valid values are exempt as well - the station's ceiling.
#'   NULL switches the station rule off.
#' @param min_valid_frac Minimum fraction of the window that must hold valid
#'   values for the window to be judged.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. The
#'   exemption threshold that applied at each station - the smaller of
#'   `sat_max` and the station's quantile - is in
#'   `qc_info$rh4_stuck_values$sat_ceiling`.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_4_stuck_values(res, window_size = "6 hours")
#' res$qc_info$rh4_stuck_values$sat_ceiling   # where each sensor saturates
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_4_stuck_values <- function(input,
                                 window_size = "6 hours",
                                 sd_tol = 0.1,
                                 sat_max = 95,
                                 sat_quantile = 0.9,
                                 min_valid_frac = 0.5,
                                 verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh4_stuck_values")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # a negative tolerance would switch the test off silently
  if (!is.numeric(sd_tol) || length(sd_tol) != 1 || sd_tol < 0)
    stop("sd_tol must be zero or a positive variation.")
  # the saturation exemption must lie inside the physical range
  if (!is.numeric(sat_max) || length(sat_max) != 1 || sat_max <= 0 || sat_max > 100)
    stop("sat_max must lie in (0, 100].")
  # the station rule needs a quantile that is a quantile; NULL means off
  if (!is.null(sat_quantile) && (!is.numeric(sat_quantile) || length(sat_quantile) != 1 ||
                                 sat_quantile <= 0 || sat_quantile >= 1))
    stop("sat_quantile must lie in (0, 1), or be NULL.")
  # the population guard must be a usable fraction
  if (!is.numeric(min_valid_frac) || length(min_valid_frac) != 1 || min_valid_frac <= 0 || min_valid_frac > 1)
    stop("min_valid_frac must be a fraction in (0, 1].")

  #-------------------------------------------------------------------------------
  # the window in points, and the population guard as a minimum count

  width <- qc_window_points(window_size, input$qc_info$dataset_humidity$time_step_sec, "window_size")
  min_valid <- max(2L, ceiling(min_valid_frac * width))

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 4

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: cells inside at least one judged window
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))
  # the exemption threshold that applied at each station
  sat_ceiling <- stats::setNames(rep(sat_max, ncol(X)), colnames(X))

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # the station's own ceiling: a window at the top of its range is
    # saturation for THIS sensor, whatever the fixed limit says
    if (!is.null(sat_quantile) && any(!is.na(v)))
      sat_ceiling[s] <- min(sat_max, stats::quantile(v, sat_quantile, na.rm = TRUE, names = FALSE))
    # the search itself lives in qc_find_stuck(), shared with T_QC_4; the
    # saturation exception is the one thing this chain adds
    st <- qc_find_stuck(v, width, sd_tol, min_valid, exempt_above = sat_ceiling[[s]])
    n_judged[s] <- sum(st$judged & !is.na(v))

    # combine the verdict with THIS station's column only
    mask <- st$hit & !is.na(v) & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the stuck stretch so later levels never see it
      X[mask, s] <- NA
      # record this level's code (4 = level 4, fixed by convention)
      previous_flag[mask, s] <- 4
      # add the number of new flags to the counters
      n_station[s] <- n_found
      n_total <- n_total + n_found
    }
  }

  # write the matrices back into the xts shells, keeping index and column names
  if (n_total > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("RH4 stuck values (%s = %d points, min valid %d, sd <= %g, exempt if median >= %g or at the station's own q%g): %d flagged",
                    if (inherits(window_size, "difftime")) format(window_size) else window_size,
                    width, min_valid, sd_tol, sat_max, if (is.null(sat_quantile)) NA else sat_quantile * 100, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                  <- x
  input$qc_data_flagged          <- flg
  input$qc_info$rh4_stuck_values <- list(n_flagged            = n_total,
                                         n_flagged_by_station = n_station,
                                         n_judged_by_station  = n_judged,
                                         window_size          = window_size,
                                         width_points         = width,
                                         min_valid            = min_valid,
                                         sd_tol               = sd_tol,
                                         sat_max              = sat_max,
                                         sat_quantile         = sat_quantile,
                                         sat_ceiling          = sat_ceiling,
                                         min_valid_frac       = min_valid_frac)
  input
}
