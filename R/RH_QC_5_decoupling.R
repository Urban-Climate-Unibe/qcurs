#' Humidity QC Level 5: decoupling from the station's own temperature
#'
#' Over a day with a normal diurnal cycle, relative humidity tracks its own
#' temperature closely and negatively. A humidity channel that stops following
#' its own temperature is broken even when its values look plausible - the
#' typical signature of a dead or detached element. This is the humidity twin
#' of temperature level 8: it needs no neighbours and therefore also covers
#' stations the spatial tests cannot reach.
#'
#' The temperature is matched cell by cell, by exact time stamp and logger
#' name (`qc_match_grid()`, the same rule as levels 1 and 7); a cell without
#' a temperature twin is not judged. Days with less than `decouple_trange` of
#' temperature amplitude are skipped, because without a diurnal cycle the
#' correlation is meaningless - fog days would otherwise be flagged. The
#' verdict is per day: all valid humidity values of a decoupled day are
#' objected to. Without any usable temperature the level says so, records
#' it, and hands the pair back unchanged.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature The temperature run: its result list (preferred; its
#'   CLEANED `qc_data` is taken) or a temperature xts with logger names.
#'   NULL (default) skips the level.
#' @param decouple_r Correlation at or above which a day counts as decoupled
#'   (healthy days are strongly negative).
#' @param decouple_trange Minimum diurnal temperature range in Kelvin.
#' @param decouple_min_n Minimum valid T/RH pairs per day.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. Stations
#'   without a temperature twin are listed in
#'   `qc_info$rh5_decoupling$skipped_stations`.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_5_decoupling(res, temperature = t_res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_5_decoupling <- function(input,
                               temperature = NULL,
                               decouple_r = -0.3,
                               decouple_trange = 3,
                               decouple_min_n = 60,
                               verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh5_decoupling")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the decoupling threshold is a correlation
  if (!is.numeric(decouple_r) || decouple_r <= -1 || decouple_r >= 1)
    stop("decouple_r must be a correlation strictly between -1 and 1.")
  # a flat day carries no signal
  if (!is.numeric(decouple_trange) || decouple_trange <= 0)
    stop("decouple_trange must be a positive range in Kelvin.")
  # a correlation from a handful of pairs is noise
  if (!is.numeric(decouple_min_n) || decouple_min_n < 3)
    stop("decouple_min_n must be at least 3.")

  #-------------------------------------------------------------------------------
  # the temperature twin. Anything unusable is not an error but a skip: say
  # why, record it, hand the pair back unchanged, let the chain go on

  skip <- function(reason) {
    if (isTRUE(verbose))
      message("RH5 decoupling: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$rh5_decoupling <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    input
  }
  if (is.null(temperature)) return(skip("no temperature supplied"))
  # the whole chain result: take its CLEANED series
  tx <- temperature
  if (is.list(tx) && !inherits(tx, "xts")) {
    if (!"qc_data" %in% names(tx)) return(skip("temperature list has no 'qc_data'"))
    tx <- tx$qc_data
  }
  # temperature on the humidity grid, by exact time stamp and logger name
  tm <- qc_match_grid(x, tx)
  if (is.null(tm)) return(skip("temperature has no logger name or time stamp in common with the humidity"))

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 5

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # the matched temperatures, NA wherever there is no twin
  TT <- tm$M
  # calendar day of every time step, in the time zone of the index (format()
  # honours it; as.Date() only does so from R 4.3 on)
  day <- format(zoo::index(x), "%Y-%m-%d")
  # the row numbers of every day, once, in the order of the record
  day_rows <- split(seq_along(day), factor(day, levels = unique(day)))
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: T/RH pairs on days that got a verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))
  # humidity stations without a temperature twin cannot be tested
  skipped_stations <- setdiff(colnames(X), tm$loggers)

  # iterate over all stations (columns) that have a temperature twin
  for (s in tm$loggers) {
    # extract humidity vector of this station and its own temperature
    v <- X[, s]; tv <- TT[, s]
    # verdict per cell, FALSE until proven otherwise
    hit <- rep(FALSE, length(v))
    # judge day by day: the diurnal cycle is the yardstick
    for (rows in day_rows) {
      # the day's valid T/RH pairs
      k <- rows[!is.na(tv[rows]) & !is.na(v[rows])]
      # too sparse: no verdict for this day
      if (length(k) < decouple_min_n) next
      # flat day: the correlation is meaningless, skip instead of flagging fog
      if (diff(range(tv[k])) < decouple_trange) next
      # the daily T-RH correlation
      r <- suppressWarnings(stats::cor(tv[k], v[k]))
      # undefined (a constant series): no verdict
      if (is.na(r)) next
      # a verdict is reached, whichever way
      n_judged[s] <- n_judged[s] + length(k)
      # healthy negative coupling: fine
      if (r < decouple_r) next
      # decoupled: mark every valid humidity value of this day
      hit[k] <- TRUE
    }

    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the decoupled days so later levels never see them
      X[mask, s] <- NA
      # record this level's code (5 = level 5, fixed by convention)
      previous_flag[mask, s] <- 5
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

  # report so a zero-hit run is visibly a run, not a skip - and name the blind spots
  if (isTRUE(verbose)) {
    message(sprintf("RH5 decoupled from own temperature (r >= %g): %d flagged, %d of %d loggers and %d of %d time steps with temperature",
                    decouple_r, n_total, length(tm$loggers), ncol(X), tm$n_time, nrow(X)))
    if (length(skipped_stations) > 0)
      message("  skipped (no temperature twin): ", paste(skipped_stations, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                <- x
  input$qc_data_flagged        <- flg
  input$qc_info$rh5_decoupling <- list(n_flagged               = n_total,
                                       n_flagged_by_station    = n_station,
                                       n_judged_by_station     = n_judged,
                                       decouple_r              = decouple_r,
                                       decouple_trange         = decouple_trange,
                                       decouple_min_n          = decouple_min_n,
                                       loggers_with_temperature = tm$loggers,
                                       n_time_with_temperature  = tm$n_time,
                                       skipped_stations        = skipped_stations)
  input
}
