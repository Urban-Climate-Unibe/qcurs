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
#' correlation is meaningless - fog days would otherwise be flagged. A day on
#' which the humidity does not move at all while the temperature swings counts
#' as decoupled.
#'
#' The correlation also breaks down on days with an air mass change: a warm
#' front with rain lifts temperature and humidity together, cold advection
#' drops both. On healthy MeteoSwiss stations the plain daily test objects
#' to 5 percent of all days, 13 percent in winter. Two rules separate weather
#' from a sensor. A day on which more than `max_network_frac` of the judged
#' stations (and at least two) look decoupled is a weather day and is judged
#' at no station; it is listed in `days_unjudged`. And a station is objected
#' to only after `min_consecutive_days` days in a row: a front lasts a
#' day, a dead element lasts weeks. A single station without a network can
#' only use the second rule, and fog that lasts gets past it: on ten years of
#' quality-checked MeteoSwiss series, each station run alone, the level
#' objected to 28 to 51 days per station, most of them from October to
#' February, and a reference sensor standing flat at 99 percent for three
#' days among low-cost elements that kept moving was objected to as well.
#' The verdict is per day: all valid humidity values of an objected day are
#' removed. Without any usable temperature the level says so, records it, and
#' hands the pair back unchanged.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature The temperature run: its result list (preferred; its
#'   CLEANED `qc_data` is taken) or a temperature xts with logger names.
#'   NULL (default) skips the level.
#' @param decouple_r Correlation at or above which a day counts as decoupled
#'   (healthy days are strongly negative).
#' @param decouple_trange Minimum diurnal temperature range in Kelvin.
#' @param decouple_min_n Minimum valid T/RH pairs per day; NULL (default)
#'   takes ten hours of values at the time step of the data (60 at 10 minutes).
#' @param min_consecutive_days Days in a row a station must look
#'   decoupled before any of them is objected to.
#' @param max_network_frac Share of the judged stations above which a day
#'   counts as a weather day and is not judged at all (needs at least two
#'   decoupled stations). 1 switches the rule off.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. Stations
#'   without a temperature twin are listed in
#'   `qc_info$rh5_decoupling$never_evaluated`, the weather days in
#'   `qc_info$rh5_decoupling$days_unjudged`.
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
                               decouple_min_n = NULL,
                               min_consecutive_days = 2,
                               max_network_frac = 0.25,
                               verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh5_decoupling")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  if (!is.numeric(decouple_r) || decouple_r <= -1 || decouple_r >= 1)
    stop("decouple_r must be a correlation strictly between -1 and 1.")
  # a flat day carries no signal
  if (!is.numeric(decouple_trange) || decouple_trange <= 0)
    stop("decouple_trange must be a positive range in Kelvin.")
  # the pairs a day needs default to ten hours of values, whatever the time step
  if (is.null(decouple_min_n)) {
    step_sec <- input$qc_info$dataset_humidity$time_step_sec
    decouple_min_n <- if (is.finite(step_sec) && step_sec > 0) max(3, round(36000 / step_sec)) else 60
  }
  # a correlation from a handful of pairs is noise
  if (!is.numeric(decouple_min_n) || decouple_min_n < 3)
    stop("decouple_min_n must be at least 3.")
  if (!is.numeric(min_consecutive_days) || length(min_consecutive_days) != 1 ||
      min_consecutive_days < 1 || min_consecutive_days != round(min_consecutive_days))
    stop("min_consecutive_days must be a whole number, at least 1.")
  # 1 switches the network rule off
  if (!is.numeric(max_network_frac) || max_network_frac <= 0 || max_network_frac > 1)
    stop("max_network_frac must lie in (0, 1].")

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
  # the whole chain result: take its cleaned series
  tx <- temperature
  if (is.list(tx) && !inherits(tx, "xts")) {
    if (!"qc_data" %in% names(tx)) return(skip("temperature list has no 'qc_data'"))
    tx <- tx$qc_data
  }
  # temperature on the humidity grid, by exact time stamp and logger name
  tm <- qc_match_grid(x, tx)
  if (is.null(tm)) return(skip("temperature has no logger name or time stamp in common with the humidity"))

  #-------------------------------------------------------------------------------
  # perform RH QC level 5

  # plain matrices of values and flags (time in rows, stations in columns)
  X <- coredata(x)
  previous_flag <- coredata(flg)
  # the matched temperatures, NA wherever there is no twin
  TT <- tm$M
  # calendar day of every time step, in the time zone of the index (format()
  # honours it; as.Date() only does so from R 4.3 on)
  day <- format(zoo::index(x), "%Y-%m-%d")
  # the row numbers of every day, once, in the order of the record
  day_rows <- split(seq_along(day), factor(day, levels = unique(day)))
  n_total <- 0
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: valid values on days that got a verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))
  # humidity stations without a temperature twin cannot be tested
  never_evaluated <- setdiff(colnames(X), tm$loggers)

  # the daily verdict per station: NA = no verdict, TRUE = looks decoupled
  V <- matrix(NA, length(day_rows), ncol(X), dimnames = list(names(day_rows), colnames(X)))
  # the valid humidity values of each judged day, for the coverage count
  NP <- matrix(0L, length(day_rows), ncol(X), dimnames = dimnames(V))

  # first pass: judge every day of every station that has a temperature twin
  for (s in tm$loggers) {
    # this station's humidity and its own temperature
    v <- X[, s]; tv <- TT[, s]
    for (d in seq_along(day_rows)) {
      rows <- day_rows[[d]]
      # the day's valid T/RH pairs
      k <- rows[!is.na(tv[rows]) & !is.na(v[rows])]
      # too sparse: no verdict for this day
      if (length(k) < decouple_min_n) next
      # flat day: the correlation is meaningless, skip instead of flagging fog
      if (diff(range(tv[k])) < decouple_trange) next
      r <- suppressWarnings(stats::cor(tv[k], v[k]))
      # a verdict is reached, whichever way. Undefined means the humidity did
      # not move at all while the temperature swung: the clearest decoupling
      V[d, s] <- is.na(r) || r >= decouple_r
      NP[d, s] <- sum(!is.na(v[rows]))
    }
  }

  # the network rule: when more than max_network_frac of the judged stations
  # (and at least two) look decoupled on the same day, it is the weather -
  # a front, not a sensor. Nobody is judged on such a day
  n_jud <- rowSums(!is.na(V))
  n_dec <- rowSums(V, na.rm = TRUE)
  weather <- n_dec >= 2 & n_dec / pmax(n_jud, 1) > max_network_frac
  V[weather, ] <- NA
  NP[weather, ] <- 0L
  # the weather days, so a zero-flag run is not a clean bill
  days_unjudged <- names(day_rows)[weather]

  # second pass: per station, only runs of decoupled days long enough count
  for (s in tm$loggers) {
    # every value on a judged day counts as judged
    n_judged[s] <- sum(NP[, s])
    # decoupled days; an unjudged day in between breaks the run
    dec <- !is.na(V[, s]) & V[, s]
    # run-length encode to find consecutive stretches
    r <- rle(dec)
    # end and start index of every run
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    # collect the day indices of every TRUE run that is long enough
    flag_days <- integer(0)
    for (k in seq_along(r$lengths))
      if (r$values[k] && r$lengths[k] >= min_consecutive_days)
        flag_days <- c(flag_days, starts[k]:ends[k])
    # nothing decoupled long enough at this station
    if (!length(flag_days)) next
    # every valid humidity value on the objected days
    v <- X[, s]
    hit <- rep(FALSE, length(v))
    hit[unlist(day_rows[flag_days])] <- TRUE
    hit <- hit & !is.na(v)

    # only cells that carry no earlier flag
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    n_found <- sum(mask)

    if (n_found > 0) {
      # blank the decoupled days so later levels never see them
      X[mask, s] <- NA
      # record this level's code (5 = level 5, fixed by convention)
      previous_flag[mask, s] <- 5
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
    message(sprintf("RH5 decoupled from own temperature (r >= %g, >= %d days in a row): %d flagged, %d weather day(s) not judged, %d of %d loggers and %d of %d time steps with temperature",
                    decouple_r, min_consecutive_days, n_total, length(days_unjudged),
                    length(tm$loggers), ncol(X), tm$n_time, nrow(X)))
    if (length(never_evaluated) > 0)
      message("  never evaluated (no temperature twin): ", paste(never_evaluated, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                <- x
  input$qc_data_flagged        <- flg
  input$qc_info$rh5_decoupling <- list(n_flagged                = n_total,
                                       n_flagged_by_station     = n_station,
                                       n_judged_by_station      = n_judged,
                                       decouple_r               = decouple_r,
                                       decouple_trange          = decouple_trange,
                                       decouple_min_n           = decouple_min_n,
                                       min_consecutive_days     = min_consecutive_days,
                                       max_network_frac         = max_network_frac,
                                       loggers_with_temperature = tm$loggers,
                                       n_time_with_temperature  = tm$n_time,
                                       never_evaluated          = never_evaluated,
                                       days_unjudged            = days_unjudged)
  input
}
