#' Temperature QC Level 10: clock offset against the network
#'
#' A logger whose clock is wrong - set to the wrong time zone, not adjusted
#' after a battery change, drifting for a year - delivers correct values
#' under wrong time stamps. Every single-value test passes; the station just
#' warms and cools an hour before or after everybody else, and level 6
#' objects to a few cells on the steep flanks without saying why.
#'
#' Per station and window of `window_days` days, the station's series is
#' cross-correlated with the network median series over lags of up to
#' `max_lag` steps either way; the lag with the highest correlation, refined
#' to a fraction of a step from the shape of the peak, is the window's
#' verdict. The lag counts the steps by which the station's stamps are
#' EARLY: a positive lag means the station shows the network's warming that
#' many steps before the network does, so its clock runs behind. When the
#' best lag rounds to at least `min_lag` steps, with the same sign, in
#' `min_consecutive_windows` judged windows in a row, every valid value of
#' those windows is objected to. The `episodes` table carries the lag in
#' minutes. The data are not shifted - that repair is for whoever knows when
#' the clock was set.
#'
#' A window is judged only when at least `min_valid_frac` of its steps have
#' a value at the station and in the network, and when the correlation at
#' the best lag reaches `min_cor`: a series that does not resemble the
#' network has no clock to speak of (a logger indoors is level 8's case).
#' Sites lag too, and by more than one would think: on the Thun network a
#' forest station runs 20 to 50 minutes behind the network (55 in one sunny
#' week), a station at the lake up to 40, the open reference site 10 to 30
#' minutes ahead - in both summers, and varying with the weather. The default
#' `min_lag` is therefore 6 steps, one hour at 10-minute data: the smallest
#' error a wrong time zone or a missed clock change makes, and above every
#' site lag seen. A planted hour at a Thun station was found in every full
#' week from its first to its last. The lag table in the record shows the
#' site lags all the same. A best lag at the edge of the search range is
#' reported as `max_lag` and means "at least that".
#'
#' The level needs no metadata, only `min_stations` stations with values on
#' the same steps. The yardstick of a station is the median of the OTHER
#' stations, so an offset station does not pull its own reference along.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param window_days Length of the windows in calendar days.
#' @param max_lag Largest lag searched, in steps either way (12 = 2 hours
#'   at 10-minute data).
#' @param min_lag Smallest best lag, in steps, that counts as a clock
#'   offset (6 = one hour at 10-minute data).
#' @param min_consecutive_windows Judged windows in a row with a lag of the
#'   same sign at or beyond `min_lag` before any of them is objected to.
#' @param min_valid_frac Share of a window's steps that need a value at the
#'   station and in the network for the window to be judged.
#' @param min_cor Correlation at the best lag below which a window gets no
#'   verdict.
#' @param min_stations Stations with a value a step needs before it can be
#'   judged; with fewer stations than that in the data the level skips
#'   itself and says why.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. The best
#'   lag of every judged station-window, in steps, is in
#'   `qc_info$t10_clock_offset$lags` (windows in rows, stations in columns,
#'   NA = not judged), the objected runs in `episodes` (station, from, to,
#'   windows, lag_min, r).
#'
#' @examples
#' \dontrun{
#' res <- T_QC_10_clock_offset(res)
#' res$qc_info$t10_clock_offset$episodes
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_10_clock_offset <- function(input,
                                 window_days = 7,
                                 max_lag = 12,
                                 min_lag = 6,
                                 min_consecutive_windows = 2,
                                 min_valid_frac = 0.5,
                                 min_cor = 0.7,
                                 min_stations = 3,
                                 verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t10_clock_offset")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # a window is a positive whole number of days
  if (!is.numeric(window_days) || length(window_days) != 1 || window_days < 1 || window_days != round(window_days))
    stop("window_days must be a positive whole number of days.")
  # the lags are whole numbers of steps, and the smallest must fit in the range
  if (!is.numeric(max_lag) || length(max_lag) != 1 || max_lag < 1 || max_lag != round(max_lag))
    stop("max_lag must be a positive whole number of steps.")
  if (!is.numeric(min_lag) || length(min_lag) != 1 || min_lag < 1 || min_lag > max_lag || min_lag != round(min_lag))
    stop("min_lag must be a whole number of steps between 1 and max_lag.")
  # at least one window in a row must be required
  if (!is.numeric(min_consecutive_windows) || min_consecutive_windows < 1)
    stop("min_consecutive_windows must be at least 1.")
  # the coverage is a fraction
  if (!is.numeric(min_valid_frac) || length(min_valid_frac) != 1 || min_valid_frac <= 0 || min_valid_frac > 1)
    stop("min_valid_frac must lie in (0, 1].")
  # the correlation floor is a correlation
  if (!is.numeric(min_cor) || length(min_cor) != 1 || min_cor < 0 || min_cor >= 1)
    stop("min_cor must lie in [0, 1).")
  # a network median of fewer than three stations is not a yardstick
  if (!is.numeric(min_stations) || min_stations < 3)
    stop("min_stations must be at least 3.")
  # too few stations, or no time step to count lags in: not an error but a
  # skip - say why, record it, hand the pair back
  step_sec <- input$qc_info$dataset_temperature$time_step_sec
  reason <- NULL
  if (ncol(input$qc_data) < min_stations)
    reason <- sprintf("only %d station(s); need at least %d for a network median", ncol(input$qc_data), min_stations)
  else if (!is.finite(step_sec) || step_sec <= 0)
    reason <- "no time step could be determined (single time stamp)"
  if (!is.null(reason)) {
    if (isTRUE(verbose))
      message("T10 clock offset: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$t10_clock_offset <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    return(input)
  }

  #-------------------------------------------------------------------------------
  # Perform QC Level 10

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # the time stamps as seconds, for the lag matching
  tt <- as.numeric(zoo::index(x))
  # calendar day of every time step, in the time zone of the index
  day <- format(zoo::index(x), "%Y-%m-%d")
  # the window of every time step: blocks of window_days calendar days
  win <- (match(day, unique(day)) - 1) %/% window_days + 1
  # the row numbers of every window, once, in the order of the record
  win_rows <- split(seq_along(win), win)
  # first and last day of every window, for the report
  win_from <- vapply(win_rows, function(rows) day[rows[1]], character(1))
  win_to   <- vapply(win_rows, function(rows) day[rows[length(rows)]], character(1))
  # the steps a full window holds
  full_n <- window_days * 86400 / step_sec
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: values in windows that got a verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))

  # a step is usable when enough stations have a value
  step_ok <- rowSums(!is.na(X)) >= min_stations
  # the lags searched, in steps
  lags <- seq(-max_lag, max_lag)
  # for every lag, the row that lies that many steps later - by time stamp,
  # not by row number, so a missing row does not bend the lag
  shifted <- lapply(lags, function(L) match(tt + L * step_sec, tt))
  # best lag and its correlation per window and station, NA = not judged
  LAG <- matrix(NA_real_, length(win_rows), ncol(X), dimnames = list(names(win_rows), colnames(X)))
  COR <- LAG

  # first pass: the best lag of every window at every station
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # the yardstick: the median of the OTHER stations at every step, NA
    # where too few stations report. Leaving the station out matters in a
    # small network, where an offset station would otherwise pull the
    # median's timing along
    net <- qc_row_median(X[, colnames(X) != s, drop = FALSE])
    net[!step_ok] <- NA_real_
    for (w in seq_along(win_rows)) {
      rows <- win_rows[[w]]
      # too sparse at the station or in the network: no verdict for this window
      if (sum(!is.na(v[rows]) & !is.na(net[rows])) < min_valid_frac * full_n) next
      # the correlation of the station with the network moved by every lag
      rr <- vapply(seq_along(lags), function(i) {
        b <- net[shifted[[i]][rows]]
        a <- v[rows]
        k <- !is.na(a) & !is.na(b)
        # too few pairs, or a flat series: no correlation
        if (sum(k) < 3) return(NA_real_)
        suppressWarnings(stats::cor(a[k], b[k]))
      }, numeric(1))
      if (all(is.na(rr))) next
      # the lag at which the station lines up best with the network
      best <- which.max(rr)
      # the station does not resemble the network: not a clock question
      if (rr[best] < min_cor) next
      # refine the peak to a fraction of a step by the parabola through its
      # three points, so a true offset of 6 steps does not read 5 in one
      # window and 6 in the next; at the edge of the range there is no parabola
      lag_est <- lags[best]
      if (best > 1 && best < length(lags) && !anyNA(rr[c(best - 1, best + 1)])) {
        curv <- rr[best - 1] - 2 * rr[best] + rr[best + 1]
        if (curv < 0) lag_est <- lag_est + max(-0.5, min(0.5, 0.5 * (rr[best - 1] - rr[best + 1]) / curv))
      }
      LAG[w, s] <- round(lag_est, 2)
      COR[w, s] <- rr[best]
    }
  }

  # the objected runs, one row each
  episodes <- list()

  # second pass: per station, only runs of offset windows long enough count
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # the judged windows of this station, in the order of the record
    jw <- which(!is.na(LAG[, s]))
    # every value in a judged window counts as judged
    n_judged[s] <- sum(!is.na(v[unlist(win_rows[jw])]))
    if (!length(jw)) next
    # the verdict per judged window: the sign of an offset that rounds to
    # min_lag or beyond, 0 otherwise - a run must keep its sign
    verdict <- ifelse(abs(round(LAG[jw, s])) >= min_lag, sign(LAG[jw, s]), 0)
    # run-length encode over the JUDGED windows only
    r <- rle(verdict)
    # end index of every run
    ends <- cumsum(r$lengths)
    # start index of every run
    starts <- ends - r$lengths + 1
    # collect the windows of every offset run that is long enough
    flag_wins <- integer(0)
    for (k in seq_along(r$lengths)) {
      if (r$values[k] == 0 || r$lengths[k] < min_consecutive_windows) next
      w_idx <- jw[starts[k]:ends[k]]
      flag_wins <- c(flag_wins, w_idx)
      # the run for the report: where, how long, how far off, how clearly
      episodes[[length(episodes) + 1]] <-
        data.frame(station = s,
                   from    = win_from[w_idx[1]],
                   to      = win_to[w_idx[length(w_idx)]],
                   windows = length(w_idx),
                   lag_min = round(stats::median(LAG[w_idx, s]) * step_sec / 60),
                   r       = round(stats::median(COR[w_idx, s]), 3),
                   stringsAsFactors = FALSE)
    }
    # no clock offset long enough at this station
    if (!length(flag_wins)) next
    # every valid value in the objected windows
    hit <- rep(FALSE, length(v))
    hit[unlist(win_rows[flag_wins])] <- TRUE
    hit <- hit & !is.na(v)

    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the offset windows so later levels never see them
      X[mask, s] <- NA
      # record this level's code (10 = level 10, fixed by convention)
      previous_flag[mask, s] <- 10
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
  # the episodes as one table, empty when nothing was objected to
  episodes <- if (length(episodes)) do.call(rbind, episodes) else
    data.frame(station = character(0), from = character(0), to = character(0),
               windows = integer(0), lag_min = numeric(0), r = numeric(0), stringsAsFactors = FALSE)

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose)) {
    message(sprintf("T10 clock offset (|lag| >= %d steps = %g min, same sign over >= %d windows of %d days): %d flagged, %d episode(s), %d of %d station-windows judged",
                    min_lag, min_lag * step_sec / 60, min_consecutive_windows, window_days,
                    n_total, nrow(episodes), sum(!is.na(LAG)), length(LAG)))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                  <- x
  input$qc_data_flagged          <- flg
  input$qc_info$t10_clock_offset <- list(n_flagged               = n_total,
                                         n_flagged_by_station    = n_station,
                                         n_judged_by_station     = n_judged,
                                         window_days             = window_days,
                                         max_lag                 = max_lag,
                                         min_lag                 = min_lag,
                                         min_consecutive_windows = min_consecutive_windows,
                                         min_valid_frac          = min_valid_frac,
                                         min_cor                 = min_cor,
                                         min_stations            = min_stations,
                                         lags                    = LAG,
                                         episodes                = episodes)
  input
}
