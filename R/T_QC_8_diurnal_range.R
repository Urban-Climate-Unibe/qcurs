#' Temperature QC Level 8: diurnal range collapse (a logger that went indoors)
#'
#' The case this level is for: a logger was taken down but not checked out
#' and keeps logging in an office or a store room. Its values pass every
#' classic test - physically plausible, smooth, varying a little - and at a
#' station without landuse-compatible neighbours they are spatially unchecked
#' too.
#'
#' Two things give such a logger away, and both are needed. First, its
#' diurnal range collapses while the network's does not: a day whose range is
#' below `range_frac` of the NETWORK's median range on the same day counts as
#' collapsed, and only `min_consecutive_days` collapsed days in a row count at
#' all. A fog week collapses every station together and stays invisible.
#' Second, over that run its DAILY MEANS no longer follow the weather: a
#' station alone under high fog has no diurnal cycle either, but it still
#' takes every cold front and every mild spell with the network (Grenchen,
#' January 2026: daily means -2.4 to +3.4 like the network, r = 0.95), while a
#' logger in a heated room sits at 21 whatever happens outside. So a run is
#' objected to only when the swing of the station's daily means stays below
#' `range_frac` of the network's swing over the same days. If the network
#' itself hardly moved during the run (less than `min_network_swing`), the
#' run cannot be decided; it is listed in `runs_undecided` and will be
#' decided once it has lasted through a change of weather - a logger indoors
#' waits that out, fog does not. A logger in an UNHEATED room follows the
#' weather in a damped way and may escape; the check-out date in the
#' metadata remains the primary source, this level is the net underneath.
#'
#' The level needs no metadata, only enough stations reporting on the same
#' day.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param range_frac Threshold as a fraction of the network: for the diurnal
#'   range of each day, and for the swing of the daily means over a run.
#' @param min_consecutive_days Collapsed days in a row before a run is judged
#'   at all (3: a swing over two days is just one difference).
#' @param min_network_swing Minimum swing of the network's daily means over a
#'   run, in Kelvin, for the run to be decidable.
#' @param min_obs_per_day Valid values a day needs to have its range judged
#'   (72 = half of a 10-minute day).
#' @param min_stations Stations with a range a day needs before it can serve
#'   as the yardstick; with fewer stations than that in the data the level
#'   skips itself and says why.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. Days
#'   without a network yardstick are in `qc_info$t8_diurnal_range$days_unjudged`,
#'   collapsed runs the weather did not allow to decide in `runs_undecided`.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_8_diurnal_range(res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_8_diurnal_range <- function(input,
                                 range_frac = 0.4,
                                 min_consecutive_days = 3,
                                 min_obs_per_day = 72,
                                 min_stations = 3,
                                 min_network_swing = 2,
                                 verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t8_diurnal_range")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the fraction must lie strictly between 0 and 1
  if (!is.numeric(range_frac) || length(range_frac) != 1 || range_frac <= 0 || range_frac >= 1)
    stop("range_frac must be a fraction between 0 and 1.")
  # at least one day in a row must be required
  if (!is.numeric(min_consecutive_days) || min_consecutive_days < 1)
    stop("min_consecutive_days must be at least 1.")
  # a range computed from too few points is meaningless
  if (!is.numeric(min_obs_per_day) || min_obs_per_day < 2)
    stop("min_obs_per_day must be at least 2.")
  # a network median of fewer than three stations is not a yardstick
  if (!is.numeric(min_stations) || min_stations < 3)
    stop("min_stations must be at least 3.")
  # a zero swing would call every calm spell decidable
  if (!is.numeric(min_network_swing) || length(min_network_swing) != 1 || min_network_swing <= 0)
    stop("min_network_swing must be a positive number of Kelvin.")
  # too few stations altogether: not an error but a skip - say why, record it, hand the pair back
  if (ncol(input$qc_data) < min_stations) {
    reason <- sprintf("only %d station(s); need at least %d for a network range", ncol(input$qc_data), min_stations)
    if (isTRUE(verbose))
      message("T8 diurnal range: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$t8_diurnal_range <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    return(input)
  }

  #-------------------------------------------------------------------------------
  # Perform QC Level 8

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # calendar day of every time step, in the time zone of the index (format()
  # honours it; as.Date() only does so from R 4.3 on)
  day <- format(zoo::index(x), "%Y-%m-%d")
  # the row numbers of every day, once, in the order of the record
  day_rows <- split(seq_along(day), factor(day, levels = unique(day)))
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: values on days that got a range verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))

  # diurnal range per day and station: max minus min, NA when the day is too sparse
  RNG <- vapply(colnames(X), function(s) vapply(day_rows, function(rows) {
    # the day's values
    z <- X[rows, s]
    # too sparse: no verdict for this day
    if (sum(!is.na(z)) < min_obs_per_day) return(NA_real_)
    # the diurnal amplitude
    diff(range(z, na.rm = TRUE))
  }, numeric(1)), numeric(length(day_rows)))
  # always a days x stations matrix, even for a single day
  RNG <- matrix(RNG, nrow = length(day_rows), dimnames = list(names(day_rows), colnames(X)))
  # daily mean per day and station, same shape, same sparsity rule: the
  # weather-coupling test of a run works on these
  DM <- vapply(colnames(X), function(s) vapply(day_rows, function(rows) {
    z <- X[rows, s]
    if (sum(!is.na(z)) < min_obs_per_day) return(NA_real_)
    mean(z, na.rm = TRUE)
  }, numeric(1)), numeric(length(day_rows)))
  DM <- matrix(DM, nrow = length(day_rows), dimnames = list(names(day_rows), colnames(X)))
  # the yardstick: the network's median range on each day, NA with too few stations
  n_rng <- rowSums(!is.na(RNG))
  ref   <- apply(RNG, 1, stats::median, na.rm = TRUE)
  ref[n_rng < min_stations | !(ref > 0)] <- NA_real_
  # the network's daily mean, same rule
  net_dm <- apply(DM, 1, stats::median, na.rm = TRUE)
  net_dm[n_rng < min_stations] <- NA_real_
  # days without a yardstick, so a zero-flag run is not a clean bill
  days_unjudged <- names(day_rows)[is.na(ref)]
  # collapsed runs the weather did not allow to decide
  runs_undecided <- character(0)

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station and its daily ranges
    v <- X[, s]; rng <- RNG[, s]
    # a day is judged when the station AND the network have a range
    judged_day <- !is.na(rng) & !is.na(ref)
    # every value on a judged day counts as judged
    n_judged[s] <- sum(!is.na(v[unlist(day_rows[judged_day])]))
    # days where this station collapsed while the network did not
    low <- judged_day & rng < range_frac * ref
    # run-length encode to find CONSECUTIVE stretches of collapsed days
    r <- rle(low)
    # end index of every run
    ends <- cumsum(r$lengths)
    # start index of every run
    starts <- ends - r$lengths + 1
    # every TRUE run that is long enough is a candidate; it is objected to
    # only when the station's daily means also stopped following the weather
    flag_days <- integer(0)
    for (k in seq_along(r$lengths)) {
      if (!r$values[k] || r$lengths[k] < min_consecutive_days) next
      d_idx <- starts[k]:ends[k]
      # how much the network moved over the run, and how much the station did
      sw_net <- suppressWarnings(diff(range(net_dm[d_idx], na.rm = TRUE)))
      sw_st  <- suppressWarnings(diff(range(DM[d_idx, s], na.rm = TRUE)))
      # a calm spell: nothing to compare against, say so and wait
      if (!is.finite(sw_net) || sw_net < min_network_swing) {
        runs_undecided <- c(runs_undecided, sprintf("%s:%s..%s(%d days, network swing %.1f K)",
                                                    s, names(day_rows)[starts[k]], names(day_rows)[ends[k]],
                                                    length(d_idx), if (is.finite(sw_net)) sw_net else NA))
        next
      }
      # the station still takes the fronts with the network: fog, not a room
      if (is.finite(sw_st) && sw_st >= range_frac * sw_net) next
      flag_days <- c(flag_days, d_idx)
    }
    # nothing collapsed long enough at this station
    if (!length(flag_days)) next
    # every valid value on the objected days
    hit <- rep(FALSE, length(v))
    hit[unlist(day_rows[flag_days])] <- TRUE
    hit <- hit & !is.na(v)

    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the collapsed days so later levels never see them
      X[mask, s] <- NA
      # record this level's code (8 = level 8, fixed by convention)
      previous_flag[mask, s] <- 8
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
  if (isTRUE(verbose)) {
    message(sprintf("T8 diurnal range collapse (< %g of network, >= %d days, daily means decoupled): %d flagged, %d of %d days judged, %d run(s) undecided",
                    range_frac, min_consecutive_days, n_total, sum(!is.na(ref)), length(ref), length(runs_undecided)))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                  <- x
  input$qc_data_flagged          <- flg
  input$qc_info$t8_diurnal_range <- list(n_flagged            = n_total,
                                         n_flagged_by_station = n_station,
                                         n_judged_by_station  = n_judged,
                                         range_frac           = range_frac,
                                         min_consecutive_days = min_consecutive_days,
                                         min_obs_per_day      = min_obs_per_day,
                                         min_stations         = min_stations,
                                         min_network_swing    = min_network_swing,
                                         days_unjudged        = days_unjudged,
                                         runs_undecided       = runs_undecided)
  input
}
