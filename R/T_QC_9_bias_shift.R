#' Temperature QC Level 9: bias shift against the network
#'
#' A sensor that reads 2 K too high after a knock, a replaced element with
#' another calibration, a shield that broke and lets the sun in: the values
#' stay plausible, smooth and within the spatial tolerance of level 6
#' (max(6 sigma, 3 K) on single values), the station simply becomes the
#' warmest one. Only a before and an after give such a sensor away.
#'
#' Per station and day, the daily mean minus the median of the OTHER
#' stations' daily means is the station's anomaly: its site climate plus its
#' calibration. The median anomaly of the station's first `ref_days` judged
#' days is its baseline. The anomaly scatters from day to day by half a
#' Kelvin even at a healthy station (a sunny day warms the open sites more
#' than the forest), so it is smoothed with a running median over
#' `min_consecutive_days` judged days before it is compared: one mild day
#' does not break a run, one front day does not start one. A run of at
#' least `min_consecutive_days` judged days on which the smoothed anomaly
#' departs from the baseline by more than `bias_tol` is objected to - every
#' valid value of those days. The level finds the CHANGE, not an absolute
#' offset: a sensor that read 2 K high from day one has a baseline of +2
#' and stays unobjected, and when a station shifts inside its baseline
#' period the level objects to the shorter side. Which side is right is a
#' decision for whoever knows the station; the `episodes` table (station,
#' from, to, days, shift) is the lead.
#'
#' Weather moves the anomalies too: a foehn or an inversion spell separates
#' the valley from the hill stations by several Kelvin. Two rules keep that
#' out. A day on which more than `max_network_frac` of the judged stations
#' (and at least two) depart from their baselines is not judged at all and
#' is listed in `days_unjudged`; and a run has to be `min_consecutive_days`
#' long - a front lasts a day, a calibration lasts until the next visit.
#' Unjudged days neither continue nor break a run, because an offset does
#' not pause for a front. On two summers of the Thun network (forest, lake
#' and hill sites among them) the smoothed anomaly of a healthy station
#' never departed by more than 1.15 K from its baseline, and a planted
#' +2 K from mid-July on was found from its first day to the end of the
#' record. A station with a climate of its own can still show seasonal
#' episodes; they stand in the table with size and sign, and a `bias_tol`
#' of 2 or 3 K is the remedy for such a network.
#'
#' The level needs no metadata, only `min_stations` stations with a daily
#' mean on the same day.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param bias_tol Departure of the smoothed anomaly from the baseline, in
#'   Kelvin, beyond which a day counts as shifted.
#' @param min_consecutive_days Shifted judged days in a row before any of
#'   them is objected to; also the width of the running median (made odd
#'   when even).
#' @param ref_days Judged days at the start of the station's record that
#'   form its baseline; a station with fewer judged days is refused.
#' @param min_obs_per_day Valid values a day needs to have a daily mean
#'   (72 = half of a 10-minute day).
#' @param min_stations Stations with a daily mean a day needs before it can
#'   be judged; with fewer stations than that in the data the level skips
#'   itself and says why.
#' @param max_network_frac Share of the judged stations above which a day
#'   counts as a weather day and is not judged at all (needs at least two
#'   shifted stations). 1 switches the rule off.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. The
#'   baselines are in `qc_info$t9_bias_shift$baseline`, the objected runs in
#'   `episodes` (one row each: station, from, to, days, shift_K), the weather
#'   days in `days_unjudged`, the stations without a baseline in `refused`.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_9_bias_shift(res)
#' res$qc_info$t9_bias_shift$episodes
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_9_bias_shift <- function(input,
                              bias_tol = 1.5,
                              min_consecutive_days = 5,
                              ref_days = 30,
                              min_obs_per_day = 72,
                              min_stations = 3,
                              max_network_frac = 0.25,
                              verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t9_bias_shift")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the tolerance is a positive number of Kelvin
  if (!is.numeric(bias_tol) || length(bias_tol) != 1 || bias_tol <= 0)
    stop("bias_tol must be a positive number of Kelvin.")
  # at least one day in a row must be required
  if (!is.numeric(min_consecutive_days) || min_consecutive_days < 1)
    stop("min_consecutive_days must be at least 1.")
  # a baseline from a handful of days is noise, not a reference
  if (!is.numeric(ref_days) || length(ref_days) != 1 || ref_days < 5)
    stop("ref_days must be at least 5.")
  # a daily mean from too few values is not a daily mean
  if (!is.numeric(min_obs_per_day) || min_obs_per_day < 2)
    stop("min_obs_per_day must be at least 2.")
  # a network median of fewer than three stations is not a yardstick
  if (!is.numeric(min_stations) || min_stations < 3)
    stop("min_stations must be at least 3.")
  # the weather share must be a usable fraction; 1 means the rule is off
  if (!is.numeric(max_network_frac) || max_network_frac <= 0 || max_network_frac > 1)
    stop("max_network_frac must lie in (0, 1].")
  # too few stations altogether: not an error but a skip - say why, record it, hand the pair back
  if (ncol(input$qc_data) < min_stations) {
    reason <- sprintf("only %d station(s); need at least %d for a network anomaly", ncol(input$qc_data), min_stations)
    if (isTRUE(verbose))
      message("T9 bias shift: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$t9_bias_shift <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    return(input)
  }

  #-------------------------------------------------------------------------------
  # Perform QC Level 9

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
  # per-station coverage: values on days that got a verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))

  # daily mean per day and station, NA when the day is too sparse
  DM <- vapply(colnames(X), function(s) vapply(day_rows, function(rows) {
    # the day's values
    z <- X[rows, s]
    # too sparse: no daily mean for this day
    if (sum(!is.na(z)) < min_obs_per_day) return(NA_real_)
    mean(z, na.rm = TRUE)
  }, numeric(1)), numeric(length(day_rows)))
  # always a days x stations matrix, even for a single day
  DM <- matrix(DM, nrow = length(day_rows), dimnames = list(names(day_rows), colnames(X)))
  # a day is usable when enough stations have a daily mean
  day_ok <- rowSums(!is.na(DM)) >= min_stations
  # the anomaly: the station's daily mean minus the median of the OTHER
  # stations. Leaving the station out matters in a small network, where a
  # shifted station would otherwise drag the yardstick along
  A <- DM
  for (j in seq_len(ncol(DM))) A[, j] <- DM[, j] - qc_row_median(DM[, -j, drop = FALSE])
  A[!day_ok, ] <- NA_real_

  # the baseline per station: the median anomaly of its first ref_days
  # judged days. Too short a record gets no baseline and is refused
  baseline <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  refused <- character(0)
  for (s in colnames(X)) {
    jd <- which(!is.na(A[, s]))
    if (length(jd) < ref_days) {
      if (length(jd) > 0 || any(!is.na(X[, s])))
        refused <- c(refused, sprintf("%s(days=%d)", s, length(jd)))
      next
    }
    baseline[s] <- stats::median(A[jd[seq_len(ref_days)], s])
  }

  # the anomaly smoothed per station with a running median over the judged
  # days, so a single mild day does not break a run and a single front day
  # does not start one. The width is min_consecutive_days, made odd
  k <- min_consecutive_days + (min_consecutive_days %% 2 == 0)
  S <- A
  for (s in colnames(X)) {
    jd <- which(!is.na(A[, s]))
    if (length(jd) > k) S[jd, s] <- stats::runmed(A[jd, s], k, endrule = "median")
  }
  # the daily verdict per station: NA = no verdict, TRUE = departed from the baseline
  V <- abs(sweep(S, 2, baseline, "-")) > bias_tol

  # the network rule: when more than max_network_frac of the judged stations
  # (and at least two) depart on the same day, it is the weather - foehn, an
  # inversion - not a sensor. Nobody is judged on such a day
  n_jud <- rowSums(!is.na(V))
  n_dep <- rowSums(V, na.rm = TRUE)
  weather <- n_dep >= 2 & n_dep / pmax(n_jud, 1) > max_network_frac
  V[weather, ] <- NA
  # the weather days, so a zero-flag run is not a clean bill
  days_unjudged <- names(day_rows)[weather]
  # the objected runs, one row each
  episodes <- list()

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # no baseline: nothing to judge against
    if (is.na(baseline[s])) next
    # extract data vector of this station
    v <- X[, s]
    # the judged days of this station, in the order of the record
    jd <- which(!is.na(V[, s]))
    # every value on a judged day counts as judged
    n_judged[s] <- sum(!is.na(v[unlist(day_rows[jd])]))
    if (!length(jd)) next
    # run-length encode over the JUDGED days only: an unjudged day in between
    # neither continues nor breaks a run, an offset does not pause for a front
    r <- rle(V[jd, s])
    # end index of every run
    ends <- cumsum(r$lengths)
    # start index of every run
    starts <- ends - r$lengths + 1
    # collect the day indices of every departed run that is long enough
    flag_days <- integer(0)
    for (k in seq_along(r$lengths)) {
      if (!r$values[k] || r$lengths[k] < min_consecutive_days) next
      d_idx <- jd[starts[k]:ends[k]]
      flag_days <- c(flag_days, d_idx)
      # the run for the report: where, how long, how far from the baseline
      episodes[[length(episodes) + 1]] <-
        data.frame(station = s,
                   from    = names(day_rows)[d_idx[1]],
                   to      = names(day_rows)[d_idx[length(d_idx)]],
                   days    = length(d_idx),
                   shift_K = round(stats::median(A[d_idx, s]) - baseline[[s]], 2),
                   stringsAsFactors = FALSE)
    }
    # nothing shifted long enough at this station
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
      # blank the shifted days so later levels never see them
      X[mask, s] <- NA
      # record this level's code (9 = level 9, fixed by convention)
      previous_flag[mask, s] <- 9
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
               days = integer(0), shift_K = numeric(0), stringsAsFactors = FALSE)

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip - and name the refusals
  if (isTRUE(verbose)) {
    message(sprintf("T9 bias shift (> %g K off the baseline of the first %d days, >= %d days): %d flagged, %d episode(s), %d weather day(s) not judged, %d of %d stations with a baseline",
                    bias_tol, ref_days, min_consecutive_days, n_total, nrow(episodes),
                    length(days_unjudged), sum(!is.na(baseline)), ncol(X)))
    if (length(refused) > 0)
      message("  refused (record shorter than ref_days): ", paste(refused, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data               <- x
  input$qc_data_flagged       <- flg
  input$qc_info$t9_bias_shift <- list(n_flagged            = n_total,
                                      n_flagged_by_station = n_station,
                                      n_judged_by_station  = n_judged,
                                      bias_tol             = bias_tol,
                                      min_consecutive_days = min_consecutive_days,
                                      ref_days             = ref_days,
                                      min_obs_per_day      = min_obs_per_day,
                                      min_stations         = min_stations,
                                      max_network_frac     = max_network_frac,
                                      baseline             = baseline,
                                      episodes             = episodes,
                                      days_unjudged        = days_unjudged,
                                      refused              = refused)
  input
}
