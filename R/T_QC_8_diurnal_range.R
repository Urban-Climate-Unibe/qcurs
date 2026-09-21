#' Temperature QC Level 8: diurnal range collapse
#'
#' A sensor reading indoors (or otherwise decoupled from
#' the atmosphere) produces values that might pass every classic test - physically
#' plausible, locally smooth, varying and at stations without
#' landuse-compatible neighbours it is also spatially unchecked.
#'
#' This level removes days whose diurnal range falls below a fraction of the
#' station's OWN median range, for a minimum number of CONSECUTIVE days. It
#' needs no neighbours and therefore also covers isolated stations
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param range_frac Threshold as a fraction of the station median range.
#' @param min_consecutive_days Days in a row below the threshold before any of
#'   them is objected to.
#' @param min_obs_per_day Valid values a day needs to have its range judged
#'   (72 = half of a 10-minute day).
#' @param min_reference_days Valid days a station needs for a robust median.
#' @param flag_code QC-code written by this level. Here, default is 8.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
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
                                 min_consecutive_days = 2,
                                 min_obs_per_day = 72,
                                 min_reference_days = 14,
                                 flag_code = 8,
                                 verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
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
  # the reference median needs a real sample
  if (!is.numeric(min_reference_days) || min_reference_days < 3)
    stop("min_reference_days must be at least 3.")

  #-------------------------------------------------------------------------------
  # Perform QC Level 8

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # calendar day of every time step
  day <- as.Date(zoo::index(x))
  # the day sequence of the record, in order
  days <- unique(day)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # diurnal range per day: max minus min, or NA when the day is too sparse
    rng <- vapply(days, function(dd) {
      # the day's values
      z <- v[day == dd]
      # too sparse: no verdict for this day
      if (sum(!is.na(z)) < min_obs_per_day) return(NA_real_)
      # the diurnal amplitude
      diff(range(z, na.rm = TRUE))
    }, numeric(1))
    # the station's own typical amplitude, robust against single odd days
    ref <- stats::median(rng, na.rm = TRUE)
    # not enough judged days for a baseline: skip the whole station
    if (sum(!is.na(rng)) < min_reference_days || !is.finite(ref) || ref <= 0) next
    # days with a collapsed amplitude
    low <- !is.na(rng) & rng < range_frac * ref
    # run-length encode to find CONSECUTIVE stretches of collapsed days
    r <- rle(low)
    # end index of every run
    ends <- cumsum(r$lengths)
    # start index of every run
    starts <- ends - r$lengths + 1
    # collect the day indices of every TRUE run that is long enough
    flag_days <- integer(0)
    for (k in seq_along(r$lengths))
      if (r$values[k] && r$lengths[k] >= min_consecutive_days)
        flag_days <- c(flag_days, starts[k]:ends[k])
    # nothing collapsed long enough at this station
    if (!length(flag_days)) next
    # every valid value on the objected days
    hit <- day %in% days[flag_days] & !is.na(v)

    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the collapsed days so later levels never see them
      X[mask, s] <- NA
      # record this level's code
      previous_flag[mask, s] <- flag_code
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
    message(sprintf("T8 diurnal range collapse (< %g of station median, >= %d days): %d flagged",
                    range_frac, min_consecutive_days, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                  <- x
  input$qc_data_flagged          <- flg
  input$qc_info$t8_diurnal_range <- list(n_flagged            = n_total,
                                         n_flagged_by_station = n_station,
                                         range_frac           = range_frac,
                                         min_consecutive_days = min_consecutive_days,
                                         min_obs_per_day      = min_obs_per_day,
                                         min_reference_days   = min_reference_days)
  input
}
