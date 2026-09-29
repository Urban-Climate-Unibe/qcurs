#' Temperature QC Level 4: stuck values
#'
#' Removes windows in which the sensor shows no variation at all. Via the
#' shared `qc_find_stuck()`: a window is judged only if it is sufficiently
#' populated - the guard whose absence in the published chain deleted whole
#' windows from two surviving identical readings - and the objected region
#' is expanded by position, never via time strings.
#'
#' The window is also the DETECTION DELAY: a stuck sensor looks fine until
#' the window is full, then the whole block is flagged retroactively. The
#' shortest safe window is the longest stretch over which genuine temperature
#' stays constant at the sensor's resolution - about 40 minutes at 0.01 K on
#' the Biel summer data, 2 hours at 0.1 K, and 10 hours at 0.5 K. The 6-hour
#' default follows the paper; shorten it per campaign, and check winter data
#' first, where fog and inversions produce far longer genuine plateaus.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param window_size Window length, either a `difftime` (e.g.
#'   `as.difftime(6, units = "hours")`) or a string "<n> <unit>" with a unit
#'   `as.difftime()` understands - secs, mins, hours, days, weeks.
#' @param na_tolerance_frac Maximum fraction of NA a window may contain and
#'   still be judged (paper: up to 50 percent).
#' @param min_non_na Minimum valid values a window needs (paper: 5).
#' @param sd_tol Variation at or below this counts as "no variation". 0 keeps
#'   exact constancy; a small value (e.g. 0.005) also catches a sensor
#'   alternating between two adjacent quantisation steps.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_4_stuck_values(res, window_size = "6 hours")
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_4_stuck_values <- function(input,
                                window_size = "6 hours",
                                na_tolerance_frac = 0.5,
                                min_non_na = 5,
                                sd_tol = 0,
                                verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t4_stuck_values")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the NA budget must be a fraction
  if (!is.numeric(na_tolerance_frac) || length(na_tolerance_frac) != 1 ||
      na_tolerance_frac < 0 || na_tolerance_frac >= 1)
    stop("na_tolerance_frac must be a fraction in [0, 1).")
  # the population floor must be at least 2, because sd() of one value is NA anyway
  if (!is.numeric(min_non_na) || length(min_non_na) != 1 || min_non_na < 2)
    stop("min_non_na must be at least 2.")
  # a negative tolerance would switch the test off silently
  if (!is.numeric(sd_tol) || length(sd_tol) != 1 || sd_tol < 0)
    stop("sd_tol must be zero or a positive variation.")

  #-------------------------------------------------------------------------------
  # the window in points, and the population guard as ONE minimum count

  width <- qc_window_points(window_size, input$qc_info$dataset_temperature$time_step_sec, "window_size")
  # "at most na_tolerance_frac NA" and "at least min_non_na valid" are both
  # lower bounds on the valid count; the stricter one applies
  min_valid <- max(min_non_na, width - floor(na_tolerance_frac * width))
  # a window that can never hold enough valid values would flag nothing, silently
  if (min_valid > width)
    stop(sprintf("window_size gives only %d points, but min_non_na = %d. Use a longer window or a smaller min_non_na.",
                 width, min_non_na))

  #-------------------------------------------------------------------------------
  # Perform QC Level 4

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

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # the search itself lives in qc_find_stuck(), shared with RH_QC_5
    st <- qc_find_stuck(v, width, sd_tol, min_valid)
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
    message(sprintf("T4 stuck values (%s = %d points, min valid %d, sd <= %g): %d flagged",
                    if (inherits(window_size, "difftime")) format(window_size) else window_size,
                    width, min_valid, sd_tol, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                 <- x
  input$qc_data_flagged         <- flg
  input$qc_info$t4_stuck_values <- list(n_flagged            = n_total,
                                        n_flagged_by_station = n_station,
                                        n_judged_by_station  = n_judged,
                                        window_size          = window_size,
                                        width_points         = width,
                                        min_valid            = min_valid,
                                        na_tolerance_frac    = na_tolerance_frac,
                                        min_non_na           = min_non_na,
                                        sd_tol               = sd_tol)
  input
}
