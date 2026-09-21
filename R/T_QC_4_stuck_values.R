#' Temperature QC Level 4: temporal persistence (stuck sensor)
#'
#' Removes windows in which the sensor shows no variation at all.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param window_size Window length, either a `difftime` (e.g.
#'   `as.difftime(6, units = "hours")`) or a string "<n> <unit>" with a unit
#'   `as.difftime()` understands - secs, mins, hours, days, weeks, abbreviated
#'   as far as it stays unambiguous ("6 h" works, "6 m" is mins).
#' @param na_tolerance_frac Maximum fraction of NA a window may contain and
#'   still be judged (paper: up to 50 percent).
#' @param min_non_na Minimum valid values a window needs (paper: 5).
#' @param sd_tol Variation at or below this counts as "no variation". 0 keeps
#'   exact constancy; a small value (e.g. 0.005) also catches a sensor
#'   alternating between two adjacent quantisation steps.
#' @param flag_code QC-code written by this level. Here, default is 4.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_4_temporal_persistence(res, window_size = "6 hours")
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
                                        flag_code = 4,
                                        verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
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
  # derive the window width in points from the resolution of the dataset

  # take the dominant time step the standard preamble determined for us
  step_min <- input$qc_info$dataset_temperature$time_step_sec / 60
  # a single-row series has no resolution; refuse instead of guessing
  if (!is.finite(step_min) || step_min <= 0)
    stop("Cannot determine the time step of the series (need at least two time stamps).")
  # the window length in minutes.
  if (inherits(window_size, "difftime")) {
    win_min <- as.numeric(window_size, units = "mins")
  } else {
    parts <- strsplit(trimws(window_size), "\\s+")[[1]]
    if (length(parts) != 2) stop("window_size must look like '6 hours' or '90 mins'.")
    unit <- match.arg(tolower(parts[2]), c("secs", "mins", "hours", "days", "weeks"))
    num  <- suppressWarnings(as.numeric(parts[1]))
    if (!is.finite(num) || num <= 0)
      stop(sprintf("'%s' is not a positive number of %s.", parts[1], unit))
    win_min <- as.numeric(as.difftime(num, units = unit), units = "mins")
  }
  # window length in points, right-aligned including the end point:
  # 6 h at 10-min data = 36 intervals = 37 points
  n_steps <- win_min / step_min
  # a window that is not a whole number of steps cannot be used --> round
  if (abs(n_steps - round(n_steps)) > 1e-9)
    warning(sprintf("window_size '%s' is %.2f time steps at a %g min resolution; using %d steps (%g min).",
                    window_size, n_steps, step_min, round(n_steps), round(n_steps) * step_min))
  width <- as.integer(round(n_steps)) + 1
  if (width < 3)
    stop(sprintf("window_size '%s' is only %d time steps at a %g min resolution; need at least 2.",
                 window_size, width - 1, step_min))
  # absolute NA budget derived from the fraction
  max_na <- floor(na_tolerance_frac * width)

  #-------------------------------------------------------------------------------
  # Perform QC Level 4

  # plain numeric matrix of the values (time in rows, stations in columns).
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station and determine its length
    v <- X[, s]; n <- length(v)
    # a station with fewer values than one window can never be judged: skip it
    if (sum(!is.na(v)) < width) next
    # create a vector with length n. All entries are FALSE
    hit <- rep(FALSE, n)

    # rolling variation of the window ENDING at each point
    sdv <- zoo::rollapply(v, width, function(z) stats::sd(z, na.rm = TRUE),
                          align = "right", fill = NA)
    # rolling NA count of the same windows
    nna <- zoo::rollapply(v, width, function(z) sum(is.na(z)),
                          align = "right", fill = NA)
    # rolling count of valid values in the same windows
    nok <- width - nna
    # window ends that show (near) zero variation, respect the NA budget, and
    # hold enough valid values for the verdict to mean anything
    ends <- which(!is.na(sdv) & sdv <= sd_tol & nna <= max_na & nok >= min_non_na)
    # nothing stuck at this station
    if (!length(ends)) next
    # expand each window end back over its width, BY POSITION (never via time
    # strings), and mark those positions in the verdict vector
    pos <- unique(unlist(lapply(ends, function(i) max(1, i - width + 1):i)))
    hit[pos] <- TRUE

    # combine the verdict with THIS station's column only: hit is a vector of
    # length n, so it must meet vectors, not the whole matrix
    mask <- hit & !is.na(v) &
      (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the stuck stretch so later levels never see it
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
    message(sprintf("T4 persistence (%s = %d points, NA budget %d, min valid %d): %d flagged",
                    if (inherits(window_size, "difftime")) format(window_size) else window_size,
                    width, max_na, min_non_na, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                         <- x
  input$qc_data_flagged                 <- flg
  input$qc_info$t4_temporal_persistence <- list(n_flagged            = n_total,
                                                n_flagged_by_station = n_station,
                                                window_size          = window_size,
                                                width_points         = width,
                                                na_tolerance_frac    = na_tolerance_frac,
                                                min_non_na           = min_non_na,
                                                sd_tol               = sd_tol)
  input
}
