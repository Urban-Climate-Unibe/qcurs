#' Humidity QC Level 5: stuck sensor, with a saturation exception
#'
#' Rolling standard deviation below `stuck_sd` over `stuck_window` marks a
#' stuck sensor - EXCEPT when the window median sits at saturation, because in
#' fog and continuous rain the humidity genuinely stops moving; without this
#' exception the test flags every winter systematically. A window must also be
#' sufficiently populated - the guard whose absence in the published
#' temperature chain deleted whole windows from two surviving identical
#' readings.
#'
#' The window width is derived from `qc_info$dataset_humidity$time_step_sec`,
#' the DOMINANT spacing of the series determined by `qc_prepare_input()`.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param stuck_window Window length, either a `difftime` or a string
#'   "<n> <unit>" with a unit `as.difftime()` understands (secs, mins, hours,
#'   days, weeks).
#' @param stuck_sd Variation below this counts as "no variation". Unlike the
#'   temperature twin the default is not 0: capacitive elements quantise more
#'   coarsely.
#' @param stuck_sat_max Windows whose median is at or above this are exempt
#'   (the saturation exception).
#' @param stuck_min_frac Minimum fraction of the window that must hold valid
#'   values for the window to be judged.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_5_persistence(res, stuck_window = "6 hours")
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_5_persistence <- function(input,
                                stuck_window = "6 hours",
                                stuck_sd = 0.1,
                                stuck_sat_max = 95,
                                stuck_min_frac = 0.5,
                                verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # a zero or negative tolerance would switch the test off silently
  if (!is.numeric(stuck_sd) || length(stuck_sd) != 1 || stuck_sd <= 0)
    stop("stuck_sd must be a positive variation.")
  # the saturation exemption must lie inside the physical range
  if (!is.numeric(stuck_sat_max) || stuck_sat_max <= 0 || stuck_sat_max > 100)
    stop("stuck_sat_max must lie in (0, 100].")
  # the population guard must be a usable fraction
  if (!is.numeric(stuck_min_frac) || stuck_min_frac <= 0 || stuck_min_frac > 1)
    stop("stuck_min_frac must be a fraction in (0, 1].")
  
  #-------------------------------------------------------------------------------
  # derive the window width in points from the resolution of the dataset
  
  # take the dominant time step the standard preamble determined for us
  step_min <- input$qc_info$dataset_humidity$time_step_sec / 60
  # a single-row series has no resolution; refuse instead of guessing
  if (!is.finite(step_min) || step_min <= 0)
    stop("Cannot determine the time step of the series (need at least two time stamps).")
  # the window length in minutes, via R's own duration type (see T_QC_4)
  if (inherits(stuck_window, "difftime")) {
    win_min <- as.numeric(stuck_window, units = "mins")
  } else {
    parts <- strsplit(trimws(stuck_window), "\\s+")[[1]]
    if (length(parts) != 2) stop("stuck_window must look like '6 hours' or '90 mins'.")
    unit <- match.arg(tolower(parts[2]), c("secs", "mins", "hours", "days", "weeks"))
    num  <- suppressWarnings(as.numeric(parts[1]))
    if (!is.finite(num) || num <= 0)
      stop(sprintf("'%s' is not a positive number of %s.", parts[1], unit))
    win_min <- as.numeric(as.difftime(num, units = unit), units = "mins")
  }
  # window length in points, right-aligned including the end point
  n_steps <- win_min / step_min
  # round and say what was used instead of truncating in silence
  if (abs(n_steps - round(n_steps)) > 1e-9)
    warning(sprintf("stuck_window '%s' is %.2f time steps at a %g min resolution; using %d steps (%g min).",
                    stuck_window, n_steps, step_min, round(n_steps), round(n_steps) * step_min))
  width <- as.integer(round(n_steps)) + 1
  if (width < 3)
    stop(sprintf("stuck_window '%s' is only %d time steps at a %g min resolution; need at least 2.",
                 stuck_window, width - 1, step_min))
  # absolute population floor derived from the fraction
  min_cnt <- ceiling(stuck_min_frac * width)
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 5
  
  # plain numeric matrix of the values (time in rows, stations in columns)
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
    # rolling level of the same windows, for the saturation exception
    med <- zoo::rollapply(v, width, function(z) stats::median(z, na.rm = TRUE),
                          align = "right", fill = NA)
    # rolling valid count of the same windows, for the population guard
    cnt <- zoo::rollapply(v, width, function(z) sum(!is.na(z)),
                          align = "right", fill = NA)
    # window ends with no variation, sufficient population, and a level BELOW
    # saturation - fog is allowed to be constant
    ends <- which(!is.na(sdv) & sdv < stuck_sd & cnt >= min_cnt &
                    (is.na(med) | med < stuck_sat_max))
    if (!length(ends)) next
    # expand each window end back over its width, BY POSITION
    pos <- unique(unlist(lapply(ends, function(i) max(1, i - width + 1):i)))
    hit[pos] <- TRUE
    
    # combine the verdict with THIS station's column only
    mask <- hit & !is.na(v) & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)
    
    # apply only if something was found
    if (n_found > 0) {
      # blank the stuck stretch so later levels never see it
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
  
  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("RH5 stuck (sd < %g over %s = %d points, skip if median >= %g): %d flagged",
                    stuck_sd,
                    if (inherits(stuck_window, "difftime")) format(stuck_window) else stuck_window,
                    width, stuck_sat_max, n_total))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                 <- x
  input$qc_data_flagged         <- flg
  input$qc_info$rh5_persistence <- list(n_flagged            = n_total,
                                        n_flagged_by_station = n_station,
                                        stuck_window         = stuck_window,
                                        width_points         = width,
                                        stuck_sd             = stuck_sd,
                                        stuck_sat_max        = stuck_sat_max,
                                        stuck_min_frac       = stuck_min_frac)
  input
}
