#' Window length in points from a duration and the series resolution (internal)
#'
#' Turns "6 hours" (or a `difftime`) into the number of points a right-aligned
#' window spans at the given time step: 6 h at 10-minute data = 36 intervals
#' = 37 points. Used by every window based level, so the parsing and its
#' guards exist once.
#'
#' `as.difftime()` knows the unit conversions and `match.arg()` does the
#' abbreviation matching ("h", "hour", "hours") plus the error message for
#' anything else. A window that is not a whole number of steps is rounded,
#' with a warning that says what was actually used.
#'
#' @param window A `difftime`, or a string "<n> <unit>" with a unit
#'   `as.difftime()` understands: secs, mins, hours, days, weeks.
#' @param step_sec The dominant time step in seconds, from
#'   `qc_info$dataset_<what>$time_step_sec`.
#' @param arg Name of the caller's argument, for the messages.
#'
#' @return Integer: the window length in points, at least 3.
#'
#' @keywords internal
#' @noRd
qc_window_points <- function(window, step_sec, arg = "window_size") {
  # a single-row series has no resolution; refuse instead of guessing
  step_min <- step_sec / 60
  if (!is.finite(step_min) || step_min <= 0)
    stop("Cannot determine the time step of the series (need at least two time stamps).")
  # the window length in minutes, via R's own duration type
  if (inherits(window, "difftime")) {
    win_min <- as.numeric(window, units = "mins")
  } else {
    parts <- strsplit(trimws(window), "\\s+")[[1]]
    if (length(parts) != 2) stop(sprintf("%s must look like '6 hours' or '90 mins'.", arg))
    unit <- match.arg(tolower(parts[2]), c("secs", "mins", "hours", "days", "weeks"))
    num  <- suppressWarnings(as.numeric(parts[1]))
    if (!is.finite(num) || num <= 0)
      stop(sprintf("'%s' is not a positive number of %s.", parts[1], unit))
    win_min <- as.numeric(as.difftime(num, units = unit), units = "mins")
  }
  label <- if (inherits(window, "difftime")) format(window) else window
  # window length in intervals; round and SAY so instead of truncating in silence
  n_steps <- win_min / step_min
  if (abs(n_steps - round(n_steps)) > 1e-9)
    warning(sprintf("%s '%s' is %.2f time steps at a %g min resolution; using %d steps (%g min).",
                    arg, label, n_steps, step_min, round(n_steps), round(n_steps) * step_min))
  # points = intervals + 1 (right-aligned, end point included)
  width <- as.integer(round(n_steps)) + 1L
  if (width < 3)
    stop(sprintf("%s '%s' is only %d time steps at a %g min resolution; need at least 2.",
                 arg, label, width - 1, step_min))
  width
}
