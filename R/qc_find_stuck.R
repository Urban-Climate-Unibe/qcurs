#' Find stuck stretches in one series (internal)
#'
#' A window ending at each point is "stuck" when its valid values show no
#' variation (standard deviation at or below `sd_tol`) and enough of the
#' window holds valid values. Every point of a stuck window is marked - the
#' region is expanded by POSITION, never by reconstructing time stamps.
#'
#' `exempt_above` is the humidity chain's saturation exception: a window whose
#' median is at or above it is not judged, because in fog and continuous rain
#' the humidity genuinely stops moving. Temperature passes NULL.
#'
#' Used by `T_QC_4_stuck_values()` and `RH_QC_5_stuck_values()`, so the
#' search exists once and a correction reaches both levels.
#'
#' @param v Numeric vector of one station, NA for gaps.
#' @param width Window length in points (see `qc_window_points()`).
#' @param sd_tol Variation at or below this counts as "no variation".
#' @param min_valid Minimum valid values a window needs to be judged.
#' @param exempt_above Windows with a median at or above this are skipped;
#'   NULL for no exemption.
#'
#' @return Logical vector of the same length as `v`: TRUE at every stuck point.
#'
#' @keywords internal
#' @noRd
qc_find_stuck <- function(v, width, sd_tol, min_valid, exempt_above = NULL) {
  hit <- rep(FALSE, length(v))
  # a station with fewer values than one window can never be judged
  if (sum(!is.na(v)) < width) return(hit)
  # rolling variation and valid count of the window ENDING at each point
  sdv <- zoo::rollapply(v, width, function(z) stats::sd(z, na.rm = TRUE), align = "right", fill = NA)
  cnt <- zoo::rollapply(v, width, function(z) sum(!is.na(z)),           align = "right", fill = NA)
  # windows with no variation and enough valid values
  stuck <- !is.na(sdv) & sdv <= sd_tol & cnt >= min_valid
  # the saturation exception: a constant window at the top is not a fault
  if (!is.null(exempt_above)) {
    med <- zoo::rollapply(v, width, function(z) stats::median(z, na.rm = TRUE), align = "right", fill = NA)
    stuck <- stuck & (is.na(med) | med < exempt_above)
  }
  ends <- which(stuck)
  if (!length(ends)) return(hit)
  # expand each window end back over its width, BY POSITION
  pos <- unique(unlist(lapply(ends, function(i) max(1, i - width + 1):i)))
  hit[pos] <- TRUE
  hit
}
