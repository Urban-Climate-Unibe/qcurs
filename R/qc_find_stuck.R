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
#' @return list(hit = logical, TRUE at every stuck point; judged = logical,
#'   TRUE at every point covered by at least one judged window). Both as long
#'   as `v`, both FALSE at gaps.
#'
#' @keywords internal
#' @noRd
qc_find_stuck <- function(v, width, sd_tol, min_valid, exempt_above = NULL) {
  n <- length(v)
  hit    <- rep(FALSE, n)
  judged <- rep(FALSE, n)
  # a series shorter than one window, or without enough values to fill one, is never judged
  if (n < width || sum(!is.na(v)) < min_valid) return(list(hit = hit, judged = judged))
  # rolling variation and valid count of the window ENDING at each point
  sdv <- zoo::rollapply(v, width, function(z) stats::sd(z, na.rm = TRUE), align = "right", fill = NA)
  cnt <- zoo::rollapply(v, width, function(z) sum(!is.na(z)),           align = "right", fill = NA)
  # windows with enough valid values get a verdict
  ok <- !is.na(cnt) & cnt >= min_valid
  # the saturation exception: a constant window at the top is not judged
  if (!is.null(exempt_above)) {
    med <- zoo::rollapply(v, width, function(z) stats::median(z, na.rm = TRUE), align = "right", fill = NA)
    ok <- ok & !is.na(med) & med < exempt_above
  }
  # the stuck windows among the judged ones
  stuck <- ok & !is.na(sdv) & sdv <= sd_tol
  # expand each window end back over its width, BY POSITION; gaps inside a
  # window are not points and stay FALSE
  expand <- function(ends) unique(unlist(lapply(ends, function(i) max(1, i - width + 1):i)))
  judged[expand(which(ok))] <- TRUE
  hit[expand(which(stuck))] <- TRUE
  list(hit = hit & !is.na(v), judged = judged & !is.na(v))
}
