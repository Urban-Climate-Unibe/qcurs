#' Find isolated spikes in one series (internal)
#'
#' The one spike test both chains use: a value is a spike only if it deviates
#' by at least `threshold` from the median of its window AND from the values
#' before it AND from the values after it. The triple condition makes the
#' test immune to genuine fronts and step changes: those move the neighbours
#' with them, so at most the single point of the jump can ever qualify.
#' Persistent offsets are invisible to this test by construction.
#'
#' Used by `T_QC_3_time_consistency()` and `RH_QC_3_spike()`, so the search
#' exists once and a correction reaches both levels.
#'
#' @param v Numeric vector of one station, NA for gaps.
#' @param dt Half window in time steps.
#' @param threshold Minimum deviation, in the unit of `v`.
#'
#' @return Logical vector of the same length as `v`: TRUE at every spike.
#'
#' @keywords internal
#' @noRd
qc_find_spikes <- function(v, dt, threshold) {
  n <- length(v)
  # every point innocent until proven guilty
  hit <- rep(FALSE, n)
  # walk every time step
  for (i in seq_len(n)) {
    # a gap cannot be judged; the first and last point have no two-sided context
    if (is.na(v[i]) || i == 1 || i == n) next
    # the window around i, clamped to the edges of the series
    nb <- max(1, i - dt):min(n, i + dt)
    # too little context to judge: refuse instead of guessing
    if (sum(!is.na(v[nb])) < (1 + dt)) next
    # condition 1: off the window median (the window includes v[i] itself,
    # which a single outlier among 2*dt + 1 values does not move)
    if (abs(v[i] - stats::median(v[nb], na.rm = TRUE)) < threshold) next
    # the values BEFORE i, gaps removed
    bef <- v[max(1, i - dt):(i - 1)]; bef <- bef[!is.na(bef)]
    # condition 2: off the past
    if (!length(bef) || abs(v[i] - stats::median(bef)) < threshold) next
    # the values AFTER i, gaps removed
    aft <- v[(i + 1):min(n, i + dt)]; aft <- aft[!is.na(aft)]
    # condition 3: off the future
    if (!length(aft) || abs(v[i] - stats::median(aft)) < threshold) next
    # all three conditions met: an isolated spike
    hit[i] <- TRUE
  }
  hit
}
