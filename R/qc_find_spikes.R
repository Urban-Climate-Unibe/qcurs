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
#' @return list(hit = logical, TRUE at every spike; judged = logical, TRUE at
#'   every point the test could actually be applied to). Both as long as `v`.
#'
#' @keywords internal
#' @noRd
qc_find_spikes <- function(v, dt, threshold) {
  n <- length(v)
  hit    <- rep(FALSE, n)
  judged <- rep(FALSE, n)
  for (i in seq_len(n)) {
    # a gap cannot be judged; the first and last point have no two-sided context
    if (is.na(v[i]) || i == 1 || i == n) next
    # the window around i, clamped to the edges of the series
    nb <- max(1, i - dt):min(n, i + dt)
    # too little context: refuse instead of guessing
    if (sum(!is.na(v[nb])) < (1 + dt)) next
    # the valid values before and after i
    bef <- v[max(1, i - dt):(i - 1)]; bef <- bef[!is.na(bef)]
    aft <- v[(i + 1):min(n, i + dt)]; aft <- aft[!is.na(aft)]
    # one side empty: no two-sided verdict possible
    if (!length(bef) || !length(aft)) next
    # from here on the test is applied, whatever it says
    judged[i] <- TRUE
    # off the window median (v[i] is inside the window; one outlier among
    # 2*dt + 1 values does not move it) AND off the past AND off the future
    hit[i] <- abs(v[i] - stats::median(v[nb], na.rm = TRUE)) >= threshold &&
              abs(v[i] - stats::median(bef)) >= threshold &&
              abs(v[i] - stats::median(aft)) >= threshold
  }
  list(hit = hit, judged = judged)
}
