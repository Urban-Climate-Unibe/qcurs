#' Linear gap filling with bookkeeping, station by station (internal)
#'
#' Fills gaps of at most `maxgap` steps by linear interpolation and writes the
#' interpolation codes of the convention: 50 for a filled plain gap, 50 + N
#' for a cell that level N had removed and is now refilled. The base 50 is
#' fixed - `RH_QC_1_inherit_temperature()` recognises a flag matrix by exactly
#' this code set. Series ends are never extrapolated. With the time stamps
#' given, a gap is measured on the full grid: rows missing from the index
#' count as steps, so a hole of six hours is not bridged as one step.
#'
#' Used by `T_QC_10_interpolate()` and `RH_QC_8_interpolate()`, so the filling
#' exists once.
#'
#' @param X Numeric value matrix, QC removals already NA.
#' @param previous_flag Flag matrix of the same shape.
#' @param maxgap Longest gap (in steps) that is filled.
#' @param refill_flagged Also fill cells removed by the QC levels.
#' @param tt Time stamps of the rows in seconds, or NULL to count rows.
#' @param step_sec Time step of the data in seconds.
#'
#' @return list(X, previous_flag = the updated matrices; n_gap = plain gaps
#'   filled; n_ref = QC removals refilled).
#'
#' @keywords internal
#' @noRd
qc_fill_gaps <- function(X, previous_flag, maxgap, refill_flagged, tt = NULL, step_sec = NA_real_) {
  n_gap <- 0L; n_ref <- 0L
  # the position of every row on the full grid; the row number where the time
  # stamps are not given or do not sit on one grid
  pos <- seq_len(nrow(X))
  if (!is.null(tt) && is.finite(step_sec) && step_sec > 0) {
    p <- round((tt - tt[1]) / step_sec) + 1
    if (!anyDuplicated(p)) pos <- p
  }
  for (s in colnames(X)) {
    v <- X[, s]; f <- previous_flag[, s]
    # everything that is currently missing (true gaps and QC removals)
    was_na <- is.na(v)
    # optionally protect QC removals from being refilled
    keep_na <- if (isTRUE(refill_flagged)) rep(FALSE, length(v)) else was_na & !is.na(f) & f > 0
    # linear interpolation with a bounded gap length; ends are never extrapolated
    full <- rep(NA_real_, max(pos)); full[pos] <- v
    vi <- zoo::na.approx(full, na.rm = FALSE, maxgap = maxgap)[pos]
    vi[keep_na] <- NA
    # cells that received a value
    filled <- was_na & !is.na(vi)
    if (!any(filled)) next
    X[filled, s] <- vi[filled]
    # fills over former QC removals get 50 + the original level number,
    # fills over plain gaps get 50
    qc_hit <- filled & !is.na(f) & f > 0
    plain  <- filled & (is.na(f) | f == 0)
    previous_flag[qc_hit, s] <- 50 + f[qc_hit]
    previous_flag[plain,  s] <- 50
    n_ref <- n_ref + sum(qc_hit)
    n_gap <- n_gap + sum(plain)
  }
  list(X = X, previous_flag = previous_flag, n_gap = n_gap, n_ref = n_ref)
}
