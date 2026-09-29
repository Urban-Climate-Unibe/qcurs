#' Linear gap filling with bookkeeping, station by station (internal)
#'
#' Fills gaps of at most `maxgap` steps by linear interpolation and writes the
#' interpolation codes of the convention: 50 for a filled plain gap, 50 + N
#' for a cell that level N had removed and is now refilled. The base 50 is
#' fixed - `RH_QC_1_inherit_temperature()` recognises a flag matrix by exactly
#' this code set. Series ends are never extrapolated.
#'
#' Used by `T_QC_9_interpolate()` and `RH_QC_8_interpolate()`, so the filling
#' exists once.
#'
#' @param X Numeric value matrix, QC removals already NA.
#' @param previous_flag Flag matrix of the same shape.
#' @param maxgap Longest gap (in steps) that is filled.
#' @param refill_flagged Also fill cells removed by the QC levels.
#'
#' @return list(X, previous_flag = the updated matrices; n_gap = plain gaps
#'   filled; n_ref = QC removals refilled).
#'
#' @keywords internal
#' @noRd
qc_fill_gaps <- function(X, previous_flag, maxgap, refill_flagged) {
  n_gap <- 0L; n_ref <- 0L
  for (s in colnames(X)) {
    v <- X[, s]; f <- previous_flag[, s]
    # everything that is currently missing (true gaps AND QC removals)
    was_na <- is.na(v)
    # optionally protect QC removals from being refilled
    keep_na <- if (isTRUE(refill_flagged)) rep(FALSE, length(v)) else was_na & !is.na(f) & f > 0
    # linear interpolation with a bounded gap length; ends are never extrapolated
    vi <- zoo::na.approx(v, na.rm = FALSE, maxgap = maxgap)
    vi[keep_na] <- NA
    # cells that actually received a value
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
