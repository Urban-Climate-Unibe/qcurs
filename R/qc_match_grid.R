#' Lay a second series onto the grid of the first, by time stamp and name (internal)
#'
#' The one alignment rule of the humidity chain: a value from the temperature
#' run belongs to a humidity cell only if it has EXACTLY the same time stamp
#' and the same logger name. Everything else stays NA. No positional
#' matching, no guessing.
#'
#' Used by `RH_QC_1_inherit_temperature()` (flags), `RH_QC_5_decoupling()`
#' and `RH_QC_7_dewpoint_consistency()` (values).
#'
#' @param x The target xts (the humidity series).
#' @param src The source xts with logger names (a temperature series or flag
#'   matrix).
#'
#' @return list(M = numeric matrix like `x`, source values where matched, NA
#'   elsewhere; loggers = matched logger names; n_time = matched time steps),
#'   or NULL when nothing matches at all.
#'
#' @keywords internal
#' @noRd
qc_match_grid <- function(x, src) {
  if (!inherits(src, "xts") || is.null(colnames(src))) return(NULL)
  # loggers present in BOTH series, by name
  loggers <- intersect(colnames(x), colnames(src))
  # for every target time stamp: the row of the SAME stamp in the source, or NA
  row_s <- match(as.numeric(zoo::index(x)), as.numeric(zoo::index(src)))
  ok <- !is.na(row_s)
  if (!length(loggers) || !any(ok)) return(NULL)
  # the source values on the target grid: NA wherever there is no (time, logger) match
  M <- matrix(NA_real_, nrow(x), ncol(x), dimnames = list(NULL, colnames(x)))
  M[ok, loggers] <- as.matrix(zoo::coredata(src))[row_s[ok], loggers]
  list(M = M, loggers = loggers, n_time = sum(ok))
}
