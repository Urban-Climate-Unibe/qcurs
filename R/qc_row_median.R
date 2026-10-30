#' Row medians of a numeric matrix, NA aware (internal)
#'
#' The median of every row, ignoring NA, by one sort of the whole matrix
#' instead of one `median()` call per row: ten times faster on a season of
#' 10-minute data, which matters when a level needs one median per station
#' (the median of all OTHER stations, as levels 9 and 10 do).
#'
#' @param X Numeric matrix, stations in columns.
#'
#' @return Numeric vector, one median per row; NA where a row has no value.
#'
#' @keywords internal
#' @noRd
qc_row_median <- function(X) {
  nr <- nrow(X); nc <- ncol(X)
  # a single column: the median is the value
  if (nc == 1) return(as.vector(X))
  # sort every row at once: by row first, by value second, NA last
  v <- X[order(row(X), X, na.last = TRUE)]
  # valid values per row; the sorted row starts at off + 1
  cnt <- rowSums(!is.na(X))
  off <- (seq_len(nr) - 1) * nc
  # the middle element(s) of the valid part; their mean is the median
  lo <- v[off + pmax(1, floor((cnt + 1) / 2))]
  hi <- v[off + pmax(1, ceiling((cnt + 1) / 2))]
  out <- (lo + hi) / 2
  out[cnt == 0] <- NA_real_
  out
}
