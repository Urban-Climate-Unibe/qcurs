#' Pairwise station distances from the metadata (internal)
#'
#' Validates the coordinates and returns the great-circle distance matrix in
#' metres, one row and column per station in the order of `md`, with NA on
#' the diagonal (a station is not its own neighbour).
#'
#' Used by `qc_neighbour_weights()` (levels T6 and RH7) and
#' `T_QC_7_spatiotemporal_consistency()`, so the geometry exists once.
#'
#' @param md Data frame with the columns ID, LAT, LON, already reduced to the
#'   stations of interest.
#'
#' @return Numeric matrix, dimnames = station IDs, NA diagonal.
#'
#' @keywords internal
#' @noRd
qc_distance_matrix <- function(md) {
  # coordinates must be numeric and present: a character column dies deep
  # inside geosphere, and an NA would silently zero the station's weights
  if (!is.numeric(md$LAT) || !is.numeric(md$LON) || anyNA(md$LAT) || anyNA(md$LON))
    stop("Metadata LAT/LON must be numeric and complete.")
  # all pairs at once, in metres
  D <- geosphere::distm(cbind(md$LON, md$LAT), fun = geosphere::distHaversine)
  dimnames(D) <- list(md$ID, md$ID)
  # a station is not its own neighbour
  diag(D) <- NA
  D
}
