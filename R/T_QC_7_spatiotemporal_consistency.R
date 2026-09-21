#' Temperature QC Level 7: spatiotemporal consistency
#'
#' Removes values that are simultaneously extreme in space (against ALL nearby
#' stations) and in time (against their own preceding and following steps).
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param metadata Data frame with columns ID, LAT, LON.
#' @param radius_m Neighbour search radius in metres (paper: 2500).
#' @param n_neighbours Required number of nearby stations (paper: 5).
#' @param threshold Quantile for both the spatial and the temporal criteria.
#' @param verbose Report the tally and the skipped stations.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. Stations
#'   without a full neighbourhood are listed in
#'   `qc_info$t7_spatiotemporal$skipped`.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_7_spatiotemporal_consistency(res, metadata = meta_biel)
#' }
#'
#' @import xts
#' @import zoo
#' @import geosphere
#' @export
T_QC_7_spatiotemporal_consistency <- function(input,
                                              metadata,
                                              radius_m = 2500,
                                              n_neighbours = 5,
                                              threshold = 0.99,
                                              verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters and the metadata

  # the quantile must be a probability above 0.5: both tails are derived from it
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0.5 || threshold >= 1)
    stop("threshold must be a probability between 0.5 and 1.")
  if (!is.numeric(radius_m) || radius_m <= 0) stop("radius_m must be a positive distance in metres.")
  if (!is.numeric(n_neighbours) || n_neighbours < 1) stop("n_neighbours must be at least 1.")
  # the metadata are required for this level
  if (missing(metadata) || is.null(metadata)) stop("This level needs the metadata (ID, LAT, LON).")
  # tolerate tibbles and friends
  md <- as.data.frame(metadata)
  # the three columns this level relies on
  if (!all(c("ID", "LAT", "LON") %in% names(md)))
    stop("Metadata must contain the columns ID, LAT, LON.")
  # duplicated IDs would make the match ambiguous
  if (anyDuplicated(md$ID)) stop("Metadata contains duplicated IDs.")
  # neighbours are searched among DATA-BEARING stations only (repair, see header)
  ids <- intersect(colnames(x), md$ID)
  # align the metadata rows to those stations
  md <- md[match(ids, md$ID), ]
  # number of stations actually usable
  ns <- length(ids)

  #-------------------------------------------------------------------------------
  # pairwise distances once, they do not change over time

  # pairwise distance matrix in metres
  D <- matrix(0, ns, ns, dimnames = list(ids, ids))
  for (i in seq_len(ns)) for (j in seq_len(ns))
    D[i, j] <- geosphere::distHaversine(c(md$LON[i], md$LAT[i]), c(md$LON[j], md$LAT[j]))
  # a station is not its own neighbour
  diag(D) <- NA

  # internal shift helpers, so this file needs no dplyr dependency
  shift_back <- function(v) c(NA, v[-length(v)])   # the previous value at each position
  shift_fwd  <- function(v) c(v[-1], NA)           # the next value at each position

  #-------------------------------------------------------------------------------
  # Perform QC Level 7

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # stations without a full neighbourhood, so a zero-flag run is not a clean bill
  skipped <- character(0)

  # iterate over all target stations that have metadata
  for (s in ids) {
    # distances from this target station to everyone else
    d <- D[s, ]
    # all neighbours within the radius, nearest first
    nb <- names(sort(d[!is.na(d) & d <= radius_m]))
    # not enough for the "extreme against ALL of them" logic: remember and move on
    if (length(nb) < n_neighbours) {
      skipped <- c(skipped, sprintf("%s(%d)", s, length(nb)))
      next
    }
    # exactly the required number, the nearest ones
    nb <- nb[seq_len(n_neighbours)]
    # extract data vector of this station (already cleaned by the earlier levels)
    v <- X[, s]

    # target minus each neighbour, one column per neighbour
    dif_nb <- sweep(X[, nb, drop = FALSE], 1, v, FUN = function(a, b) b - a)
    # per-neighbour upper bound of that difference distribution
    up <- apply(dif_nb, 2, function(z) stats::quantile(z, threshold, na.rm = TRUE))
    # per-neighbour lower bound
    lo <- apply(dif_nb, 2, function(z) stats::quantile(z, 1 - threshold, na.rm = TRUE))
    # extreme against that neighbour, in either tail
    ex <- sweep(dif_nb, 2, up, ">=") | sweep(dif_nb, 2, lo, "<=")
    # how many neighbour comparisons exist per time step
    n_avail <- rowSums(!is.na(dif_nb))
    # how many of them are extreme
    n_ex <- rowSums(ex, na.rm = TRUE)
    # spatial criterion: the value exists, the FULL neighbourhood is available,
    # and the value is extreme against ALL of it
    spatial <- !is.na(v) & n_avail >= n_neighbours & n_ex >= n_avail

    # change from the previous step
    d_pre <- v - shift_back(v)
    # change to the next step
    d_nex <- shift_fwd(v) - v
    # symmetric threshold on the ABSOLUTE change, past branch (see header)
    q_pre <- stats::quantile(abs(d_pre), threshold, na.rm = TRUE)
    # symmetric threshold, future branch
    q_nex <- stats::quantile(abs(d_nex), threshold, na.rm = TRUE)
    # temporal criterion: both changes exist and both are extreme
    temporal <- !is.na(d_pre) & !is.na(d_nex) & abs(d_pre) > q_pre & abs(d_nex) > q_nex

    # extreme in space AND in time: the objection condition
    hit <- spatial & temporal
    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the values so later levels never see them
      X[mask, s] <- NA
      # record this level's code (7 = level 7, fixed by convention)
      previous_flag[mask, s] <- 7
      # add the number of new flags to the counters
      n_station[s] <- n_found
      n_total <- n_total + n_found
    }
  }

  # write the matrices back into the xts shells, keeping index and column names
  if (n_total > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip - one summary line, no spam
  if (isTRUE(verbose)) {
    message(sprintf("T7 spatiotemporal (q=%g, %d neighbours <= %g m): %d flagged",
                    threshold, n_neighbours, radius_m, n_total))
    if (length(skipped) > 0)
      message("  skipped (too few neighbours in radius): ", paste(skipped, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                   <- x
  input$qc_data_flagged           <- flg
  input$qc_info$t7_spatiotemporal <- list(n_flagged            = n_total,
                                          n_flagged_by_station = n_station,
                                          radius_m             = radius_m,
                                          n_neighbours         = n_neighbours,
                                          threshold            = threshold,
                                          skipped              = skipped)
  input
}
