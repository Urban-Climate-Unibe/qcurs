#' Temperature QC Level 6: spatial consistency (weighted neighbour consensus)
#'
#' Implements equations (1)-(4) of Amini et al. (2026): a value is removed when
#' it deviates from the Gaussian-distance- and landuse-weighted consensus of
#' its neighbours by more than max(k * sigma', delta). The network geometry
#' comes from the shared `qc_neighbour_weights()` (k-nearest cap AFTER the
#' landuse masking, both forest spellings) and the test itself from the shared
#' `qc_find_spatial_outliers()` (NO sigma gate), so both exist once for the
#' temperature and the dewpoint level.
#'
#' Stations that never reach the minimum number of compatible valid
#' neighbours are listed in `qc_info$t6_spatial_consistency$never_evaluated`,
#' and `n_judged_by_station` counts the cells the test actually reached.
#' A station this test cannot see is a coverage statement, not a clean bill.
#' With fewer than `min_neighbours + 1` stations in data AND metadata the
#' level skips itself, records why, and hands the pair back unchanged.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param metadata Data frame with columns ID, LAT, LON, Landuse. Stations are
#'   matched to data columns by ID.
#' @param k_sigma Multiplier k on the weighted neighbour spread (paper: 6).
#' @param abs_floor Absolute floor delta in Kelvin (paper: 3).
#' @param radius_m Maximum neighbour distance in metres (paper: 3000).
#' @param k_neighbours At most this many nearest COMPATIBLE neighbours (paper: 5).
#' @param min_neighbours Minimum valid neighbours to judge at all (paper: 2).
#' @param landuse_mode "graded" (paper: same class 1.0, vegetated/forest pairs
#'   0.4, else 0) or "strict" (same class only).
#' @param verbose Report the tally and the never-evaluated stations.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_6_spatial_consistency(res, metadata = meta_biel)
#' }
#'
#' @import xts
#' @import zoo
#' @import geosphere
#' @export
T_QC_6_spatial_consistency <- function(input,
                                       metadata,
                                       k_sigma = 6,
                                       abs_floor = 3,
                                       radius_m = 3000,
                                       k_neighbours = 5,
                                       min_neighbours = 2,
                                       landuse_mode = c("graded", "strict"),
                                       verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t6_spatial_consistency")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  landuse_mode <- match.arg(landuse_mode)

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  if (!is.numeric(k_sigma) || k_sigma <= 0) stop("k_sigma must be a positive number.")
  if (!is.numeric(abs_floor) || abs_floor <= 0) stop("abs_floor must be a positive number of Kelvin.")
  if (!is.numeric(radius_m) || radius_m <= 0) stop("radius_m must be a positive distance in metres.")
  if (!is.numeric(min_neighbours) || min_neighbours < 1) stop("min_neighbours must be at least 1.")
  if (!is.numeric(k_neighbours) || k_neighbours < min_neighbours)
    stop("k_neighbours must be at least min_neighbours.")
  if (missing(metadata) || is.null(metadata))
    stop("This level needs the metadata (ID, LAT, LON, Landuse).")

  #-------------------------------------------------------------------------------
  # the network geometry, once (metadata validation lives in the helper)

  net <- qc_neighbour_weights(metadata, colnames(x), radius_m, k_neighbours, landuse_mode)
  # without at least min_neighbours + 1 stations the condition can never be
  # met: not an error but a skip - say why, record it, hand the pair back
  if (length(net$ids) < min_neighbours + 1) {
    reason <- sprintf("only %d station(s) match data and metadata; need at least %d",
                      length(net$ids), min_neighbours + 1)
    if (isTRUE(verbose))
      message("T6 spatial consistency: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$t6_spatial_consistency <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    return(input)
  }

  #-------------------------------------------------------------------------------
  # Perform QC Level 6

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # the test itself lives in qc_find_spatial_outliers(), shared with RH_QC_7
  res <- qc_find_spatial_outliers(X[, net$ids, drop = FALSE], net$W, k_sigma, abs_floor, min_neighbours)
  # the verdict on the full grid: stations outside the network are never hit
  hit <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
  hit[, net$ids] <- res$hit

  # only cells that carry no earlier flag
  mask <- hit & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))
  # per-station coverage: cells the consensus test was applied to (0 outside the network)
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))
  n_judged[net$ids] <- colSums(res$judged)

  # apply only if something was found
  if (n_total > 0) {
    # blank the objected values so later levels never see them
    X[mask] <- NA
    # record this level's code (6 = level 6, fixed by convention)
    previous_flag[mask] <- 6
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # the blind spots: stations the test never reached, including those outside the network
  never_evaluated <- names(n_judged)[n_judged == 0]
  if (isTRUE(verbose)) {
    message(sprintf("T6 spatial consistency (max(%g*sigma, %g K), k=%d after landuse): %d flagged",
                    k_sigma, abs_floor, k_neighbours, n_total))
    if (length(never_evaluated))
      message("  never evaluated (insufficient compatible neighbours): ",
              paste(never_evaluated, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                        <- x
  input$qc_data_flagged                <- flg
  input$qc_info$t6_spatial_consistency <- list(n_flagged            = n_total,
                                               n_flagged_by_station = n_station,
                                               n_judged_by_station  = n_judged,
                                               sigma_d              = net$sigma_d,
                                               k_sigma              = k_sigma,
                                               abs_floor            = abs_floor,
                                               radius_m             = radius_m,
                                               k_neighbours         = k_neighbours,
                                               min_neighbours       = min_neighbours,
                                               landuse_mode         = landuse_mode,
                                               never_evaluated      = never_evaluated)
  input
}
