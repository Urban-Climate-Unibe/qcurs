#' Humidity QC Level 7: dewpoint consistency with the neighbouring stations
#'
#' The spatial comparison runs on the DEWPOINT, not on relative humidity: RH
#' is spatially noisy because it carries the temperature field with it, while
#' the dewpoint is a smooth tracer of the air mass. Dewpoint after Magnus with
#' the Alduchov & Eskridge (1996) coefficients. Geometry and test are the
#' shared `qc_neighbour_weights()` and `qc_find_spatial_outliers()` of
#' temperature level 6, with all its fixes. A "Td <= T" check is deliberately
#' absent: with Td derived from T and RH it holds automatically for RH <= 100.
#'
#' The temperature is matched cell by cell, by exact time stamp and logger
#' name (`qc_match_grid()`); a cell without a temperature twin has no dewpoint
#' and is not judged. Objections land on the HUMIDITY cells: the dewpoint is
#' only the yardstick. Without a usable temperature, without metadata, or
#' with fewer than three matched stations the level says so, records it, and
#' hands the pair back unchanged so the chain can go on.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature The temperature run: its result list (preferred; its
#'   CLEANED `qc_data` is taken) or a temperature xts with logger names.
#'   NULL (default) skips the level.
#' @param metadata Data frame with columns ID, LAT, LON, Landuse. NULL
#'   (default) skips the level.
#' @param k_sigma Multiplier k on the weighted neighbour spread of the dewpoint.
#' @param abs_floor Absolute floor in Kelvin of dewpoint deviation.
#' @param radius_m Maximum neighbour distance in metres.
#' @param k_neighbours At most this many nearest COMPATIBLE neighbours.
#' @param min_neighbours Minimum valid neighbours to judge at all.
#' @param landuse_mode "graded" (same class 1.0, vegetated/forest pairs 0.4,
#'   else 0) or "strict" (same class only).
#' @param verbose Report the tally and the never-evaluated stations.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_7_dewpoint_consistency(res, temperature = t_res, metadata = meta)
#' }
#'
#' @import xts
#' @import zoo
#' @import geosphere
#' @export
RH_QC_7_dewpoint_consistency <- function(input,
                                         temperature = NULL,
                                         metadata = NULL,
                                         k_sigma = 6,
                                         abs_floor = 3,
                                         radius_m = 3000,
                                         k_neighbours = 5,
                                         min_neighbours = 2,
                                         landuse_mode = c("graded", "strict"),
                                         verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh7_dewpoint_consistency")
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

  #-------------------------------------------------------------------------------
  # what this level needs from outside: a temperature twin and the geometry.
  # Anything unusable is not an error but a skip: say why, record it, hand the
  # pair back unchanged, let the chain go on

  skip <- function(reason) {
    if (isTRUE(verbose))
      message("RH7 dewpoint consistency: ", reason, " - level skipped, continuing with the next level.")
    input$qc_info$rh7_dewpoint_consistency <- list(n_flagged = 0L, skipped = TRUE, reason = reason)
    input
  }
  if (is.null(temperature)) return(skip("no temperature supplied"))
  if (is.null(metadata))    return(skip("no metadata supplied"))
  # the whole chain result: take its CLEANED series
  tx <- temperature
  if (is.list(tx) && !inherits(tx, "xts")) {
    if (!"qc_data" %in% names(tx)) return(skip("temperature list has no 'qc_data'"))
    tx <- tx$qc_data
  }
  # temperature on the humidity grid, by exact time stamp and logger name
  tm <- qc_match_grid(x, tx)
  if (is.null(tm)) return(skip("temperature has no logger name or time stamp in common with the humidity"))
  # the network: stations with humidity, temperature AND metadata
  net <- qc_neighbour_weights(metadata, tm$loggers, radius_m, k_neighbours, landuse_mode)
  # a consensus over fewer than 3 stations is not a network statement
  if (length(net$ids) < 3)
    return(skip(sprintf("only %d station(s) have humidity, temperature and metadata; need at least 3",
                        length(net$ids))))

  #-------------------------------------------------------------------------------
  # the dewpoint field the test runs on

  # Magnus formula with the Alduchov & Eskridge (1996) coefficients
  dewpoint <- function(t_c, rh_pct) {
    a <- 17.625; b <- 243.04
    # clamp into the defined domain of the logarithm
    r <- pmin(pmax(rh_pct, 0.1), 100)
    g <- log(r / 100) + (a * t_c) / (b + t_c)
    (b * g) / (a - g)
  }

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 7

  # plain numeric matrix of the humidity values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # dewpoint per cell of the network, from the CURRENT (cleaned) humidity and
  # the matched temperature; NA wherever either is missing
  TD <- dewpoint(tm$M[, net$ids, drop = FALSE], X[, net$ids, drop = FALSE])
  # the test itself lives in qc_find_spatial_outliers(), shared with T_QC_6
  res <- qc_find_spatial_outliers(TD, net$W, k_sigma, abs_floor, min_neighbours)
  # the verdict on the full grid: stations outside the network are never hit
  hit <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
  hit[, net$ids] <- res$hit

  # only cells that carry no earlier flag - the objection lands on the HUMIDITY
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
    # blank the objected humidity values so later levels never see them
    X[mask] <- NA
    # record this level's code (7 = level 7, fixed by convention)
    previous_flag[mask] <- 7
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # the blind spots: stations the test never reached, including those outside the network
  never_evaluated <- names(n_judged)[n_judged == 0]
  if (isTRUE(verbose)) {
    message(sprintf("RH7 dewpoint consistency (max(%g*sigma, %g K), k=%d after landuse): %d flagged, %d of %d loggers and %d of %d time steps with temperature",
                    k_sigma, abs_floor, k_neighbours, n_total,
                    length(tm$loggers), ncol(X), tm$n_time, nrow(X)))
    if (length(never_evaluated))
      message("  never evaluated (no temperature, no metadata, or insufficient compatible neighbours): ",
              paste(never_evaluated, collapse = ", "))
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                          <- x
  input$qc_data_flagged                  <- flg
  input$qc_info$rh7_dewpoint_consistency <- list(n_flagged            = n_total,
                                                 n_flagged_by_station = n_station,
                                                 n_judged_by_station  = n_judged,
                                                 sigma_d              = net$sigma_d,
                                                 k_sigma              = k_sigma,
                                                 abs_floor            = abs_floor,
                                                 radius_m             = radius_m,
                                                 k_neighbours         = k_neighbours,
                                                 min_neighbours       = min_neighbours,
                                                 landuse_mode         = landuse_mode,
                                                 loggers_with_temperature = tm$loggers,
                                                 n_time_with_temperature  = tm$n_time,
                                                 never_evaluated      = never_evaluated)
  input
}
