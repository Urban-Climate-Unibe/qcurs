#' Humidity QC Level 6: saturation drift
#'
#' The classic ageing failure of capacitive elements: fog or continuous rain
#' drives the whole network to saturation, and a drifted sensor stays well
#' below it while its values remain perfectly plausible. No other test sees
#' this. The criterion runs in segments of `window_days`: judged over the
#' whole record, the healthy first half of a campaign would acquit a sensor
#' that drifted in the second half. Objections land on the saturated time
#' steps of the offending station - where the evidence is - not across the
#' whole segment.
#'
#' Two properties to keep in mind. A station that reaches `station_max` even
#' once during a segment's saturated steps is acquitted for that whole
#' segment; if a drift starts mid-segment, shorten `window_days`. And the
#' network median includes the drifting stations themselves, so when a large
#' share of a small network drifts at once, the network no longer counts as
#' saturated and nothing is caught.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param network_min Network median above which a step counts as saturated.
#' @param station_max Ceiling the station must stay below, across ALL
#'   saturated steps of the segment, to count as a drifter.
#' @param min_hits Minimum saturated observations per segment to judge.
#' @param min_stations Minimum stations reporting for a network median.
#' @param window_days Segment length in days.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_6_saturation_drift(res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_6_saturation_drift <- function(input,
                                     network_min = 95,
                                     station_max = 90,
                                     min_hits = 12,
                                     min_stations = 3,
                                     window_days = 30,
                                     verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity", level = "rh6_saturation_drift")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # the station ceiling must lie below the network criterion, or the logic inverts
  if (!is.numeric(network_min) || !is.numeric(station_max) ||
      station_max >= network_min)
    stop("station_max must lie below network_min.")
  # a median of one station is not a network statement
  if (!is.numeric(min_stations) || min_stations < 2)
    stop("min_stations must be at least 2.")
  # at least one saturated observation is needed to judge at all
  if (!is.numeric(min_hits) || min_hits < 1)
    stop("min_hits must be at least 1.")
  # segments must have a positive length
  if (!is.numeric(window_days) || window_days <= 0)
    stop("window_days must be a positive number of days.")
  
  #-------------------------------------------------------------------------------
  # find the saturated time steps of the network, once
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # how many stations report at each time step
  valid_per_step <- rowSums(!is.na(X))
  # the network median at each time step
  net_med <- apply(X, 1, stats::median, na.rm = TRUE)
  # a step is "saturated" when enough stations report AND their median is at the top
  is_sat <- valid_per_step >= min_stations & !is.na(net_med) & net_med >= network_min
  # segment id of every time step (0, 1, 2, ... in window_days blocks)
  seg <- as.integer(floor(as.numeric(difftime(zoo::index(x), zoo::index(x)[1],
                                              units = "days")) / window_days))
  # the saturated time steps of every segment, once; segments without any are dropped
  sat_by_seg <- split(which(is_sat), seg[is_sat])

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 6

  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: values at saturated steps of station-segments with a verdict
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # verdict per cell, FALSE until proven otherwise
    hit <- rep(FALSE, length(v))
    # judge segment by segment: a drift is a property of a period, not of the record
    for (sat_steps in sat_by_seg) {
      # the station's values during the saturation episodes
      vs <- v[sat_steps]
      # the evidence cells: observed values at saturated steps
      obs <- sat_steps[!is.na(vs)]
      # too few saturated observations: no verdict for this station-segment
      if (length(obs) < min_hits) next
      # a verdict is reached, whichever way
      n_judged[s] <- n_judged[s] + length(obs)
      # the station does reach saturation at least once: not a drifter
      if (max(vs, na.rm = TRUE) >= station_max) next
      # a drifter: mark its observed values at the evidence steps
      hit[obs] <- TRUE
    }

    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)
    
    # apply only if something was found
    if (n_found > 0) {
      # blank the drifted values so later levels never see them
      X[mask, s] <- NA
      # record this level's code (6 = level 6, fixed by convention)
      previous_flag[mask, s] <- 6
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
  
  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("RH6 saturation drift (network >= %g / station < %g, %g-day segments): %d flagged",
                    network_min, station_max, window_days, n_total))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                      <- x
  input$qc_data_flagged              <- flg
  input$qc_info$rh6_saturation_drift <- list(n_flagged            = n_total,
                                             n_flagged_by_station = n_station,
                                             n_judged_by_station  = n_judged,
                                             network_min          = network_min,
                                             station_max          = station_max,
                                             min_hits             = min_hits,
                                             min_stations         = min_stations,
                                             window_days          = window_days)
  input
}
