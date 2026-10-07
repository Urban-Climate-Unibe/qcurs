#' Humidity QC Level 3: isolated spikes
#'
#' Same triple condition as temperature level 3, via the shared
#' `qc_find_spikes()`: a value is objected to only if it deviates by at least
#' `threshold` from the median of the surrounding window AND from the values
#' before it AND from the values after it.
#' Thunderstorm outflow moves the humidity by 30 percent within half an hour -
#' but as a step, not as a single isolated point, so the triple condition
#' preserves such fronts.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param dt Half window in time steps.
#' @param threshold Minimum deviation in percent relative humidity.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_3_spike(res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_3_spike <- function(input,
                          dt = 3,
                          threshold = 20,
                          verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity", level = "rh3_spike")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # the half window must be a positive whole number of steps
  if (!is.numeric(dt) || length(dt) != 1 || dt < 1 || dt != round(dt))
    stop("dt must be a positive whole number of time steps.")
  # the threshold must be a positive humidity difference
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0)
    stop("threshold must be a positive deviation in percent.")
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 3
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: how many cells the test was applied to
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # the spike search itself lives in qc_find_spikes(), shared with T_QC_3
    sp <- qc_find_spikes(v, dt, threshold)
    n_judged[s] <- sum(sp$judged)

    # combine the verdict with THIS station's column only
    mask <- sp$hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)
    
    # apply only if something was found
    if (n_found > 0) {
      # blank the spikes so later levels never see them
      X[mask, s] <- NA
      # record this level's code (3 = level 3, fixed by convention)
      previous_flag[mask, s] <- 3
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
    message(sprintf("RH3 spike (>= %g%% in +/-%d steps): %d flagged",
                    threshold, dt, n_total))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data           <- x
  input$qc_data_flagged   <- flg
  input$qc_info$rh3_spike <- list(n_flagged            = n_total,
                                  n_flagged_by_station = n_station,
                                  n_judged_by_station  = n_judged,
                                  dt                   = dt,
                                  threshold            = threshold)
  input
}
