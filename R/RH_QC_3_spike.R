#' Humidity QC Level 3: isolated spikes
#'
#' Same triple condition as temperature level 3: a value is objected to only
#' if it deviates by at least `spike_diff` from the median of the surrounding
#' window AND from the values before it AND from the values after it.
#' Thunderstorm outflow moves the humidity by 30 percent within half an hour -
#' but as a step, not as a single isolated point, so the triple condition
#' preserves such fronts.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param spike_dt Half window in time steps.
#' @param spike_diff Threshold in percent relative humidity.
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
                          spike_dt = 3,
                          spike_diff = 20,
                          verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # the half window must be a positive whole number of steps
  if (!is.numeric(spike_dt) || length(spike_dt) != 1 || spike_dt < 1 || spike_dt != round(spike_dt))
    stop("spike_dt must be a positive whole number of time steps.")
  # the threshold must be a positive humidity difference
  if (!is.numeric(spike_diff) || length(spike_diff) != 1 || spike_diff <= 0)
    stop("spike_diff must be a positive threshold in percent.")
  
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
  
  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station and determine its length
    v <- X[, s]; n <- length(v)
    # create a vector with length n. All entries are FALSE
    hit <- rep(FALSE, n)
    # walk every time step
    for (i in seq_len(n)) {
      # skip if there is a NA or at the ends of the vector
      if (is.na(v[i]) || i == 1 || i == n) next
      # create the window around i and ensure it respects the edges
      nb <- max(1, i - spike_dt):min(n, i + spike_dt)
      # too little context to judge: refuse instead of guessing
      if (sum(!is.na(v[nb])) < (1 + spike_dt)) next
      # condition 1: deviates at least spike_diff from the window median
      if (abs(v[i] - stats::median(v[nb], na.rm = TRUE)) < spike_diff) next
      # extract all values before i and remove NAs
      bef <- v[max(1, i - spike_dt):(i - 1)]; bef <- bef[!is.na(bef)]
      # condition 2: deviates at least spike_diff from previous values
      if (!length(bef) || abs(v[i] - stats::median(bef)) < spike_diff) next
      # extract all values after i and remove NAs
      aft <- v[(i + 1):min(n, i + spike_dt)]; aft <- aft[!is.na(aft)]
      # condition 3: deviates at least spike_diff from future values
      if (!length(aft) || abs(v[i] - stats::median(aft)) < spike_diff) next
      # all three conditions met: an isolated spike has been found
      hit[i] <- TRUE
    }
    
    # combine the verdict with THIS station's column only
    mask <- hit & !is.na(v) & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
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
                    spike_diff, spike_dt, n_total))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data           <- x
  input$qc_data_flagged   <- flg
  input$qc_info$rh3_spike <- list(n_flagged            = n_total,
                                  n_flagged_by_station = n_station,
                                  spike_dt             = spike_dt,
                                  spike_diff           = spike_diff)
  input
}
