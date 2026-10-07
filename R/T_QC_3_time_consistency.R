#' Temperature QC Level 3: time consistency (isolated spikes)
#'
#' A value is removed only if it deviates by at least `diff` from the median of
#' the surrounding window AND from the values before it AND from the values
#' after it.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param dt Half window in time steps (default is 3 which means +/-30 min at
#'   10-minute data).
#' @param threshold Minimum deviation in Kelvin (default is 6 K).
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_3_time_consistency(res)   # chained after level 2
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_3_time_consistency <- function(input,
                                    dt = 3,
                                    threshold = 6,
                                    verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t3_time_consistency")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the half window must be a positive whole number of steps
  if (!is.numeric(dt) || length(dt) != 1 || dt < 1 || dt != round(dt))
    stop("dt must be a positive whole number of time steps.")
  # the threshold must be a positive temperature difference
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0)
    stop("threshold must be a positive deviation in Kelvin.")

  #-------------------------------------------------------------------------------
  # Perform QC Level 3

  # plain numeric matrix of the values (time in rows, stations in columns).
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
    # the spike search itself lives in qc_find_spikes(), shared with RH_QC_3
    sp <- qc_find_spikes(v, dt, threshold)
    n_judged[s] <- sum(sp$judged)

    # combine the verdict with THIS station's column only
    mask <- sp$hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station.
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
    message(sprintf("T3 time consistency (spike >= %g K, +/-%d steps): %d flagged",
                    threshold, dt, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                     <- x
  input$qc_data_flagged             <- flg
  input$qc_info$t3_time_consistency <- list(n_flagged            = n_total,
                                            n_flagged_by_station = n_station,
                                            n_judged_by_station  = n_judged,
                                            dt                   = dt,
                                            threshold            = threshold)
  input
}
