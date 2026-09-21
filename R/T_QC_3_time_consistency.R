#' Temperature QC Level 3: time consistency (isolated spikes)
#'
#' A value is removed only if it deviates by at least `diff` from the median of
#' the surrounding window AND from the values before it AND from the values
#' after it.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param dt Half window in time steps (default is 3 which means +/-30 min at
#'   10-minute data).
#' @param diff Threshold in Kelvin (default is 6 K).
#' @param flag_code QC-code written by this level. Here, default is 3.
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
                                    diff = 6,
                                    flag_code = 3,
                                    verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the half window must be a positive whole number of steps
  if (!is.numeric(dt) || length(dt) != 1 || dt < 1 || dt != round(dt))
    stop("dt must be a positive whole number of time steps.")
  # the threshold must be a positive temperature difference
  if (!is.numeric(diff) || length(diff) != 1 || diff <= 0)
    stop("diff must be a positive threshold in Kelvin.")

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

  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station and determine its length
    v <- X[, s]; n <- length(v)
    # create a vector with length n. All entries are FALSE
    hit <- rep(FALSE, n)
    # generate a sequence-vector of length n and iterate over it
    for (i in seq_len(n)) {
      # skip if there is a NA or at the ends of the vector
      if (is.na(v[i]) || i == 1 || i == n) next
      # create the window around i and ensure it respects the edges
      nb <- max(1, i - dt):min(n, i + dt)
      # if the window holds fewer than 1 + dt values, there is too little
      # context to judge --> refuse
      if (sum(!is.na(v[nb])) < (1 + dt)) next
      # condition 1: deviates at least diff from the window median
      if (abs(v[i] - stats::median(v[nb], na.rm = TRUE)) < diff) next
      # extract all values before i and remove NAs
      bef <- v[max(1, i - dt):(i - 1)]; bef <- bef[!is.na(bef)]
      # condition 2: deviates at least diff from previous values
      if (!length(bef) || abs(v[i] - stats::median(bef)) < diff) next
      # extract all values after i and remove NAs
      aft <- v[(i + 1):min(n, i + dt)]; aft <- aft[!is.na(aft)]
      # condition 3: deviates at least diff from future values
      if (!length(aft) || abs(v[i] - stats::median(aft)) < diff) next
      # if all three conditions met: an isolated spike has been found
      hit[i] <- TRUE
    }

    # creatte a mask based on the conditions
    mask <- hit & !is.na(v) & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station.
    n_found <- sum(mask)

    # apply only if something was found
    if (n_found > 0) {
      # blank the spikes so later levels never see them
      X[mask, s] <- NA
      # record this level's code
      previous_flag[mask, s] <- flag_code
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
                    diff, dt, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                     <- x
  input$qc_data_flagged             <- flg
  input$qc_info$t3_time_consistency <- list(n_flagged            = n_total,
                                            n_flagged_by_station = n_station,
                                            dt                   = dt,
                                            diff                 = diff)
  input
}
