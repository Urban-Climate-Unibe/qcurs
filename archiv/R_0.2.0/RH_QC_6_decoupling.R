#' Humidity QC Level 6: decoupling from the station's own temperature
#'
#' Over a day with a normal diurnal cycle, relative humidity tracks its own
#' temperature closely and negatively. A humidity channel that stops following
#' its own temperature is broken even when its values look plausible - the
#' typical signature of a dead or detached element. This is the humidity twin
#' of temperature level 8: it needs no neighbours and therefore also covers
#' stations the spatial tests cannot reach.
#'
#' Days with less than `decouple_trange` of temperature amplitude are skipped,
#' because without a diurnal cycle the correlation is meaningless - fog days
#' would otherwise be flagged. The verdict is per day: all valid humidity
#' values of a decoupled day are objected to.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature The station temperatures: an xts with the same time
#'   index, or the whole result list of the temperature chain (its cleaned
#'   `qc_data` is taken). Stations are matched by name; humidity stations
#'   without a temperature twin are skipped and reported.
#' @param decouple_r Correlation at or above which a day counts as decoupled
#'   (healthy days are strongly negative).
#' @param decouple_trange Minimum diurnal temperature range in Kelvin.
#' @param decouple_min_n Minimum valid T/RH pairs per day.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. Stations
#'   without a temperature twin are listed in
#'   `qc_info$rh6_decoupling$skipped_stations`.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_6_decoupling(res, temperature = t_res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_6_decoupling <- function(input,
                               temperature,
                               decouple_r = -0.3,
                               decouple_trange = 3,
                               decouple_min_n = 60,
                               verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters and the temperature series
  
  # the decoupling threshold is a correlation
  if (!is.numeric(decouple_r) || decouple_r <= -1 || decouple_r >= 1)
    stop("decouple_r must be a correlation strictly between -1 and 1.")
  # a flat day carries no signal
  if (!is.numeric(decouple_trange) || decouple_trange <= 0)
    stop("decouple_trange must be a positive range in Kelvin.")
  # a correlation from a handful of pairs is noise
  if (!is.numeric(decouple_min_n) || decouple_min_n < 3)
    stop("decouple_min_n must be at least 3.")
  # the temperature series is required for this level
  if (missing(temperature) || is.null(temperature))
    stop("This level needs the temperature series.")
  # accept the whole temperature chain result: its CLEANED data are the best twin
  tx <- temperature
  if (is.list(tx) && !is.data.frame(tx) && !inherits(tx, "xts")) {
    if (!"qc_data" %in% names(tx)) stop("temperature list must contain 'qc_data'.")
    tx <- tx$qc_data
  }
  if (is.null(dim(tx))) stop("temperature must be an xts, a matrix, or a chain result list.")
  # the same dataset or nothing
  if (nrow(tx) != nrow(x))
    stop(sprintf("Row mismatch: humidity has %d time steps, temperature has %d.", nrow(x), nrow(tx)))
  if (inherits(tx, "xts") &&
      !isTRUE(all.equal(as.numeric(zoo::index(tx)), as.numeric(zoo::index(x)))))
    stop("Time index of the temperature series differs from the humidity series.")
  if (is.null(colnames(tx))) stop("temperature needs station names to be matched to the humidity columns.")
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 6
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # plain numeric matrix of the temperatures
  TT <- as.matrix(coredata(tx)); storage.mode(TT) <- "double"
  # calendar day of every time step
  day <- as.Date(zoo::index(x))
  # the day sequence of the record
  days <- unique(day)
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # humidity stations without a temperature twin cannot be tested - collect them
  skipped_stations <- setdiff(colnames(X), colnames(TT))
  
  # iterate over all stations (columns) that have a temperature twin
  for (s in intersect(colnames(X), colnames(TT))) {
    # extract humidity vector of this station and determine its length
    v <- X[, s]; n <- length(v)
    # its own temperature, matched by name
    tv <- TT[, s]
    # create a vector with length n. All entries are FALSE
    hit <- rep(FALSE, n)
    # judge day by day: the diurnal cycle is the yardstick
    for (dd in days) {
      # the day's valid T/RH pairs
      k <- which(day == dd & !is.na(tv) & !is.na(v))
      # too sparse: no verdict for this day
      if (length(k) < decouple_min_n) next
      # flat day: the correlation is meaningless, skip instead of flagging fog
      if (diff(range(tv[k])) < decouple_trange) next
      # the daily T-RH correlation
      r <- suppressWarnings(stats::cor(tv[k], v[k]))
      # healthy negative coupling (or an undefined correlation): fine
      if (is.na(r) || r < decouple_r) next
      # decoupled: mark every valid humidity value of this day
      hit[k] <- TRUE
    }
    
    # combine the verdict with THIS station's column only
    mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
    # how many cells this level objects to at this station
    n_found <- sum(mask)
    
    # apply only if something was found
    if (n_found > 0) {
      # blank the decoupled days so later levels never see them
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
  
  # report so a zero-hit run is visibly a run, not a skip - and name the blind spots
  if (isTRUE(verbose)) {
    message(sprintf("RH6 decoupled from own temperature (r >= %g): %d flagged", decouple_r, n_total))
    if (length(skipped_stations) > 0)
      message("  skipped (no temperature twin): ", paste(skipped_stations, collapse = ", "))
  }
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                <- x
  input$qc_data_flagged        <- flg
  input$qc_info$rh6_decoupling <- list(n_flagged            = n_total,
                                       n_flagged_by_station = n_station,
                                       decouple_r           = decouple_r,
                                       decouple_trange      = decouple_trange,
                                       decouple_min_n       = decouple_min_n,
                                       skipped_stations     = skipped_stations)
  input
}
