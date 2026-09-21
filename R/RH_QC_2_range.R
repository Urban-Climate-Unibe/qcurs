#' Humidity QC Level 2: physical range
#'
#' Readings marginally above 100 percent are common in condensation because of
#' calibration tolerance and are physically meaningful, so values between 100
#' and `rh_max` are CLAMPED to 100 rather than discarded. Below `rh_min` and
#' above `rh_max` is impossible and objected to. Set `rh_max = 100` for a hard
#' cut without the tolerance zone.
#'
#' The clamping is an alteration of the data and is therefore reported and
#' recorded (`n_clamped`) with the same visibility as the objections.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param rh_min Lower physical bound in percent.
#' @param rh_max Upper tolerance bound in percent; values in (100, rh_max]
#'   are clamped to 100, values above are objected to.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_2_range(res)
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_2_range <- function(input,
                          rh_min = 0,
                          rh_max = 105,
                          verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # a ceiling below 100 would flag genuine saturation
  if (!is.numeric(rh_max) || rh_max < 100)
    stop("rh_max below 100 would flag genuine saturation; use rh_max >= 100.")
  # the floor must lie below saturation
  if (!is.numeric(rh_min) || rh_min < 0 || rh_min >= 100)
    stop("rh_min must lie in [0, 100).")
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 2
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  
  # the condensation-tolerance zone: clamp BEFORE flagging, so these survive as 100
  clamp <- !is.na(X) & X > 100 & X <= rh_max
  n_clamped <- sum(clamp)
  if (n_clamped > 0) X[clamp] <- 100
  
  # what remains impossible after clamping. Elementwise over the whole matrix -
  # there is no per-station computation, so no station loop is needed
  hit <- !is.na(X) & (X < rh_min | X > 100)
  # only cells that carry no earlier flag
  mask <- hit & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))
  
  # apply if anything was clamped OR objected to
  if (n_total > 0 || n_clamped > 0) {
    # blank the impossible values so later levels never see them
    X[mask] <- NA
    # record this level's code (2 = level 2, fixed by convention)
    previous_flag[mask] <- 2
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }
  
  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level
  
  # report BOTH actions: the clamping alters data and must be as visible as the flags
  if (isTRUE(verbose))
    message(sprintf("RH2 range %g-%g: %d flagged, %d clamped to 100",
                    rh_min, rh_max, n_total, n_clamped))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data           <- x
  input$qc_data_flagged   <- flg
  input$qc_info$rh2_range <- list(n_flagged            = n_total,
                                  n_flagged_by_station = n_station,
                                  n_clamped            = n_clamped,
                                  rh_min               = rh_min,
                                  rh_max               = rh_max)
  input
}
