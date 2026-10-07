#' Humidity QC Level 2: physical range
#'
#' Objects to values outside `[rh_min, rh_max]`. Readings a few percent above
#' 100 are common in condensation because of calibration tolerance and are
#' physically meaningful, so the default ceiling is 105 rather than 100. They
#' are kept as measured and unflagged - this level never clamps. Whoever wants
#' a hard cap at 100 for display or publication does that in the analysis
#' script.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param rh_min Lower physical bound in percent.
#' @param rh_max Upper bound in percent, including the condensation tolerance.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_2_range(res)
#' res <- RH_QC_2_range(res, rh_max = 100)   # no tolerance zone
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

  input <- qc_prepare_input(input, what = "humidity", level = "rh2_range")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # a ceiling below 100 would flag genuine saturation
  if (!is.numeric(rh_max) || length(rh_max) != 1 || rh_max < 100)
    stop("rh_max below 100 would flag genuine saturation; use rh_max >= 100.")
  # the floor must lie below saturation
  if (!is.numeric(rh_min) || length(rh_min) != 1 || rh_min < 0 || rh_min >= 100)
    stop("rh_min must lie in [0, 100).")

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 2

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)

  # cells that hold a value outside the physical range. Elementwise over the
  # whole matrix - there is no per-station computation, so no station loop
  hit <- !is.na(X) & (X < rh_min | X > rh_max)
  # only cells that carry no earlier flag
  mask <- hit & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))
  # per-station coverage: every cell that held a value was judged
  n_judged <- stats::setNames(colSums(!is.na(X)), colnames(X))

  # apply only if something was found
  if (n_total > 0) {
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

  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("RH2 range %g-%g: %d flagged", rh_min, rh_max, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data           <- x
  input$qc_data_flagged   <- flg
  input$qc_info$rh2_range <- list(n_flagged            = n_total,
                                  n_flagged_by_station = n_station,
                                  n_judged_by_station  = n_judged,
                                  rh_min               = rh_min,
                                  rh_max               = rh_max)
  input
}
