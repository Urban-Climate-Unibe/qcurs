#' Temperature QC Level 1: gross errors
#'
#' Removes physically impossible readings and logger fault codes (e.g. -45 or
#' 999). Flag convention for the whole chain: 0 = checked and unobjected,
#' NA = no observation, N = objected by level N.
#'
#' @param input xts of temperature, or list(qc_data, qc_data_flagged, qc_info)
#'   from a previous level.
#' @param threshold_max Upper bound in deg C --> values ABOVE it are flagged.
#' @param threshold_min Lower bound in deg C --> values AT or BELOW it are
#'   flagged (inclusive, so the fault code -40.0 itself is caught).
#' @param flag_code Code written by this level. Default 1.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_1_gross_error(biel_xts)
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_1_gross_error <- function(input, threshold_max = 60, threshold_min = -40,
                               flag_code = 1, verbose = TRUE) {

#-------------------------------------------------------------------------------
# normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

#-------------------------------------------------------------------------------
# validate the thresholds, because the caller may set them freely

  # the lower bound must lie below the upper bound
  if (threshold_min >= threshold_max) stop("threshold_min must be below threshold_max.")

#-------------------------------------------------------------------------------
# work on the plain matrices: no xts recycling surprises, one pass over the data

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  F <- coredata(flg)
  # cells that hold a value, violate the bounds, and carry no earlier flag
  mask <- !is.na(X) & (X > threshold_max | X <= threshold_min) & (is.na(F) | F == 0)
  # how many cells this level objects to
  n <- sum(mask)
  if (n > 0) {
    # blank the objected values so later levels never see them
    X[mask] <- NA
    # record this level's code
    F[mask] <- flag_code
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- F
  }

#-------------------------------------------------------------------------------
# report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("T1 gross error (>%g or <=%g): %d flagged", threshold_max, threshold_min, n))

  # write the updated matrices back and append this level under its own name
  input$qc_data               <- x
  input$qc_data_flagged       <- flg
  input$qc_info$t1_gross_error <- list(n_flagged     = n,
                                       threshold_max = threshold_max,
                                       threshold_min = threshold_min)
  input
}
