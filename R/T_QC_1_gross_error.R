#' Temperature QC Level 1: gross errors
#'
#' Removes physically impossible readings and logger fault codes (e.g. -45 or 999).
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param threshold_max Upper bound in deg C --> values ABOVE it are flagged.
#' @param threshold_min Lower bound in deg C --> values AT or BELOW it are
#'   flagged (inclusive, so the fault code -40.0 itself is caught).
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_1_gross_error(x)
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_1_gross_error <- function(input,
                               threshold_max = 60,
                               threshold_min = -40,
                               verbose = TRUE) {

  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the thresholds, because the caller may set them freely

  # both bounds must be single numbers
  if (!is.numeric(threshold_max) || length(threshold_max) != 1 ||
      !is.numeric(threshold_min) || length(threshold_min) != 1)
    stop("threshold_max and threshold_min must be single numbers.")
  # the lower bound must lie below the upper bound
  if (threshold_min >= threshold_max)
    stop("threshold_min must be below threshold_max.")

  #-------------------------------------------------------------------------------
  # Perform QC Level 1

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)

  # look for cells that hold a value and violate the bounds
  hit <- !is.na(X) & (X > threshold_max | X <= threshold_min)
  # only cells that carry no earlier flag
  mask <- hit & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))

  # apply only if something was found
  if (n_total > 0) {
    # blank the objected values so later levels never see them
    X[mask] <- NA
    # record this level's code (1 = level 1, fixed by convention)
    previous_flag[mask] <- 1
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("T1 gross error (> %g or <= %g): %d flagged",
                    threshold_max, threshold_min, n_total))

  # write the updated matrices back and append this level under its own name
  input$qc_data                <- x
  input$qc_data_flagged        <- flg
  input$qc_info$t1_gross_error <- list(n_flagged            = n_total,
                                       n_flagged_by_station = n_station,
                                       threshold_max        = threshold_max,
                                       threshold_min        = threshold_min)
  input
}
