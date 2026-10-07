#' Humidity QC Level 8 (optional): linear gap interpolation
#'
#' NOT a quality test - this level ADDS values. It linearly interpolates gaps
#' of at most `maxgap` steps, via the shared `qc_fill_gaps()`. Keep it
#' optional and OFF by default in any runner, because the underlying method
#' paper states that QC only flags and never alters the original data;
#' whoever publishes the interpolated series must say so and must document
#' the codes.
#'
#' Flag codes written here, fixed by the convention: 50 for a filled plain
#' gap, 50 + N for a value that level N had removed and this level refilled -
#' 51 therefore means "inherited from the temperature chain, then
#' refilled", 53 "spike, refilled", and so on.
#'
#' @param input xts of humidity, or list from a previous QC level. A bare
#'   xts gets a fresh flag matrix, so every fill is then a plain gap fill;
#'   the refill codes only appear when the chain ran before.
#' @param maxgap Longest gap (in steps) that is filled (5 = 50 min at
#'   10-minute data).
#' @param refill_flagged Also fill cells removed by the QC levels. If FALSE,
#'   QC removals stay NA and keep their original codes.
#' @param verbose Report what was ADDED.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_8_interpolate(res, maxgap = 5)   # only on explicit request
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_8_interpolate <- function(input,
                                maxgap = 5,
                                refill_flagged = TRUE,
                                verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh8_interpolate")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the gap bound must be a positive whole number of steps
  if (!is.numeric(maxgap) || length(maxgap) != 1 || maxgap < 1 || maxgap != round(maxgap))
    stop("maxgap must be a positive whole number of steps.")

  #-------------------------------------------------------------------------------
  # Perform RH8 - this level creates values, it does not test

  # plain numeric matrix of the values (QC removals are NA here)
  X <- coredata(x)
  # plain numeric matrix of the accumulated flags
  previous_flag <- coredata(flg)
  # the filling itself lives in qc_fill_gaps(), shared with T_QC_9_interpolate
  res <- qc_fill_gaps(X, previous_flag, maxgap, refill_flagged)

  # write the matrices back into the xts shells, keeping index and column names
  if (res$n_gap + res$n_ref > 0) {
    x[]   <- res$X
    flg[] <- res$previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # honest report of what was ADDED
  if (isTRUE(verbose))
    message(sprintf("RH8 interpolation (maxgap %d): %d gaps filled (code 50), %d QC removals refilled (codes 50+level)",
                    maxgap, res$n_gap, res$n_ref))

  # write the updated matrices back and append this level under its own name
  input$qc_data                 <- x
  input$qc_data_flagged         <- flg
  input$qc_info$rh8_interpolate <- list(n_gap_filled   = res$n_gap,
                                      n_refilled     = res$n_ref,
                                      maxgap         = maxgap,
                                      refill_flagged = refill_flagged)
  input
}
