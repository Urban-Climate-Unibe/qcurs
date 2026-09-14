#' Humidity QC Level 8 (optional): linear gap interpolation
#'
#' NOT a quality test - this level ADDS values. It linearly interpolates gaps
#' of at most `maxgap` steps. Keep it optional and OFF by default in any
#' runner, because the underlying method paper states that QC only flags and
#' never alters the original data; whoever publishes the interpolated series
#' must say so and must document the codes.
#'
#' Flag codes written here: `base_code` (50) for a filled plain gap, and
#' `base_code + N` for a value that RH level N had removed and this level
#' refilled - 51 therefore means "inherited from the temperature chain, then
#' refilled", 53 "spike, refilled", and so on.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#'   A bare xts gets a fresh flag matrix, so every fill is then a plain gap
#'   fill; the refill codes only appear when the chain ran before.
#' @param maxgap Longest gap (in steps) that is filled (5 = 50 min at
#'   10-minute data).
#' @param refill_flagged Also fill cells removed by the QC levels. If FALSE,
#'   QC removals stay NA and keep their original codes.
#' @param base_code Code for a plain gap fill; refills get base_code + level.
#'   Must stay clear of the QC level codes (default 50).
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
                                base_code = 50,
                                verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely
  
  # the gap bound must be a positive number of steps
  if (!is.numeric(maxgap) || maxgap < 1) stop("maxgap must be at least 1 step.")
  # the base code must sit above every QC level code, or refills collide with tests
  if (!is.numeric(base_code) || base_code <= 10)
    stop("base_code must be larger than 10 to stay clear of the QC level codes.")
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 8
  
  # plain numeric matrix of the values (QC removals are NA here)
  X <- coredata(x)
  # plain numeric matrix of the accumulated flags
  previous_flag <- coredata(flg)
  # tallies: filled plain gaps and refilled QC removals
  n_gap <- 0
  n_ref <- 0
  
  # iterate over all stations (columns)
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # this station's flags
    f <- previous_flag[, s]
    # everything that is currently missing (true gaps AND QC removals)
    was_na <- is.na(v)
    if (!refill_flagged) {
      # optionally protect QC removals from being refilled
      keep_na <- was_na & !is.na(f) & f > 0
    } else {
      # or allow everything to be considered
      keep_na <- rep(FALSE, length(v))
    }
    # linear interpolation with a bounded gap length; ends are never extrapolated
    vi <- zoo::na.approx(v, na.rm = FALSE, maxgap = maxgap)
    # undo fills on protected cells
    vi[keep_na] <- NA
    # cells that actually received a value
    filled <- was_na & !is.na(vi)
    # nothing filled at this station
    if (!any(filled)) next
    # write the interpolated values
    X[filled, s] <- vi[filled]
    # fills over former QC removals get base_code + the original level number
    qc_hit <- filled & !is.na(f) & f > 0
    # fills over plain gaps get the base code
    plain <- filled & (is.na(f) | f == 0)
    previous_flag[qc_hit, s] <- base_code + f[qc_hit]
    previous_flag[plain, s]  <- base_code
    # add to the tallies
    n_ref <- n_ref + sum(qc_hit)
    n_gap <- n_gap + sum(plain)
  }
  
  # write the matrices back into the xts shells, keeping index and column names
  if (n_gap + n_ref > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }
  
  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level
  
  # honest report of what was ADDED - this level creates values, it does not test
  if (isTRUE(verbose))
    message(sprintf("RH8 interpolation (maxgap %d): %d gaps filled (code %d), %d QC removals refilled (codes %d+level)",
                    maxgap, n_gap, base_code, n_ref, base_code))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                 <- x
  input$qc_data_flagged         <- flg
  input$qc_info$rh8_interpolate <- list(n_gap_filled   = n_gap,
                                        n_refilled     = n_ref,
                                        maxgap         = maxgap,
                                        refill_flagged = refill_flagged,
                                        base_code      = base_code)
  input
}
