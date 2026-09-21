#' Temperature QC Level 9 (optional): linear gap interpolation
#'
#' this level ADDS values. It linearly interpolates gaps
#' of at most `maxgap` steps.
#'
#' Flag codes written here: `base_code` (50) for a filled plain gap, and
#' `base_code + N` for a value that level N had removed and this level
#' refilled. (The published pipeline said "Flag 11-17" in its message while
#' writing 81-87; this one writes what it says.)
#'
#' @param input xts of temperature, or list(qc_data, qc_data_flagged, qc_info).
#'   A bare xts gets a fresh flag matrix, so every fill is then a plain gap
#'   fill; the refill codes only appear when the chain ran before.
#' @param maxgap Longest gap (in steps) that is filled (5 = 50 min at 10-min
#'   data).
#' @param refill_flagged Also fill cells removed by the QC levels. If FALSE,
#'   QC removals stay NA and keep their original codes.
#' @param base_code Code for a plain gap fill; refills get base_code + level.
#'   Must leave room above the highest QC code (default 50 > 8).
#' @param verbose Report what was ADDED.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_9_interpolate(res, maxgap = 5)   # only on explicit request
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_9_interpolate <- function(input, maxgap = 5, refill_flagged = TRUE,
                               base_code = 50, verbose = TRUE) {

#-------------------------------------------------------------------------------
# normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
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
# work on the plain matrices, station by station

  # plain numeric matrix of the values (QC removals are NA here)
  X <- coredata(x)
  # plain numeric matrix of the accumulated flags
  F <- coredata(flg)
  # tallies: filled plain gaps and refilled QC removals
  n_gap <- 0
  n_ref <- 0

  for (s in colnames(X)) {
    # this station's values as a plain vector
    v <- X[, s]
    # this station's flags
    f <- F[, s]
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
    F[qc_hit, s] <- base_code + f[qc_hit]
    F[plain, s]  <- base_code
    # add to the tallies
    n_ref <- n_ref + sum(qc_hit)
    n_gap <- n_gap + sum(plain)
  }
  if (n_gap + n_ref > 0) {
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- F
  }

#-------------------------------------------------------------------------------
# report and hand the pair on to the next level

  # honest report of what was ADDED - this level creates values, it does not test
  if (isTRUE(verbose))
    message(sprintf("T9 interpolation (maxgap %d): %d gaps filled (code %d), %d QC removals refilled (codes %d+level)",
                    maxgap, n_gap, base_code, n_ref, base_code))

  # write the updated matrices back and append this level under its own name
  input$qc_data                 <- x
  input$qc_data_flagged         <- flg
  input$qc_info$t9_interpolate  <- list(n_gap_filled   = n_gap,
                                        n_refilled     = n_ref,
                                        maxgap         = maxgap,
                                        refill_flagged = refill_flagged,
                                        base_code      = base_code)
  input
}
