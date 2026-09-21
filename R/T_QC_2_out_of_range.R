#' Temperature QC Level 2: seasonal out-of-range
#'
#' Flags values outside seasonal plausibility bounds, evaluated per time step
#' against the bound of the season that time step falls into.
#'
#' @param input xts of temperature, or list(qc_data, qc_data_flagged).
#' @param season_thresholds Named list winter/spring/summer/autumn, each
#'   list(min_val, max_val), in deg C. Values BELOW min_val or ABOVE max_val
#'   are flagged.
#' @param summer_min Convenience override for the most problematic bound;
#'   overwrites `season_thresholds$summer$min_val` when not NULL.
#' @param verbose Report the per-season tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_2_out_of_range(res)
#' res <- T_QC_2_out_of_range(res, summer_min = -10)   # colder site
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_2_out_of_range <- function(input,
                                season_thresholds = list(
                                  winter = list(min_val = -24.44, max_val = 23.54),
                                  spring = list(min_val = -21.05, max_val = 38.08),
                                  summer = list(min_val =  -5.85, max_val = 46.26),
                                  autumn = list(min_val =  -8.42, max_val = 40.87)),
                                summer_min = NULL,
                                verbose = TRUE) {

#-------------------------------------------------------------------------------
# normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

#-------------------------------------------------------------------------------
# validate the thresholds, because the caller may replace them freely

  seasons <- c("winter", "spring", "summer", "autumn")
  # ... all four seasons must be present
  if (!all(seasons %in% names(season_thresholds)))
    stop("season_thresholds must contain winter, spring, summer and autumn.")
  # ... the convenience override wins over the list entry
  if (!is.null(summer_min)) {
    if (!is.numeric(summer_min) || length(summer_min) != 1)
      stop("summer_min must be a single number.")
    season_thresholds$summer$min_val <- summer_min
  }
  for (s in seasons) {
    b <- season_thresholds[[s]]
    if (!all(c("min_val", "max_val") %in% names(b)))
      stop(sprintf("season_thresholds$%s must contain min_val and max_val.", s))
    if (!is.numeric(b$min_val) || !is.numeric(b$max_val))
      stop(sprintf("season_thresholds$%s bounds must be numeric.", s))
    if (b$min_val >= b$max_val)
      stop(sprintf("season_thresholds$%s: min_val must be below max_val.", s))
  }

#-------------------------------------------------------------------------------
# assign every time step to its season and build the bound vectors

  # extract the month number of every time step ...
  mon <- as.numeric(format(zoo::index(x), "%m"))
  # ... and assign it to its corresponding season (DJF, MAM, JJA, SON)
  season <- ifelse(mon %in% c(12, 1, 2), "winter",
            ifelse(mon %in% c(3, 4, 5), "spring",
            ifelse(mon %in% c(6, 7, 8), "summer", "autumn")))
  # lower bound valid at each time step (length = number of time steps)
  lo <- vapply(season, function(s) season_thresholds[[s]]$min_val, numeric(1), USE.NAMES = FALSE)
  # upper bound valid at each time step
  hi <- vapply(season, function(s) season_thresholds[[s]]$max_val, numeric(1), USE.NAMES = FALSE)

#-------------------------------------------------------------------------------
# work on the plain matrices: no xts recycling surprises, one pass over the data

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # the bounds repeated across all stations, so the comparison is cell by cell
  LO <- matrix(lo, nrow(X), ncol(X))
  HI <- matrix(hi, nrow(X), ncol(X))
  # cells that hold a value and violate their season's bound
  hit <- !is.na(X) & (X < LO | X > HI)
  # only cells that carry no earlier flag
  mask <- hit & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))
  # per-season tally: row(mask)[mask] gives the time index of every flagged cell
  n_season <- c(winter = 0L, spring = 0L, summer = 0L, autumn = 0L)

  # apply only if something was found
  if (n_total > 0) {
    tab <- table(factor(season[row(mask)[mask]], levels = seasons))
    n_season[seasons] <- as.integer(tab[seasons])
    # blank the objected values so later levels never see them
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
    message(sprintf("T2 out of range: %d flagged (winter %d, spring %d, summer %d, autumn %d)",
                    n_total, n_season["winter"], n_season["spring"],
                    n_season["summer"], n_season["autumn"]))

  # write the updated matrices back and append this level under its own name
  input$qc_data                <- x
  input$qc_data_flagged        <- flg
  input$qc_info$t2_out_of_range <- list(n_flagged            = n_total,
                                        n_flagged_by_station = n_station,
                                        n_flagged_by_season  = n_season,
                                        season_thresholds    = season_thresholds)
  input
}
