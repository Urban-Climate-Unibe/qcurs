#' Temperature QC Level 5: climatic outliers
#'
#' Removes values beyond Q1 - m*IQR / Q3 + m*IQR of the station's own
#' distribution per calendar month, all years pooled.
#'
#' Two guards apply: a real minimum sample size, and - by default -
#' the requirement that a month's sample spans more than one calendar year.
#' set `require_multi_year = FALSE` to override.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param ext_lim_factor Multiplier m on the IQR (paper: 4, empirical).
#' @param min_n Minimum valid values a station-month needs to be judged.
#' @param require_multi_year Only judge a station-month whose sample spans
#'   more than one calendar year.
#' @param flag_code QC-code written by this level. Here, default is 5.
#' @param verbose Report the tally and the refusals.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. The
#'   refused station-months are listed in
#'   `qc_info$t5_climatic_outliers$refused`.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_5_climatic_outliers(res)                             # honest self-skip
#' res <- T_QC_5_climatic_outliers(res, require_multi_year = FALSE) # force it anyway
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_5_climatic_outliers <- function(input,
                                     ext_lim_factor = 4,
                                     min_n = 300,
                                     require_multi_year = TRUE,
                                     flag_code = 5,
                                     verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the IQR multiplier must be positive
  if (!is.numeric(ext_lim_factor) || length(ext_lim_factor) != 1 || ext_lim_factor <= 0)
    stop("ext_lim_factor must be a positive number.")
  # a quartile estimate from fewer than ~30 values is noise, not a climatology
  if (!is.numeric(min_n) || length(min_n) != 1 || min_n < 30)
    stop("min_n must be at least 30.")

  #-------------------------------------------------------------------------------
  # Perform QC Level 5

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # calendar month of every time step
  mon <- as.numeric(format(zoo::index(x), "%m"))
  # calendar year of every time step, for the multi-year requirement
  yr <- as.numeric(format(zoo::index(x), "%Y"))
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # which station-months were refused, so a zero-flag run is not mistaken for a clean one
  refused <- character(0)

  # iterate over all stations (columns): a climatology is site-specific
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # every calendar month present in the record
    for (mm in sort(unique(mon))) {
      # this station's valid values in this month, all years pooled
      sel <- mon == mm & !is.na(v)
      # sample size of the pooled month
      nv <- sum(sel)
      # over how many distinct years the sample spreads
      years_span <- length(unique(yr[sel]))
      # guards: enough data AND (if required) a genuine multi-year climatology
      if (nv < min_n || (require_multi_year && years_span < 2)) {
        # remember the refusal instead of silently doing nothing
        if (nv > 0) refused <- c(refused, sprintf("%s/%02d(n=%d,y=%d)", s, mm, nv, years_span))
        next
      }
      # quartiles of the pooled monthly sample
      q <- stats::quantile(v[sel], c(0.25, 0.75))
      # interquartile range
      iqr <- q[2] - q[1]
      # lower and upper climatic bound
      lo <- q[1] - ext_lim_factor * iqr
      hi <- q[2] + ext_lim_factor * iqr
      # strictly outside ("exceed", per the paper wording), this month only
      hit <- mon == mm & !is.na(v) & (v < lo | v > hi)

      # combine the verdict with THIS station's column only
      mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
      # how many cells this level objects to in this station-month
      n_found <- sum(mask)

      # apply only if something was found
      if (n_found > 0) {
        # blank the outliers so later levels never see them
        X[mask, s] <- NA
        # record this level's code
        previous_flag[mask, s] <- flag_code
        # add the number of new flags to the counters
        n_station[s] <- n_station[s] + n_found
        n_total <- n_total + n_found
      }
    }
  }

  # write the matrices back into the xts shells, keeping index and column names
  if (n_total > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip - and name the refusals
  if (isTRUE(verbose)) {
    message(sprintf("T5 climatic outliers (m=%g): %d flagged, %d station-months refused",
                    ext_lim_factor, n_total, length(refused)))
    if (length(refused) > 0 && length(refused) <= 12)
      message("  refused: ", paste(refused, collapse = ", "))
    else if (length(refused) > 12)
      message("  refused (first 12): ", paste(utils::head(refused, 12), collapse = ", "), ", ...")
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                      <- x
  input$qc_data_flagged              <- flg
  input$qc_info$t5_climatic_outliers <- list(n_flagged            = n_total,
                                             n_flagged_by_station = n_station,
                                             ext_lim_factor       = ext_lim_factor,
                                             min_n                = min_n,
                                             require_multi_year   = require_multi_year,
                                             refused              = refused)
  input
}
